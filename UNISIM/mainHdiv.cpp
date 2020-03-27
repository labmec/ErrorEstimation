#include "Tools.h"

#include <Material/REAL/mixedpoisson.h>
#include <Material/TPZVecL2.h>
#include <Material/pzbndcond.h>
#include <Matrix/pzstepsolver.h>
#include <Mesh/TPZCompMeshTools.h>
#include <Mesh/TPZMultiphysicsCompMesh.h>
#include <Post/TPZVTKGeoMesh.h>
#include <Pre/TPZGmshReader.h>
#include <Pre/TPZHybridizeHDiv.h>
#include <StrMatrix/pzstrmatrix.h>
#include <libInterpolate/Interpolate.hpp>

#include <Analysis/pzanalysis.h>
#include <ProblemConfig.h>
#include <StrMatrix/TPZSSpStructMatrix.h>
#include <iostream>
#include <map>
#include <vector>

TPZGeoMesh *CreateFlatGeoMesh(std::string &geometry_file2D);

void ModifyZCoordinates(TPZGeoMesh *gmesh, std::string &filename);

void ReadReservoirGeometryData(const std::string &name, std::vector<double> &x,
                               std::vector<double> &y, std::vector<double> &z);

void PrintGeometry(TPZGeoMesh *gmesh, const std::string &file_name);

void UNISIMHDiv();

TPZMultiphysicsCompMesh *CreateMixedMesh(const ProblemConfig &problem);

TPZCompMesh *CreateFluxMesh(const ProblemConfig &problem);

TPZCompMesh *CreatePressureL2Mesh(const ProblemConfig &problem);

void hAdaptivity(TPZCompMesh *postProcessMesh, TPZGeoMesh *gmeshToRefine);

void SolveMixedHybridProblem(TPZCompMesh *Hybridmesh, int InterfaceMatId,
                        const ProblemConfig &problem);
int main() {
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    UNISIMHDiv();

    return 0;
}

void UNISIMHDiv() {

    std::string geometry_file2D = "InputData/UNISIMFlatMesh.msh";
    TPZGeoMesh *gmesh = CreateFlatGeoMesh(geometry_file2D);

    std::string name = "InitialGeoMesh";
    PrintGeometry(gmesh, name);

    ProblemConfig config;
    config.porder = 1;
    config.hdivmais = 1;
    config.dimension = 2;
    config.makepressurecontinuous = true;

    TLaplaceExample1 example;
    config.exact.fExact = example.EConst;
    config.problemname = "SinSin";
    config.dir_name = "TesteUNISIM";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    config.materialids.insert(1);
    config.materialids.insert(2);
    config.bcmaterialids.insert(3);
    config.bcmaterialids.insert(4);
    config.bcmaterialids.insert(5);

    config.gmesh = gmesh;

    TPZMultiphysicsCompMesh *cmesh_HDiv = CreateMixedMesh(config);
    cmesh_HDiv->InitializeBlock();

    // Hybridizes mesh
    TPZHybridizeHDiv hybrid;
    auto HybridMesh = hybrid.Hybridize(cmesh_HDiv);
    HybridMesh->CleanUpUnconnectedNodes();
    HybridMesh->AdjustBoundaryElements();

    cmesh_HDiv = HybridMesh;

    SolveMixedHybridProblem(cmesh_HDiv, hybrid.fInterfaceMatid, config);

    //    {
    //
    //        if(RunMark){
    //            TPZHDivErrorEstimatorH1 HDivEstimate(*cmesh_HDiv);
    //            HDivEstimate.fProblemConfig = config;
    //            HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
    //            HDivEstimate.SetAnalyticSolution(config.exact);
    //            HDivEstimate.fperformUplift = true;
    //            HDivEstimate.fUpliftOrder = config.hdivmais;
    //
    //            HDivEstimate.PotentialReconstruction();
    //
    //            TPZManVector<REAL> elementerrors;
    //            HDivEstimate.ComputeErrors(elementerrors);
    //
    //        }
    //
    //        else{
    //
    //            TPZHybridHDivErrorEstimator HDivEstimate(*cmesh_HDiv);
    //            HDivEstimate.fProblemConfig = config;
    //            HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
    //            HDivEstimate.SetAnalyticSolution(config.exact);
    //            HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
    //
    //            HDivEstimate.fPostProcesswithHDiv = false;
    //
    //            HDivEstimate.PotentialReconstruction();
    //
    //            TPZManVector<REAL> elementerrors;
    //            HDivEstimate.ComputeErrors(elementerrors);
    //            hAdaptivity(&HDivEstimate.fPostProcMesh, hybridEstimatorMesh);
    //
    //        }
}

TPZGeoMesh *CreateFlatGeoMesh(std::string &geometry_file2D) {

    TPZGmshReader Geometry;
    TPZGeoMesh *gmesh2d;
    REAL l = 1.0;
    Geometry.SetCharacteristiclength(l);
    Geometry.SetFormatVersion("4.1");
    gmesh2d = Geometry.GeometricGmshMesh(geometry_file2D);
    Geometry.PrintPartitionSummary(std::cout);

    std::string filename1 = "InputData/UNISIMPointCloud.txt";
    //ModifyZCoordinates(gmesh2d, filename1);
    return gmesh2d;
}

void ModifyZCoordinates(TPZGeoMesh *gmesh, std::string &filename) {

    std::vector<double> x, y, z;
    ReadReservoirGeometryData(filename, x, y, z);

    _2D::ThinPlateSplineInterpolator<double> interp;
    interp.setData(x, y, z);

    int64_t nCoordinates = gmesh->NodeVec().NElements();
    double sum = 0.0;
    for (auto val : z) {
        sum += val;
    }
    double val_storage = sum / z.size();

    for (int icoord = 0; icoord < nCoordinates; icoord++) {
        TPZGeoNode node = gmesh->NodeVec()[icoord];
        TPZVec<REAL> co(3);
        node.GetCoordinates(co);
        double val_interp = interp(co[0], co[1]);

        if (val_interp == 0.0) {
            co[2] = val_storage;
        }
        if (val_interp > 1.0) {
            co[2] = val_interp;
        }

        gmesh->NodeVec()[icoord].SetCoord(co);
    }
}

void ReadReservoirGeometryData(const std::string &name, std::vector<double> &x,
                               std::vector<double> &y, std::vector<double> &z) {
    std::ifstream file;
    file.open(name);

    std::string line;
    int i = 1;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        char l = line[0];
        if (l != '/') {
            i = i + 1;
            int val = i % 15;
            if (val == 0) {
                double a, b, c;
                if (iss >> a >> b >> c)
                    ;
                x.push_back(a);
                y.push_back(b);
                z.push_back(c);
            }
        }
    }

    if (x.empty()) {
        std::cout << "No data read.\n";
        DebugStop();
    }

    std::cout << "File successfully read!\n";
}

void PrintGeometry(TPZGeoMesh *gmesh, const std::string &file_name) {
    std::stringstream txt_name;
    std::stringstream vtk_name;
    txt_name << file_name << ".txt";
    vtk_name << file_name << ".vtk";
    std::ofstream textfile(txt_name.str().c_str());
    gmesh->Print(textfile);
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
}

TPZCompMesh *CreatePressureL2Mesh(const ProblemConfig &problem) {
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial *mat = 0;
    for (auto matid : problem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }

    cmesh->SetDefaultOrder(problem.porder + problem.hdivmais);
    // cmesh->SetDefaultOrder(problem.porder);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->AutoBuild();
    int64_t n_connects = cmesh->NConnects();
    for (int64_t i = 0; i < n_connects; ++i) {
        cmesh->ConnectVec()[i].SetLagrangeMultiplier(1);
    }

    return cmesh;
}

TPZCompMesh *CreateFluxMesh(const ProblemConfig &problem) {
    int dim = problem.gmesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial *mat = NULL;
    problem.gmesh->ResetReference();
    for (auto matid : problem.materialids) {
        TPZVecL2 *mix = new TPZVecL2(matid);
        mix->SetDimension(dim);
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 1.);
        int bctype;
        if (matid == 3 || matid == 4 || matid == 5) {
            bctype = 0;
        } else {
            bctype = 1;
        }
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(problem.exact.Exact());
        cmesh->InsertMaterialObject(bc);
    }
    cmesh->SetDefaultOrder(problem.porder);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmesh->AutoBuild();

    cmesh->InitializeBlock();
    return cmesh;
}

TPZMultiphysicsCompMesh *CreateMixedMesh(const ProblemConfig &problem) {

    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(problem.gmesh);
    TPZFMatrix<REAL> K(3, 3, 0), invK(3, 3, 0);
    TPZMaterial *mat = NULL;
    K.Identity();
    invK.Identity();

    for (auto matid : problem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        mix->SetForcingFunction(problem.exact.ForcingFunction());
        mix->SetForcingFunctionExact(problem.exact.Exact());
        mix->SetPermeabilityTensor(K, invK);

        if (!mat) mat = mix;

        cmesh->InsertMaterialObject(mix);
    }
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype;
        if (matid == 3 || matid == 4 || matid == 5) {
            bctype = 0;
        } else {
            bctype = 1;
        }
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(problem.exact.Exact());
        cmesh->InsertMaterialObject(bc);
    }
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

    TPZManVector<int> active(2, 1);
    TPZManVector<TPZCompMesh *> meshvector(2, 0);

    meshvector[0] = CreateFluxMesh(problem);
    meshvector[1] = CreatePressureL2Mesh(problem);

    TPZCompMeshTools::AdjustFluxPolynomialOrders(meshvector[0],
                                                 problem.hdivmais);
    TPZCompMeshTools::SetPressureOrders(meshvector[0], meshvector[1]);

    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    bool keepmatrix = false;
    bool keeponelagrangian = true;
    TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian,
                                               keepmatrix);

    return cmesh;
}

void SolveMixedHybridProblem(TPZCompMesh *Hybridmesh, int InterfaceMatId,
                        const ProblemConfig &problem) {

    TPZAnalysis an(Hybridmesh);

#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(Hybridmesh);
    strmat.SetNumThreads(4);
#else
    TPZSkylineStructMatrix strmat(Hybridmesh);
    strmat.SetNumThreads(0);
#endif

    std::set<int> matIds;

    for (auto matid : problem.materialids) {

        matIds.insert(matid);
    }

    for (auto matidbc : problem.bcmaterialids) {

        matIds.insert(matidbc);
    }

    matIds.insert(InterfaceMatId);

    strmat.SetMaterialIds(matIds);

    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();

    std::cout << "Writing output files...\n";
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("ExactPressure");
    scalnames.Push("Pressure");
    vecnames.Push("ExactFlux");
    vecnames.Push("Flux");

    std::stringstream sout;
    sout << problem.dir_name << "/"
         << "OriginalHybrid_Order_" << problem.porder << "Nref_"
         << problem.ndivisions << "NAdapStep_" << problem.adaptivityStep
         << ".vtk";
    an.DefineGraphMesh(2, scalnames, vecnames, sout.str());
    int resolution = 2;
    an.PostProcess(resolution, Hybridmesh->Dimension());
}

void hAdaptivity(TPZCompMesh *postProcessMesh, TPZGeoMesh *gmeshToRefine) {

    // Column of the flux error estimate on the element solution matrix
    const int fluxErrorEstimateCol = 3;

    int64_t nelem = postProcessMesh->ElementSolution().Rows();

    // Iterates through element errors to get the maximum value
    REAL maxError = 0.;
    for (int64_t iel = 0; iel < nelem; iel++) {
        TPZCompEl *cel = postProcessMesh->ElementVec()[iel];
        if (!cel) continue;
        if (cel->Dimension() != postProcessMesh->Dimension()) continue;
        REAL elementError =
            postProcessMesh->ElementSolution()(iel, fluxErrorEstimateCol);

        if (elementError > maxError) {
            maxError = elementError;
        }
    }

    // Refines elements which error are bigger than 30% of the maximum error
    REAL threshold = 0.3 * maxError;

    for (int64_t iel = 0; iel < nelem; iel++) {
        TPZCompEl *cel = postProcessMesh->ElementVec()[iel];
        if (!cel) continue;
        if (cel->Dimension() != postProcessMesh->Dimension()) continue;

        REAL elementError =
            postProcessMesh->ElementSolution()(iel, fluxErrorEstimateCol);
        if (elementError > threshold) {
            TPZGeoEl *gel = cel->Reference();
            TPZVec<TPZGeoEl *> sons;
            TPZGeoEl *gelToRefine = gmeshToRefine->ElementVec()[gel->Id()];
            if (gelToRefine && !gelToRefine->HasSubElement()) {
                gelToRefine->Divide(sons);
            }
        }
    }
    DivideLowerDimensionalElements(gmeshToRefine);
}
