#include "Tools.h"

#include <Material/REAL/mixedpoisson.h>
#include <Material/TPZNullMaterial.h>
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
#include <string>
#include <vector>

TPZGeoMesh *CreateFlatGeoMesh(std::string &gmshFile);

void ModifyZCoordinates(TPZGeoMesh *gmesh, std::string &filename);

void ReadReservoirGeometryData(const std::string &name, std::vector<double> &x,
                               std::vector<double> &y, std::vector<double> &z);

void PrintGeometry(TPZGeoMesh *gmesh, const std::string &file_name);

void UNISIMHDiv();

TPZMultiphysicsCompMesh *CreateMixedCMesh(const ProblemConfig &problem);

TPZCompMesh *CreateFluxCMesh(const ProblemConfig &problem);

TPZCompMesh *CreatePressureCMesh(const ProblemConfig &problem);

void hAdaptivity(TPZCompMesh *postProcessMesh, TPZGeoMesh *gmeshToRefine);

void SolveMixedHybridProblem(TPZCompMesh *Hybridmesh,
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
    config.dir_name = "TesteUNISIM";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.bcmaterialids.insert(-2);
    config.bcmaterialids.insert(-3);

    config.gmesh = gmesh;

    TPZMultiphysicsCompMesh *cmesh_HDiv = CreateMixedCMesh(config);
    cmesh_HDiv->InitializeBlock();

    // Hybridizes mesh
    TPZHybridizeHDiv hybrid;
    auto HybridMesh = hybrid.Hybridize(cmesh_HDiv);
    HybridMesh->CleanUpUnconnectedNodes();
    HybridMesh->AdjustBoundaryElements();

    cmesh_HDiv = HybridMesh;

    SolveMixedHybridProblem(cmesh_HDiv, config);

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

TPZGeoMesh *CreateFlatGeoMesh(std::string &gmshFile) {

    TPZGmshReader gmeshReader;
    TPZGeoMesh *gmesh;

    TPZManVector<std::map<std::string, int>, 4> meshMaterialData(3);
    // B.C. materials
    meshMaterialData[1].insert({"ZeroFlux", -1});
    meshMaterialData[1].insert({"Productors", -2});
    meshMaterialData[1].insert({"Injectors", -3});
    // Domain materials
    meshMaterialData[2].insert({"RockMatrix", 1});
    meshMaterialData[2].insert({"RockMatrix2", 1});
    gmeshReader.SetDimNamePhysical(meshMaterialData);

    gmeshReader.SetFormatVersion("4.1");
    gmesh = gmeshReader.GeometricGmshMesh(gmshFile);

    // std::string filename1 = "InputData/UNISIMPointCloud.txt";
    // ModifyZCoordinates(gmesh, filename1);
    return gmesh;
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

TPZCompMesh *CreatePressureCMesh(const ProblemConfig &problem) {
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);

    TPZNullMaterial *pressureMat = new TPZNullMaterial(1);
    cmesh->InsertMaterialObject(pressureMat);

    cmesh->SetDefaultOrder(problem.porder + problem.hdivmais);

    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    cmesh->AutoBuild();
    // This puts pressure equations after flux equations. It is necessary so we
    // don't end up with '0' values in the diagonal of the global matrix.
    int64_t n_connects = cmesh->NConnects();
    for (int64_t i = 0; i < n_connects; ++i) {
        cmesh->ConnectVec()[i].SetLagrangeMultiplier(1);
    }

    return cmesh;
}

TPZCompMesh *CreateFluxCMesh(const ProblemConfig &problem) {
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);

    TPZVecL2 *fluxMat = new TPZVecL2(1);
    fluxMat->SetDimension(cmesh->Dimension());
    cmesh->InsertMaterialObject(fluxMat);

    for (auto bcID : {-1, -2, -3}) {
        // BC values defined here won't be used
        TPZFNMatrix<1, REAL> val1, val2(1, 1, 0.);
        int bctype = -1;
        TPZBndCond *bc = fluxMat->CreateBC(fluxMat, bcID, bctype, val1, val2);
        cmesh->InsertMaterialObject(bc);
    }

    cmesh->SetDefaultOrder(problem.porder);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(cmesh->Dimension());
    cmesh->AutoBuild();

    cmesh->InitializeBlock();
    return cmesh;
}

TPZMultiphysicsCompMesh *CreateMixedCMesh(const ProblemConfig &problem) {
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(problem.gmesh);

    TPZMixedPoisson *mix = new TPZMixedPoisson(1, cmesh->Dimension());

    TPZFMatrix<REAL> K(3, 3, 0), invK(3, 3, 0);
    K.Identity();
    invK.Identity();
    mix->SetPermeabilityTensor(K, invK);

    cmesh->InsertMaterialObject(mix);

    // Insert boundary conditions
    TPZFNMatrix<1, REAL> val1(1, 1, 1.e12); // Not used by the material
    TPZFNMatrix<1, REAL> val2(1, 1, 0.);
    const int dirichlet = 0;
    const int neumann = 1;

    // Zero flux (reservoir boundary)
    TPZBndCond *zeroFlux = mix->CreateBC(mix, -1, neumann, val1, val2);
    // Productors
    val2(0, 0) = -10.;
    TPZBndCond *productors = mix->CreateBC(mix, -2, dirichlet, val1, val2);
    // Injectors
    val2(0, 0) = 20.;
    TPZBndCond *injectors = mix->CreateBC(mix, -3, dirichlet, val1, val2);

    cmesh->InsertMaterialObject(zeroFlux);
    cmesh->InsertMaterialObject(productors);
    cmesh->InsertMaterialObject(injectors);

    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

    TPZManVector<int> active(2, 1);
    TPZManVector<TPZCompMesh *> meshVector(2, 0);

    meshVector[0] = CreateFluxCMesh(problem);
    meshVector[1] = CreatePressureCMesh(problem);

    TPZCompMeshTools::AdjustFluxPolynomialOrders(meshVector[0],
                                                 problem.hdivmais);
    TPZCompMeshTools::SetPressureOrders(meshVector[0], meshVector[1]);

    cmesh->BuildMultiphysicsSpace(active, meshVector);
    cmesh->LoadReferences();
    bool keepmatrix = false;
    bool keeponelagrangian = true;
    TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian,
                                               keepmatrix);

    return cmesh;
}

void SolveMixedHybridProblem(TPZCompMesh *Hybridmesh,
                             const ProblemConfig &problem) {

    TPZAnalysis an(Hybridmesh);

#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(Hybridmesh);
    strmat.SetNumThreads(4);
#else
    TPZSkylineStructMatrix strmat(Hybridmesh);
    strmat.SetNumThreads(4);
#endif

    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;

    an.Assemble();
    an.Solve();

    std::cout << "Writing output files...\n";

    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    vecnames.Push("Flux");

    std::stringstream sout;
    sout << problem.dir_name + "/Hdiv-Order" << problem.porder << ".vtk";
    an.DefineGraphMesh(2, scalnames, vecnames, sout.str());
    an.PostProcess(1, Hybridmesh->Dimension());
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
