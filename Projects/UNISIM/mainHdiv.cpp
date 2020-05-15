#include "Tools.h"

#include <Analysis/pzanalysis.h>
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
#include <ProblemConfig.h>
#include <Refine/TPZRefPatternTools.h>
#include <StrMatrix/TPZSSpStructMatrix.h>
#include <StrMatrix/pzstrmatrix.h>

#include <libInterpolate/Interpolate.hpp>

#include <Analysis/pzanalysis.h>
#include <ProblemConfig.h>
#include <StrMatrix/TPZSSpStructMatrix.h>
#include "pzcheckgeom.h"
#include <iostream>
#include <map>
#include <string>
#include <vector>

TPZGeoMesh *CreateFlatGeoMesh();

void ModifyZCoordinates(TPZGeoMesh *gmesh, std::string &filename);

void ReadReservoirGeometryData(const std::string &name, std::vector<double> &x, std::vector<double> &y,
                               std::vector<double> &z);

void PrintGeometry(TPZGeoMesh *gmesh, const std::string &file_name, bool printTXT = true, bool printVTK = true);

void UNISIMHDiv(TPZGeoMesh *gmesh);

TPZMultiphysicsCompMesh *CreateMixedCMesh(const ProblemConfig &problem);

TPZCompMesh *CreateFluxCMesh(const ProblemConfig &problem);

TPZCompMesh *CreatePressureCMesh(const ProblemConfig &problem);

void SolveMixedHybridProblem(TPZCompMesh *Hybridmesh, const ProblemConfig &problem);

void hAdaptivity(TPZHybridHDivErrorEstimator &estimator, REAL thresholdRatio);

void ApplyDirectionalRefinement(TPZGeoMesh *gmesh, int nRef);

int main() {
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    gRefDBase.InitializeRefPatterns(2);

    TPZGeoMesh *gmesh = CreateFlatGeoMesh();
    int nDirectionalRefinements = 1;
    ApplyDirectionalRefinement(gmesh, nDirectionalRefinements);


    TPZManVector<REAL,3> nod0(3,0.);
    gmesh->NodeVec()[0].GetCoordinates(nod0);
    int64_t nnodes = gmesh->NNodes();
    for (int64_t no = 0; no<nnodes; no++) {
        TPZManVector<REAL, 3> co(3);
        gmesh->NodeVec()[no].GetCoordinates(co);
        for(int ic=0; ic<3; ic++) co[ic] -= nod0[ic];
        gmesh->NodeVec()[no].SetCoord(co);
    }
//    {
//        std::ofstream out("unisim.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
//        std::ofstream outgeo("gmesh.txt");
//        gmesh->Print(outgeo);
//        TPZCheckGeom check(gmesh);
//        check.UniformRefine(1);
//        std::ofstream out2("unisim2.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
//    }
//    return 0;
    int nSteps = 4;
    for (int i = 0; i < nSteps; i++) {
        UNISIMHDiv(gmesh);
    }

    return 0;
}

void ApplyDirectionalRefinement(TPZGeoMesh *gmesh, int nRef) {
    // Mat IDs of productors and injectors BCs
    set<int> matids{-2, -3};

    for (auto i = 0; i < nRef; i++) {
        int nelements = gmesh->NElements();
        cout << "Refinement step: " << i << "\nNumber of elements = " << nelements << '\n';
        for (auto el = 0; el < nelements; el++) {
            TPZGeoEl *element = gmesh->ElementVec()[el];
            if (!element) continue;
            TPZRefPatternTools::RefineDirectional(element, matids);
        }
        stringstream meshfilename;
        meshfilename << "DirectionalRefinementGMesh" << i;
        PrintGeometry(gmesh, meshfilename.str(), false, true);
    }
}

void UNISIMHDiv(TPZGeoMesh *gmesh) {

    static int adaptivityStep = 0;

    ProblemConfig config;
    config.porder = 1;
    config.hdivmais = 1;
    config.dimension = 2;
    config.adaptivityStep = adaptivityStep;
    config.makepressurecontinuous = true;

    config.dir_name = "TesteUNISIM";
    config.problemname = "UNISIM_Errors";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    std::stringstream gmeshFileName;
    gmeshFileName << config.dir_name << "/GeoMesh" << adaptivityStep;
    PrintGeometry(gmesh, gmeshFileName.str(), false, true);

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

    // Solves FEM problem
    SolveMixedHybridProblem(cmesh_HDiv, config);

    // Estimates error
    TPZHybridHDivErrorEstimator HDivEstimate(*cmesh_HDiv);
    HDivEstimate.fProblemConfig = config;
    HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;

    HDivEstimate.fPostProcesswithHDiv = false;

    HDivEstimate.PotentialReconstruction();

    TPZManVector<REAL> elementerrors;
    HDivEstimate.ComputeErrors(elementerrors);

    // h-refinement
    hAdaptivity(HDivEstimate, 0.5);
    adaptivityStep++;
}

TPZGeoMesh *CreateFlatGeoMesh() {

#ifdef MACOSX
    std::string gmshFile = "../InputData/UNISIMFlatMesh.msh";
#else
    std::string gmshFile = "InputData/UNISIMFlatMesh.msh";
#endif
    TPZGmshReader gmeshReader;
    TPZGeoMesh *gmesh;

    TPZManVector<std::map<std::string, int>, 4> meshMaterialData(3);
    // BC materials
    meshMaterialData[1].insert({"ZeroFlux", -1});
    meshMaterialData[1].insert({"Productors", -2});
    meshMaterialData[1].insert({"Injectors", -3});
    meshMaterialData[1].insert({"Faults", 99}); // Won't be used
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

void ReadReservoirGeometryData(const std::string &name, std::vector<double> &x, std::vector<double> &y,
                               std::vector<double> &z) {
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

void PrintGeometry(TPZGeoMesh *gmesh, const std::string &file_name, bool printTXT, bool printVTK) {
    if (printTXT) {
        std::stringstream txt_name;
        txt_name << file_name << ".txt";
        std::ofstream textfile(txt_name.str().c_str());
        gmesh->Print(textfile);
    }
    if (printVTK) {
        std::stringstream vtk_name;
        vtk_name << file_name << ".vtk";
        std::ofstream vtkfile(vtk_name.str().c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
    }
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
        // The information here is not important.
        // The materials are needed only so the mesh creates BC elements.
        TPZFNMatrix<1, REAL> val1, val2;
        int bctype = 999;
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

    // Zero flux (reservoir boundary)
    TPZBndCond *zeroFlux = mix->CreateBC(mix, -1, dirichlet, val1, val2);
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

    TPZCompMeshTools::AdjustFluxPolynomialOrders(meshVector[0], problem.hdivmais);
    TPZCompMeshTools::SetPressureOrders(meshVector[0], meshVector[1]);

    cmesh->BuildMultiphysicsSpace(active, meshVector);
    cmesh->LoadReferences();
    bool keepmatrix = false;
    bool keeponelagrangian = true;
    TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian, keepmatrix);

    return cmesh;
}

void SolveMixedHybridProblem(TPZCompMesh *Hybridmesh, const ProblemConfig &problem) {

    TPZAnalysis an(Hybridmesh);

#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(Hybridmesh);
    strmat.SetNumThreads(0);
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
    sout << problem.dir_name + "/Hdiv-Order" << problem.porder << "-AdaptivityStep" << problem.adaptivityStep << ".vtk";
    an.DefineGraphMesh(2, scalnames, vecnames, sout.str());
    an.PostProcess(2, Hybridmesh->Dimension());
}

void hAdaptivity(TPZHybridHDivErrorEstimator &estimator, REAL thresholdRatio) {
    // Column of the flux error estimate on the element solution matrix
    const int fluxErrorCol = 3;

    TPZCompMesh *postProcessMesh = &estimator.fPostProcMesh;
    TPZGeoMesh *gmesh = postProcessMesh->Reference();
    postProcessMesh->LoadReferences();

    // Cleans geometric elements of interfaces and other auxiliary materials
    // and iterates through element errors to get the maximum value
    std::set<int> matIdsToDelete;
    matIdsToDelete.insert(estimator.fHybridizer.fLagrangeInterface);
    matIdsToDelete.insert(estimator.fHybridizer.fHDivWrapMatid);
    matIdsToDelete.insert(estimator.fHybridizer.fInterfaceMatid);

    REAL maxError = 0.;
    int64_t nElemsToBeKept = 0;

    int64_t nelem = gmesh->NElements();
    for (int64_t iel = 0; iel < nelem; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if (!gel) continue;
        int elMatId = gel->MaterialId();
        if (matIdsToDelete.find(elMatId) != matIdsToDelete.end()) {
            if (nElemsToBeKept == 0) nElemsToBeKept = iel;
            gmesh->DeleteElement(gel, iel);
            continue;
        }

        if (gel->Dimension() != postProcessMesh->Dimension()) continue;

        TPZCompEl *cel = gel->Reference();
        if (!cel) continue;
        int64_t celId = cel->Index();

        REAL elemError = postProcessMesh->ElementSolution()(celId, fluxErrorCol);
        if (elemError > maxError) {
            maxError = elemError;
        }
    }

    REAL threshold = thresholdRatio * maxError;

    nelem = postProcessMesh->NElements();
    for (int64_t iel = 0; iel < nelem; iel++) {
        TPZCompEl *cel = postProcessMesh->ElementVec()[iel];
        if (!cel) continue;
        if (cel->Dimension() != postProcessMesh->Dimension()) continue;

        TPZGeoEl *gel = cel->Reference();
        if (!gel) DebugStop();

        REAL elementError = postProcessMesh->ElementSolution()(iel, fluxErrorCol);
        if (elementError > threshold) {
            TPZVec<TPZGeoEl *> sons;
            if (!gel->HasSubElement()) {
                gel->Divide(sons);
            }
        }
    }
    DivideLowerDimensionalElements(gmesh);
}
