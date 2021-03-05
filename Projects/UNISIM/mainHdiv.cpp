#include "Tools.h"
#include "ToolsUNISIM.h"

#include <Analysis/pzanalysis.h>
#include <Geom/pzgeoquad.h>
#include <Material/REAL/mixedpoisson.h>
#include <Material/TPZNullMaterial.h>
#include <Material/TPZVecL2.h>
#include <Material/pzbndcond.h>
#include <Matrix/pzstepsolver.h>
#include <Mesh/TPZCompMeshTools.h>
#include <Mesh/TPZMultiphysicsCompMesh.h>
#include <Mesh/tpzgeoelrefpattern.h>
#include <Pre/TPZHybridizeHDiv.h>
#include <ProblemConfig.h>
#include <StrMatrix/TPZSSpStructMatrix.h>
#include <StrMatrix/pzstrmatrix.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

void UNISIMMHM(TPZGeoMesh *gmesh, std::vector<std::pair<REAL, int64_t>> &results);

TPZMultiphysicsCompMesh *CreateMixedCMesh(const ProblemConfig &problem);

TPZCompMesh *CreateFluxCMesh(const ProblemConfig &problem);

TPZCompMesh *CreatePressureCMesh(const ProblemConfig &problem);

void SolveMHMProblem(TPZCompMesh *Hybridmesh, const ProblemConfig &problem);

int main() {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    gRefDBase.InitializeRefPatterns(2);
    gRefDBase.InitializeAllUniformRefPatterns();

    bool modifyZCoords = false;
    TPZGeoMesh *gmesh = CreateUNISIMSurfaceGeoMesh(false);
    std::string meshFileName{"UNISIMMesh"};
    int nDirectionalRefinements = 0;
    PrintGeometry(gmesh, meshFileName, false, true);
    ApplyDirectionalRefinement(gmesh, nDirectionalRefinements);
    meshFileName.append("AfterDirectionalRef");
    PrintGeometry(gmesh, meshFileName, false, true);

    int nSteps = 1;
    std::vector<std::pair<REAL, int64_t>> results; // Stores error and nDOF
    for (int i = 0; i < nSteps; i++) {
        UNISIMMHM(gmesh, results);
    }
    for (size_t i = 0; i < results.size(); i++) {
        std::cout << "Step: " << i << " Energy Error: " << results[i].first << " nDOF: " << results[i].second << '\n';
    }
    delete gmesh;
    return 0;
}

void UNISIMMHM(TPZGeoMesh *gmesh, std::vector<std::pair<REAL, int64_t>> &results) {

    static int adaptivityStep = 0;

    ProblemConfig config;
    config.porder = 1;
    config.hdivmais = 1;
    config.dimension = 2;
    config.adaptivityStep = adaptivityStep;
    config.makepressurecontinuous = true;

#ifdef DEBUGTEST
    config.dir_name = "DebugTest";
    {
        config.exact = new TLaplaceExample1;
        config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
    }
#else
    config.dir_name = "UNISIM_Flat_AdaptivityMore";
#endif
    config.problemname = "UNISIM_HDIV";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    std::stringstream gmeshFileName;
    gmeshFileName << config.dir_name << "/GeoMesh" << adaptivityStep;
    PrintGeometry(gmesh, gmeshFileName.str(), false, true);

    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.bcmaterialids.insert(-2);
    config.bcmaterialids.insert(-3);

    config.gmesh = new TPZGeoMesh(*gmesh);

    TPZMultiphysicsCompMesh *cmesh_HDiv = CreateMixedCMesh(config);
    cmesh_HDiv->InitializeBlock();

    // Solves FEM problem
    SolveMHMProblem(cmesh_HDiv, config);
    int64_t neq = cmesh_HDiv->NEquations();
    std::cout << "Finished simulation!\n";
    std::cout << "Starting error estimation procedure...\n";

    TPZManVector<REAL, 6> errorVec;
    TPZManVector<REAL> elementerrors;
    {
        // Estimates error
        bool postProcWithHDiv = false;
        TPZHybridHDivErrorEstimator HDivEstimate(*cmesh_HDiv, postProcWithHDiv);
        HDivEstimate.SetProblemConfig(config);
        HDivEstimate.SetPostProcUpliftOrder(config.hdivmais);

#ifdef DEBUGTEST
        HDivEstimate.SetAnalyticSolution(config.exact);
#endif
        std::cout << "Reconstructing potential...\n";
        HDivEstimate.PotentialReconstruction();

        std::cout << "Computing errors...\n";
        std::string vtkPath = config.dir_name + "/UNISIM_error_results" + std::to_string(config.adaptivityStep) + ".vtk";
        HDivEstimate.ComputeErrors(errorVec, elementerrors, vtkPath);
    }
    delete config.gmesh;
    // h-refinement on elements with bigger errors
    // hAdaptivity(gmesh, elementerrors, 0.3);
    // UniformRefinement(1, gmesh);
    results.emplace_back(errorVec[1], neq);
    adaptivityStep++;
}

TPZCompMesh *CreatePressureCMesh(const ProblemConfig &problem) {
    auto *cmesh = new TPZCompMesh(problem.gmesh);

    auto *pressureMat = new TPZNullMaterial(1, problem.gmesh->Dimension(), 1);
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
    auto *cmesh = new TPZCompMesh(problem.gmesh);

    auto *fluxMat = new TPZVecL2(1);
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
    auto *cmesh = new TPZMultiphysicsCompMesh(problem.gmesh);

    auto *mix = new TPZMixedPoisson(1, cmesh->Dimension());

#ifdef DEBUGTEST
    {
        mix->SetForcingFunction(problem.exact.operator*().ForcingFunction());
        mix->SetForcingFunctionExact(problem.exact.operator*().Exact());
    }
#endif
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
#ifdef DEBUGTEST
    productors->TPZMaterial::SetForcingFunction(problem.exact.operator*().Exact());
    injectors->TPZMaterial::SetForcingFunction(problem.exact.operator*().Exact());
#endif

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

void SolveMHMProblem(TPZCompMesh *Hybridmesh, const ProblemConfig &problem) {

    TPZAnalysis an(Hybridmesh);

#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(Hybridmesh);
    strmat.SetNumThreads(0);
#else
    TPZSkylineStructMatrix strmat(Hybridmesh);
    strmat.SetNumThreads(4);
#endif

    an.SetStructuralMatrix(strmat);

    auto *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;

    an.Assemble();
    an.Solve();
}
