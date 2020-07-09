//
// Created by Gustavo A. Batistela on 06/07/2020.
//
// This file contains the numerical tests to be shown in the CILAMCE 2020 article.
//

#include <Mesh/pzgmesh.h>
#include <Pre/TPZMHMixedMeshControl.h>
#include <TPZMFSolutionTransfer.h>
#include <Tools.h>
#include <ToolsMHM.h>
#include <Util/pzlog.h>

void RunSinSinProblem();
void RunOscillatoryProblem();
void Run3DProblem();
void RunSingularPoblem();

TPZGeoMesh *CreateQuadGeoMesh(int nCoarseDiv, int nInternalRef);
TPZGeoMesh *CreateCubeGeoMesh();
TPZGeoMesh *CreateLGeoMesh();

TPZMultiphysicsCompMesh *CreateMixedCompMesh(const ProblemConfig &problem);
TPZCompMesh *CreateFluxCompMesh(const ProblemConfig &problem);
TPZCompMesh *CreatePressureCompMesh(const ProblemConfig &problem);
void SubstructureCompMesh(TPZMHMixedMeshControl *mhm, const ProblemConfig &config, int nInternalRef);

void SolveMHMProblem(TPZMHMixedMeshControl *mhm, const ProblemConfig &config);

void EstimateError(const ProblemConfig &config, const TPZMHMixedMeshControl *mhm);

int main() {
    InitializePZLOG();

    RunSinSinProblem();
    RunOscillatoryProblem();

    return 0;
}

void RunSinSinProblem() {

    ProblemConfig config;

    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
    config.problemname = "SinSin";
    config.dir_name = "CILAMCE";
    config.porder = 1;
    config.hdivmais = 1;
    config.ndivisions = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);

    int nCoarseDiv = 4;
    int nInternalRef = 4;
    config.gmesh = CreateQuadGeoMesh(nCoarseDiv, nInternalRef);

    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    TPZMultiphysicsCompMesh *multiphysicsMesh = CreateMixedCompMesh(config);

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    SubstructureCompMesh(mhm, config, nInternalRef);
    SolveMHMProblem(mhm, config);
    EstimateError(config, mhm);
}

void RunOscillatoryProblem() {
    ProblemConfig config;

    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::EBubble;
    config.problemname = "Bubble";
    config.dir_name = "CILAMCE";
    config.porder = 1;
    config.hdivmais = 1;
    config.ndivisions = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);

    int nCoarseDiv = 4;
    int nInternalRef = 4;
    config.gmesh = CreateQuadGeoMesh(nCoarseDiv, nInternalRef);

    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    TPZMultiphysicsCompMesh *multiphysicsMesh = CreateMixedCompMesh(config);

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    SubstructureCompMesh(mhm, config, nInternalRef);
    SolveMHMProblem(mhm, config);
    EstimateError(config, mhm);

}

void Run3DProblem() {}

void RunSingularPoblem() {}

TPZGeoMesh *CreateQuadGeoMesh(int nCoarseDiv, int nInternalRef) {

    TPZManVector<int, 4> bcIDs(4, -1);
    TPZGeoMesh *gmesh = CreateGeoMesh(nCoarseDiv, bcIDs);

    UniformRefinement(nInternalRef, gmesh);
    DivideLowerDimensionalElements(gmesh);

    return gmesh;
}

TPZMultiphysicsCompMesh *CreateMixedCompMesh(const ProblemConfig &problem) {

    auto *cmesh = new TPZMultiphysicsCompMesh(problem.gmesh);
    TPZMaterial *mat = nullptr;
    TPZFMatrix<REAL> K(3, 3, 0), invK(3, 3, 0);
    K.Identity();
    invK.Identity();

    STATE Km = problem.Km;

    if (problem.TensorNonConst && problem.gmesh->Dimension() == 3) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (i == j) {
                    K(i, j) = 2.;
                    invK(i, j) = 3. / 4.;
                } else {
                    K(i, j) = 1.;
                    invK(i, j) = (-1.) / 4.;
                }
            }
        }
    }

    //    K.Print(std::cout);
    //    invK.Print(std::cout);

    for (auto matid : problem.materialids) {
        auto *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        mix->SetForcingFunction(problem.exact.operator*().ForcingFunction());
        mix->SetForcingFunctionExact(problem.exact.operator*().Exact());
        mix->SetPermeabilityTensor(K, invK);

        if (!mat) mat = mix;

        cmesh->InsertMaterialObject(mix);
    }

    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);

        int bctype;
        switch (matid) {
        case -1: {
            bctype = 0;
            break;
        }
        case -2: {
            bctype = 1;
            break;
        }
        default: {
            bctype = -99;
            DebugStop();
        }
        }

        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(problem.exact.operator*().Exact());
        cmesh->InsertMaterialObject(bc);
    }
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

    TPZManVector<int> active(2, 1);
    TPZManVector<TPZCompMesh *> meshvector(2, 0);

    meshvector[0] = CreateFluxCompMesh(problem);
    meshvector[1] = CreatePressureCompMesh(problem);

    TPZCompMeshTools::AdjustFluxPolynomialOrders(meshvector[0], problem.hdivmais);
    TPZCompMeshTools::SetPressureOrders(meshvector[0], meshvector[1]);

    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    bool keepmatrix = false;
    bool keeponelagrangian = true;
    TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian, keepmatrix);

    cmesh->InitializeBlock();
    return cmesh;
}

TPZCompMesh *CreateFluxCompMesh(const ProblemConfig &problem) {
    int dim = problem.gmesh->Dimension();
    auto *cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial *mat = nullptr;
    problem.gmesh->ResetReference();
    for (auto matid : problem.materialids) {
        auto *mix = new TPZVecL2(matid);
        mix->SetDimension(dim);
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 1.);
        int bctype;
        if (matid == -1 || matid == 2) {
            bctype = 0;
        } else {
            bctype = 1;
        }
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(problem.exact.operator*().Exact());
        cmesh->InsertMaterialObject(bc);
    }
    cmesh->SetDefaultOrder(problem.porder);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmesh->AutoBuild();

    cmesh->InitializeBlock();
    return cmesh;
}

TPZCompMesh *CreatePressureCompMesh(const ProblemConfig &problem) {
    auto *cmesh = new TPZCompMesh(problem.gmesh);

    for (const auto matid : problem.materialids) {
        auto *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        cmesh->InsertMaterialObject(mix);
    }

    cmesh->SetDefaultOrder(problem.porder + problem.hdivmais);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    cmesh->AutoBuild();

    for (int64_t i = 0; i < cmesh->NConnects(); ++i) {
        cmesh->ConnectVec()[i].SetLagrangeMultiplier(1);
    }

    return cmesh;
}

void SubstructureCompMesh(TPZMHMixedMeshControl *mhm, const ProblemConfig &config, int nInternalRef) {

    TPZAutoPointer<TPZGeoMesh> gmeshauto(config.gmesh);

    TPZVec<int64_t> elementIndexes;
    ComputeCoarseIndices(gmeshauto.operator->(), elementIndexes);

    mhm->DefinePartitionbyCoarseIndices(elementIndexes);

    // Indicate material indices to the MHM control structure
    mhm->fMaterialIds = config.materialids;
    mhm->fMaterialBCIds = config.bcmaterialids;

    // Insert the material objects in the multiphysics mesh
    InsertMaterialObjects(*mhm, config);

    // General approximation order settings
    mhm->SetInternalPOrder(config.porder);
    mhm->SetSkeletonPOrder(config.porder);
    mhm->SetHdivmaismaisPOrder(config.hdivmais);

    // Refine skeleton elements
    mhm->DivideSkeletonElements(nInternalRef);
    mhm->DivideBoundarySkeletonElements();

    // Creates MHM mesh
    bool substructure = true;
    mhm->BuildComputationalMesh(substructure);
}

void SolveMHMProblem(TPZMHMixedMeshControl *mhm, const ProblemConfig &config) {

    TPZAutoPointer<TPZCompMesh> cmesh = mhm->CMesh();

    bool shouldrenumber = true;
    TPZAnalysis an(cmesh, shouldrenumber);

#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(0 /*config.n_threads*/);
#else
    TPZSkylineStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(config.n_threads);
#endif

    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "Assembling\n";
    an.Assemble();

    std::cout << "Solving\n";
    an.Solve();
    std::cout << "Finished\n";
    an.LoadSolution(); // compute internal dofs

    TPZMFSolutionTransfer transfer;
    transfer.BuildTransferData(cmesh.operator->());
    transfer.TransferFromMultiphysics();

    TPZStack<std::string> scalnames, vecnames;
    TPZMaterial *mat = cmesh->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }

    TLaplaceExample1 analytic = config.exact.operator*();
    an.SetExact(analytic.ExactSolution());

    scalnames.Push("Pressure");
    scalnames.Push("ExactPressure");
    vecnames.Push("Flux");
    vecnames.Push("ExactFlux");

    int resolution = 0;
    std::string plotname = "MHMResults.vtk";
    an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotname);
    an.PostProcess(resolution, cmesh->Dimension());

    TPZManVector<REAL> errors(4, 0.);
    an.SetThreadsForError(2);
    an.SetExact(analytic.ExactSolution());
    an.PostProcessError(errors, false);
}

void EstimateError(ProblemConfig &config, TPZMHMixedMeshControl *mhm) {
    cout << "Error Estimation processing for MHM-Hdiv problem " << endl;

    // Error estimation
    TPZMultiphysicsCompMesh *InputMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm->CMesh().operator->());
    if (!InputMesh) DebugStop();

    TPZMHMHDivErrorEstimator ErrorEstimator(*InputMesh, mhm);
    ErrorEstimator.fOriginalIsHybridized = false;
    ErrorEstimator.SetAnalyticSolution(config.exact);

    ErrorEstimator.fPostProcesswithHDiv = false;

    ErrorEstimator.fProblemConfig = config;

    ErrorEstimator.PotentialReconstruction();

    {
        string command = "mkdir " + config.dir_name;
        system(command.c_str());

        TPZManVector<REAL, 6> errors;
        TPZManVector<REAL> elementerrors;
        bool store_errors = true;
        ErrorEstimator.ComputeErrors(errors, elementerrors, store_errors);
    }
}
