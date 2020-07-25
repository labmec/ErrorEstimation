//
// Created by Gustavo A. Batistela on 06/07/2020.
//
// This file contains the numerical tests to be shown in the CILAMCE 2020 article.
//

#include <Mesh/pzgmesh.h>
#include <Pre/TPZGenGrid3D.h>
#include <Pre/TPZMHMixedMeshControl.h>
#include <TPZMFSolutionTransfer.h>
#include <Tools.h>
#include <ToolsMHM.h>
#include <Util/pzlog.h>

void RunCosCosProblem();
void RunOscillatoryProblem();
void RunNonConvexProblem();
void Run3DProblem();
void RunSingularProblem();

TPZGeoMesh *CreateQuadGeoMesh(int nCoarseDiv, int nInternalRef);
TPZGeoMesh *CreateCubeGeoMesh(int nCoarseRef, int nInternalRef, TPZStack<int64_t> &coarseIndexes);
TPZGeoMesh *CreateLShapeGeoMesh(int nCoarseRef, int nInternalRef, TPZStack<int64_t> &mhmIndexes);

void InsertMaterialsInMHMMesh(TPZMHMixedMeshControl &control, const ProblemConfig &config);
void InsertUNISIMMaterialsInMHMMesh(TPZMHMixedMeshControl &control, const ProblemConfig &config);
void CreateMHMCompMesh(TPZMHMixedMeshControl *mhm, const ProblemConfig &config, int nInternalRef,
                       bool definePartitionByCoarseIndex, TPZManVector<int64_t> mhmIndexes);

void SolveMHMProblem(TPZMHMixedMeshControl *mhm, const ProblemConfig &config);

void EstimateError(ProblemConfig &config, TPZMHMixedMeshControl *mhm);

int main() {
    InitializePZLOG();
    gRefDBase.InitializeAllUniformRefPatterns();

    //RunCosCosProblem();
    //RunOscillatoryProblem();
    RunNonConvexProblem();
    //Run3DProblem();
    //RunSingularProblem();

    return 0;
}

void RunCosCosProblem() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ECosCos;
    config.problemname = "CosCos";
    config.dir_name = "CILAMCE";
    config.porder = 1;
    config.hdivmais = 1;
    config.ndivisions = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nCoarseDiv = 2;
    int nInternalRef = 3;
    config.gmesh = CreateQuadGeoMesh(nCoarseDiv, nInternalRef);

    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    TPZManVector<int64_t> coarseIndexes;
    ComputeCoarseIndices(config.gmesh, coarseIndexes);
    bool definePartitionByCoarseIndexes = true;
    CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, coarseIndexes);

    SolveMHMProblem(mhm, config);
    EstimateError(config, mhm);
}

void RunOscillatoryProblem() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::EArcTan;
    config.problemname = "Oscillatory";
    config.dir_name = "CILAMCE";
    config.porder = 1;
    config.hdivmais = 1;
    config.ndivisions = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nCoarseDiv = 4;
    int nInternalRef = 3;
    config.gmesh = CreateQuadGeoMesh(nCoarseDiv, nInternalRef);

    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    TPZManVector<int64_t> coarseIndexes;
    ComputeCoarseIndices(config.gmesh, coarseIndexes);
    bool definePartitionByCoarseIndexes = true;
    CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, coarseIndexes);

    SolveMHMProblem(mhm, config);
    EstimateError(config, mhm);
}

void RunNonConvexProblem() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
    config.problemname = "NonConvex";
    config.dir_name = "CILAMCE";
    config.porder = 1;
    config.hdivmais = 1;
    config.ndivisions = 0;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nDiv = 1;
    int nInternalRef = 0;

    TPZVec<int64_t> coarseIndexes;
    config.gmesh = CreateLMHMMesh(nDiv, coarseIndexes);

    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    bool definePartitionByCoarseIndexes = false;
    CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, coarseIndexes);

    SolveMHMProblem(mhm, config);
    EstimateError(config, mhm);
}

void Run3DProblem() {
    ProblemConfig config;
    config.dimension = 3;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
    config.problemname = "SinSinCube";
    config.dir_name = "CILAMCE";
    config.porder = 1;
    config.hdivmais = 1;
    config.ndivisions = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nCoarseDiv = 3;
    int nInternalRef = 1;

    TPZStack<int64_t> mhmIndexes;
    config.gmesh = CreateCubeGeoMesh(nCoarseDiv, nInternalRef, mhmIndexes);

    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    bool definePartitionByCoarseIndexes = true;
    CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, mhmIndexes);

    SolveMHMProblem(mhm, config);
    EstimateError(config, mhm);
}

void RunSingularProblem() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinMark;
    config.problemname = "SinMarkLShape";
    config.dir_name = "CILAMCE";
    config.porder = 1;
    config.hdivmais = 1;
    config.ndivisions = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nCoarseRef = 1;
    int nInternalRef = 0;

    TPZStack<int64_t> mhmIndexes;
    config.gmesh = CreateLShapeGeoMesh(nCoarseRef, nInternalRef, mhmIndexes);

    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    bool definePartitionByCoarseIndexes = true;
    CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, mhmIndexes);

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    SolveMHMProblem(mhm, config);
    EstimateError(config, mhm);
}

TPZGeoMesh *CreateQuadGeoMesh(int nCoarseDiv, int nInternalRef) {

    TPZManVector<int, 4> bcIDs(4, -1);
    TPZGeoMesh *gmesh = CreateGeoMesh(nCoarseDiv, bcIDs);
    gmesh->SetDimension(2);

    UniformRefinement(nInternalRef, gmesh);
    DivideLowerDimensionalElements(gmesh);

    return gmesh;
}

TPZGeoMesh *CreateCubeGeoMesh(int nCoarseRef, int nInternalRef, TPZStack<int64_t> &coarseIndexes) {

    auto *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(3);
    int matID = 1;

    // Creates matrix with node coordinates
    const int NodeNumber = 8;
    REAL coordinates[NodeNumber][3] = {
        {0., 0., 0.},
        {1., 0., 0.},
        {1., 1., 0.},
        {0., 1., 0.},
        {0., 0., 1.},
        {1., 0., 1.},
        {1., 1., 1.},
        {0., 1., 1.}
    };

    // Inserts coordinates in the TPZGeoMesh object
    for (int i = 0; i < NodeNumber; i++) {
        int64_t nodeID = gmesh->NodeVec().AllocateNewElement();

        TPZVec<REAL> nodeCoord(3);
        nodeCoord[0] = coordinates[i][0];
        nodeCoord[1] = coordinates[i][1];
        nodeCoord[2] = coordinates[i][2];

        gmesh->NodeVec()[nodeID] = TPZGeoNode(i, nodeCoord, *gmesh);
    }

    // Creates cube element
    TPZManVector<int64_t> nodeIDs(8);
    nodeIDs[0] = 0;
    nodeIDs[1] = 1;
    nodeIDs[2] = 2;
    nodeIDs[3] = 3;
    nodeIDs[4] = 4;
    nodeIDs[5] = 5;
    nodeIDs[6] = 6;
    nodeIDs[7] = 7;
    new TPZGeoElRefPattern<pzgeom::TPZGeoCube>(nodeIDs, matID, *gmesh);

    // Creates boundary faces
    nodeIDs.Resize(4);
    matID = -1;
    nodeIDs[0] = 0;
    nodeIDs[1] = 1;
    nodeIDs[2] = 2;
    nodeIDs[3] = 3;
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeIDs, matID, *gmesh);
    nodeIDs[0] = 4;
    nodeIDs[1] = 5;
    nodeIDs[2] = 6;
    nodeIDs[3] = 7;
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeIDs, matID, *gmesh);

    for (int i = 0; i < 4; i++) {
        nodeIDs[0] = (0 + i) % 4;
        nodeIDs[1] = (1 + i) % 4;
        nodeIDs[2] = nodeIDs[1] + 4;
        nodeIDs[3] = nodeIDs[0] + 4;
        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeIDs, matID, *gmesh);
    }

    gmesh->BuildConnectivity();

    UniformRefinement(nCoarseRef, gmesh);
    DivideLowerDimensionalElements(gmesh);

    int64_t nElem = gmesh->NElements();

    coarseIndexes.clear();
    for (int64_t i = 0; i < nElem; i++) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (gel->Dimension() != gmesh->Dimension() || gel->NSubElements() > 0) continue;
        coarseIndexes.Push(i);
    }

    UniformRefinement(nInternalRef, gmesh);

    return gmesh;
}

TPZGeoMesh *CreateLShapeGeoMesh(int nCoarseRef, int nInternalRef, TPZStack<int64_t> &mhmIndexes) {

    TPZVec<int> bcIDs(8, -1);
    TPZGeoMesh *gmesh = CreateQuadLShapeMesh(bcIDs);
    gmesh->SetDimension(2);
    gmesh->BuildConnectivity();

    UniformRefinement(nCoarseRef, gmesh);
    DivideLowerDimensionalElements(gmesh);

    int64_t nElem = gmesh->NElements();
    for (int64_t i = 0; i < nElem; i++) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (gel->Dimension() != gmesh->Dimension() || gel->NSubElements() > 0) continue;
        mhmIndexes.Push(i);
    }

    UniformRefinement(nInternalRef, gmesh);
    DivideLowerDimensionalElements(gmesh);

    for (int64_t i = 0; i < mhmIndexes.size(); i++) {
        std::cout << mhmIndexes[i] << '\n';
    }
    std::cout << '\n';

    return gmesh;
}

void CreateMHMCompMesh(TPZMHMixedMeshControl *mhm, const ProblemConfig &config, int nInternalRef,
                       bool definePartitionByCoarseIndex, TPZManVector<int64_t> mhmIndexes) {

    if (definePartitionByCoarseIndex) {
        mhm->DefinePartitionbyCoarseIndices(mhmIndexes);
    } else {
        mhm->DefinePartition(mhmIndexes);
    }

    // Indicate material indices to the MHM control structure
    mhm->fMaterialIds = config.materialids;
    mhm->fMaterialBCIds = config.bcmaterialids;

    // Insert the material objects in the multiphysics mesh
    InsertMaterialsInMHMMesh(*mhm, config);

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
    {
        std::string fileName = "CompMesh.txt";
        std::ofstream file(fileName);
        mhm->CMesh()->Print(file);
    }
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

    TLaplaceExample1 *analytic = &config.exact.operator*();
    an.SetExact(analytic->ExactSolution());

    scalnames.Push("Pressure");
    scalnames.Push("ExactPressure");
    vecnames.Push("Flux");
    vecnames.Push("ExactFlux");

    int resolution = 0;
    std::string plotname = config.dir_name + "/" + config.problemname + "Results.vtk";
    an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotname);
    an.PostProcess(resolution, cmesh->Dimension());

    //TPZManVector<REAL> errors(4, 0.);
    //an.SetThreadsForError(2);
    //an.SetExact(analytic->ExactSolution());
    //an.PostProcessError(errors, false);
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

void InsertMaterialsInMHMMesh(TPZMHMixedMeshControl &control, const ProblemConfig &config) {
    TPZCompMesh &cmesh = control.CMesh();

    int dim = control.GMesh()->Dimension();
    cmesh.SetDimModel(dim);

    TPZMixedPoisson *mat = new TPZMixedPoisson(1, dim);

    TPZFMatrix<REAL> K(3, 3, 0), invK(3, 3, 0);
    K.Identity();
    invK.Identity();

    mat->SetForcingFunctionExact(config.exact.operator*().Exact());
    mat->SetForcingFunction(config.exact.operator*().ForcingFunction());
    mat->SetPermeabilityTensor(K, invK);

    cmesh.InsertMaterialObject(mat);

    for (auto matid : config.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(config.exact.operator*().Exact());
        cmesh.InsertMaterialObject(bc);
    }
}

void InsertUNISIMMaterialsInMHMMesh(TPZMHMixedMeshControl &control, const ProblemConfig &config) {
    TPZCompMesh &cmesh = control.CMesh();

    int dim = control.GMesh()->Dimension();
    cmesh.SetDimModel(dim);

    TPZMixedPoisson *mat = new TPZMixedPoisson(1, dim);

    TPZFMatrix<REAL> K(3, 3, 0), invK(3, 3, 0);
    K.Identity();
    invK.Identity();

    mat->SetPermeabilityTensor(K, invK);

    cmesh.InsertMaterialObject(mat);

    int bctype = 0;
    TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);

    val1(0, 0) = 10;
    TPZBndCond *bc = mat->CreateBC(mat, -1, bctype, val1, val2);
    cmesh.InsertMaterialObject(bc);

    val1(0, 0) = -20;
    bc = mat->CreateBC(mat, -2, bctype, val1, val2);
    cmesh.InsertMaterialObject(bc);
}
