//
// Created by Gustavo Batistela on 3/31/21.
//

#include "pzgmesh.h"
#include <Pre/TPZGenGrid3D.h>
#include <Pre/TPZMHMixedMeshControl.h>
#include <TPZMFSolutionTransfer.h>
//#include <Tools.h>
#include <ToolsMHM.h>
#include <Util/pzlog.h>
#include "TPZLinearAnalysis.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "Elasticity/TPZMixedElasticityND.h"
#include "TPZBndCondT.h"
#include "Pre/TPZMHMeshControl.h"
#include "TPZElasticityMHMHDivErrorEstimator.h"

void RunSmoothProblem(int nCoarseDiv, int nInternalRef);
void RunElasticityProblem(int nCoarseDiv, int nInternalRef);
void RunNonConvexProblem();
void RunHighGradientProblem(int nCoarseDiv, int nInternalRef);
void RunInnerSingularityProblem(int nCoarseDiv, int nInternalRef);
void RunPeriodicPermProblem(int nCoarseDiv, int nInternalRef);

TPZGeoMesh *CreateQuadGeoMesh(int nCoarseDiv, int nInternalRef);
TPZGeoMesh *CreateLShapeGeoMesh(int nCoarseRef, int nInternalRef, TPZStack<int64_t> &mhmIndexes);

void InsertMaterialsInMHMMesh(TPZMHMixedMeshControl &control, const ProblemConfig &config);
void CreateMHMCompMesh(TPZMHMixedMeshControl *mhm, const ProblemConfig &config, int nInternalRef,
                       bool definePartitionByCoarseIndex, TPZManVector<int64_t>& mhmIndexes);
void CreateMHMCompMeshHeteroPerm(TPZMHMixedMeshControl *mhm, const ProblemConfig &config, int nInternalRef,
                                 bool definePartitionByCoarseIndex, TPZManVector<int64_t>& mhmIndexes);

void SolveMHMProblem(TPZMHMixedMeshControl *mhm, const ProblemConfig &config);

void EstimateError(ProblemConfig &config, TPZMHMixedMeshControl *mhm);
void EstimateErrorElasticity(ProblemConfig &config, TPZMHMixedMeshControl *mhm);

void MHMAdaptivity(TPZMHMixedMeshControl *mhm, TPZGeoMesh* gmeshToRefine, ProblemConfig& config);

void CreateMHMCompMeshPermFunction(TPZMHMixedMeshControl &mhm);
STATE PeriodicPermeabilityFunction(const TPZVec<REAL> &coord);
void PeriodicProblemForcingFunction(const TPZVec <REAL> &pt, TPZVec <STATE> &result);

int main() {
    TPZLogger::InitializePZLOG();
    gRefDBase.InitializeAllUniformRefPatterns();

    // const std::set<int> nCoarseDiv = {3, 4, 5, 6};
    const std::set<int> nCoarseDiv = {3};
    const std::set<int> nInternalRef = {0};
    // const std::set<int> nInternalRef = {0, 1, 2, 3};
    for (const auto coarse_div : nCoarseDiv) {
        for (const auto internal_ref : nInternalRef) {
            // RunSmoothProblem(coarse_div, internal_ref);
            RunElasticityProblem(coarse_div, internal_ref);
            // RunInnerSingularityProblem(coarse_div, internal_ref);
            //RunHighGradientProblem(coarse_div, internal_ref);
            //RunPeriodicPermProblem(coarse_div, internal_ref);
        }
    }
    // RunNonConvexProblem();

    return 0;
}

void RunSmoothProblem(const int nCoarseDiv, const int nInternalRef) {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
    config.problemname = "Smooth";
    config.dir_name = "Journal";
    config.porder = 1;
    config.hdivmais = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    config.ndivisions = nCoarseDiv;
    config.ninternalref = nInternalRef;
    config.gmesh = CreateQuadGeoMesh(nCoarseDiv, nInternalRef);

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "-GeoMesh.vtk";
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

void RunElasticityProblem(const int nCoarseDiv, const int nInternalRef) {
    ProblemConfig config; 
    config.dimension = 2;
    config.exactElast = new TElasticity2DAnalytic;
    config.exactElast.operator*().fProblemType = TElasticity2DAnalytic::EDispx;
    config.problemname = "Elasticity";
    config.dir_name = "Journal";
    config.porder = 1;
    config.hdivmais = 1;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;
    config.problemtype = ProblemConfig::TProbType::EElasticity;

    config.ndivisions = nCoarseDiv;
    config.ninternalref = nInternalRef;
    config.gmesh = CreateQuadGeoMesh(nCoarseDiv, nInternalRef);

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "-GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    TPZManVector<int64_t> coarseIndexes;
    ComputeCoarseIndices(config.gmesh, coarseIndexes);
    bool definePartitionByCoarseIndexes = true;
    mhm->SetProblemType(TPZMHMeshControl::MProblemType::EElasticity2D);
    CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, coarseIndexes);

    SolveMHMProblem(mhm, config);
    EstimateErrorElasticity(config, mhm);
}

void RunHighGradientProblem(const int nCoarseDiv, const int nInternalRef) {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::EBoundaryLayer;
    config.problemname = "HighGradient";
    config.dir_name = "Journal";
    config.porder = 1;
    config.hdivmais = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    config.ndivisions = nCoarseDiv;
    config.ninternalref = nInternalRef;
    config.gmesh = CreateQuadGeoMesh(nCoarseDiv, nInternalRef);

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "-GeoMesh.vtk";
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
    config.dir_name = "Journal";
    config.porder = 1;
    config.hdivmais = 3;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nDiv = 1;
    int nInternalRef = 0;

    TPZManVector<int64_t> coarseIndexes;
    config.ndivisions = nDiv;
    config.gmesh = CreateLMHMMesh(nDiv, coarseIndexes);

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "-GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    bool definePartitionByCoarseIndexes = false;
    CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, coarseIndexes);

    SolveMHMProblem(mhm, config);
    EstimateError(config, mhm);
}

void RunInnerSingularityProblem(const int nCoarseDiv, const int nInternalRef) {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESteklovNonConst;
    config.problemname = "InnerSingularity";
    config.dir_name = "Journal";
    config.porder = 1;
    config.hdivmais = 2;
    config.materialids = {1, 2};
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    config.ndivisions = nCoarseDiv;
    config.ninternalref = nInternalRef;

    TPZManVector<REAL> x0(3, -1.), x1(3, 1.);
    x1[2] = 0.;
    TPZGenGrid2D gen(TPZManVector<int,2>(2, nCoarseDiv), x0, x1, 1, 0);

    gen.SetRefpatternElements(true);
    config.gmesh = new TPZGeoMesh;
    gen.Read(config.gmesh);

    TPZManVector<int, 4> bcIDs(4, -1);
    gen.SetBC(config.gmesh, 4, bcIDs[0]);
    gen.SetBC(config.gmesh, 5, bcIDs[1]);
    gen.SetBC(config.gmesh, 6, bcIDs[2]);
    gen.SetBC(config.gmesh, 7, bcIDs[3]);

    config.gmesh->SetDimension(2);

    Tools::UniformRefinement(nInternalRef, config.gmesh);
    Tools::DivideLowerDimensionalElements(config.gmesh);

    for (int i = 0; i < config.gmesh->NElements(); i++) {
        TPZGeoEl * gel = config.gmesh->Element(i);
        if (gel->HasSubElement()) continue;
        if (gel->Dimension() != 2) continue;

        TPZManVector<REAL,2> qsi(2,0.);
        TPZManVector<REAL,3> result(3,0.);
        gel->X(qsi, result);

        if (result[0] * result[1] < 0) {
            gel->SetMaterialId(2);
        }
    }

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "-GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    TPZManVector<int64_t> coarseIndexes;
    ComputeCoarseIndices(config.gmesh, coarseIndexes);
    bool definePartitionByCoarseIndexes = true;
    CreateMHMCompMeshHeteroPerm(mhm, config, nInternalRef, definePartitionByCoarseIndexes, coarseIndexes);

    SolveMHMProblem(mhm, config);
    EstimateError(config, mhm);
}

void RunPeriodicPermProblem(const int nCoarseDiv, const int nInternalRef) {
    ProblemConfig config;
    config.dimension = 2;
    // TODO: do we know the exact solution?
    //config.exact = new TLaplaceExample1;
    //config.exact.operator*().fExact = TLaplaceExample1::EBoundaryLayer;
    config.problemname = "PeriodicPerm";
    config.dir_name = "Journal";
    config.porder = 1;
    config.hdivmais = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    config.ndivisions = nCoarseDiv;
    config.ninternalref = nInternalRef;
    config.gmesh = CreateQuadGeoMesh(nCoarseDiv, nInternalRef);

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "-GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    CreateMHMCompMeshPermFunction(*mhm);

    SolveMHMProblem(mhm, config);
    EstimateError(config, mhm);
}

TPZGeoMesh *CreateQuadGeoMesh(int nCoarseDiv, int nInternalRef) {

    TPZManVector<int, 4> bcIDs(4, -1);
    TPZGeoMesh *gmesh = Tools::CreateGeoMesh(nCoarseDiv, bcIDs);
    gmesh->SetDimension(2);

    Tools::UniformRefinement(nInternalRef, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

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

    Tools::UniformRefinement(nCoarseRef, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    int64_t nElem = gmesh->NElements();

    coarseIndexes.clear();
    for (int64_t i = 0; i < nElem; i++) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (gel->Dimension() != gmesh->Dimension() || gel->NSubElements() > 0) continue;
        coarseIndexes.Push(i);
    }

    Tools::UniformRefinement(nInternalRef, gmesh);

    return gmesh;
}

TPZGeoMesh *CreateLShapeGeoMesh(int nCoarseRef, int nInternalRef, TPZStack<int64_t> &mhmIndexes) {

    TPZVec<int> bcIDs(8, -1);
    TPZGeoMesh *gmesh = Tools::CreateQuadLShapeMesh(bcIDs);
    gmesh->SetDimension(2);
    gmesh->BuildConnectivity();

    Tools::UniformRefinement(nCoarseRef, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    int64_t nElem = gmesh->NElements();
    for (int64_t i = 0; i < nElem; i++) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (gel->Dimension() != gmesh->Dimension() || gel->NSubElements() > 0) continue;
        mhmIndexes.Push(i);
    }

    Tools::UniformRefinement(nInternalRef, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    for (int64_t i = 0; i < mhmIndexes.size(); i++) {
        std::cout << mhmIndexes[i] << '\n';
    }
    std::cout << '\n';

    return gmesh;
}

void CreateMHMCompMesh(TPZMHMixedMeshControl *mhm, const ProblemConfig &config, int nInternalRef,
                       bool definePartitionByCoarseIndex, TPZManVector<int64_t>& mhmIndexes) {

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
    mhm->DivideSkeletonElements(0);
    mhm->DivideBoundarySkeletonElements();
    // Creates MHM mesh
    bool substructure = false;
    mhm->BuildComputationalMesh(substructure);
}

void SolveMHMProblem(TPZMHMixedMeshControl *mhm, const ProblemConfig &config) {

    TPZAutoPointer<TPZCompMesh> cmesh = mhm->CMesh();

    TPZLinearAnalysis an(cmesh, RenumType::ESloan);

#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<> strmat(cmesh.operator->());
    strmat.SetNumThreads(8);
#else
    TPZSkylineStructMatrix<STATE> strmat(cmesh.operator->());
    strmat.SetNumThreads(0);
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

    // TODO: Phil, this is the part that needs a fix!
    TPZMFSolutionTransfer transfer;
    transfer.BuildTransferData(cmesh.operator->());
    transfer.TransferFromMultiphysics();

    TPZStack<std::string> scalnames, vecnames;
    TPZMaterial *mat = cmesh->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }

    if (config.problemtype == ProblemConfig::TProbType::EDarcy){
        scalnames.Push("Pressure");
        scalnames.Push("Permeability");
        vecnames.Push("Flux");

        if (config.exact) {
            scalnames.Push("ExactPressure");
            vecnames.Push("ExactFlux");
        }
    } else if (config.problemtype == ProblemConfig::TProbType::EElasticity){
        scalnames.Push("SigmaX");
        scalnames.Push("SigmaY");
        scalnames.Push("TauXY");
        vecnames.Push("Displacement");

        if (config.exactElast) {
            vecnames.Push("ExactDisplacement");
        }
    }
    
    std::cout << "Post Processing...\n";

    int resolution = 2;
    std::stringstream plotname;
    plotname << config.dir_name << "/" << config.problemname << "-" << config.ndivisions << "-" << config.ninternalref
             << "-Results.vtk";
    an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotname.str());
    an.PostProcess(resolution, cmesh->Dimension());
}

void EstimateError(ProblemConfig &config, TPZMHMixedMeshControl *mhm) {

    std::cout << "\nError Estimation processing for MHM-Hdiv problem " << std::endl;

    // Error estimation
    TPZMultiphysicsCompMesh *originalMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm->CMesh().operator->());
    if (!originalMesh) DebugStop();

    bool postProcWithHDiv = false;
    TPZDarcyMHMHDivErrorEstimator ErrorEstimator(*originalMesh, mhm, postProcWithHDiv);
    ErrorEstimator.SetAnalyticSolution(config.exact);
    
    ErrorEstimator.PrimalReconstruction();

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    TPZManVector<REAL, 6> errors;
    TPZManVector<REAL> elementerrors;
    std::stringstream outVTK;
    outVTK << config.dir_name << "/" << config.problemname << "-" << config.ndivisions << "-" << config.ninternalref
           << "-Errors.vtk";
    std::string outVTKstring = outVTK.str();
    ErrorEstimator.ComputeErrors(errors, elementerrors, outVTKstring);

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "-GlobalErrors.txt";
        std::ofstream file(fileName, std::ios::app);
        Tools::PrintErrors(file, config, errors);
    }
}


void EstimateErrorElasticity(ProblemConfig &config, TPZMHMixedMeshControl *mhm) {

    std::cout << "\nError Estimation processing for MHM-Hdiv problem " << std::endl;

    // Error estimation
    TPZMultiphysicsCompMesh *originalMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm->CMesh().operator->());
    if (!originalMesh) DebugStop();

    bool postProcWithHDiv = false;
    TPZElasticityMHMHDivErrorEstimator ErrorEstimator(*originalMesh, mhm, postProcWithHDiv);
    ErrorEstimator.SetAnalyticSolution(config.exactElast);
    
    ErrorEstimator.PrimalReconstruction();

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    TPZManVector<REAL, 6> errors;
    TPZManVector<REAL, 6> elementerrors;
    std::stringstream outVTK;
    outVTK << config.dir_name << "/" << config.problemname << "-" << config.ndivisions << "-" << config.ninternalref
           << "-Errors.vtk";
    std::string outVTKstring = outVTK.str();
    ErrorEstimator.ComputeErrors(errors, elementerrors, outVTKstring);

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "-GlobalErrors.txt";
        std::ofstream file(fileName, std::ios::app);
        Tools::PrintElasticityErrors(file, config, errors);
    }
}

void InsertMaterialsInMHMMesh(TPZMHMixedMeshControl &control, const ProblemConfig &config) {
    TPZCompMesh &cmesh = control.CMesh();

    int dim = control.GMesh()->Dimension();
    cmesh.SetDimModel(dim);

    TPZMaterial *mat;
    if (config.problemtype == ProblemConfig::TProbType::EElasticity){
        double E= 1000.;
        double nu = 0.3;
        auto matelast = new TPZMixedElasticityND(1,E,nu,0,0,1,dim);
        matelast->SetExactSol(config.exactElast->ExactSolution(),5);
        matelast->SetForcingFunction(config.exactElast->ForceFunc(),5);
        mat = matelast;
    } else {
        auto matdarcy = new TPZMixedDarcyFlow(1, dim);
        matdarcy->SetExactSol(config.exact->ExactSolution(),5);
        matdarcy->SetForcingFunction(config.exact->ForceFunc(),5);
        matdarcy->SetConstantPermeability(1.);
        mat=matdarcy;
    }

    cmesh.InsertMaterialObject(mat);

    for (auto matid : config.bcmaterialids) {
        
        auto *matdarcy = dynamic_cast<TPZMixedDarcyFlow *> (mat);
        auto *matelast = dynamic_cast<TPZMixedElasticityND *> (mat);
        if (matdarcy){
            TPZFNMatrix<1, REAL> val1(1, 1, 0.);
            TPZManVector<REAL> val2(1, 0.);
            int bctype = 0;
            TPZBndCondT<STATE> *bc = matdarcy->CreateBC(mat, matid, bctype, val1, val2);
            bc->SetForcingFunctionBC(config.exact->ExactSolution(),4);
            cmesh.InsertMaterialObject(bc);
        }
        if (matelast){
            TPZFNMatrix<1, REAL> val1(dim, dim, 0.);
            TPZManVector<REAL> val2(dim, 0.);
            int bctype = 0;
            TPZBndCondT<STATE> *bc = matelast->CreateBC(mat, matid, bctype, val1, val2);
            bc->SetForcingFunctionBC(config.exactElast->ExactSolution(),4);
            cmesh.InsertMaterialObject(bc);
        }
    }
}

void CreateMHMCompMeshHeteroPerm(TPZMHMixedMeshControl *mhm, const ProblemConfig &config, int nInternalRef,
                       bool definePartitionByCoarseIndex, TPZManVector<int64_t>& mhmIndexes) {

    if (definePartitionByCoarseIndex) {
        mhm->DefinePartitionbyCoarseIndices(mhmIndexes);
    } else {
        mhm->DefinePartition(mhmIndexes);
    }

    // Indicate material indices to the MHM control structure
    mhm->fMaterialIds = config.materialids;
    mhm->fMaterialBCIds = config.bcmaterialids;

    // Insert the material objects in the multiphysics mesh
    TPZCompMesh &cmesh = mhm->CMesh();

    int dim = mhm->GMesh()->Dimension();
    cmesh.SetDimModel(dim);

    auto *mat = new TPZMixedDarcyFlow(1, dim);

\
    mat->SetExactSol(config.exact->ExactSolution(),3);
    
    mat->SetForcingFunction(config.exact->ForceFunc(),3);
    mat->SetConstantPermeability(15.);

    cmesh.InsertMaterialObject(mat);

    auto *mat2 = new TPZMixedDarcyFlow(2, dim);


    mat2->SetExactSol(config.exact->ExactSolution(),1);
    mat2->SetForcingFunction(config.exact->ForceFunc(),2);
    mat2->SetConstantPermeability(1.);

    cmesh.InsertMaterialObject(mat2);

    for (auto matid : config.bcmaterialids) {
        TPZFNMatrix<1, STATE> val1(1, 1, 0.);
        TPZManVector<STATE,1> val2(1, 0.);
        int bctype = 0;
        TPZBndCondT<STATE> *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->SetForcingFunctionBC(config.exact->ExactSolution(),4);
        cmesh.InsertMaterialObject(bc);
    }

    // General approximation order settings
    mhm->SetInternalPOrder(config.porder);
    mhm->SetSkeletonPOrder(config.porder);
    mhm->SetHdivmaismaisPOrder(config.hdivmais);

    // Refine skeleton elements
    mhm->DivideSkeletonElements(0);
    mhm->DivideBoundarySkeletonElements();
    // Creates MHM mesh
    bool substructure = false;
    mhm->BuildComputationalMesh(substructure);
}

void CreateMHMCompMeshPermFunction(TPZMHMixedMeshControl &mhm) {

    TPZGeoMesh *gmesh = mhm.GMesh().operator->();
    TPZManVector<int64_t, 22 * 6> coarse_indexes;
    ComputeCoarseIndices(gmesh, coarse_indexes);

    int nInternalRef = 0;
    Tools::UniformRefinement(nInternalRef, 2, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    mhm.DefinePartitionbyCoarseIndices(coarse_indexes);

    // Indicate material indices to the MHM control structure
    mhm.fMaterialIds = {1};
    mhm.fMaterialBCIds = {-1};

    // Insert the material objects in the multiphysics mesh
    TPZCompMesh *cmesh = mhm.CMesh().operator->();
    auto *mix = new TPZMixedDarcyFlow(1, cmesh->Dimension());
    PermeabilityFunctionType permfunc;
    permfunc = PeriodicPermeabilityFunction;
    mix->SetPermeabilityFunction(permfunc);
    mix->SetForcingFunction(PeriodicProblemForcingFunction, 3);

    TPZFNMatrix<1, REAL> val1(1, 1, 0.);
    TPZManVector<REAL,1> val2(1, 0.);
    constexpr int dirichlet_bc = 0;
    TPZBndCondT<STATE> *pressure_left = mix->CreateBC(mix, -1, dirichlet_bc, val1, val2);

    cmesh->InsertMaterialObject(mix);
    cmesh->InsertMaterialObject(pressure_left);

    // General approximation order settings
    mhm.SetInternalPOrder(2);
    mhm.SetSkeletonPOrder(2);
    mhm.SetHdivmaismaisPOrder(2);

    // Refine skeleton elements
    mhm.DivideSkeletonElements(0);
    mhm.DivideBoundarySkeletonElements();

    // Creates MHM mesh
    bool substructure = true;
    mhm.BuildComputationalMesh(substructure);
}

STATE PeriodicPermeabilityFunction(const TPZVec<REAL> &coord) {

    constexpr auto epsilon = 0.04;
    constexpr auto P = 1.8;
    const auto x = coord[0];
    const auto y = coord[1];

    REAL term_1 = 2 + P * cos(2 * M_PI * (x - 0.5) / epsilon);
    REAL term_2 = 2 + P * cos(2 * M_PI * (y - 0.5) / epsilon);

    auto perm = 1 / (term_1 * term_2);
    return perm;
}

void PeriodicProblemForcingFunction(const TPZVec<REAL> &pt, TPZVec<STATE> &result) { result[0] = -1; }
