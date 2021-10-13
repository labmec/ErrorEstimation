//
// Created by Gustavo Batistela on 3/31/21.
//
#include "pzgmesh.h"
#include <Pre/TPZGenGrid3D.h>
#include <Pre/TPZMHMixedMeshControl.h>
#include <TPZMFSolutionTransfer.h>
//#include <Tools.h>
#include <Material/DarcyFlow/TPZMixedDarcyFlow.h>
#include <ToolsMHM.h>
#include <Util/pzlog.h>

struct ErrorResult {
    ErrorResult(const std::string &s, const int nIntRef, const int nCoarseDiv, const int kOrder, const int nOrder,
                const REAL estimatedError, const REAL exactError, const REAL ieff) {
        problem_name = s;
        n_internal_ref = nIntRef;
        n_coarse_div = nCoarseDiv;
        k_order = kOrder;
        n_order = nOrder;
        estimated_error = estimatedError;
        exact_error = exactError;
        effectivity_index = ieff;
    }
    std::string problem_name = "NULL";
    int n_internal_ref = -1;
    int n_coarse_div = -1;
    int k_order = -1;
    int n_order = -1;
    REAL estimated_error = -1.;
    REAL exact_error = -1.;
    REAL effectivity_index = -1.;
};

std::vector<ErrorResult> result_vec;

const int porder = 1;
const int hdivmais = 4;
std::string dir_name = "Journal_k" + std::to_string(porder) + "n" + std::to_string(hdivmais);

void RunSmoothProblem(int nCoarseDiv, int nInternalRef, int k, int n);
void RunNonConvexProblem();
void RunHighGradientProblem(int nCoarseDiv, int nInternalRef);
void RunInnerSingularityProblem(int nCoarseDiv, int nInternalRef, int k, int n);
void RunPeriodicPermProblem(int nCoarseDiv, int nInternalRef);
void RunHighGradientAdaptivityProblem(int n_divisions);

TPZGeoMesh *CreateQuadGeoMesh(int nCoarseDiv, int nInternalRef);
TPZGeoMesh *CreateLShapeGeoMesh(int nCoarseRef, int nInternalRef, TPZStack<int64_t> &mhmIndexes);

void InsertMaterialsInMHMMesh(TPZMHMixedMeshControl &control, const ProblemConfig &config);
void CreateMHMCompMesh(TPZMHMixedMeshControl *mhm, const ProblemConfig &config, bool definePartitionByCoarseIndex,
                       TPZManVector<int64_t> &mhmIndexes);

void CreateMHMCompMesh(TPZMHMixedMeshControl *mhm, const ProblemConfig &config, bool definePartitionByCoarseIndex,
                       TPZManVector<int64_t> &mhmIndexes, const std::vector<int> &refLevelPerSubdomain);

void CreateMHMCompMeshHeteroPerm(TPZMHMixedMeshControl *mhm, const ProblemConfig &config, int nInternalRef,
                                 bool definePartitionByCoarseIndex, TPZManVector<int64_t>& mhmIndexes);

void SolveMHMProblem(TPZMHMixedMeshControl *mhm, const ProblemConfig &config);

void EstimateError(ProblemConfig &config, TPZMHMixedMeshControl *mhm);
void EstimateError(ProblemConfig &config, TPZMHMixedMeshControl *mhm, TPZMHMHDivErrorEstimator &estimator);

void MHMAdaptivity(TPZMHMixedMeshControl *mhm, ProblemConfig &config, TPZCompMesh *postProcMesh,
                   std::vector<int> &refLevelPerSubdomain);

void CreateMHMCompMeshPermFunction(TPZMHMixedMeshControl &mhm);
void PeriodicProblemForcingFunction(const TPZVec <REAL> &pt, TPZVec <STATE> &result);

void PrintLatexGraphs(std::ostream & out);

void RunSmoothProblemSuite(const std::set<int> &nCoarseDiv, const std::set<int> &nInternalRef,
                           const std::set<int> &kOrder, const std::set<int> &nOrder);
void RunInnerSingularityProblemSuite(const std::set<int> &nCoarseDiv, const std::set<int> &nInternalRef,
                                     const std::set<int> &kOrder, const std::set<int> &nOrder);



int main() {
    TPZLogger::InitializePZLOG();
    gRefDBase.InitializeAllUniformRefPatterns();

    const std::set<int> nCoarseDiv = {8};
    const std::set<int> nInternalRef = {1};
    const std::set<int> kOrder = {1};
    const std::set<int> nOrder = {2};

    //RunSmoothProblemSuite(nCoarseDiv, nInternalRef, kOrder, nOrder);

    //RunInnerSingularityProblemSuite(nCoarseDiv, nInternalRef, kOrder, nOrder);

    RunHighGradientAdaptivityProblem(4);

    //RunNonConvexProblem();
    //std::ofstream out_latex("LatexGraphs.txt", std::ios::trunc);
    //PrintLatexGraphs(out_latex);

    return 0;
}

void RunInnerSingularityProblemSuite(const std::set<int> &nCoarseDiv, const std::set<int> &nInternalRef,
                                     const std::set<int> &kOrder, const std::set<int> &nOrder) {
    for (const auto k : kOrder) {
        for (const auto n : nOrder) {
            for (const auto internal_ref : nInternalRef) {
                for (const auto coarse_div : nCoarseDiv) {
                    RunInnerSingularityProblem(coarse_div, internal_ref, k, n);
                }
            }
        }
    }
}
void RunSmoothProblemSuite(const std::set<int> &nCoarseDiv, const std::set<int> &nInternalRef,
                           const std::set<int> &kOrder, const std::set<int> &nOrder) {
    for (const auto k : kOrder) {
        for (const auto n : nOrder) {
            for (const auto internal_ref : nInternalRef) {
                for (const auto coarse_div : nCoarseDiv) {
                    RunSmoothProblem(coarse_div, internal_ref, k, n);
                }
            }
        }
    }
}

void RunSmoothProblem(const int nCoarseDiv, const int nInternalRef, const int k, const int n) {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
    config.problemname = "Smooth";
    config.dir_name = dir_name;
    config.porder = k;
    config.hdivmais = n;
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
    CreateMHMCompMesh(mhm, config, definePartitionByCoarseIndexes, coarseIndexes);

    SolveMHMProblem(mhm, config);
    TPZCompMesh *cmesh = nullptr;
    EstimateError(config, mhm);
}

void RunHighGradientProblem(const int nCoarseDiv, const int nInternalRef) {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::EBoundaryLayer;
    config.problemname = "HighGradient";
    config.dir_name = dir_name;
    config.porder = porder;
    config.hdivmais = hdivmais;
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
    CreateMHMCompMesh(mhm, config, definePartitionByCoarseIndexes, coarseIndexes);

    SolveMHMProblem(mhm, config);
    EstimateError(config, mhm);
}

void RunNonConvexProblem() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
    config.problemname = "NonConvex";
    config.dir_name = dir_name;
    config.porder = porder;
    config.hdivmais = hdivmais;
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
    CreateMHMCompMesh(mhm, config, definePartitionByCoarseIndexes, coarseIndexes);

    SolveMHMProblem(mhm, config);
    TPZCompMesh *cmesh = nullptr;
    EstimateError(config, mhm);
}

void RunInnerSingularityProblem(const int nCoarseDiv, const int nInternalRef, int k, int n) {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESteklovNonConst;
    config.problemname = "InnerSingularity";
    config.dir_name = dir_name;
    config.porder = k;
    config.hdivmais = n;
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
    TPZCompMesh *cmesh = nullptr;
    EstimateError(config, mhm);

}

void RunPeriodicPermProblem(const int nCoarseDiv, const int nInternalRef) {
    ProblemConfig config;
    config.dimension = 2;
    // TODO: do we know the exact solution?
    //config.exact = new TLaplaceExample1;
    //config.exact.operator*().fExact = TLaplaceExample1::EBoundaryLayer;
    config.problemname = "PeriodicPerm";
    config.dir_name = dir_name;
    config.porder = porder;
    config.hdivmais = hdivmais;
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

void CreateMHMCompMesh(TPZMHMixedMeshControl *mhm, const ProblemConfig &config, bool definePartitionByCoarseIndex,
                       TPZManVector<int64_t> &mhmIndexes, const std::vector<int> &refLevelPerSubdomain) {

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

    const auto interfaces = mhm->GetInterfaces();
    if (!refLevelPerSubdomain.empty()) {
        for (auto interface : interfaces) {
            if (interface.first == interface.second.first || interface.first == interface.second.second) continue;
            auto right_ref_level = refLevelPerSubdomain[interface.second.first];
            auto left_ref_level = refLevelPerSubdomain[interface.second.second];
            const auto n_divisions = std::max(right_ref_level, left_ref_level);
            mhm->DivideSkeletonElement(interface.first, n_divisions);
        }
    }

    mhm->DivideBoundarySkeletonElements();
    // Creates MHM mesh
    bool substructure = true;
    mhm->BuildComputationalMesh(substructure);
}

void CreateMHMCompMesh(TPZMHMixedMeshControl *mhm, const ProblemConfig &config,
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
    bool substructure = true;
    mhm->BuildComputationalMesh(substructure);
}

void SolveMHMProblem(TPZMHMixedMeshControl *mhm, const ProblemConfig &config) {

    TPZAutoPointer<TPZCompMesh> cmesh = mhm->CMesh();

    bool shouldrenumber = true;
    TPZLinearAnalysis an(cmesh, shouldrenumber);

#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> strmat(cmesh.operator->());
    strmat.SetNumThreads(4);
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

    TPZMFSolutionTransfer transfer;
    transfer.BuildTransferData(cmesh.operator->());
    transfer.TransferFromMultiphysics();

    TPZStack<std::string> scalnames, vecnames;
    TPZMaterial *mat = cmesh->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }

    scalnames.Push("Pressure");
    scalnames.Push("Permeability");
    vecnames.Push("Flux");

    if (config.exact) {
        scalnames.Push("ExactPressure");
        vecnames.Push("ExactFlux");
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
    TPZMHMHDivErrorEstimator ErrorEstimator(*originalMesh, mhm, postProcWithHDiv);
    ErrorEstimator.SetAnalyticSolution(config.exact);
    ErrorEstimator.PotentialReconstruction();

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

    result_vec.emplace_back(config.problemname,
                            config.ninternalref,
                            config.ndivisions,
                            config.porder,
                            config.hdivmais,
                            errors[3],
                            errors[2],
                            sqrt(errors[4] + errors[3]) / sqrt(errors[2])
                            );
}

void InsertMaterialsInMHMMesh(TPZMHMixedMeshControl &control, const ProblemConfig &config) {
    TPZCompMesh &cmesh = control.CMesh();

    int dim = control.GMesh()->Dimension();
    cmesh.SetDimModel(dim);

    auto *mat = new TPZMixedDarcyFlow(1, dim);

    auto ff_lambda = [config](const TPZVec<REAL> &loc, TPZVec<STATE> &result) {
        config.exact.operator*().ForcingFunction()->Execute(loc, result);
    };
    auto exact_sol_lambda = [config](const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv) {
        config.exact.operator*().Exact()->Execute(loc, result, deriv);
    };

    mat->SetForcingFunction(ff_lambda, 5);
    mat->SetExactSol(exact_sol_lambda, 5);
    mat->SetConstantPermeability(1);

    cmesh.InsertMaterialObject(mat);

    for (auto matid : config.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.);
        TPZManVector<REAL, 1> val2(1, 0.);
        int bctype = 0;
        auto *bc = mat->CreateBC(mat, matid, bctype, val1, val2);

        auto ff_bc_lambda = [config](const TPZVec<REAL> &loc,
                                     TPZVec<STATE> &rhsVal,
                                     TPZFMatrix<STATE> &matVal) {
            config.exact.operator*().Exact()->Execute(loc, rhsVal, matVal);
        };
        bc->SetForcingFunctionBC(ff_bc_lambda);
        cmesh.InsertMaterialObject(bc);
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

    auto ff_lambda = [config](const TPZVec<REAL> &loc, TPZVec<STATE> &result) {
        config.exact.operator*().ForcingFunction()->Execute(loc, result);
    };
    auto exact_sol_lambda = [config](const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv) {
        config.exact.operator*().Exact()->Execute(loc, result, deriv);
    };

    mat->SetForcingFunction(ff_lambda, 5);
    mat->SetExactSol(exact_sol_lambda, 5);
    mat->SetConstantPermeability(5);

    cmesh.InsertMaterialObject(mat);

    auto *mat2 = new TPZMixedDarcyFlow(2, dim);
    mat2->SetForcingFunction(ff_lambda, 5);
    mat2->SetExactSol(exact_sol_lambda, 5);
    mat2->SetConstantPermeability(1);

    cmesh.InsertMaterialObject(mat2);

    for (auto matid : config.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.);
        TPZManVector<REAL, 1> val2(1, 0.);
        int bctype = 0;
        auto *bc = mat->CreateBC(mat, matid, bctype, val1, val2);

        auto ff_bc_lambda = [config](const TPZVec<REAL> &loc,
                                     TPZVec<STATE> &rhsVal,
                                     TPZFMatrix<STATE> &matVal) {
            config.exact.operator*().Exact()->Execute(loc, rhsVal, matVal);
        };
        bc->SetForcingFunctionBC(ff_bc_lambda);
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
    bool substructure = true;
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

    std::function<STATE(const TPZVec<REAL> &coord)> perm_function = [](const TPZVec<REAL> &coord) {
        constexpr auto epsilon = 0.04;
        constexpr auto P = 1.8;
        const auto x = coord[0];
        const auto y = coord[1];

        REAL term_1 = 2 + P * cos(2 * M_PI * (x - 0.5) / epsilon);
        REAL term_2 = 2 + P * cos(2 * M_PI * (y - 0.5) / epsilon);

        auto perm = 1 / (term_1 * term_2);
        return perm;
    };

    mix->SetPermeabilityFunction(perm_function);
    mix->SetForcingFunction(PeriodicProblemForcingFunction, 1);

    TPZFNMatrix<1, REAL> val1(1, 1, 0.);
    TPZManVector<REAL, 1> val2(1, 0.);
    constexpr int dirichlet_bc = 0;
    TPZBndCond *pressure_left = mix->CreateBC(mix, -1, dirichlet_bc, val1, val2);

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

void PeriodicProblemForcingFunction(const TPZVec<REAL> &pt, TPZVec<STATE> &result) { result[0] = -1; }

void PrintLatexGraphs(std::ostream & out) {

    std::stringstream latex_text;
    if (result_vec.empty()) DebugStop();

    auto it = result_vec.begin();
    while (it != result_vec.end()) {
        std::string current_problem = it->problem_name;
        latex_text << "\\section{" << it->problem_name << "}\n";
        while (it->problem_name == current_problem) {

            int previous_k = it->k_order;
            while (it->k_order == previous_k) {
                int previous_n = it->n_order;
                while (it->n_order == previous_n) {
                    latex_text << "  \\subsection{$k = " << it->k_order << ", n = " << it->n_order << "$}\n";

                    std::stringstream error_graph, ieff_graph, estimated_error_pts, exact_error_pts,
                        estimated_error_legend, exact_error_legend, ieff_error_legend;

                    error_graph << "    \\begin{figure}[ht!]\n"
                                   "      \\centering\n"
                                   "      \\begin{tikzpicture}\n"
                                   "      \\begin{loglogaxis}[\n"
                                   "        width=0.6\\textwidth,\n"
                                   "        height=0.45\\textwidth,\n"
                                   "        xlabel={$h_{sk}$},\n"
                                   "        xtick scale label code/.code={},\n"
                                   "        legend pos=outer north east,\n"
                                   "        grid style=dashed,\n"
                                   "        cycle list name=mycycle,\n"
                                   "      ]\n";

                    ieff_graph << "    \\begin{figure}[ht!]\n"
                                  "      \\centering\n"
                                  "      \\begin{tikzpicture}\n"
                                  "      \\begin{semilogxaxis}[\n"
                                  "        width=0.6\\textwidth,\n"
                                  "        height=0.4\\textwidth,\n"
                                  "        xlabel={$h_{sk}$},\n"
                                  "        %ymin=1.0, ymax=1.1,\n"
                                  "        xtick scale label code/.code={},\n"
                                  "        legend pos=outer north east,\n"
                                  "        grid style=dashed,\n"
                                  "        cycle list name=mycycle,\n"
                                  "        y tick label style={\n"
                                  "        /pgf/number format/.cd,\n"
                                  "        fixed,\n"
                                  "        fixed zerofill,\n"
                                  "        precision=2,\n"
                                  "        /tikz/.cd\n"
                                  "      },\n"
                                  "      ]\n";
                    int current_int_ref = it->n_internal_ref;
                    bool is_first = true;
                    while (it->n_internal_ref == current_int_ref) {
                        if (is_first) {
                            estimated_error_pts << "      \\addplot\n      coordinates{\n";
                            exact_error_pts << "      \\addplot\n      coordinates{\n";
                            ieff_graph << "      \\addplot\n      coordinates{\n";
                            is_first = false;
                            estimated_error_legend << "      $h_{in} = h_{sk}/"
                                                   << std::to_string((int)pow(2, it->n_internal_ref))
                                                   << "$ (estimated),\n";
                            exact_error_legend << "      $h_{in} = h_{sk}/"
                                               << std::to_string((int)pow(2, it->n_internal_ref)) << "$ (exact),\n";
                            ieff_error_legend << "      $h_{in} = h_{sk}/"
                                              << std::to_string((int)pow(2, it->n_internal_ref)) << "$,\n";
                        }

                        estimated_error_pts << "        (1/" << it->n_coarse_div << ", " << it->estimated_error
                                            << ")\n";
                        exact_error_pts << "        (1/" << it->n_coarse_div << ", " << it->exact_error << ")\n";
                        ieff_graph << "        (1/" << it->n_coarse_div << ", " << it->effectivity_index << ")\n";

                        previous_n = it->n_order;
                        previous_k = it->k_order;
                        current_int_ref = it->n_internal_ref;
                        current_problem = it->problem_name;
                        if (it != result_vec.end()) {
                            it++;
                            if (it->n_order == previous_n && it->k_order == previous_k &&
                                it->n_internal_ref != current_int_ref) {
                                is_first = true;
                                current_int_ref = it->n_internal_ref;
                                estimated_error_pts << "      };\n";
                                exact_error_pts << "      };\n";
                                ieff_graph << "      };\n";
                            }
                        }
                    }
                    estimated_error_pts << "      };\n";
                    exact_error_pts << "      };\n";
                    ieff_graph << "      };\n";

                    error_graph << estimated_error_pts.str();
                    error_graph << exact_error_pts.str();
                    error_graph << "      \\legend{\n"
                                << estimated_error_legend.str() << exact_error_legend.str()
                                << "      }\n"
                                   "      \\end{loglogaxis}\n"
                                   "      \\end{tikzpicture}\n"
                                   "      \\caption{"
                                << "Estimated and exact errors, $k = " << previous_k << "$, $n = " << previous_n
                                << "$}\n    \\end{figure}\n\n";

                    ieff_graph << "      \\legend{\n"
                               << ieff_error_legend.str()
                               << "      }\n"
                                  "      \\end{semilogxaxis}\n"
                                  "      \\end{tikzpicture}\n"
                                  "      \\caption{"
                                  "Effectivity index, $k = " << previous_k << "$, $n = " << previous_n
                               << "$}\n    \\end{figure}\n\n";

                    latex_text << "    % Error graph\n" << error_graph.str();
                    latex_text << "    % Ieff graph\n" << ieff_graph.str();
                    latex_text << "\\clearpage\n";
                }
            }
        }
    }
    out << latex_text.str();
    std::cout << latex_text.str();
}

void RunHighGradientAdaptivityProblem(const int n_divisions){

    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::EArcTan;
    config.problemname = "EArcTan";
    config.dir_name = "MHMAdaptivity";
    config.porder = 1;
    config.hdivmais = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nCoarseRef = 7;
    int nInternalRef = n_divisions;

    config.ndivisions = nCoarseRef;

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    std::vector<int> refLevelPerSubdomain;
    int max_refinement = 0;
    int adaptivity_step = 0;
    while (max_refinement <= n_divisions) {
        config.gmesh = CreateQuadGeoMesh(nCoarseRef, nInternalRef);

        auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
        TPZManVector<int64_t> coarseIndexes;
        ComputeCoarseIndices(config.gmesh, coarseIndexes);
        bool definePartitionByCoarseIndexes = true;
        CreateMHMCompMesh(mhm, config, definePartitionByCoarseIndexes, coarseIndexes, refLevelPerSubdomain);
        if (max_refinement == 0) {
            const auto n_subdomains = mhm->Coarse_to_Submesh().size();
            refLevelPerSubdomain.resize(n_subdomains, 0);
        }

        {
            std::string fileName = config.dir_name + "/" + config.problemname + "GMesh" + std::to_string(nCoarseRef) +
                                   "x" + std::to_string(nCoarseRef) + std::to_string(adaptivity_step) + ".vtk";
            std::ofstream file(fileName);
            TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
        }

        SolveMHMProblem(mhm, config);

        bool postProcWithHDiv = false;
        TPZMultiphysicsCompMesh *originalMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm->CMesh().operator->());
        if (!originalMesh) DebugStop();
        TPZMHMHDivErrorEstimator estimator(*originalMesh, mhm, postProcWithHDiv);
        estimator.SetAdaptivityStep(adaptivity_step);
        EstimateError(config, mhm, estimator);

        auto *postprocmesh = estimator.PostProcMesh();
        if (!postprocmesh) DebugStop();

        MHMAdaptivity(mhm, config, postprocmesh, refLevelPerSubdomain);
        adaptivity_step++;

        max_refinement = *std::max_element(refLevelPerSubdomain.begin(), refLevelPerSubdomain.end());
        std::cout << "Max ref level: " << max_refinement << '\n';
    }
}

void MHMAdaptivity(TPZMHMixedMeshControl *mhm, ProblemConfig &config, TPZCompMesh *postProcMesh,
                   std::vector<int> &refLevelPerSubdomain) {

    // Column of the flux error estimate on the element solution matrix
    const int fluxErrorEstimateCol = 3;

    TPZMultiphysicsCompMesh *cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm->CMesh().operator->());
    if (!cmesh) DebugStop();

    TPZFMatrix<STATE> &sol = postProcMesh->ElementSolution();
    TPZSolutionMatrix &elsol = postProcMesh->ElementSolution();
    int64_t nelem = elsol.Rows();

    // Iterates through element errors to get the maximum value
    STATE maxError = 0.;
    for (int64_t iel = 0; iel < nelem; iel++) {
        auto * submesh = dynamic_cast<TPZSubCompMesh*>(postProcMesh->ElementVec()[iel]);
        if (!submesh) continue;

        STATE submeshError = sol(iel, fluxErrorEstimateCol);
        if (submeshError > maxError) {
            maxError = submeshError;
        }
    }

    std::cout << "Max error: " << maxError << "\n";

    // Refines elements which error are bigger than 30% of the maximum error
    REAL threshold = 0.5 * maxError;
    const auto geoToMHM = mhm->GetGeoToMHMDomain();
    const auto interfaces = mhm->GetInterfaces();

    for (int64_t iel = 0; iel < nelem; iel++) {
        auto * submesh = dynamic_cast<TPZSubCompMesh*>(postProcMesh->ElementVec()[iel]);
        if (!submesh) continue;

        STATE submeshError = sol(iel, fluxErrorEstimateCol);
        if (submeshError > threshold) {
            TPZGeoEl * gel = submesh->Element(0)->Reference();
            const auto submesh_id = geoToMHM[gel->Index()];
            refLevelPerSubdomain[submesh_id]++;
            std::cout << "Refining submesh " << submesh_id << " which error is " << submeshError << ".\n";
        }
    }

    for (auto i = 0; i < refLevelPerSubdomain.size(); i++) {
        if (refLevelPerSubdomain[i] != 0) {
            std::cout << i << ": " << refLevelPerSubdomain[i] << '\n';
        }
    }
    std::cout << "Done!\n";
}

void EstimateError(ProblemConfig &config, TPZMHMixedMeshControl *mhm, TPZMHMHDivErrorEstimator& estimator) {

    std::cout << "\nError Estimation processing for MHM-Hdiv problem " << std::endl;

    // Error estimation
    TPZMultiphysicsCompMesh *originalMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm->CMesh().operator->());
    if (!originalMesh) DebugStop();

    estimator.SetAnalyticSolution(config.exact);
    estimator.PotentialReconstruction();

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    TPZManVector<REAL, 6> errors;
    TPZManVector<REAL> elementerrors;
    std::stringstream outVTK;
    outVTK << config.dir_name << "/" << config.problemname << "-" << config.ndivisions << "-" << config.ninternalref
           << "-Errors.vtk";
    std::string outVTKstring = outVTK.str();
    estimator.ComputeErrors(errors, elementerrors, outVTKstring);

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "-GlobalErrors.txt";
        std::ofstream file(fileName, std::ios::app);
        Tools::PrintErrors(file, config, errors);
    }

    result_vec.emplace_back(config.problemname,
                            config.ninternalref,
                            config.ndivisions,
                            config.porder,
                            config.hdivmais,
                            errors[3],
                            errors[2],
                            sqrt(errors[4] + errors[3]) / sqrt(errors[2])
    );
}
