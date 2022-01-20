//
// Created by Gustavo Batistela on 4/5/21.
//

#include "Tools.h"
#include "ToolsSPE10.h"
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <TPZGenGrid2D.h>
#include <TPZGenGrid3D.h>
#include <TPZMFSolutionTransfer.h>
#include <TPZMHMHDivErrorEstimator.h>
#include <TPZMHMixedMeshControl.h>
#include <ToolsMHM.h>
#include <iostream>
#include <memory>
#include <pzgmesh.h>

// Global variables
constexpr int nx = 220;
constexpr int ny = 60;
constexpr int n_cells = nx * ny;
TPZManVector<REAL, n_cells> perm_vec(n_cells, 1);

std::map<int, REAL> neq_error;

// Function declarations
void ReadSPE10CellPermeabilities(TPZVec<REAL>*perm_vec, int layer);
TPZGeoMesh *CreateSPE10CoarseGeoMesh();
STATE PermeabilityFunction(const TPZVec<REAL> &x);
void InsertMaterials(TPZCompMesh *cmesh);
void CreateSPE10MHMCompMesh(TPZMHMixedMeshControl &mhm, const std::vector<int64_t> &skelsToDivide);
void SolveMHMProblem(TPZMHMixedMeshControl &mhm, int adaptivity_step);
void EstimateError(TPZMHMixedMeshControl &mhm, TPZMHMHDivErrorEstimator &estimator);
void MHMAdaptivity(TPZMHMixedMeshControl *mhm, TPZCompMesh *postProcMesh, std::vector<int64_t> &skelsToRefine);

constexpr int adaptivity_steps = 5;
int main() {

    gRefDBase.InitializeAllUniformRefPatterns();

    constexpr int layer = 36;

    ReadSPE10CellPermeabilities(&perm_vec, layer);

    std::vector<int64_t> skelsToRefine;
    for (int step = 0; step < adaptivity_steps; step++) {

        TPZGeoMesh *gmesh = CreateSPE10CoarseGeoMesh();
        std::cout << "SPE10 initial grid created. NElem: " << gmesh->NElements() << "\n";

        TPZMHMixedMeshControl mhm(gmesh);
        CreateSPE10MHMCompMesh(mhm, skelsToRefine);
        std::stringstream filename;
        filename << "FineCroppedSPE10GeoMesh-Step" << step;
        Tools::PrintGeometry(gmesh, filename.str(), false, true);

        SolveMHMProblem(mhm, step);

        bool postProcWithHDiv = false;
        TPZMultiphysicsCompMesh *originalMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm.CMesh().operator->());
        if (!originalMesh) DebugStop();
        TPZMHMHDivErrorEstimator estimator(*originalMesh, &mhm, postProcWithHDiv);
        estimator.SetAdaptivityStep(step);
        EstimateError(mhm, estimator);

        auto *postprocmesh = estimator.PostProcMesh();
        if (!postprocmesh) DebugStop();

        MHMAdaptivity(&mhm, postprocmesh, skelsToRefine);
    }
    for (auto it: neq_error) {
        std::cout << "(" << it.first << ", " << it.second << ")\n";
    }

    return 0;
}

TPZGeoMesh *CreateSPE10CoarseGeoMesh() {
    std::cout << "Creating SPE10 initial grid...\n";

    const TPZManVector<REAL, 3> x0 = {0, 0, 0};
    const TPZManVector<REAL, 3> x1 = {208., 48., 0.};
    const TPZManVector<int, 3> ndiv = {13, 3, 0};

    TPZGenGrid2D gen(ndiv, x0, x1);

    gen.SetRefpatternElements(true);
    auto gmesh = new TPZGeoMesh;
    gen.Read(gmesh);

    gen.SetBC(gmesh, 5, -1);
    gen.SetBC(gmesh, 6, -2);
    gen.SetBC(gmesh, 7, -2);
    //gen.SetBC(gmesh, 4, -2);

    std::cout << "SPE10 initial grid created. NElem: " << gmesh->NElements() << "\n";

    return gmesh;
}

void ReadSPE10CellPermeabilities(TPZVec<REAL> *perm_vec, const int layer) {

    std::cout << "Reading permeability data...\n";

    std::ifstream perm_file("InputData/spe_perm.dat", std::ios::in);
    if (!perm_file) {
        std::cerr << "Unable to open input file\n";
        DebugStop();
    }

    int cell_id = 0;
    const auto n_cells = perm_vec->size();
    const auto start_line = 1 + n_cells * (layer - 1) / 6;

    int line_num = 0;
    int line_num2 = 0;
    while (perm_file) {
        line_num++;
        line_num2++;
        std::string line;
        std::getline(perm_file, line, '\n');

        if (line_num < start_line) continue;

        std::stringstream stream(line);
        for (int i = 0; i < 6; i++) {
            stream >> perm_vec->operator[](cell_id);
            cell_id++;
        }
        if (cell_id == n_cells) break;
    }
    std::cout << "Finished reading permeability data from input file!\n";
}

STATE PermeabilityFunction(const TPZVec<REAL> &x) {
    auto rounded_x = static_cast<int>(x[0]);
    auto rounded_y = static_cast<int>(x[1]);
    if (rounded_x == 220) rounded_x = 219;
    if (rounded_y == 60) rounded_y = 59;
    return perm_vec[rounded_x * 60 + rounded_y];
}

void CreateSPE10MHMCompMesh(TPZMHMixedMeshControl &mhm, const std::vector<int64_t> &skelsToDivide) {

    TPZGeoMesh *gmesh = mhm.GMesh().operator->();
    TPZManVector<int64_t, 22 * 6> coarse_indexes;
    ComputeCoarseIndices(gmesh, coarse_indexes);

    int nInternalRef = adaptivity_steps - 1;
    Tools::UniformRefinement(nInternalRef, 2, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    mhm.DefinePartitionbyCoarseIndices(coarse_indexes);

    // Indicate material indices to the MHM control structure
    mhm.fMaterialIds = {1};
    mhm.fMaterialBCIds = {-1, -2, -3};

    // Insert the material objects in the multiphysics mesh
    TPZCompMesh *cmesh = mhm.CMesh().operator->();
    InsertMaterials(cmesh);

    // General approximation order settings
    mhm.SetInternalPOrder(1);
    mhm.SetSkeletonPOrder(1);
    mhm.SetHdivmaismaisPOrder(3);

    // Refine skeleton elements
    for (auto skelid : skelsToDivide) {
        mhm.DivideSkeletonElement(skelid);
    }
    mhm.DivideBoundarySkeletonElements();

    // Creates MHM mesh
    bool substructure = true;
    mhm.BuildComputationalMesh(substructure);
}

void SolveMHMProblem(TPZMHMixedMeshControl &mhm, const int adaptivity_step) {

    TPZAutoPointer<TPZCompMesh> cmesh = mhm.CMesh();
    std::cout << "*** NEQ: ," << mhm.CMesh()->NEquations() << '\n';
    bool should_renumber = true;
    TPZLinearAnalysis an(cmesh, should_renumber);

#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> strmat(cmesh.operator->());
    strmat.SetNumThreads(0);
#else
    TPZSkylineStructMatrix<STATE> strmat(cmesh.operator->());
    strmat.SetNumThreads(8);
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
    an.LoadSolution();

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

    int resolution = 0;
    std::stringstream plotname;
    plotname << "FineCroppedSPE10-Results-Step" << adaptivity_step << ".vtk";
    an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotname.str());
    an.PostProcess(resolution, cmesh->Dimension());
}

void InsertMaterials(TPZCompMesh *cmesh) {

    auto *mix = new TPZMixedDarcyFlow(1, cmesh->Dimension());
    std::function<STATE(const TPZVec<REAL> &coord)> func = PermeabilityFunction;
    mix->SetPermeabilityFunction(func);

    TPZFNMatrix<1, REAL> val1(1, 1, 0.);
    TPZManVector<REAL, 1> val2(1, 0.);
    constexpr int dirichlet_bc = 0;
    constexpr int neumann_bc = 1;

    // Pressure at reservoir boundary
    val2[0] = 10;
    TPZBndCond *pressure_left = mix->CreateBC(mix, -1, dirichlet_bc, val1, val2);
    val2[0] = 0;
    TPZBndCond *pressure_general = mix->CreateBC(mix, -2, neumann_bc, val1, val2);

    cmesh->InsertMaterialObject(mix);
    cmesh->InsertMaterialObject(pressure_left);
    cmesh->InsertMaterialObject(pressure_general);
}

void EstimateError(TPZMHMixedMeshControl &mhm, TPZMHMHDivErrorEstimator &estimator) {

    std::cout << "\nError Estimation processing for MHM-Hdiv problem " << std::endl;

    // Error estimation
    TPZMultiphysicsCompMesh *originalMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm.CMesh().operator->());
    if (!originalMesh) DebugStop();

    estimator.PotentialReconstruction();

    std::string command = "mkdir SPE10";
    system(command.c_str());

    TPZManVector<REAL, 6> errors;
    TPZManVector<REAL> elementerrors;
    std::stringstream plotname;
    plotname << "FineCroppedSPE10-Errors-Step" << estimator.AdaptivityStep() << ".vtk";
    auto plotname_str = plotname.str();
    estimator.ComputeErrors(errors, elementerrors, plotname_str);
    Tools::PrintErrors(std::cout, errors);

    neq_error.insert({mhm.CMesh()->NEquations(), errors[3]});

    std::cout << "Finished computing errors!\n";
}

void MHMAdaptivity(TPZMHMixedMeshControl *mhm, TPZCompMesh *postProcMesh, std::vector<int64_t> &skelsToRefine) {

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


    // Refines elements which error are bigger than 30% of the maximum error
    REAL threshold = 0.5 * maxError;
    const auto geoToMHM = mhm->GetGeoToMHMDomain();
    const auto interfaces = mhm->GetInterfaces();

    std::set<int64_t> interfacesToRefine;
    for (int64_t iel = 0; iel < nelem; iel++) {
        auto * submesh = dynamic_cast<TPZSubCompMesh*>(postProcMesh->ElementVec()[iel]);
        if (!submesh) continue;

        STATE submeshError = sol(iel, fluxErrorEstimateCol);
        if (submeshError > threshold) {
        //if (true) {
            TPZGeoEl * gel = submesh->Element(0)->Reference();
            const auto submesh_id = geoToMHM[gel->Index()];
            //std::cout << "Refining submesh " << submesh_id << " which error is " << submeshError << ".\n";
            for (auto interface : interfaces) {
                if (interface.second.first == submesh_id || interface.second.second == submesh_id) {
                    interfacesToRefine.insert(interface.first);
                }
            }
        }
    }

    for (auto it : interfacesToRefine) {
        skelsToRefine.push_back(it);
    }
}