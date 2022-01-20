//
// Created by Gustavo Batistela on 10/26/21.
//

#include "Tools.h"
#include "ToolsMHM.h"
#include "ToolsSPE10.h"
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <TPZGenGrid2D.h>
#include <TPZMFSolutionTransfer.h>
#include <TPZSSpStructMatrix.h>
#include <pzcmesh.h>
#include <pzgeoquad.h>
#include <pzstepsolver.h>
#include <tpzgeoelrefpattern.h>

[[maybe_unused]] TPZGeoMesh *SPE10::CreateFineGridGeoMesh() {

    TPZManVector<int, 4> bcIDs = {-2, -1, -2, -2};
    TPZManVector<int, 2> n_elements = {220, 60};
    TPZManVector<REAL, 3> x0 = {0, 0, 0};
    TPZManVector<REAL, 3> x1 = {220, 60, 0};

    TPZGenGrid2D gen(n_elements, x0, x1, 1, 0);
    gen.SetRefpatternElements(true);

    auto *gmesh = new TPZGeoMesh;
    gen.Read(gmesh);

    gen.SetBC(gmesh, 4, bcIDs[0]);
    gen.SetBC(gmesh, 5, bcIDs[1]);
    gen.SetBC(gmesh, 6, bcIDs[2]);
    gen.SetBC(gmesh, 7, bcIDs[3]);

    gmesh->SetDimension(2);

    return gmesh;
}

[[maybe_unused]] TPZGeoMesh *SPE10::CreateMHMGeoMesh() {

    auto *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);

    auto *gmesh2x1 = CreateRefinementGeoMesh(2, 1);
    TPZRefPattern ref_pat2x1(*gmesh2x1);
    TPZAutoPointer<TPZRefPattern> ref2x1(&ref_pat2x1);

    auto *gmesh1x2 = CreateRefinementGeoMesh(1, 2);
    TPZRefPattern ref_pat1x2(*gmesh1x2);
    TPZAutoPointer<TPZRefPattern> ref1x2(&ref_pat1x2);

    auto *gmesh1x1 = CreateRefinementGeoMesh(1, 1);
    TPZRefPattern ref_pat1x1(*gmesh1x1);
    TPZAutoPointer<TPZRefPattern> ref1x1(&ref_pat1x1);

    TPZManVector<REAL, 3> coord(3, 0.);
    for (int y = 0; y <= 8; y++) {
        for (int x = 0; x <= 28; x++) {
            coord = {8.0 * x, 8.0 * y, 0};
            if (x == 28) coord[0] = 220.;
            if (y == 8) coord[1] = 60.;

            // Create new node
            const auto newID = gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[newID].Initialize(coord, *gmesh);
        }
    }

    constexpr int porousMediaMatId = 1;
    constexpr int bcMatId1 = -1;
    constexpr int bcMatId2 = -2;

    // Inserts quad elements
    TPZManVector<int64_t, 4> nodesIdVec(4);
    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 28; x++) {
            nodesIdVec[0] = 29 * y + x;
            nodesIdVec[1] = 29 * y + x + 1;
            nodesIdVec[2] = 29 * (y + 1) + x + 1;
            nodesIdVec[3] = 29 * (y + 1) + x + 0;
            auto *gel = new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec, porousMediaMatId, *gmesh);
            if (x == 27 && y != 7) gel->SetRefPattern(ref1x2);
            if (x != 27 && y == 7) gel->SetRefPattern(ref2x1);
            if (x == 27 && y == 7) gel->SetRefPattern(ref1x1);
        }
    }

    auto *gmesh1x1line = CreateLineRefinementGeoMesh(1);
    TPZRefPattern ref_pat1x1line(*gmesh1x1line);
    TPZAutoPointer<TPZRefPattern> ref1x1line(&ref_pat1x1line);

    nodesIdVec.resize(2);
    for (int x = 0; x < 28; x++) {
        nodesIdVec[0] = x;
        nodesIdVec[1] = x + 1;
        auto *bc1 = new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, bcMatId1, *gmesh);
        nodesIdVec[0] = x + 8 * 29;
        nodesIdVec[1] = x + 8 * 29 + 1;
        auto *bc2 = new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, bcMatId1, *gmesh);
        if (x == 27) {
            bc1->SetRefPattern(ref1x1line);
            bc2->SetRefPattern(ref1x1line);
        }
    }

    for (int y = 0; y < 8; y++) {
        nodesIdVec[0] = y * 29;
        nodesIdVec[1] = (y + 1) * 29;
        auto *bc1 = new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, bcMatId1, *gmesh);
        nodesIdVec[0] = y * 29 + 28;
        nodesIdVec[1] = (y + 1) * 29 + 28;
        auto *bc2 = new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, bcMatId2, *gmesh);
        if (y == 7) {
            bc1->SetRefPattern(ref1x1line);
            bc2->SetRefPattern(ref1x1line);
        }
    }

    gmesh->BuildConnectivity();

    TPZManVector<TPZGeoEl *, 4> sons;
    for (int div = 0; div < 3; div++) {
        auto nelem = gmesh->NElements();
        for (int64_t i = 0; i < nelem; i++) {
            auto *gel = gmesh->ElementVec()[i];
            const int has_sub = gel->HasSubElement();
            if (has_sub == 0) {
                gel->Divide(sons);
            }
        }
    }
    return gmesh;
}

[[maybe_unused]] TPZGeoMesh *SPE10::CreateRefinementGeoMesh(const int nx, const int ny) {
    auto *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    TPZManVector<REAL, 3> coord(3, 0.);

    for (int y = 0; y <= ny; y++) {
        for (int x = 0; x <= nx; x++) {
            coord = {static_cast<REAL>(x), static_cast<REAL>(y), 0.};
            // Create new node
            const auto newID = gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[newID].Initialize(coord, *gmesh);
        }
    }

    constexpr int matId = 1;

    // Inserts quad elements
    TPZManVector<int64_t, 4> nodesIdVec(4, 0);
    nodesIdVec[1] = nx;
    nodesIdVec[2] = (ny + 1) * (nx + 1) - 1;
    nodesIdVec[3] = ny * (nx + 1);
    auto *father_gel = new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec, matId, *gmesh);

    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            nodesIdVec[0] = (nx + 1) * y + x;
            nodesIdVec[1] = (nx + 1) * y + x + 1;
            nodesIdVec[2] = (nx + 1) * (y + 1) + x + 1;
            nodesIdVec[3] = (nx + 1) * (y + 1) + x + 0;
            auto *gel = new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec, matId, *gmesh);
            gel->SetFather(father_gel);
        }
    }

    gmesh->BuildConnectivity();
    return gmesh;
}

[[maybe_unused]] TPZGeoMesh *SPE10::CreateLineRefinementGeoMesh(const int nx) {
    auto *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(1);

    TPZManVector<REAL, 3> coord(3, 0.);
    for (int x = 0; x <= nx; x++) {
        coord = {static_cast<REAL>(x), 0, 0.};
        // Create new node
        const auto newID = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[newID].Initialize(coord, *gmesh);
    }

    constexpr int matId = 1;

    // Inserts line elements
    TPZManVector<int64_t, 4> nodesIdVec(2, 0);
    nodesIdVec[1] = nx;
    auto *father_gel = new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, matId, *gmesh);

    for (int x = 0; x < nx; x++) {
        nodesIdVec[0] = x;
        nodesIdVec[1] = x + 1;
        auto *gel = new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, matId, *gmesh);
        gel->SetFather(father_gel);
    }

    gmesh->BuildConnectivity();
    return gmesh;
}

STATE SPE10::PermeabilityFunction(const TPZVec<REAL> &x) {
    auto rounded_x = static_cast<int>(x[0]);
    auto rounded_y = static_cast<int>(x[1]);
    if (rounded_x == 220) rounded_x = 219;
    if (rounded_y == 60) rounded_y = 59;
    return perm_vec->operator[](rounded_x * 60 + rounded_y);
}

void SPE10::ReadSPE10CellPermeabilities() {

    std::cout << "Reading permeability data...\n";

    std::ifstream perm_file("InputData/spe_perm.dat", std::ios::in);
    if (!perm_file) {
        std::cerr << "Unable to open input file\n";
        DebugStop();
    }

    perm_vec = new TPZManVector<REAL, n_cells>(n_cells, 0);

    int cell_id = 0;
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

void SPE10::InsertMaterials(TPZCompMesh *cmesh) {

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

void SPE10::EstimateError(TPZMHMixedMeshControl &mhm, TPZMHMHDivErrorEstimator &estimator) {

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
    plotname << "CroppedSPE10-Errors-Step" << estimator.AdaptivityStep() << ".vtk";
    auto plotname_str = plotname.str();
    estimator.ComputeErrors(errors, elementerrors, plotname_str);
    std::cout << "Finished computing errors!\n";
}

void SPE10::EstimateError(TPZMHMixedMeshControl &mhm, TPZMHMHDivErrorEstimator &estimator, TPZMultiphysicsCompMesh * ref_sol) {

    std::cout << "\nError Estimation processing for MHM-Hdiv problem " << std::endl;

    // Error estimation
    TPZMultiphysicsCompMesh *originalMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm.CMesh().operator->());
    if (!originalMesh) DebugStop();

    estimator.SetReferenceSolution(true);
    estimator.SetReferenceSolutionMeshes(ref_sol->MeshVector()[1], ref_sol->MeshVector()[1]);

    estimator.PotentialReconstruction();


    std::string command = "mkdir SPE10";
    system(command.c_str());

    TPZManVector<REAL, 6> errors;
    TPZManVector<REAL> elementerrors;
    std::stringstream plotname;
    plotname << "CroppedSPE10-Errors-Step" << estimator.AdaptivityStep() << ".vtk";
    auto plotname_str = plotname.str();
    estimator.ComputeErrors(errors, elementerrors, plotname_str);
    std::cout << "Finished computing errors!\n";
}

void SPE10::SolveMHMProblem(TPZMHMixedMeshControl &mhm, const int adaptivity_step) {

    TPZAutoPointer<TPZCompMesh> cmesh = mhm.CMesh();

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
    plotname << "CroppedSPE10-Results-Step" << adaptivity_step << ".vtk";
    an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotname.str());
    an.PostProcess(resolution, cmesh->Dimension());
}

TPZGeoMesh *SPE10::CreateSPE10CoarseGeoMesh() {
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
    // gen.SetBC(gmesh, 4, -2);

    std::cout << "SPE10 initial grid created. NElem: " << gmesh->NElements() << "\n";

    return gmesh;
}

void SPE10::CreateSPE10MHMCompMesh(TPZMHMixedMeshControl &mhm, const std::vector<int64_t> &skelsToDivide, const int nInternalRef) {

    TPZGeoMesh *gmesh = mhm.GMesh().operator->();
    TPZManVector<int64_t, 22 * 6> coarse_indexes;
    ComputeCoarseIndices(gmesh, coarse_indexes);

    Tools::UniformRefinement(nInternalRef, 2, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    mhm.DefinePartitionbyCoarseIndices(coarse_indexes);

    // Indicate material indices to the MHM control structure
    mhm.fMaterialIds = {1};
    mhm.fMaterialBCIds = {-1, -2, -3};

    // Insert the material objects in the multiphysics mesh
    TPZCompMesh *cmesh = mhm.CMesh().operator->();
    SPE10::InsertMaterials(cmesh);

    // General approximation order settings
    mhm.SetInternalPOrder(1);
    mhm.SetSkeletonPOrder(1);
    mhm.SetHdivmaismaisPOrder(2);

    // Refine skeleton elements
    for (auto skelid : skelsToDivide) {
        mhm.DivideSkeletonElement(skelid);
    }
    mhm.DivideBoundarySkeletonElements();

    // Creates MHM mesh
    bool substructure = true;
    mhm.BuildComputationalMesh(substructure);
}
