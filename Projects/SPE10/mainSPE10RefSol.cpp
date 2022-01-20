//
// Created by Gustavo Batistela on 4/5/21.
//

#include "Tools.h"
#include "ToolsSPE10.h"
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <TPZNullMaterial.h>
#include <TPZNullMaterialCS.h>
#include <TPZGenGrid3D.h>
#include <TPZMFSolutionTransfer.h>
#include <TPZMHMHDivErrorEstimator.h>
#include <TPZMHMixedMeshControl.h>
#include <iostream>
#include <memory>
#include <pzgmesh.h>

constexpr int porder = 1;
constexpr int korder = 2;
constexpr int nInternalRef = 5;

TPZMultiphysicsCompMesh RunFineScaleProblem();

int main() {

    gRefDBase.InitializeAllUniformRefPatterns();

    SPE10::ReadSPE10CellPermeabilities();

    auto fine_cmesh = RunFineScaleProblem();
    // Create geomesh
    TPZGeoMesh *gmesh = SPE10::CreateSPE10CoarseGeoMesh();
    std::cout << "SPE10 initial grid created. NElem: " << gmesh->NElements() << "\n";

    // Create compmesh
    TPZMHMixedMeshControl mhm(gmesh);

    std::vector<int64_t> skelsToRefine;
    SPE10::CreateSPE10MHMCompMesh(mhm, skelsToRefine, nInternalRef);
    constexpr int step = 0;
    std::stringstream filename;
    filename << "CroppedSPE10GeoMesh-Step" << step;
    Tools::PrintGeometry(gmesh, filename.str(), false, true);

    SPE10::SolveMHMProblem(mhm, step);

    // Estimate error
    bool postProcWithHDiv = false;
    TPZMultiphysicsCompMesh *originalMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm.CMesh().operator->());
    if (!originalMesh) DebugStop();

    TPZMHMHDivErrorEstimator estimator(*originalMesh, &mhm, postProcWithHDiv);
    estimator.SetAdaptivityStep(step);
    SPE10::EstimateError(mhm, estimator, &fine_cmesh);

    return 0;
}

TPZMultiphysicsCompMesh RunFineScaleProblem() {
    // Create geomesh
    TPZGeoMesh *gmesh = SPE10::CreateSPE10CoarseGeoMesh();
    std::cout << "SPE10 initial grid created. NElem: " << gmesh->NElements() << "\n";

    Tools::UniformRefinement(nInternalRef, 2, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);
    const int dim = gmesh->Dimension();

    TPZMultiphysicsCompMesh mixed_cmesh(gmesh);

    SPE10::InsertMaterials(&mixed_cmesh);
    mixed_cmesh.ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

    TPZManVector<int> active(2, 1);
    TPZManVector<TPZCompMesh*> meshvector(2, 0);

    auto* cmesh_flux = new TPZCompMesh(gmesh);
    gmesh->ResetReference();
    auto *null_mat = new TPZNullMaterial<STATE>(1);
    null_mat->SetDimension(dim);
    cmesh_flux->InsertMaterialObject(null_mat);

    for (auto matid : {-1, -2}) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.);
        TPZManVector<REAL, 1> val2(1, 1.);
        const int bctype = 0;
        auto * bc = null_mat->CreateBC(null_mat, matid, bctype, val1, val2);
        cmesh_flux->InsertMaterialObject(bc);
    }
    cmesh_flux->SetDefaultOrder(porder);
    cmesh_flux->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmesh_flux->AutoBuild();

    cmesh_flux->InitializeBlock();

    auto *cmesh_pressure = new TPZCompMesh(gmesh);
    auto *null_mat_p = new TPZNullMaterial<STATE>(1);
    null_mat_p->SetDimension(cmesh_pressure->Dimension());
    cmesh_pressure->InsertMaterialObject(null_mat_p);

    cmesh_pressure->SetDefaultOrder(porder + korder);
    // cmesh_pressure->SetDefaultOrder(problem.porder);
    cmesh_pressure->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh_pressure->ApproxSpace().CreateDisconnectedElements(true);
    cmesh_pressure->AutoBuild();
    int64_t n_connects = cmesh_pressure->NConnects();
    for (int64_t i = 0; i < n_connects; ++i) {
        cmesh_pressure->ConnectVec()[i].SetLagrangeMultiplier(1);
    }

    meshvector[0] = cmesh_flux;
    meshvector[1] = cmesh_pressure;

    TPZCompMeshTools::AdjustFluxPolynomialOrders(meshvector[0], korder);
    TPZCompMeshTools::SetPressureOrders(meshvector[0], meshvector[1]);

    mixed_cmesh.BuildMultiphysicsSpace(active, meshvector);
    mixed_cmesh.LoadReferences();
    bool keepmatrix = false;
    bool keeponelagrangian = true;
    TPZCompMeshTools::CreatedCondensedElements(&mixed_cmesh, keeponelagrangian, keepmatrix);
    mixed_cmesh.InitializeBlock();

    TPZHybridizeHDiv hybridizer;
    TPZMultiphysicsCompMesh* hybridMesh = hybridizer.Hybridize(&mixed_cmesh);
    hybridMesh->CleanUpUnconnectedNodes();
    hybridMesh->AdjustBoundaryElements();

    TPZLinearAnalysis an(hybridMesh);

#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> strmat(hybridMesh);
    strmat.SetNumThreads(0);
#else
    TPZSkylineStructMatrix<STATE> strmat(Hybridmesh);
    strmat.SetNumThreads(0);
#endif

    std::set<int> matIds = {1, -1, -2, hybridizer.fInterfaceMatid.first, hybridizer.fInterfaceMatid.second};
    strmat.SetMaterialIds(matIds);

    an.SetStructuralMatrix(strmat);

    auto *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = nullptr;
    an.Assemble();
    an.Solve();

    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("ExactPressure");
    scalnames.Push("Pressure");
    scalnames.Push("Permeability");
    vecnames.Push("ExactFlux");
    vecnames.Push("Flux");

    std::stringstream sout;
    sout << "SPE10FINE-Results.vtk";
    an.DefineGraphMesh(2, scalnames, vecnames, sout.str());
    constexpr int resolution = 0;
    an.PostProcess(resolution, hybridMesh->Dimension());
    return mixed_cmesh;
}
