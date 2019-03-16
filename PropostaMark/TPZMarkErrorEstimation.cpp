//
// Created by gus on 14/03/19.
//

#include "TPZMarkErrorEstimation.h"


TPZMarkErrorEstimation::TPZMarkErrorEstimation(ProblemConfig problem) {
    fProblem = problem;
}

void TPZMarkErrorEstimation::CreateMultiphysicsMesh() {

    TPZCompMesh *cmesh = new TPZCompMesh(fProblem.gmesh);
    TPZMaterial *mat = nullptr;

    for (auto matid : fProblem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        // mix->SetForcingFunction(problem.exact.ForcingFunction());
        TPZAutoPointer<TPZFunction<STATE>> force1 = new TPZDummyFunction<STATE>(Forcing, 0);
        mix->SetForcingFunction(force1);
        //mix->SetForcingFunction(problem.forcingCte.);
        mix->SetInternalFlux(1);
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    for (auto matid : fProblem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        if (matid == -2) {
            bctype = 0;
            val2.Zero();
        }
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        //   bc->TPZMaterial::SetForcingFunction(problem.exact.Exact());
        cmesh->InsertMaterialObject(bc);
    }

    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    cmesh->AutoBuild();

    CreateSolutionPressureMesh();
    CreateSolutionFluxMesh();

    TPZManVector<TPZCompMesh *,2> meshvector;
    meshvector[0] = fMeshVector[0];
    meshvector[1] = fMeshVector[1];

    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh);
    cmesh->LoadReferences();
    bool keepmatrix = false;
    TPZCompMeshTools::CreatedCondensedElements(cmesh, true, keepmatrix);

    cmesh->InitializeBlock();
    fMeshVector[2] = cmesh;
}

void TPZMarkErrorEstimation::CreateSolutionFluxMesh() {

    int dim = fProblem.gmesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(fProblem.gmesh);
    TPZMaterial *mat = NULL;
    fProblem.gmesh->ResetReference();
    for (auto matid : fProblem.materialids) {
        TPZVecL2 *mix = new TPZVecL2(matid);
        mix->SetDimension(dim);
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    for (auto matid : fProblem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 1.);
        TPZBndCond *bc = mat->CreateBC(mat, matid, 0, val1, val2);
        cmesh->InsertMaterialObject(bc);
    }
    cmesh->SetDefaultOrder(fProblem.porder);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmesh->AutoBuild();
    if (fProblem.hdivmais) {
        int64_t nel = cmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if (gel->Dimension() == dim) {
                int side = gel->NSides() - 1;
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
                intel->SetSideOrder(side, fProblem.porder + fProblem.hdivmais);//seta ordem +hdivmais
                intel->SetPreferredOrder(fProblem.porder + fProblem.hdivmais);
            }
        }
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if (gel->Dimension() == dim - 1) {
                int side = gel->NSides() - 1;
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
                intel->SetSideOrder(side, fProblem.porder + fProblem.hdivmais);
                intel->SetPreferredOrder(fProblem.porder + fProblem.hdivmais);
            }
        }
    }

    cmesh->InitializeBlock();
    fMeshVector[0] = cmesh;
}

void TPZMarkErrorEstimation::CreateSolutionPressureMesh() {

    TPZCompMesh *cmesh = new TPZCompMesh(fProblem.gmesh);
    TPZMaterial *mat = 0;
    for (auto matid : fProblem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    cmesh->SetDefaultOrder(fProblem.porder + fProblem.hdivmais);//ordem + hdivmais
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->AutoBuild();
    int64_t n_connects = cmesh->NConnects();
    for (int64_t i = 0; i < n_connects; ++i) {
        cmesh->ConnectVec()[i].SetLagrangeMultiplier(1);
    }

    fMeshVector[1] = cmesh;
}

void TPZMarkErrorEstimation::CreatePostProcessMeshes() {

}

void TPZMarkErrorEstimation::SolveMixedProblem() {

    TPZAnalysis an(fMeshVector[2], false);
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(fMeshVector[2]);
    strmat.SetNumThreads(0);
//        strmat.SetDecomposeType(ELDLt);
    an.SetStructuralMatrix(strmat);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmesh_HDiv);
    strmat.SetNumThreads(0);
    //        TPZSkylineStructMatrix strmat3(cmesh_HDiv);
    //        strmat3.SetNumThreads(8);
#endif

    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();//resolve o problema misto ate aqui

    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    vecnames.Push("Flux");
    an.DefineGraphMesh(2, scalnames, vecnames, "Original.vtk");

    // Post processing
    an.PostProcess(0, 2);
}

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp) {
    disp[0] = 1.;
    return;
}
