//
// Created by victor on 16/02/23.
//

#include "TPZHybridH1CreateRecMeshes.h"
#include "TPZAnalysis.h"
#include "TPZCompMeshTools.h"
#include "TPZCreateMultiphysicsSpace.h"
#include "TPZElementMatrixT.h"
#include "TPZHybridH1ErrorEstimator.h"
#include "TPZHybridH1HdivFluxRecMaterial.h"
#include "TPZInterfaceEl.h"
#include "TPZMatLaplacianHybrid.h"
#include "TPZMixedHdivErrorEstimate.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"
#include "pzintel.h"
#include "pzsubcmesh.h"
#include <TPZVTKGeoMesh.h>

// a method for generating the HDiv mesh
TPZCompMesh *TPZHybridH1CreateRecMeshes::CreateFluxReconstructionHDivMesh()
{
    TPZCompMesh *HDivAtomicMesh = fHybridH1EE->fOriginal->MeshVector()[0]->Clone(); // HDIV-BOUND elements: clone might be unnecessary

#define ERRORESTIMATION_DEBUG33
#ifdef ERRORESTIMATION_DEBUG33
    std::string command = "mkdir -p " + fHybridH1EE->fDebugDirName;
    std::string dirPath = fHybridH1EE->fDebugDirName + "/";
    system(command.c_str());
    {
        std::ofstream outCon(dirPath + "OriginalFluxConnects.txt");
        if(HDivAtomicMesh->Dimension() == 2) {
            //TPZCompMeshTools::PrintConnectInfoByGeoElement(cmeshHdiv, outCon, {}, false, true);
        }
        std::ofstream outOriginalP(dirPath + "OriginalFlux.txt");
        HDivAtomicMesh->Print(outOriginalP);
        std::ofstream foutOriginalP(dirPath + "fOriginalFlux.txt");
        fHybridH1EE->fOriginal->MeshVector()[0]->Print(foutOriginalP);
        std::ofstream outGOriginalVTK(dirPath + "gFlux.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(HDivAtomicMesh, outGOriginalVTK);
    }
#endif

    // Verify neighbouring information of the HDiv-bound mesh (works only for 2D meshes)
    VerifyBoundaryFluxConsistency(HDivAtomicMesh);
    std::ofstream myoutput("chekckingfOriginalMesh.txt"); fHybridH1EE->fOriginal->Print(myoutput); myoutput.flush();
    int meshdim = HDivAtomicMesh->Dimension();
    for (auto mat : fHybridH1EE->fOriginal->MaterialVec()) {
        if (!dynamic_cast<TPZBndCondT<STATE> *>(mat.second)) {
            if (mat.second->Dimension() == meshdim) {
                auto mymat = new TPZNullMaterial(mat.first, mat.second->Dimension());
                HDivAtomicMesh->InsertMaterialObject(mymat);
            }
        }
    }

    HDivAtomicMesh->SetDefaultOrder(fHybridH1EE->fProblemConfig.k);

    HDivAtomicMesh->AutoBuild();
    HDivAtomicMesh->InitializeBlock();

    TPZCompMeshTools::AdjustFluxPolynomialOrders(HDivAtomicMesh, fHybridH1EE->fProblemConfig.n); //Increases internal flux order by "hdivmais"

#ifdef ERRORESTIMATION_DEBUG33
    {
        std::ofstream outCon(dirPath + "HdivFluxConnects.txt");
        //TPZCompMeshTools::PrintConnectInfoByGeoElement(cmeshHdiv, outCon, {}, false, true);
        std::ofstream outOriginalP(dirPath + "HdivFlux.txt");
        HDivAtomicMesh->Print(outOriginalP);
        std::ofstream outGOriginalVTK(dirPath + "gHdivFlux.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(HDivAtomicMesh, outGOriginalVTK);
    }
#endif

    return HDivAtomicMesh;
}

TPZCompMesh *TPZHybridH1CreateRecMeshes::CreateFluxReconstructionL2Mesh(){

#ifdef ERRORESTIMATION_DEBUG
    {
        std::ofstream outOriginalP(dirPath + "HdivCondFlux.txt");
        cmeshHdiv->Print(outOriginalP);
    }
#endif

    TPZCompMesh *L2AtomicMesh = new TPZCompMesh(fHybridH1EE->fProblemConfig.gmesh);

    int dimMesh = fHybridH1EE->fProblemConfig.gmesh->Dimension();

    int potential_order = fHybridH1EE->fProblemConfig.k+fHybridH1EE->fProblemConfig.n;
    L2AtomicMesh->SetDefaultOrder(potential_order);
    L2AtomicMesh->SetDimModel(dimMesh);

    L2AtomicMesh->SetAllCreateFunctionsContinuous(); //H1 functions
    L2AtomicMesh->ApproxSpace().CreateDisconnectedElements(true);

    for(auto matid:fHybridH1EE->fProblemConfig.materialids){
        TPZNullMaterial<> *material = new TPZNullMaterial<>(matid); material->SetDimension(dimMesh);
        L2AtomicMesh->InsertMaterialObject(material);
    }
    L2AtomicMesh->AutoBuild();
    L2AtomicMesh->ExpandSolution();

    return L2AtomicMesh;
}

TPZCompMesh *TPZHybridH1CreateRecMeshes::CreateFluxReconstructionConstantMesh(){
    TPZCompMesh *constant = new TPZCompMesh(fHybridH1EE->fProblemConfig.gmesh);
    {
        for (auto matid:fHybridH1EE->fProblemConfig.materialids) {
            TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(matid);
            nullmat->SetDimension(fHybridH1EE->fProblemConfig.gmesh->Dimension());
            nullmat->SetNStateVariables(1);
            constant->InsertMaterialObject(nullmat);
        }
        constant->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
        constant->SetDefaultOrder(0);
        constant->AutoBuild();
    }
    return constant;
}

TPZMultiphysicsCompMesh *TPZHybridH1CreateRecMeshes::CreateFluxReconstructionMesh()
{
    TPZCompMesh *HDivAtomicMesh = CreateFluxReconstructionHDivMesh();
    TPZCompMesh *L2AtomicMesh = CreateFluxReconstructionL2Mesh();
    TPZCompMesh *gspace = CreateFluxReconstructionConstantMesh();

    auto HdivRecMesh = new TPZMultiphysicsCompMesh(fHybridH1EE->fOriginal->Reference());

    TPZManVector<TPZCompMesh *> mesh_vectors(5, 0);
    TPZManVector<int> active(5, 0);

    mesh_vectors[0] = HDivAtomicMesh;
    mesh_vectors[1] = L2AtomicMesh;
    mesh_vectors[2] = gspace;
    mesh_vectors[3] = fHybridH1EE->fOriginal->MeshVector()[3]->Clone(); // avg-space
    mesh_vectors[4] = fHybridH1EE->fOriginal->MeshVector()[1];

    active[0] = 1;
    active[1] = 1;
    active[2] = 1;



    // Not good
    // Insert materials into Multiphysics mesh
    // The Wrap and interface material does not need to be created
    for(auto mat: fHybridH1EE->fOriginal->MaterialVec()) {
        TPZMatLaplacianHybrid *matlaplacian = dynamic_cast<TPZMatLaplacianHybrid *>(mat.second);
        // Casting lagrange material to error estimate material for bc and lagrange coefficient objects
        if (matlaplacian) {
            TPZHybridH1HdivFluxRecMaterial *EEMat = new TPZHybridH1HdivFluxRecMaterial(*matlaplacian);
            HdivRecMesh->InsertMaterialObject(EEMat);

            for(auto trybcmat: fHybridH1EE->fOriginal->MaterialVec()) {
                auto *bc = dynamic_cast<TPZBndCondT<STATE> *>(trybcmat.second);
                if (bc) {
                    // add bc mtf;
                    if(bc->Material()->Id() != EEMat->Id()){
                        DebugStop(); // this should not be evoked unless multiple volumetric ids are invoked.
                    }
                    bc->SetMaterial(EEMat);
                }
            }
        }
        if(mat.first == fHybridH1EE->fHybridizer.fFluxMatId){
            // add lagrange material
            HdivRecMesh->InsertMaterialObject(mat.second);
        }
    }

    HdivRecMesh->SetAllCreateFunctionsMultiphysicElem();

    HdivRecMesh->BuildMultiphysicsSpace(active,mesh_vectors);

    bool keeponelagrangian = true, keepmatrix = false;
    TPZCompMeshTools::CreatedCondensedElements(HdivRecMesh, keeponelagrangian, keepmatrix);

#define ERRORESTIMATION_DEBUG37
#ifdef ERRORESTIMATION_DEBUG37
    {
        std::ofstream outMultF("HdivMultMesh.txt");
        HdivRecMesh->Print(outMultF);

        std::ofstream outMultOriginal("HdivMultOriginal.txt");
        fHybridH1EE->fOriginal->Print(outMultOriginal);

        std::ofstream outL2("L2Mesh.txt");
        HdivRecMesh->MeshVector()[1]->Print(outL2);

        std::ofstream outH1("H1Mesh.txt");
        HdivRecMesh->MeshVector()[4]->Print(outH1);
    }
#endif

    ComputeElementStiffnesses(*HdivRecMesh);

    HdivRecMesh->LoadSolution(HdivRecMesh->Solution());

    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mesh_vectors, HdivRecMesh);

#ifdef ERRORESTIMATION_DEBUG
    {
        std::ofstream outFCon(dirPath + "HdivFluxAfterLoadSolConnects.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(cmeshHdiv, outFCon, {}, false, true);
        std::ofstream outF(dirPath + "HdivFluxAfterLoadSol.txt");
        cmeshHdiv->Print(outF);
        std::ofstream outMultF(dirPath + "HdivMultFluxAfterLoadSol.txt");
        HdivRecMesh.Print(outMultF);
        std::ofstream outvtk(dirPath + "fluxmesh.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmeshHdiv, outvtk);
    }
#endif
#ifdef ERRORESTIMATION_DEBUG
    TPZLinearAnalysis an(&HdivRecMesh, false);
    TPZStack<std::string> scalnames, vecnames;
    int dim = HdivRecMesh.Reference()->Dimension();

    vecnames.Push("FluxReconstructed");
    vecnames.Push("FluxExact");
    std::stringstream out;
    out << dirPath << "fluxmesh.vtk";

    std::set<int> matids ={1};
    an.DefineGraphMesh(dim,matids, scalnames, vecnames, out.str());
    an.PostProcess(2, dim);
#endif
    std::ofstream fluxrectxt("myfluxreconstructedafternewclass.txt");
    HdivRecMesh->MeshVector()[0]->Print(fluxrectxt);
    std::ofstream outH1after("H1MeshafterLoadSol.txt");
    HdivRecMesh->MeshVector()[4]->Print(outH1after);

    return HdivRecMesh;
}

void TPZHybridH1CreateRecMeshes::VerifyBoundaryFluxConsistency(TPZCompMesh* fluxmesh){
    TPZGeoMesh *gmesh = fluxmesh->Reference();
    gmesh->ResetReference();
    fluxmesh->LoadReferences();
    int meshdim = gmesh->Dimension();
    int nel = fluxmesh->NElements(), nEffectiveFluxes = 0;
    for(int iel = 0; iel < nel ;iel++) {
        TPZCompEl *cel = fluxmesh->Element(iel);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        int matID = gel->MaterialId();
        if(meshdim != fluxmesh->Dimension()) DebugStop();
        bool IsBC = (fHybridH1EE->fProblemConfig.bcmaterialids.find(matID) != fHybridH1EE->fProblemConfig.bcmaterialids.end());
        if(!IsBC && matID != fHybridH1EE->fHybridizer.fFluxMatId) {
            DebugStop();
        }
        nEffectiveFluxes++;
        int sideminus = gel->NSides(meshdim-2);
        // Check if every corner
        for(int iside = 0 ; iside < sideminus ; iside++){
            TPZGeoElSide minusSide(gel,iside);
            if (minusSide.Dimension() != 0) DebugStop();
            TPZStack<TPZCompElSide> minusNeigh;
            minusSide.ConnectedCompElementList(minusNeigh, 1, 0);
            if(minusNeigh.size() == 0) DebugStop();
        }
    }
    if(nEffectiveFluxes == 0) DebugStop();
}

/// computing the element stifnesses will "automatically" compute the condensed form of the matrices
void TPZHybridH1CreateRecMeshes::ComputeElementStiffnesses() {
    std::cout << "Solving local Dirichlet problem " << std::endl;
#ifdef ERRORESTIMATION_DEBUG2

    {
        std::ofstream out("MeshToComputeStiff.txt");
        fPostProcMesh.Print(out);
    }
#endif
    for (auto cel: fHybridH1EE->fPostProcMesh.ElementVec()) {
        if (!cel) continue;
        TPZCondensedCompEl *condense = dynamic_cast<TPZCondensedCompEl *>(cel);
        if (condense) {
            // for Mark proposal ek correspond to local dirichlet problem
            condense->Assemble();
        }
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (subcmesh) {
            subcmesh->Assemble();
        }
#ifdef ERRORESTIMATION_DEBUG
        if(subcmesh && condense)
        {
            DebugStop();
        }
#endif
    }
}

/// computing the element stifnesses will "automatically" compute the condensed form of the matrices
void TPZHybridH1CreateRecMeshes::ComputeElementStiffnesses(TPZCompMesh &cmesh) {
    std::cout << "Computing flux stiffness matrix " << std::endl;

    for (auto cel:cmesh.ElementVec()) {
        if (!cel) continue;
        TPZCondensedCompEl *condense = dynamic_cast<TPZCondensedCompEl *>(cel);
        if (condense) {
            // for Mark proposal ek correspond to local dirichlet problem
            condense->Assemble();
        }
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (subcmesh) {
            subcmesh->Assemble();
        }
#ifdef ERRORESTIMATION_DEBUG
        if(subcmesh && condense)
        {
            DebugStop();
        }
#endif
    }
}

void TPZHybridH1CreateRecMeshes::PostProcess(TPZMultiphysicsCompMesh *postProcMesh){
    TPZLinearAnalysis an(postProcMesh, false);

    if (fHybridH1EE->fExact) {
        an.SetExact(fHybridH1EE->fExact->ExactSolution());
        //an.SetExact(fHybridH1EE->loadCases.SingularityEx(),5); // DELETE-ME
    }


    auto errorVec = new TPZVec<REAL>;
    int64_t nErrorCols = 4;
    errorVec->resize(nErrorCols);
    errorVec->Fill(0);
    for (int64_t i = 0; i < nErrorCols; i++) {
        (*errorVec)[i] = 0;
    }

    int64_t nelem = postProcMesh->NElements();
    postProcMesh->LoadSolution(postProcMesh->Solution());
    postProcMesh->ExpandSolution();
    postProcMesh->ElementSolution().Redim(nelem, nErrorCols-1);
    for(int64_t el = 0; el<nelem; el++)
    {
        TPZCompEl *cel = postProcMesh->Element(el);
        TPZSubCompMesh *subc = dynamic_cast<TPZSubCompMesh *>(cel);
        if(subc)
        {
            int64_t nelsub = subc->NElements();
            subc->ElementSolution().Redim(nelsub, 6);
        }
    }

    std::ofstream mamesh("CheckingPostProcessMeshBeforePrinting.txt");
    postProcMesh->Print(mamesh);
    bool store=true;
    an.PostProcessError(*errorVec, store);//calculo do erro com sol exata e aprox e armazena no elementsolution

    std::cout << "\n\n############\n\n";

    TPZCompMeshTools::UnCondensedElements(postProcMesh);
    TPZCompMeshTools::UnGroupElements(postProcMesh);

    //Erro global
    std::ofstream myfile;
    myfile.open("HDivReconstructionErrors.txt", std::ios::app);
    myfile << "\n\n Estimator errors for HDiv reconstruction " << fHybridH1EE->fProblemConfig.problemname;
    myfile << "\n-------------------------------------------------- \n";
    myfile << "Ndiv = " << fHybridH1EE->fProblemConfig.ndivisions <<" Order k= " << fHybridH1EE->fProblemConfig.k << " Order n= "<< fHybridH1EE->fProblemConfig.n<<"\n";
    myfile << "e_{ex}: ||K^{0.5}.grad(u_h-u)|| = " << (*errorVec)[0] << "\n";
    myfile << "n_{F} : ||K^{0.5}.[grad(u_h)-invK.T_h]|| = " << (*errorVec)[1] << "\n";
    myfile << "||div(HDiv)-f||= "<< (*errorVec)[2] << "\n";

    myfile.close();
    DebugStop();
    PrintSolutionVTK(an);

}

void TPZHybridH1CreateRecMeshes::PrintSolutionVTK(TPZAnalysis &an){

    TPZMaterial *mat = fHybridH1EE->fPostProcMesh.FindMaterial(*fHybridH1EE->fProblemConfig.materialids.begin());
    int varindex = -1;
    if (mat) varindex = mat->VariableIndex("PressureFem");
    if (varindex != -1) {
        TPZStack<std::string> scalnames, vecnames;
        if (fHybridH1EE->fExact) {
            scalnames.Push("EnergyErrorExact");
            vecnames.Push("FluxExact");
        }
        scalnames.Push("PressureFEM");
        vecnames.Push("FluxFem");
        vecnames.Push("FluxSigmaReconstructed");
        scalnames.Push("POrder");


        int dim = fHybridH1EE->fPostProcMesh.Reference()->Dimension();

        std::stringstream out;
        out << fHybridH1EE->fProblemConfig.dir_name << "/" << fHybridH1EE->fProblemConfig.problemname
            << "_k_" << fHybridH1EE->fProblemConfig.k << "_n_"
            << fHybridH1EE->fProblemConfig.n;
        if (fHybridH1EE->fProblemConfig.ndivisions != -1) {
            out << "_Ndiv_" << fHybridH1EE->fProblemConfig.ndivisions;
        }
        if (fHybridH1EE->fProblemConfig.adaptivityStep != -1) {
            out << "_AdaptivityStep_" << fHybridH1EE->fProblemConfig.adaptivityStep;
        }
        out << ".vtk";

        int res =2;
        if(fHybridH1EE->fOriginal->NEquations()<100){
            res=6;
        }
        an.DefineGraphMesh(dim, scalnames, vecnames, out.str());
        an.PostProcess(res, dim);
    }
    else {
        std::cout << __PRETTY_FUNCTION__ << "\nPost Processing variable not found!\n";
    }
}