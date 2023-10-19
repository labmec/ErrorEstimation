//
// Created by victor on 16/02/23.
//

#include "TPZAnalysis.h"
#include "TPZCompMeshTools.h"
#include "TPZCreateMultiphysicsSpace.h"
#include "TPZElementMatrixT.h"
#include "TPZHybridH1CreateHDivReconstruction.h"
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
#include "TPZHybridH1ReconstructionBase.h"


TPZHybridH1CreateHDivReconstruction::~TPZHybridH1CreateHDivReconstruction()
{
    // prevent the base class from deleting a mesh component of the original mesh
    fMultiphysicsReconstructionMesh->MeshVector()[4] = 0;
}

// a method for generating the HDiv mesh
TPZCompMesh *TPZHybridH1CreateHDivReconstruction::CreateFluxReconstructionHDivMesh()
{
    TPZCompMesh *HDivAtomicMesh = fOriginal->MeshVector()[0]->Clone(); // HDIV-BOUND elements: clone might be unnecessary
    HDivAtomicMesh->SetName("Flux Reconstruction");
#ifdef ERRORESTIMATION_DEBUG
    std::string command = "mkdir -p " + fFolderOutput;
    std::string dirPath = fFolderOutput + "/";
    system(command.c_str());
    {
        std::ofstream outCon(dirPath + "OriginalFluxConnects.txt");
        if(HDivAtomicMesh->Dimension() == 2) {
            //TPZCompMeshTools::PrintConnectInfoByGeoElement(cmeshHdiv, outCon, {}, false, true);
        }
        std::ofstream outOriginalP(dirPath + "OriginalFlux.txt");
        HDivAtomicMesh->Print(outOriginalP);
        std::ofstream foutOriginalP(dirPath + "fOriginalFlux.txt");
        fOriginal->MeshVector()[0]->Print(foutOriginalP);
        std::ofstream outGOriginalVTK(dirPath + "gFlux.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(HDivAtomicMesh, outGOriginalVTK);
    }
#endif

    // Verify neighbouring information of the HDiv-bound mesh (works only for 2D meshes)
    VerifyBoundaryFluxConsistency(HDivAtomicMesh);
    {
        std::ofstream myoutput("checkingfOriginalMesh.txt");
        fOriginal->Print(myoutput);
        myoutput.flush();
    }
    int meshdim = HDivAtomicMesh->Dimension();
    for (auto mat : fOriginal->MaterialVec()) {
        if (!dynamic_cast<TPZBndCondT<STATE> *>(mat.second)) {
            if (mat.second->Dimension() == meshdim) {
                auto mymat = new TPZNullMaterial(mat.first, mat.second->Dimension());
                HDivAtomicMesh->InsertMaterialObject(mymat);
            }
        }
    }

    HDivAtomicMesh->SetDefaultOrder(forderFEM_k);

    HDivAtomicMesh->AutoBuild();
    HDivAtomicMesh->InitializeBlock();

    TPZCompMeshTools::AdjustFluxPolynomialOrders(HDivAtomicMesh, forderFEM_n); //Increases internal flux order by "hdivmais"

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

TPZCompMesh *TPZHybridH1CreateHDivReconstruction::CreateFluxReconstructionL2Mesh(){

#ifdef ERRORESTIMATION_DEBUG
    {
        // std::ofstream outOriginalP("HdivCondFlux.txt");
        // cmeshHdiv->Print(outOriginalP);
    }
#endif

    TPZCompMesh *L2AtomicMesh = new TPZCompMesh(fOriginal->Reference());
    L2AtomicMesh->SetName("L2 mesh for flux reconstruction");
    int dimMesh = fOriginal->Reference()->Dimension();

    int potential_order = forderFEM_k + forderFEM_n;
    L2AtomicMesh->SetDefaultOrder(potential_order);
    L2AtomicMesh->SetDimModel(dimMesh);

    L2AtomicMesh->SetAllCreateFunctionsContinuous(); //H1 functions
    L2AtomicMesh->ApproxSpace().CreateDisconnectedElements(true);

    for(auto matid:(fmaterialids)){
        TPZNullMaterial<> *material = new TPZNullMaterial<>(matid); material->SetDimension(dimMesh);
        L2AtomicMesh->InsertMaterialObject(material);
    }
    L2AtomicMesh->AutoBuild();
    L2AtomicMesh->ExpandSolution();

    return L2AtomicMesh;
}

TPZCompMesh *TPZHybridH1CreateHDivReconstruction::CreateFluxReconstructionConstantMesh(){
    TPZCompMesh *constant = new TPZCompMesh(fOriginal->Reference());
    constant->SetName("Flux reconstruction constant");
    {
        for (auto matid:fmaterialids) {
            TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(matid);
            nullmat->SetDimension(fOriginal->Reference()->Dimension());
            nullmat->SetNStateVariables(1);
            constant->InsertMaterialObject(nullmat);
        }
        constant->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
        constant->SetDefaultOrder(0);
        constant->AutoBuild();
    }
    return constant;
}

TPZMultiphysicsCompMesh *TPZHybridH1CreateHDivReconstruction::CreateFluxReconstructionMesh()
{
    fHDivReconstructionAtomicMesh = CreateFluxReconstructionHDivMesh();
    TPZCompMesh *L2AtomicMesh = CreateFluxReconstructionL2Mesh();
    TPZCompMesh *gspace = CreateFluxReconstructionConstantMesh();

    TPZManVector<TPZCompMesh *> mesh_vectors(5, 0);
    TPZManVector<int> active(5, 0);

    mesh_vectors[0] = fHDivReconstructionAtomicMesh;
    mesh_vectors[1] = L2AtomicMesh;
    mesh_vectors[2] = gspace;
    mesh_vectors[3] = fOriginal->MeshVector()[3]->Clone(); // avg-space
    mesh_vectors[3]->SetName("Flux reconstruction avg pressure");
    mesh_vectors[4] = fOriginal->MeshVector()[1];

    active[0] = 1;
    active[1] = 1;
    active[2] = 1;

    // Insert materials into Multiphysics mesh
    // The Wrap and interface material does not need to be created
    for(auto mat: fOriginal->MaterialVec()) {
        TPZMatLaplacianHybrid *matlaplacian = dynamic_cast<TPZMatLaplacianHybrid *>(mat.second);
        // Casting lagrange material to error estimate material for bc and lagrange coefficient objects
        if (matlaplacian) {
            TPZHybridH1HdivFluxRecMaterial *EEMat = new TPZHybridH1HdivFluxRecMaterial(*matlaplacian);
            fMultiphysicsReconstructionMesh->InsertMaterialObject(EEMat);
            EEMat->SetForcingFunction(matlaplacian->ForcingFunction(), matlaplacian->ForcingFunctionPOrder());

            for(auto trybcmat: fOriginal->MaterialVec()) {
                auto *bc = dynamic_cast<TPZBndCondT<STATE> *>(trybcmat.second);
                if (bc) {
                    TPZMaterial *matclone = trybcmat.second->NewMaterial();
                    bc = dynamic_cast<TPZBndCondT<STATE> *>(matclone);
                    // add bc mtf;
                    if(bc->Material()->Id() != EEMat->Id()){
                        continue;
                    }
                    bc->SetMaterial(EEMat);
                    fMultiphysicsReconstructionMesh->InsertMaterialObject(bc);
                }
            }
        }
        if(mat.first == fLagrangeMatId){
            // add lagrange material
            fMultiphysicsReconstructionMesh->InsertMaterialObject(mat.second->NewMaterial());
        }
    }

    fMultiphysicsReconstructionMesh->SetAllCreateFunctionsMultiphysicElem();

    fMultiphysicsReconstructionMesh->BuildMultiphysicsSpace(active,mesh_vectors);

    // copy the integration rule order
    {
        auto gmesh = fOriginal->Reference();
        gmesh->ResetReference();
        fOriginal->LoadReferences();
        int64_t nel = fMultiphysicsReconstructionMesh->NElements();
        for(int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = fMultiphysicsReconstructionMesh->Element(el);
            auto mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
            auto gel = cel->Reference();
            auto mfcelorig = dynamic_cast<TPZMultiphysicsElement *>(gel->Reference());
            TPZVec<int> order(gel->Dimension());
            mfcelorig->GetIntegrationRule().GetOrder(order);
            mfcel->SetIntegrationRule(order[0]);
        }
    }
    bool keeponelagrangian = true, keepmatrix = false;
    TPZCompMeshTools::CreatedCondensedElements(fMultiphysicsReconstructionMesh, keeponelagrangian, keepmatrix);

#ifdef ERRORESTIMATION_DEBUG
    {
        // std::ofstream outMultF("HdivMultMesh.txt");
        // HdivRecMesh->Print(outMultF);

        std::ofstream outMultOriginal("HdivMultOriginal.txt");
        fOriginal->Print(outMultOriginal);

        // std::ofstream outL2("L2Mesh.txt");
        // HdivRecMesh->MeshVector()[1]->Print(outL2);

        // std::ofstream outH1("H1Mesh.txt");
        // HdivRecMesh->MeshVector()[4]->Print(outH1);
    }
#endif

    ComputeElementStiffnesses();

    fMultiphysicsReconstructionMesh->LoadSolution(fMultiphysicsReconstructionMesh->Solution());

    //fMultiphysicsReconstructionMesh->TransferMultiphysicsSolution();
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mesh_vectors, fMultiphysicsReconstructionMesh);

#define ERRORESTIMATION_DEBUG
#ifdef ERRORESTIMATION_DEBUG
    {
        //std::ofstream outFCon(fFolderOutput + "HdivFluxAfterLoadSolConnects.txt");
        //TPZCompMeshTools::PrintConnectInfoByGeoElement(cmeshHdiv, outFCon, {}, false, true);
        std::ofstream outF(fFolderOutput + "HdivRecAtomic.txt");
        fMultiphysicsReconstructionMesh->MeshVector()[0]->Print(outF);
//        std::ofstream outMultF(fFolderOutput + "HDivRecMulti.txt");
//        fMultiphysicsReconstructionMesh->Print(outMultF);
//        std::ofstream outvtk(fFolderOutput + "HDivRecAtomicGeoMesh.vtk");
//        TPZVTKGeoMesh::PrintCMeshVTK(fMultiphysicsReconstructionMesh->MeshVector()[0], outvtk);
    }
#endif

    return fMultiphysicsReconstructionMesh;
}

void TPZHybridH1CreateHDivReconstruction::VerifyBoundaryFluxConsistency(TPZCompMesh* fluxmesh){
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
        bool IsBC = (fbcmaterialids.find(matID) != fbcmaterialids.end());
        if(!IsBC && matID != fLagrangeMatId) {
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

TPZVec<REAL>  TPZHybridH1CreateHDivReconstruction::PostProcess(){
    
    TPZLinearAnalysis an(fMultiphysicsReconstructionMesh,RenumType::ENone);

    // The solution is expanded to store errors,
    // Therefore it is required to account for the original solution and the errors.
    int numErrors = 3;
    numErrors++;
    TPZVec<REAL> errorVec = ComputeErrors(&an,numErrors);

    std::cout << "\n############\n";
    std::cout << "Computing Error HDiv reconstruction\n";
    std::cout << "||K^0.5(Grad(u_h)-Grad(u))||:  \t" << (errorVec)[0] <<
               "\n||K^-0.5(KGrad(u_h)+t_h)||:    \t" << (errorVec)[1]<<
               "\n||div(t_h)-f||:                \t"<< (errorVec)[2]<<"\n";

    TPZCompMeshTools::UnCondensedElements(fMultiphysicsReconstructionMesh);
    TPZCompMeshTools::UnGroupElements(fMultiphysicsReconstructionMesh);

    //Erro global
    std::ofstream myfile;
    myfile.open("HDivReconstructionErrors.txt", std::ios::app);
    myfile << "\n\n Estimator errors for HDiv reconstruction " << *fproblemname;
    myfile << "\n-------------------------------------------------- \n";
    myfile << "Ndiv = " << fnDivisions <<" Order k= " << forderFEM_k << " Order n= "<< forderFEM_n<<"\n";
    myfile << "e_{ex}: ||K^{0.5}.grad(u_h-u)|| = " << (errorVec)[0] << "\n";
    myfile << "n_{F} : ||K^{0.5}.[grad(u_h)-invK.T_h]|| = " << (errorVec)[1] << "\n";
    myfile << "||div(HDiv)-f||= "<< (errorVec)[2] << "\n";

    myfile.close();

    PrintSolutionVTK(an);

    return errorVec;
}

void TPZHybridH1CreateHDivReconstruction::FillVTKoutputVariables(TPZStack<std::string> &scalnames,TPZStack<std::string> &vecnames){
    if (fExact) {
        scalnames.Push("KGradUh_minus_KGradU");
        vecnames.Push("minus_KGradU");
    }
    scalnames.Push("uh");
    scalnames.Push("residual");
    scalnames.Push("th_plus_KGradUh");
    vecnames.Push("minus_KGradUh");
    vecnames.Push("th");
    scalnames.Push("POrder");
}

