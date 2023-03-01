//
// Created by victor on 22/02/23.
//

#include "TPZHybridH1CreateH1Reconstruction.h"
#include "TPZBndCondT.h"
#include "TPZCompMeshTools.h"
#include "TPZGeoElSideAncestors.h"
#include "TPZGeoElSidePartition.h"
#include "TPZHybridH1ErrorEstimator.h"
#include "TPZMaterialT.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"
#include "pzintel.h"
#include "pzsubcmesh.h"
#include <TPZElementMatrixT.h>
#include <TPZVTKGeoMesh.h>
#include <pzbuildmultiphysicsmesh.h>
#include "DarcyFlow/TPZDarcyFlow.h"
#include "TPZHybridH1ErrorEstimator.h"

#include "TPZPressureProjection.h"

// Delete me eventually
#include "TPZHybridH1ErrorEstimateMaterial.h"
#include "TPZHybridH1PressureRecMaterial.h"

TPZHybridH1CreateH1Reconstruction::TPZHybridH1CreateH1Reconstruction(TPZHybridH1ErrorEstimator *pEstimator){ 
       fHybridH1EE = pEstimator;
       fOriginal = pEstimator->fOriginal;
       fPressureMesh = pEstimator->fOriginal->MeshVector()[1]->Clone();
       fMultiphysicsH1reconstructionMesh = new TPZMultiphysicsCompMesh(fOriginal->Reference());

        auto &matvec = fPressureMesh->MaterialVec();
        for(auto mat : matvec){
            int matdim = mat.second->Dimension();
            if(matdim == fPressureMesh->Dimension()){
                auto isbc = dynamic_cast<TPZBndCond *> (mat.second);
                if (isbc) continue;
                auto singleSpaceMat = new TPZHybridH1PressureSingleSpace(mat.first,matdim);
                matvec[mat.first] = singleSpaceMat;
            }
        }

       std::string foldername = fPressureReconstructionFolderOutput;
       foldername.pop_back();
       std::string command = "mkdir -p " + foldername;
       system(command.c_str());
}

TPZCompMesh *TPZHybridH1CreateH1Reconstruction::CreateH1ReconstructionMesh(){

    // Delete wrap comp-elements, and associate bc and create skeleton elements to H1 reconstruction atomic mesh.
    PrepareGeometricElements();

    // Build the multiphysics cmesh on which the pressure is to be reconstructed.
    BuildMultiphysicsSpace();

    {
        std::ofstream filecmeshTXT(fPressureReconstructionFolderOutput + "toto.txt");
        fMultiphysicsH1reconstructionMesh->Print(filecmeshTXT);
    }
    CreateGroupedAndCondensedElements();

    //Compute continuos pressure on the skeleton;
    MakeSkeletonContinuous();

    VerifySkeletonContinuity();

#ifdef ERRORESTIMATION_DEBUG
    {
        std::ofstream out("MeshWithSmoothPressure.txt");
        fMultiphysicsH1reconstructionMesh->Print(out);
        std::ofstream out2("PressureMeshSmooth.txt");
        fMultiphysicsH1reconstructionMesh->MeshVector()[1]->Print(out2);
    }
#endif

    ComputeElementStiffnesses();

#ifdef ERRORESTIMATION_DEBUG
    {
        std::string command = "mkdir -p " + fDebugDirName +  "/DebuggingLoadSol";
        system(command.c_str());
        std::string dirPath = fDebugDirName + "/DebuggingLoadSol/";
        std::ofstream out(dirPath + "MeshBeforeLoadSol.txt");
        fMultiphysicsH1reconstructionMesh->Print(out);
        std::ofstream out2(dirPath + "SolBeforeLoadSolution.nb");
        fMultiphysicsH1reconstructionMesh->Solution().Print("SolBeforeLoadSolution=",out2,EMathematicaInput);
        std::ofstream out3(dirPath + "PressureWithAverage.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(fMultiphysicsH1reconstructionMesh->MeshVector()[1], out3, {1,2,3}, false, true);
        std::ofstream out4(dirPath + "PressureWithAverage.txt");
    }
#endif
     fMultiphysicsH1reconstructionMesh->LoadSolution(fMultiphysicsH1reconstructionMesh->Solution());

#ifdef ERRORESTIMATION_DEBUG
    {
        std::string dirPath = fDebugDirName + "/DebuggingLoadSol/";
        std::ofstream out(dirPath + "MeshAfterLoadSol.txt");
        fMultiphysicsH1reconstructionMesh->Print(out);
        std::ofstream outP(dirPath + "PotentialAfterLoadSol.txt");
        fMultiphysicsH1reconstructionMesh->MeshVector()[1]->Print(outP);
        //fMultiphysicsH1reconstructionMesh->Solution().Print("SolAfterLoadSolution");
        std::ofstream out2(dirPath + "SolAfterLoadSolution.nb");
        fMultiphysicsH1reconstructionMesh->Solution().Print("SolAfterLoadSolution=",out2,EMathematicaInput);
        std::ofstream out3(dirPath + "PressureAfterLoadSolution.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(fMultiphysicsH1reconstructionMesh->MeshVector()[1], out3, {1,2,3}, false, true);
    }
#endif

    {
        fMultiphysicsH1reconstructionMesh->LoadSolutionFromMultiPhysics();

#ifdef ERRORESTIMATION_DEBUG
        {
            std::string dirPath = fDebugDirName + "/";
            std::ofstream outCon(dirPath + "PressureConnectsAfterLoadSolution.txt");
            TPZCompMeshTools::PrintConnectInfoByGeoElement(fMultiphysicsH1reconstructionMesh->MeshVector()[1], outCon, {1,2,3}, false, true);
            //fMultiphysicsH1reconstructionMesh->MeshVector()[0]->Print(out2);

        }
        VerifySolutionConsistency(PressureMesh());
#endif
    }
    return fPressureMesh;
}

void TPZHybridH1CreateH1Reconstruction::BuildMultiphysicsSpace(){

#ifdef ERRORESTIMATION_DEBUG
    {
        std::ofstream out("OriginalFlux.txt");
        fOriginal->MeshVector()[0]->Print(out);
        std::ofstream out2("OriginalPotential.txt");
        fOriginal->MeshVector()[1]->Print(out2);
        std::ofstream out3("OriginalMeshHybrid.txt");
        fMultiphysicsH1reconstructionMesh->Print(out3);
    }
#endif

    int dim = fOriginal->Dimension();
    fOriginal->CopyMaterials(*fMultiphysicsH1reconstructionMesh);

    TPZManVector<TPZCompMesh *> mesh_vectors(3, 0);
    TPZManVector<int> active(3, 0);

    mesh_vectors[0] = fOriginal->MeshVector()[0]; // Lagrange coeff. : used in material contributeBC
    mesh_vectors[1] = fPressureMesh; //sh
    mesh_vectors[2] = fOriginal->MeshVector()[1];// potential

    active[1] = 1;

    for (auto mat : fMultiphysicsH1reconstructionMesh->MaterialVec()) {
        TPZMatLaplacianHybrid *matlaplacian=
                dynamic_cast<TPZMatLaplacianHybrid *>(mat.second);
        if (matlaplacian) {
            TPZHybridH1PressureRecMaterial *newmat = new TPZHybridH1PressureRecMaterial(*matlaplacian);

            for (auto bcmat : fMultiphysicsH1reconstructionMesh->MaterialVec()) {
                TPZBndCondT<STATE> *bc = dynamic_cast<TPZBndCondT<STATE> *>(bcmat.second);
                if (bc) {
                    bc->SetMaterial(newmat);
                }
            }
            fMultiphysicsH1reconstructionMesh->DeleteMaterial(mat.first);
            fMultiphysicsH1reconstructionMesh->InsertMaterialObject(newmat);
        }
    }

    fMultiphysicsH1reconstructionMesh->SetAllCreateFunctionsMultiphysicElem();

    fMultiphysicsH1reconstructionMesh->BuildMultiphysicsSpace(active, mesh_vectors);
}

///  Insert BC material into the pressure mesh material vector,
///  Create computational element on BC elements
void TPZHybridH1CreateH1Reconstruction::AddBC2PressureMesh(){
    TPZCompMesh *mult = fOriginal;
    TPZGeoMesh *gmesh = fPressureMesh->Reference();

    // Insert BC materials in pressure reconstruction mesh
    std::set<int> bcMatIDs = fHybridH1EE->fProblemConfig.bcmaterialids;
    for (auto bcID : bcMatIDs) {
        TPZMaterial *mat = mult->FindMaterial(bcID);
        TPZBndCondT<STATE> *bc = dynamic_cast<TPZBndCondT<STATE> *>(mat);
        if (!bc) DebugStop();

        int volumetricMatId = bc->Material()->Id();
        TPZMaterial *pressuremat = fPressureMesh->FindMaterial(volumetricMatId);
        if (!pressuremat) DebugStop();
        TPZHybridH1PressureSingleSpace *press = dynamic_cast<TPZHybridH1PressureSingleSpace *>(pressuremat);
        TPZBndCondT<STATE> *newbc = press->CreateBC(pressuremat, bc->Id(), bc->Type(), bc->Val1(), bc->Val2());
        if (bc->HasForcingFunctionBC()) {
            newbc->SetForcingFunctionBC(bc->ForcingFunctionBC(),5);
        } else DebugStop();
        fPressureMesh->InsertMaterialObject(newbc);
    }

    // This map stores the index of the BC geometric element and a pointer
    // to the computational element of its volumetric neighbour
    std::map<int64_t, TPZCompElSide> bcGeoElToNeighCompEl;
    for (int64_t el = 0; el < gmesh->NElements(); el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel) continue;
        if(gel->HasSubElement()) continue;

        // Filters BC elements
        int matID = gel->MaterialId();
        if (bcMatIDs.find(matID) != bcMatIDs.end()) {

            TPZGeoElSide bcSide(gel, gel->NSides() - 1);
            TPZStack<TPZCompElSide> compNeighSides;
            bcSide.EqualLevelCompElementList(compNeighSides, 1, 0);
            if (compNeighSides.size() != 1) {
                DebugStop();
            }

            TPZCompEl *cel = compNeighSides[0].Element();
            if (!cel) DebugStop();

            bcGeoElToNeighCompEl.insert({gel->Index(), compNeighSides[0]});
        }
    }

    // Create BC elements. We reset references at each step to allow for a discontinuous mesh globally
    gmesh->ResetReference();
    for (const auto &it : bcGeoElToNeighCompEl) {
        int64_t bcElID = it.first;
        TPZCompElSide neighSide = it.second;
        TPZCompEl *neighCel = neighSide.Element();
        TPZGeoEl *bcGeoEl = gmesh->Element(bcElID);
        // We load the volumetric element reference so it shares its connects with the BC element to be created
        neighCel->LoadElementReference();
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(neighCel);
        if (!intel) DebugStop();
        // Get order of the neighbour side of the volumetric element
        int iside = neighSide.Side();
        int order = intel->Connect(iside).Order();
        fPressureMesh->SetDefaultOrder(order);

        TPZCompEl *bcCel = fPressureMesh->CreateCompEl(bcGeoEl);
        // Reset references so that future elements will not share connects with these two elements
        bcCel->Reference()->ResetReference();
        neighCel->Reference()->ResetReference();
    }
}

void TPZHybridH1CreateH1Reconstruction::CreateSkeletonElements() {

    TPZCompMesh* cmesh = fOriginal;
    TPZGeoMesh* gmesh = fOriginal->Reference();
    gmesh->ResetReference();
    cmesh->LoadReferences();
    int dim = gmesh->Dimension();

#ifdef ERRORESTIMATION_DEBUG
    {
        std::string dirPath = fDebugDirName + "/";
        std::ofstream fileVTK(dirPath + "GeoMeshBeforePressureSkeleton.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, fileVTK);
        std::ofstream fileTXT(dirPath + "GeoMeshBeforePressureSkeleton.txt");
        gmesh->Print(fileTXT);
    }
#endif

    // Creation of geometric elements
    const int nel = gmesh->NElements();
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->Element(iel);
        TPZCompEl* cel = gel->Reference();
        if (!cel) continue;
        if (gel->Dimension() != dim) continue;

        // Iterates through the sides of the element
        int nsides = gel->NSides();
        for (int iside = 0; iside < nsides; iside++) {
            TPZGeoElSide gelside(gel, iside);

            // Filters boundary sides
            if (gelside.Dimension() != dim - 1) continue;

            //Create Geometric element if there is no boundary neighbour and no skeleton neighbour were created.
            std::set<int> matIDs = fHybridH1EE->fProblemConfig.bcmaterialids;
            matIDs.insert(fHybridH1EE->fPressureSkeletonMatId);
            TPZGeoElSide neighSide = gelside.HasNeighbour(matIDs);
            if(!neighSide.Exists())
            {
                TPZGeoElBC gbc(gelside, fHybridH1EE->fPressureSkeletonMatId);
            }
        }
    }

#ifdef ERRORESTIMATION_DEBUG
    {
        std::ofstream fileVTK("GeoMeshAfterPressureSkeleton.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, fileVTK);
        std::ofstream fileTXT("GeoMeshAfterPressureSkeleton.txt");
        gmesh->Print(fileTXT);
    }
#endif

    // Create skeleton elements in pressure mesh
    TPZNullMaterial<> *skeletonMat = new TPZNullMaterial<>(fHybridH1EE->fPressureSkeletonMatId);
    TPZNullMaterialCS<> *skeletonMatCS = new TPZNullMaterialCS<>(fHybridH1EE->fPressureSkeletonMatId);

    skeletonMat->SetDimension(dim - 1);
    fPressureMesh->InsertMaterialObject(skeletonMat);
    fMultiphysicsH1reconstructionMesh->InsertMaterialObject(skeletonMatCS->NewMaterial());

    std::set<int> matIdSkeleton = { fHybridH1EE->fPressureSkeletonMatId };
    gmesh->ResetReference();

    fPressureMesh->ApproxSpace().CreateDisconnectedElements(true);
    fPressureMesh->AutoBuild(matIdSkeleton);
    fPressureMesh->ExpandSolution();

    // increase the order of the dim-1 elements to the maximum of both neighbouring elements
    IncreasePressureSideOrders();//malha da pressao

    // Restrain the spaces of smaller to larger skeleton elements on hanging nodes
    RestrainSkeletonSides();
}

void TPZHybridH1CreateH1Reconstruction::IncreasePressureSideOrders() {

    TPZCompMesh *cmesh = fPressureMesh;
    TPZGeoMesh *gmesh = cmesh->Reference();

    gmesh->ResetReference();
    cmesh->LoadReferences();

#ifdef ERRORESTIMATION_DEBUG
    std::string dirPath = fDebugDirName + "/";
    std::set<int> matIDs = fProblemConfig.materialids;
    matIDs.insert(fProblemConfig.bcmaterialids.begin(),fProblemConfig.bcmaterialids.end());
    matIDs.insert(fPressureSkeletonMatId);
    {
        std::ofstream outCon(dirPath + "PressureConnectsB4IncreaseSideOrder.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(cmesh, outCon, matIDs, false, true);
    }
#endif

    int OrigOrder = cmesh->GetDefaultOrder();
    int64_t nel = cmesh->NElements();
    int dim = cmesh->Dimension();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != dim - 1) {
            continue;
        }
        TPZMaterial *mat=cel->Material();

        //   std::cout<<"material "<<mat->Id()<<std::endl;
        TPZGeoElSide gelside(gel, gel->NSides() - 1);
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 1, 0);


        int nneigh = celstack.NElements();
        if (nneigh == 1) {
            TPZCompElSide celside = gelside.LowerLevelCompElementList2(1);
            if (!celside) {
                continue;//DebugStop();/// para nao incremenentar ordem na condicao de contorno
            }
            celstack.Push(celside);
            nneigh++;
        } else if (nneigh != 2) DebugStop();

        int maxOrder = 0;

        for (int ineigh = 0; ineigh < nneigh; ineigh++) {
            TPZInterpolatedElement *intelS = dynamic_cast<TPZInterpolatedElement *>(celstack[ineigh].Element());
            int orderEl = intelS->GetPreferredOrder();

            //   std::cout<<"ordem El "<<orderEl<< std::endl;

            maxOrder = (orderEl > maxOrder) ? orderEl : maxOrder;
        }

        // std::cout<<"max order "<<maxOrder<< std::endl;

        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        int nsides = gel->NSides();
        int ncorner = gel->NCornerNodes();

        for (int side = ncorner; side < nsides; side++) {
            if (intel->NSideConnects(side)) {
                intel->SetSideOrder(side, maxOrder);
            }
        }
        //        intel->Print();
    }
    cmesh->InitializeBlock();

#ifdef ERRORESTIMATION_DEBUG
    {
        std::ofstream outCon(dirPath + "PressureConnectsAFTERIncreaseSideOrder.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(cmesh, outCon, matIDs, false, true);
    }
#endif
}

void TPZHybridH1CreateH1Reconstruction::RestrainSkeletonSides() {

    TPZGeoMesh *gmesh = fPressureMesh->Reference();
    gmesh->ResetReference();
    fPressureMesh->LoadReferences();

#ifdef ERRORESTIMATION_DEBUG
    {
        std::string dirPath = fDebugDirName + '/';
        std::ofstream out(dirPath + "MeshBeforeRestrainSkeleton.txt");
        pressure_mesh->Print(out);
    }
#endif

    int64_t nel = gmesh->NElements();
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->Element(iel);
        TPZCompEl *cel = gel->Reference();
        if (!cel) continue;
        if (gel->MaterialId() != fHybridH1EE->fPressureSkeletonMatId) continue;

        // If the element is a small skeleton, restrain its highest dimension side and then its subsides
        int nsides = gel->NSides();
        TPZGeoElSide small(gel, nsides - 1);
        TPZGeoElSideAncestors ancestors(small);
        TPZGeoElSide largerNeigh = ancestors.HasLarger(fHybridH1EE->fPressureSkeletonMatId);
        if (!largerNeigh) continue;

        TPZCompEl *smallCel = small.Element()->Reference();
        if (!smallCel) DebugStop();
        TPZInterpolatedElement *smallIntel = dynamic_cast<TPZInterpolatedElement *>(smallCel);
        if (!smallIntel) DebugStop();

        TPZCompEl *largeCel = largerNeigh.Element()->Reference();
        if (!largeCel) DebugStop();
        TPZInterpolatedElement *largeIntel = dynamic_cast<TPZInterpolatedElement *>(largeCel);
        if (!largeIntel) DebugStop();

        TPZManVector<REAL, 3> xicenter(small.Element()->Dimension(), 0.);
        gel->CenterPoint(nsides - 1, xicenter);
        TPZManVector<REAL> xcenter(3, 0.);
        gel->X(xicenter, xcenter);
        //std::cout << "Restriction @ [" << xcenter << "]:"
        //          << "  Small El: " << small.Element()->Index() << ", Side: " << small.Side()
        //          << "  Large El: " << largerNeigh.Element()->Index() << ", Side: " << largerNeigh.Side() << "\n";
        smallIntel->RestrainSide(small.Side(), largeIntel, largerNeigh.Side());

        // Restrain subsides
        for (int iside = 0; iside < nsides - 1; iside++) {
            TPZGeoElSide subsmall(gel, iside);
            TPZGeoElSideAncestors subancestors(subsmall);
            TPZGeoElSide subLargerNeigh = subancestors.LargeSide(largerNeigh.Element());
            if (!subLargerNeigh) DebugStop();
            if (subLargerNeigh.Element() == subsmall.Element()) DebugStop();
            gel->CenterPoint(iside, xicenter);
            gel->X(xicenter, xcenter);
            //std::cout << "SubRestriction @ [" << xcenter << "]:"
            //          << "  Small El: " << small.Element()->Index() << ", Side: " << subsmall.Side()
            //         << "  Large El: " << largerNeigh.Element()->Index() << ", Side: " << subLargerNeigh.Side() << "\n";
            smallIntel->RestrainSide(subsmall.Side(), largeIntel, subLargerNeigh.Side());
        }
    }

    fPressureMesh->CleanUpUnconnectedNodes();
#ifdef ERRORESTIMATION_DEBUG
    {
        std::string dirPath = fDebugDirName + '/';
        std::ofstream out(dirPath + "MeshAfterRestrainSkeleton.txt");
        pressure_mesh->Print(out);
    }
#endif
}

/// create dim-2 skeleton mesh based on the dim-1 faces
// will do nothing if the dimension of the mesh == 2
void TPZHybridH1CreateH1Reconstruction::CreateEdgeSkeletonMesh() {

    if (fPressureMesh->MaterialVec().find(fHybridH1EE->fPressureSkeletonMatId) != fPressureMesh->MaterialVec().end()) {
        DebugStop();
    }
    TPZNullMaterialCS<> *nullmat = new TPZNullMaterialCS<>(fHybridH1EE->fPressureSkeletonMatId);
    fPressureMesh->InsertMaterialObject(nullmat);
    int dim = fMultiphysicsH1reconstructionMesh->Dimension();
    int64_t nel = fPressureMesh->NElements();
    std::map<int64_t, int> gelpressures;
    // create the geometrical elements
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fPressureMesh->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = cel->Reference();
        if (!gel) DebugStop();
        //        int matid = gel->MaterialId();
        if (gel->Dimension() != dim - 1) continue;
        int ncorner = gel->NCornerNodes();
        int nsides = gel->NSides();
        int polynomialorder = intel->Connect(nsides - 1).Order();
        for (int side = ncorner; side < nsides - 1; side++) {
            if (gel->SideDimension(side) != dim - 2) DebugStop();
            TPZGeoElSide gelside(gel, side);
            TPZGeoElSide hasneigh = HasNeighbour(gelside, fHybridH1EE->fPressureSkeletonMatId);
            if (!hasneigh) {
                TPZGeoElBC gbc(gelside, fHybridH1EE->fPressureSkeletonMatId);
                TPZGeoEl *createdelement = gbc.CreatedElement();
                hasneigh = TPZGeoElSide(createdelement, createdelement->NSides() - 1);
                gelpressures[createdelement->Index()] = polynomialorder;
            } else {
                int64_t gelindex = hasneigh.Element()->Index();
#ifdef ERRORESTIMATION_DEBUG
                if (gelpressures.find(gelindex) == gelpressures.end()) {
                    DebugStop();
                }
#endif
                int polorder = gelpressures[gelindex];
                if (polorder != polynomialorder) {
                    polorder = std::max(polorder, polynomialorder);
                    gelpressures[gelindex] = polorder;
                }
            }
        }
    }
    // create the pressure computational elements
    // we assume there are no pressure elements
    fPressureMesh->Reference()->ResetReference();
    TPZGeoMesh *gmesh = fPressureMesh->Reference();
    nel = gmesh->NElements();
    for (auto indexpair : gelpressures) {
        int64_t index = indexpair.first;
        int polynomialorder = indexpair.second;
        TPZGeoEl *gel = gmesh->Element(index);
        if (!gel) DebugStop();
        TPZCompEl *cel = 0;
        fPressureMesh->SetDefaultOrder(polynomialorder);
        cel = fPressureMesh->ApproxSpace().CreateCompEl(gel, *fPressureMesh);
#ifdef ERRORESTIMATION_DEBUG
        {
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            if (!intel) DebugStop();
            int porder = intel->GetPreferredOrder();
            if (porder != polynomialorder) DebugStop();
        }
#endif
        gel->ResetReference();
    }
    AdjustNeighbourPolynomialOrders();
    fPressureMesh->ExpandSolution();
    RestrainSmallEdges();
}

/// searches for a neighbour whose element has the proper dimension and materialid
TPZGeoElSide TPZHybridH1CreateH1Reconstruction::HasNeighbour(const TPZGeoElSide &gelside, int matid) {
    TPZGeoElSide neighbour = gelside.Neighbour();
    int dim = gelside.Dimension();
    while (neighbour != gelside) {
        if (neighbour.Element()->Dimension() == dim && neighbour.Element()->MaterialId() == matid) {
            return neighbour;
        }
        neighbour = neighbour.Neighbour();
    }
    return TPZGeoElSide();
}

/// restrain the edge elements that have larger elements as neighbours
void TPZHybridH1CreateH1Reconstruction::RestrainSmallEdges() {
    //    TPZCompMesh *fPressureMesh = fMultiphysicsH1reconstructionMesh->MeshVector()[1];
    TPZGeoMesh *gmesh = fPressureMesh->Reference();
    gmesh->ResetReference();
    int dim = fMultiphysicsH1reconstructionMesh->Dimension();
    int64_t nel = fPressureMesh->NElements();
    // load the face and edge elements
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fPressureMesh->Element(el);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (!gel) DebugStop();
        int geldim = gel->Dimension();
        if (geldim == dim - 2) {
            gel->SetReference(cel);
        }
    }
    // look for elements that neighbour a larger element
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fPressureMesh->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = cel->Reference();
        if (!gel) DebugStop();
        int geldim = gel->Dimension();
        if (geldim != dim - 2) {
            continue;
        }
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel, nsides - 1);
        bool onlyinterpolated = true;
        TPZCompElSide large_celside = gelside.LowerLevelCompElementList2(onlyinterpolated);
        if (large_celside) {
            TPZInterpolatedElement *largeintel = dynamic_cast<TPZInterpolatedElement *>(large_celside.Element());
            if (!largeintel) DebugStop();
            int largeside = large_celside.Side();
            intel->RestrainSide(nsides - 1, largeintel, largeside);
            // restrain the corner nodes
            for (int side = 0; side < nsides - 1; side++) {
                TPZGeoElSide gelside_small(gel, side);
                TPZCompElSide celside_restraint = gelside_small.LowerLevelCompElementList2(onlyinterpolated);
                if (celside_restraint) {
                    TPZInterpolatedElement *largeintel = dynamic_cast<TPZInterpolatedElement *>(celside_restraint.Element());
                    if (!largeintel) DebugStop();
                    int largeside = large_celside.Side();
                    intel->RestrainSide(side, largeintel, largeside);
                }
            }
        }
    }
}

/// adjust the interpolation orders so as to create an H1/2 boundary mesh
// this method is called by the CreateEdgeSkeletonMesh method
void TPZHybridH1CreateH1Reconstruction::AdjustNeighbourPolynomialOrders() {
    //    TPZCompMesh *pressure_mesh = fMultiphysicsH1reconstructionMesh->MeshVector()[1];
    TPZGeoMesh *gmesh = fPressureMesh->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    // load the elements of lower dimension than dim
    int64_t nel = fPressureMesh->NElements();
    std::map<std::pair<int64_t, int>, int> polynomialorders;
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fPressureMesh->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() >= dim) continue;
        gel->SetReference(cel);
    }
    bool changed = true;
    while (changed) {
        changed = false;
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = fPressureMesh->Element(el);
            if (!cel) continue;
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            if (!intel) DebugStop();
            TPZGeoEl *gel = intel->Reference();
            if (gel->Dimension() >= dim) continue;
            int nsides = gel->NSides();
            int adjustporder = 0;
            for (int side = 0; side < nsides; side++) {
                if (gel->SideDimension(side) < 1) continue;
                int nconn = intel->NSideConnects(side);
                int porder = intel->SideConnect(nconn - 1, side).Order();
                std::pair<int64_t, int> elside(el, side);
                if (polynomialorders.find(elside) != polynomialorders.end()) {
                    adjustporder = polynomialorders[elside];
                    porder = adjustporder;
                }
                int maxorder = porder;
                // verify if any neighbour has a different polynomial order
                int onlyinterpolated = true;
                int removeduplicates = false;
                TPZStack<TPZCompElSide> celstack;
                TPZGeoElSide gelside(gel, side);
                gelside.EqualLevelCompElementList(celstack, onlyinterpolated, removeduplicates);
                TPZCompElSide large = gelside.LowerLevelCompElementList2(onlyinterpolated);
                if (large) {
                    celstack.Push(large);
                }
                celstack.Push(gelside.Reference());
                int nequal = celstack.size();
                // compute the maximum polynomial order of all neighbours
                for (int ieq = 0; ieq < nequal; ieq++) {
                    TPZCompEl *celneigh = celstack[ieq].Element();
                    int celside = celstack[ieq].Side();
                    TPZInterpolatedElement *intelneigh = dynamic_cast<TPZInterpolatedElement *>(celneigh);
                    if (!intelneigh) DebugStop();
                    int nneighconnects = intelneigh->NSideConnects(celside);
                    int neighporder = intelneigh->SideConnect(nneighconnects - 1, celside).Order();
                    std::pair<int64_t, int> neighsideorder(intelneigh->Index(), celside);
                    if (polynomialorders.find(neighsideorder) != polynomialorders.end()) {
                        adjustporder = polynomialorders[neighsideorder];
                        neighporder = adjustporder;
                    }
                    if (neighporder > maxorder) maxorder = neighporder;
                }
                // verify if the polynomial order of a side needs to be adjusted
                for (int ieq = 0; ieq < nequal; ieq++) {
                    TPZCompEl *celneigh = celstack[ieq].Element();
                    int celside = celstack[ieq].Side();
                    TPZInterpolatedElement *intelneigh = dynamic_cast<TPZInterpolatedElement *>(celneigh);
                    if (!intelneigh) DebugStop();
                    int nneighconnects = intelneigh->NSideConnects(celside);
                    int neighporder = intelneigh->SideConnect(nneighconnects - 1, celside).Order();
                    std::pair<int64_t, int> neighsideorder(intelneigh->Index(), celside);
                    if (polynomialorders.find(neighsideorder) != polynomialorders.end()) {
                        adjustporder = polynomialorders[neighsideorder];
                        neighporder = adjustporder;
                    }
                    if (neighporder != maxorder) {
                        std::pair<int64_t, int> neighside(intelneigh->Index(), celside);
                        polynomialorders[neighside] = maxorder;
                        changed = true;
                    }
                }
            }
        }
    }
    gmesh->ResetReference();
    for (auto it : polynomialorders) {
        int64_t index = it.first.first;
        int side = it.first.second;
        int porder = it.second;
        TPZCompEl *cel = fPressureMesh->Element(index);
        if (!cel) DebugStop();
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        intel->SetSideOrder(side, porder);
    }
}

/// Create and delete geometric elements
void TPZHybridH1CreateH1Reconstruction::PrepareGeometricElements(){
    TPZGeoMesh *gmesh = fPressureMesh->Reference();
#define ERRORESTIMATION_DEBUG55
#ifdef ERRORESTIMATION_DEBUG55
    {
        std::ofstream outTXT(fPressureReconstructionFolderOutput + "OriginalPressureMesh.txt");
        std::ofstream outVTK(fPressureReconstructionFolderOutput + "OriginalPressureMesh.vtk");
        fPressureMesh->Print(outTXT);
        TPZVTKGeoMesh::PrintCMeshVTK(fPressureMesh, outVTK);
    }
#endif

    gmesh->ResetReference();
    int dim = gmesh->Dimension();

    // Delete compels of dimension dim - 1
    // Case HybridH1 : Delete all MatWrapId computational elements
    for (int64_t el = 0; el < fPressureMesh->NElements(); el++) {
        TPZCompEl *cel = fPressureMesh->Element(el);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (!gel) DebugStop();
        if (gel->Dimension() == dim - 1) {
            delete cel;
        }
    }
    //RemoveNullCompEl(pressureMesh);

    fPressureMesh->ComputeNodElCon();
    fPressureMesh->CleanUpUnconnectedNodes();
    fPressureMesh->LoadReferences();
    fPressureMesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    fPressureMesh->ApproxSpace().CreateDisconnectedElements(true);

    //  Insert BC material into the pressure mesh material vector,
    //  Create computational element on BC elements
    AddBC2PressureMesh();

    // Create internal skeleton elements for pressure reconstruction
    CreateSkeletonElements();

    // TODO : Make it work after 2D works
    if (dim == 3) {
        CreateEdgeSkeletonMesh();
    }

#ifdef ERRORESTIMATION_DEBUG55
    {
        std::ofstream outTXT(fPressureReconstructionFolderOutput + "PressureMeshAfterPrepareGeometricElems.txt");
        std::ofstream outVTK(fPressureReconstructionFolderOutput + "PressureMeshAfterPrepareGeometricElems.vtk");
        fPressureMesh->Print(outTXT);
        TPZVTKGeoMesh::PrintCMeshVTK(fPressureMesh, outVTK);
    }
#endif
}

void TPZHybridH1CreateH1Reconstruction::MakeSkeletonContinuous(){

    {
        std::ofstream out(fPressureReconstructionFolderOutput+"mypmesh.txt");
        fPressureMesh->Print(out);
    }

    //Compute weights for bc and volumetric elements,
    //if bc, uses bignumber for Dirichlet and 0 for Neumann;
    //if Volumetric, uses the material permeability function/constant.
    ComputePressureWeights();

    {
        std::ofstream out(fPressureReconstructionFolderOutput+"mypmeshComputePressureWeights.txt");
        fPressureMesh->Print(out);
    }

    int target_dim = 1;
    // Update the mesh solution at dirichlet boundary locations
    // with exact values
    ComputeBoundaryL2Projection(fPressureMesh, target_dim);

     {
        std::ofstream out(fPressureReconstructionFolderOutput+"mypmeshComputeBoundaryL2Projection.txt");
        fPressureMesh->Print(out);
    }

    // Calculates average pressure on interface edges and vertices
    int dim = fPressureMesh->Dimension();
    ComputeAveragePressures(dim - 1);
    // in three dimensions make the one-d polynoms compatible
    if (dim == 3) {
        ComputeAveragePressures(1);
    }

    {
        std::ofstream out(fPressureReconstructionFolderOutput+"mypmeshComputeAveragePressures.txt");
        fPressureMesh->Print(out);
    }


    {
        std::string dirPath = fHybridH1EE->fDebugDirName + "/";
        std::ofstream out(dirPath+"PressureAverageMesh.txt");
        fPressureMesh->Print(out);
        //PlotLagrangeMultiplier(dirPath+"BeforeNodalAverage");
        dirPath = "DebuggingConsistency/";
        std::ofstream outCon(dirPath + "AverageB4NodalAverage.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(fPressureMesh, outCon, {fHybridH1EE->fPressureSkeletonMatId}, false, true);
    }

    ComputeNodalAverages();

    {
        std::string dirPath = "DebuggingConsistency/";
        std::ofstream outCon(dirPath + "AverageAfterNodalAverage.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(fPressureMesh, outCon, {fHybridH1EE->fPressureSkeletonMatId}, false, true);
        dirPath = fHybridH1EE->fDebugDirName + "/";
        std::ofstream out(dirPath + "PressureNodalMesh.txt");
        fPressureMesh->Print(out);
        //PlotLagrangeMultiplier("AfterNodalAverage");
        std::ofstream outCon2(dirPath + "AverageAfterLoadSolAverage.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(fPressureMesh, outCon2, {1,2,3,fHybridH1EE->fPressureSkeletonMatId}, false, true);
    }

    CopySolutionFromSkeleton();

    // transfer the continuous pressures to the multiphysics space
    {
        //TPZManVector<TPZCompMesh *, 1> meshvec(1);
        //meshvec[0] = fMultiphysicsH1reconstructionMesh->MeshVector()[1];//pressure
        fMultiphysicsH1reconstructionMesh->CleanUpUnconnectedNodes();


        std::string command = "mkdir -p " + fHybridH1EE->fDebugDirName +  "/DebuggingTransfer";
        system(command.c_str());

        {
            std::string dirPath = fHybridH1EE->fDebugDirName +  "/DebuggingTransfer/";
            std::ofstream out(dirPath + "PressureBeforeTransferFromMeshes.txt");
            TPZCompMeshTools::PrintConnectInfoByGeoElement(fPressureMesh, out);

            std::ofstream outMultiphysics(dirPath + "MultiphysicsBeforeTransferFromMeshes.txt");
            std::set<int> matIDs;
            GetPressureMatIDs(matIDs);
            TPZCompMeshTools::PrintConnectInfoByGeoElement(fMultiphysicsH1reconstructionMesh, outMultiphysics,matIDs);
        }
        fMultiphysicsH1reconstructionMesh->LoadSolutionFromMeshes();
        {
            std::string dirPath = fHybridH1EE->fDebugDirName +  "/DebuggingTransfer/";
            std::ofstream out(dirPath + "PressureAfterTransferFromMeshes.txt");
            TPZCompMeshTools::PrintConnectInfoByGeoElement(fMultiphysicsH1reconstructionMesh->MeshVector()[1], out);
            std::ofstream outMultiphysics(dirPath + "MultiphysicsAfterTransferFromMeshes.txt");
            std::set<int> matIDs;
            GetPressureMatIDs(matIDs);
            TPZCompMeshTools::PrintConnectInfoByGeoElement(fMultiphysicsH1reconstructionMesh, outMultiphysics,matIDs);
        }
    }
}

void TPZHybridH1CreateH1Reconstruction::CreateGroupedAndCondensedElements() {

    // This vector stores the connects from elements which have a neighbour of
    // an internal boundary material. We don't want to condense these connects,
    // so we are later incrementing the number of elements connected to them.
    // Then we compute the stiffness matrix and load the solution of the
    // internal degrees of freedom.
    TPZManVector<int64_t> connectsToIncrement(fMultiphysicsH1reconstructionMesh->NConnects(), -1);
    fMultiphysicsH1reconstructionMesh->ComputeNodElCon();

    fPressureMesh->LoadReferences();

#ifdef ERRORESTIMATION_DEBUG
    {
        std::ofstream txtPostProcMesh(fPressureReconstructionFolderOutput + "PostProcMeshB4IncrementingConnects.txt");
        fMultiphysicsH1reconstructionMesh->ShortPrint(txtPostProcMesh);
        std::ofstream txtPressureMesh(dirPath + "PressureMeshB4IncrementingConnects.txt");
        pressureMesh->Print(txtPressureMesh);
    }
#endif

    // Stock the index of the connects belonging to the side of volumetric elements
    // facing skeleton objects.
    for (int64_t el = 0; el < fPressureMesh->NElements(); el++) {
        TPZCompEl *cel = fPressureMesh->Element(el);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (!gel) continue;

        if (gel->MaterialId() != fHybridH1EE->fPressureSkeletonMatId) continue;

        TPZGeoElSide skelSide(gel, gel->NSides() - 1);
        //if (skelSide.HasNeighbour(fProblemConfig.bcmaterialids)) continue;

        TPZStack<TPZCompElSide> compNeighSides;
        skelSide.EqualLevelCompElementList(compNeighSides, 1, 0);
        if (compNeighSides.size() ==0) DebugStop();

        for (int i = 0; i < compNeighSides.size(); i++) {
            TPZCompEl *neighCel = compNeighSides[i].Element();
            TPZInterpolatedElement *neighIntEl = dynamic_cast<TPZInterpolatedElement *>(neighCel);
            if (!neighIntEl) DebugStop();

            int sideNum = compNeighSides[i].Side();
            int nCon = neighIntEl->NSideConnects(sideNum);

            for (int iCon = 0; iCon < nCon; iCon++) {
                int64_t conindex = neighIntEl->SideConnectIndex(iCon, sideNum);
                connectsToIncrement[conindex] = 1;
            }
        }
    }

    // Create TPZElementGroup, grouping volumetric and boundry elements
    fMultiphysicsH1reconstructionMesh->LoadReferences();
    int64_t nel = fMultiphysicsH1reconstructionMesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fMultiphysicsH1reconstructionMesh->Element(el);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (!gel) DebugStop();

        if (gel->Dimension() != fMultiphysicsH1reconstructionMesh->Dimension()) continue;
        TPZMultiphysicsElement *mcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mcel) DebugStop();

        // Insert volumetric element index and associated bc element index into elementsToGroup set.
        std::set<int64_t> elementsToGroup;
        for (int iSide = 0; iSide < gel->NSides(); iSide++) {
            TPZGeoElSide side(gel, iSide);
            if (side.Dimension() != gel->Dimension() - 1) continue;

            TPZStack<TPZCompElSide> compNeighSides;
            side.EqualLevelCompElementList3(compNeighSides, 1, 0);

            for (int i = 0; i < compNeighSides.size(); i++) {
                TPZCompElSide compSide = compNeighSides[i];
                TPZCompEl *neighCel = compSide.Element();
                TPZMaterial *mat = neighCel->Material();
                TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
                if (bc) {
                    elementsToGroup.insert(mcel->Index());
                    elementsToGroup.insert(neighCel->Index());
                }
            }
        }

        if (elementsToGroup.size()) {
            TPZElementGroup *elGroup = new TPZElementGroup(*fMultiphysicsH1reconstructionMesh);
            for (const auto &it : elementsToGroup) {
                elGroup->AddElement(fMultiphysicsH1reconstructionMesh->Element(it));
            }
        }
    }

    fMultiphysicsH1reconstructionMesh->ComputeNodElCon();
#define ERRORESTIMATION_DEBUG33
#ifdef ERRORESTIMATION_DEBUG33
    {
        std::ofstream fileVTK(fPressureReconstructionFolderOutput + "GeoMeshBeforeCondensedCompel.vtk");
        TPZGeoMesh *gmesh = fMultiphysicsH1reconstructionMesh->Reference();
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, fileVTK);
        std::ofstream fileTXT(fPressureReconstructionFolderOutput + "GeoMeshBeforeCondensedCompel.txt");
        gmesh->Print(fileTXT);
        std::ofstream filecmeshTXT(fPressureReconstructionFolderOutput + "CompMeshBeforeCondensedCompel.txt");
        fMultiphysicsH1reconstructionMesh->Print(filecmeshTXT);
    }
#endif

    for (int64_t el = 0; el < fMultiphysicsH1reconstructionMesh->NElements(); el++) {
        TPZCompEl *cel = fMultiphysicsH1reconstructionMesh->Element(el);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            TPZElementGroup *group = dynamic_cast<TPZElementGroup *>(cel);
            if (!group) DebugStop();
        }
        if (gel && gel->Dimension() != fMultiphysicsH1reconstructionMesh->Dimension()) continue;
        TPZCondensedCompEl *condense = new TPZCondensedCompEl(cel, false);
    }

    for (auto matit : fMultiphysicsH1reconstructionMesh->MaterialVec()) {
        TPZMaterial *mat = matit.second;
        TPZHybridH1ErrorEstimateMaterial *errormat = dynamic_cast<TPZHybridH1ErrorEstimateMaterial *>(mat);
        if (errormat) {
            errormat->fNeumannLocalProblem = false;
        }
    }

    fMultiphysicsH1reconstructionMesh->CleanUpUnconnectedNodes();
    fPressureMesh->CleanUpUnconnectedNodes();

#ifdef ERRORESTIMATION_DEBUG33
    std::ofstream fileVTK(fPressureReconstructionFolderOutput + "GeoMeshAfterCondensedCompel.vtk");
    TPZGeoMesh *gmesh = fMultiphysicsH1reconstructionMesh->Reference();
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, fileVTK);
    std::ofstream fileTXT(fPressureReconstructionFolderOutput + "GeoMeshAfterCondensedCompel.txt");
    gmesh->Print(fileTXT);
    std::ofstream filecmeshTXT(fPressureReconstructionFolderOutput + "CompMeshAfterCondensedCompel.txt");
    fMultiphysicsH1reconstructionMesh->Print(filecmeshTXT);
#endif
}

void TPZHybridH1CreateH1Reconstruction::VerifySkeletonContinuity(){
    TPZGeoMesh *gmesh = fPressureMesh->Reference();
    gmesh->ResetReference();
    int nel = fPressureMesh->NElements();
    for(int iel = 0; iel < nel ;iel++) {
        TPZCompEl *cel = fPressureMesh->Element(iel);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if (gel->MaterialId() == fHybridH1EE->fPressureSkeletonMatId) {
            gel->SetReference(cel);
        }
    }

    TPZBlock &block =  fPressureMesh->Block();
    TPZFMatrix<STATE> &sol = fPressureMesh->Solution();
    for(int iel = 0; iel < nel ; iel++) {
        TPZCompEl *cel = fPressureMesh->Element(iel);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if (gel->MaterialId() != fHybridH1EE->fPressureSkeletonMatId) continue;
        for (int iside = 0; iside < gel->NSides(); iside++) {
            TPZGeoElSide gelside(gel, iside);
            if (gelside.Dimension() > 1) DebugStop();  // 3D not supported yet
            if (gelside.Dimension() == 0) {
                TPZCompElSide celside(cel, gelside.Side());
                TPZConnect &elcon = fPressureMesh->ConnectVec()[celside.ConnectIndex()];
                if (elcon.NShape() != 1 || elcon.NState() != 1) DebugStop();
                STATE elsol = sol.at(block.at(elcon.SequenceNumber(), 0, 0, 0));
                TPZStack<TPZCompElSide> neighSides;
                gelside.ConnectedCompElementList(neighSides, 1, 0);

                for (int ind = 0; ind < neighSides.size(); ind++) {
                    TPZCompElSide neighCelSide = neighSides[ind];
                    if (!(neighCelSide.Element())) continue;
                    TPZGeoElSide neighGelSide = neighCelSide.Reference();
                    if (!(neighGelSide.Element()) || neighGelSide.Element()->MaterialId() != fHybridH1EE->fPressureSkeletonMatId)
                        DebugStop();
                    if (neighGelSide.Dimension() != 0) continue;

                    TPZConnect &c1 = fPressureMesh->ConnectVec()[neighCelSide.ConnectIndex()];
                    if (c1.NShape() != 1 || c1.NState() != 1) DebugStop();
                    STATE c1sol = sol.at(block.at(c1.SequenceNumber(), 0, 0, 0));

#ifdef LOG4CXX
                    if (logger->isDebugEnabled()) {
                        if (!IsZero(elsol - c1sol)) {
                            std::stringstream sout;

                            TPZManVector<REAL> x(3);
                            TPZManVector<REAL> pt(gelside.Dimension(), 0);
                            gelside.X(pt, x);

                            TPZManVector<REAL> xcenter(3, 0);
                            TPZGeoElSide gelHigherSide(gelside.Element());
                            TPZManVector<REAL> ptcenter(gelHigherSide.Dimension(), 0);
                            gelHigherSide.X(ptcenter, xcenter);

                            TPZManVector<REAL> xneigh(3);
                            TPZManVector<REAL> ptneigh(neighGelSide.Dimension(), 0);
                            neighGelSide.X(ptneigh, xneigh);

                            TPZManVector<REAL> xneighcenter(3, 0);
                            TPZGeoElSide neighHigherSide(neighGelSide.Element());
                            TPZManVector<REAL> ptneighcenter(neighHigherSide.Dimension(), 0);
                            neighHigherSide.X(ptneighcenter, xneighcenter);

                            sout << "\ngel/side =  " << gelside.Id() << "/" << gelside.Side();
                            sout << "; gel center coordinates: [ " << xcenter[0] << ", " << xcenter[1] << ", "
                                 << xcenter[2] << "]\n";
                            sout << "neigh/side =  " << neighGelSide.Id() << "/" << neighGelSide.Side();
                            sout << "; neigh center coordinates: [ " << xneighcenter[0] << ", " << xneighcenter[1]
                                 << ", "
                                 << xneighcenter[2] << "]\n";
                            sout << "Side solution =  " << elsol << "\n";
                            sout << "Neigh solution = " << c1sol << "\n";
                            sout << "Diff = " << elsol - c1sol << "\n";
                            sout << "Side coord:  [" << x[0] << ", " << x[1] << ", " << x[2] << "]\n";
                            sout << "Neigh coord: [" << xneigh[0] << ", " << xneigh[1] << ", " << xneigh[2] << "]\n";

                            c1.Print(*pressuremesh, sout);

                            //std::cout << sout.str(); // TODO remove
                            LOGPZ_DEBUG(logger, sout.str())
                        }
                    }
#endif
                }
            }
            if (gelside.Dimension() == 1) {
                TPZStack<TPZCompElSide> celstack;
                gelside.EqualLevelCompElementList(celstack, 1, 0);

                TPZCompElSide large = gelside.LowerLevelCompElementList2(1);
                if (large) celstack.Push(large);

                if (celstack.size() == 0) continue;

                int intOrder = 3;

                TPZIntPoints *intRule = gelside.CreateIntegrationRule(intOrder);

                // Iterates through the comp sides connected to the reference gelside
                int nstack = celstack.size();
                for (int ist = 0; ist < nstack; ist++) {
                    TPZCompElSide cneighbour = celstack[ist];
                    if (!cneighbour) continue;
                    TPZGeoElSide neighbour = cneighbour.Reference();
                    if (!(neighbour.Element()) || neighbour.Element()->MaterialId() != fHybridH1EE->fPressureSkeletonMatId)
                        DebugStop();

                    // Filters comp sides in elements of highest dimension (2 or 3)
                    if (neighbour.Element()->Dimension() > 1) DebugStop(); // 3D not supported
                    if (neighbour.Element()->Dimension() == 0) DebugStop();

                    // Verifies if coordinates on neighbours are the same
                    TPZTransform<REAL> transform(gelside.Dimension());
                    gelside.SideTransform3(neighbour, transform);

                    TPZManVector<REAL> pt0(gelside.Dimension(), 0);
                    TPZManVector<REAL> pt1(neighbour.Dimension(), 0);

                    int npoints = intRule->NPoints();
                    for (int ipt = 0; ipt < npoints; ipt++) {
                        REAL weight;
                        // Gets point in side parametric space from integration rule
                        intRule->Point(ipt, pt0, weight);
                        // Gets point in neighbour parametric space
                        transform.Apply(pt0, pt1);

                        // Transform from parametric to global coordinates
                        TPZManVector<REAL> x0(3);
                        TPZManVector<REAL> x1(3);

                        gelside.X(pt0, x0);
                        neighbour.X(pt1, x1);

                        TPZManVector<STATE> sol0(1), sol1(1);
                        cel->Solution(pt0, 0, sol0);
                        cneighbour.Element()->Solution(pt1, 0, sol1);

#ifdef LOG4CXX
                        if (logger->isDebugEnabled()) {
                            if (!IsZero(sol1[0] - sol0[0])) {
                                std::stringstream sout;
                                sout << "\ngel/side =  " << gelside.Id() << "/" << gelside.Side() << "\n";
                                sout << "neigh/side =  " << neighbour.Id() << "/" << neighbour.Side() << "\n";
                                sout << "Side solution =  " << sol0[0] << "\n";
                                sout << "Neigh solution = " << sol1[0] << "\n";
                                sout << "Diff = " << sol1[0] - sol0[0] << "\n";
                                sout << "Side coord:  [" << x0[0] << ", " << x0[1] << ", " << x0[2] << "]\n";
                                sout << "Neigh coord: [" << x1[0] << ", " << x1[1] << ", " << x1[2] << "]\n";

                                //std::cout << sout.str(); // TODO remove
                                LOGPZ_DEBUG(logger, sout.str())
                            }
                        }
#endif
                    }
                }
                delete intRule;
            }
        }
    }
}

/// computing the element stifnesses will "automatically" compute the condensed form of the matrices
void TPZHybridH1CreateH1Reconstruction::ComputeElementStiffnesses() {
    std::cout << "Solving local Dirichlet problem " << std::endl;
#ifdef ERRORESTIMATION_DEBUG2

    {
        std::ofstream out("MeshToComputeStiff.txt");
        fMultiphysicsH1reconstructionMesh->Print(out);
    }
#endif
    for (auto cel: fMultiphysicsH1reconstructionMesh->ElementVec()) {
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

void TPZHybridH1CreateH1Reconstruction::PlotLagrangeMultiplier(const std::string &filename, bool reconstructed) {

    TPZCompMesh *pressure = nullptr;

    if (!reconstructed) {
        pressure = fOriginal->MeshVector()[1];
    } else {
        pressure = fPressureMesh;
    }

    std::ofstream out2("PressuretoStateGraph.txt");
    pressure->Print(out2);

    {
        TPZLinearAnalysis an(pressure, false);
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("State");
        scalnames.Push("Pressure");

        int dim = pressure->Reference()->Dimension() - 1;
        std::string plotname;
        {
            std::stringstream out;
            out << filename << ".vtk";
            plotname = out.str();
        }
        std::set<int> matids ={fHybridH1EE->fPressureSkeletonMatId};
        an.DefineGraphMesh(dim,matids, scalnames, vecnames, plotname);
        an.PostProcess(2, dim);
    }
}

void TPZHybridH1CreateH1Reconstruction::GetPressureMatIDs(std::set<int> &matIDs){
    TPZCompMesh *fPressureMesh = fMultiphysicsH1reconstructionMesh->MeshVector()[1];
    int elMatID;
    for(int iel = 0 ; iel < fPressureMesh->NElements() ; iel++){
        TPZCompEl *cel = fPressureMesh->Element(iel);
        if(!cel) continue;
        elMatID = cel->Reference()->MaterialId();
        if (matIDs.find(elMatID) == matIDs.end()){
            matIDs.insert(elMatID);
        }
    }
}

/// compute the average pressures of across edges of the H(div) mesh
void TPZHybridH1CreateH1Reconstruction::ComputeAveragePressures(int target_dim) {

    TPZCompMesh *pressure_mesh = fPressureMesh;
    TPZGeoMesh *gmesh = pressure_mesh->Reference();



    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    int64_t nel = pressure_mesh->NElements();
    // load the pressure elements of dimension target_dim+1
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() != target_dim + 1 && gel->MaterialId() != fHybridH1EE->fPressureSkeletonMatId) continue;
        gel->SetReference(cel);
    }

    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        if(!cel) {
            continue;
        }

        TPZGeoEl *gel = cel->Reference();

        if (!gel || gel->Dimension() != target_dim) {
            continue;
        }

        // Skip calculation if the element is a boundary condition
        int matid = gel->MaterialId();
        TPZMaterial *mat = pressure_mesh->FindMaterial(matid);
        // TODO change this. Look for matIDs in bcMatIds instead. Only cast in debug mode for further checks

        TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
        if (bc) continue;

        // Skip calculation if the element is a small skeleton
        bool largeSideExists = false;
        if (cel->Connect(0).HasDependency()) largeSideExists = true;

#ifdef ERRORESTIMATION_DEBUG
        int nsides = gel->NSides();
        TPZGeoElSide side(gel, nsides - 1);
        TPZGeoElSideAncestors ancestors(side);
        TPZGeoElSide largerNeigh = ancestors.HasLarger(fPressureSkeletonMatId);
        if (largeSideExists && !largerNeigh) DebugStop();
#endif
        if (largeSideExists) continue;

        ComputeAverage(pressure_mesh, el);
    } 
#ifdef ERRORESTIMATION_DEBUG 
        std::ofstream outCon(fPressureReconstructionFolderOutput + "AverageB4LoadSolution.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(pressure_mesh, outCon, {1,2,3,fHybridH1EE->fPressureSkeletonMatId}, false, true);
#endif

    // Loads solution into the connects of the smaller skeletons
    pressure_mesh->LoadSolution(pressure_mesh->Solution());

#ifdef ERRORESTIMATION_DEBUG 
        std::ofstream outCon(fPressureReconstructionFolderOutput + "AverageAfterLoadSolution.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(pressure_mesh, outCon, {1,2,3,fHybridH1EE->fPressureSkeletonMatId}, false, true);
#endif

    // apply the restraints to the edge connects
    if (target_dim == dim - 2) {
        pressure_mesh->LoadSolution(pressure_mesh->Solution());
        TransferEdgeSolution();
    }
}

/// transfer the solution of the edge functions to the face functions
void TPZHybridH1CreateH1Reconstruction::TransferEdgeSolution() {
    // copy the solution associated with one-d edge connect to the corresponding side connect of the face mesh
    DebugStop();
    /*TPZCompMesh *pressure_mesh = fPressureMesh();
    TPZGeoMesh *gmesh = pressure_mesh->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    if (dim != 3) {
        std::cout << __PRETTY_FUNCTION__ << " should not be called for mesh dimension " << dim << std::endl;
        return;
    }
    int lagrangematid = fHybridizer.fLagrangeInterface;
    TPZMaterial *mat = pressure_mesh->FindMaterial(lagrangematid);
    if (!mat) DebugStop();
    int nstate = mat->NStateVariables();
    int64_t nel = pressure_mesh->NElements();
    // load the pressure elements of dimension 1 and 2
    pressure_mesh->Reference()->ResetReference();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() > 2) continue;
        gel->SetReference(cel);
    }

    // loop over the edge elements
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = cel->Reference();
        if (!gel) DebugStop();
        // if the dimension is not 1 continue
        if (gel->Dimension() != dim - 2) continue;
        int nsides = gel->NSides();
        for (int side = 0; side < nsides; side++) {
            // transfer the information of the internal connect only
            // this excludes the corner connects
            TPZConnect &edge_connect = intel->Connect(nsides - 1);
            /// transfer the solution of corner connects only if the are dependent
            // if the side is not the last (i.e. it has dimension 0) and doesnt have dependency continue
            // copying the solution of the constrained connects will ensure continuity of these sides
            if (side != nsides - 1 && !edge_connect.HasDependency()) {
                continue;
            }
            TPZGeoElSide gelside(gel, nsides - 1);
            TPZStack<TPZCompElSide> equal;
            int onlyinterpolated = 1;
            int removeduplicated = 0;
            // transfer the connect information to all connected elements
            gelside.EqualLevelCompElementList(equal, onlyinterpolated, removeduplicated);
            int nequal = equal.size();
            if (nequal == 0) {
                DebugStop();
            }
            int64_t edge_seqnum = edge_connect.SequenceNumber();
            int nshape_edge = pressure_mesh->Block().Size(edge_seqnum);
            for (int ieq = 0; ieq < nequal; ieq++) {
                TPZCompEl *celneigh = equal[ieq].Element();
                TPZInterpolatedElement *intelneigh = dynamic_cast<TPZInterpolatedElement *>(celneigh);
                TPZGeoEl *neighgel = intelneigh->Reference();
                if (!intelneigh) DebugStop();
                if (neighgel->MaterialId() != fHybridizer.fLagrangeInterface) {
                    DebugStop();
                }
                if (neighgel->Dimension() != 2) {
                    DebugStop();
                }
                int neighside = equal[ieq].Side();
                int nsideconnects = intelneigh->NSideConnects(neighside);
                TPZConnect &neigh_connect = intelneigh->SideConnect(nsideconnects - 1, neighside);
                int64_t neighblock = neigh_connect.SequenceNumber();
                int nshape_neigh = pressure_mesh->Block().Size(neighblock);
                if (nshape_edge != nshape_neigh) DebugStop();
                for (int i = 0; i < nshape_neigh; i++) {
                    pressure_mesh->Block()(neighblock, 0, i, 0) = pressure_mesh->Block()(edge_seqnum, 0, i, 0);
                }
            }
        }
    }*/
}


/// set the cornernode values equal to the averages
void TPZHybridH1CreateH1Reconstruction::ComputeNodalAverages() {
    TPZCompMesh *pressure_mesh = fPressureMesh;
    TPZGeoMesh *gmesh = pressure_mesh->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    pressure_mesh->LoadReferences();
    TPZMaterial *mat = pressure_mesh->FindMaterial(fHybridH1EE->fPressureSkeletonMatId);
    if (!mat) DebugStop();
    int nstate = mat->NStateVariables();
    int64_t nel = pressure_mesh->NElements();
    // load the pressure elements of dimension dim-1
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
    }
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();

        if (gel->Dimension() != dim - 1 || gel->MaterialId() != fHybridH1EE->fPressureSkeletonMatId) {
            // MatId of element not important for nodal average computation
            continue;
        }

        //percorre cada no do elemento de interface
        int ncorners = gel->NCornerNodes();
        for (int side = 0; side < ncorners; side++) {
            TPZGeoElSide gelside(gel, side);
            TPZCompElSide celside(intel,side);

            int64_t conindex = intel->ConnectIndex(side);
            TPZConnect &c = pressure_mesh->ConnectVec()[conindex];
            if (!c.HasDependency()) {
                if(LiesOnHangingSide(celside)) continue;
                ComputeNodalAverage(celside);
            }
        }
    }

    // Load the solution into connects with dependency
    pressure_mesh->LoadSolution(pressure_mesh->Solution());

    // Impose the solution on nodes which lies on hanging sides
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();

        if (gel->Dimension() != dim - 1 || gel->MaterialId() != fHybridH1EE->fPressureSkeletonMatId) {
            // MatId of element not important for nodal average computation
            continue;
        }

        //percorre cada no do elemento de interface
        int ncorners = gel->NCornerNodes();
        for (int side = 0; side < ncorners; side++) {
            TPZGeoElSide gelside(gel, side);
            TPZCompElSide celside(intel,side);

            int64_t conindex = intel->ConnectIndex(side);
            TPZConnect &c = pressure_mesh->ConnectVec()[conindex];
            if (c.HasDependency()) continue;
            if(LiesOnHangingSide(celside)){
                ImposeHangingNodeSolution(celside);
            }
        }
    }
}

// If all nodal (0-dimensional) connects linked to this node have dependency, it lies* on a hanging side and
// instead of calculating the average we should impose the solution of the other connects in it.
// *We are not sure if this logic will work for every case/ref. pattern, specially in 3D. So I'm leaving this
// disclaimer. (Gustavo Batistela, 14/10/2020)
bool TPZHybridH1CreateH1Reconstruction::LiesOnHangingSide(TPZCompElSide &node_celside){

    TPZCompMesh *pressure_mesh = fPressureMesh;
    int dim = pressure_mesh->Dimension();
    TPZGeoElSide node_gelside(node_celside.Reference());

    // celstack will contain all zero dimensional sides connected to the side
    TPZStack<TPZCompElSide> celstack;
    int onlyinterpolated = 1;
    int removeduplicates = 0;

    node_gelside.ConnectedCompElementList(celstack, onlyinterpolated, removeduplicates);

    bool allNodeConnectsHaveDependency = true;
    for (int elc = 0; elc < celstack.size(); elc++) {
        TPZCompElSide celside = celstack[elc];
        TPZGeoElSide gelside = celside.Reference();
        if (gelside.Element()->Dimension() != dim - 1) {
            continue;
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(celside.Element());
        if (!intel) DebugStop();

        int64_t conindex = intel->ConnectIndex(celside.Side());
        TPZConnect &c = pressure_mesh->ConnectVec()[conindex];
        if (gelside.Dimension() == 0 && !c.HasDependency()) allNodeConnectsHaveDependency = false;
    }
    return allNodeConnectsHaveDependency;
}

void TPZHybridH1CreateH1Reconstruction::ImposeHangingNodeSolution(TPZCompElSide &node_celside){

    int skeletonMatId = fHybridH1EE->fPressureSkeletonMatId;
    TPZCompMesh *pressure_mesh = fPressureMesh;
    int dim = pressure_mesh->Dimension();

    TPZMaterial *mat = pressure_mesh->FindMaterial(skeletonMatId);
    if (!mat) DebugStop();
    int nstate = mat->NStateVariables();

    TPZGeoElSide node_gelside(node_celside.Reference());
    TPZGeoEl *gel = node_gelside.Element();
    int side = node_gelside.Side();

    // celstack will contain all zero dimensional sides connected to the side
    TPZStack<TPZCompElSide> celstack;
    int onlyinterpolated = 1;
    int removeduplicates = 0;

    node_gelside.ConnectedCompElementList(celstack, onlyinterpolated, removeduplicates);

    // Impose solution of the first node connect in the celstack we can find

    celstack.Push(node_celside);

    TPZBlock &block =  pressure_mesh->Block();
    TPZFMatrix<STATE> &sol = pressure_mesh->Solution();
    for (int elc = 0; elc < celstack.size(); elc++) {
        TPZCompElSide celside = celstack[elc];
        if (celside.Reference().Dimension() != 0) continue;

        // Get solution of the neighbour
        TPZInterpolatedElement *neigh_intel = dynamic_cast<TPZInterpolatedElement *>(celside.Element());
        if (!neigh_intel) DebugStop();

        int64_t neigh_conindex = neigh_intel->ConnectIndex(celside.Side());
        TPZConnect &neigh_c = pressure_mesh->ConnectVec()[neigh_conindex];

        int64_t neigh_seqnum = neigh_c.SequenceNumber();
        if (neigh_c.NState() != nstate || neigh_c.NShape() != 1) DebugStop();
        TPZManVector<STATE, 3> neigh_sol(nstate, 0.);
        for (int istate = 0; istate < nstate; istate++) {
            neigh_sol[istate] = sol.at(block.at(neigh_seqnum, 0, istate, 0));
        }

        // Set solution to given connect
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(node_celside.Element());
        if (!intel) DebugStop();

        int64_t conindex = intel->ConnectIndex(side);
        TPZConnect &c = pressure_mesh->ConnectVec()[conindex];

        int64_t seqnum = c.SequenceNumber();
        if (c.NState() != nstate || c.NShape() != 1) DebugStop();
        for (int istate = 0; istate < nstate; istate++) {
            sol.at(block.at(seqnum, 0, istate, 0)) = neigh_sol[istate];
        }

        // TODO remove this
        {
            TPZManVector<REAL, 3> xicenter(gel->Dimension(), 0.);
            TPZManVector<REAL, 3> x_side(3, 0.);
            gel->CenterPoint(side, xicenter);
            gel->X(xicenter, x_side);
            //std::cout << __PRETTY_FUNCTION__ << "\nAll " << celstack.size() << " Connects Have Dependency!\nGel "
            //         << gel->Index() << ", Side " << side << " @ [" << x_side[0] << ", " << x_side[1] << ", "
            //         << x_side[2] << "]\n";
        }
        break;
    }
    return;
}

/// compute the nodal average of all elements that share a point
void TPZHybridH1CreateH1Reconstruction::ComputeNodalAverage(TPZCompElSide &node_celside)
{
    int skeletonMatId = fHybridH1EE->fPressureSkeletonMatId;
    TPZCompMesh *pressure_mesh = fPressureMesh;
    int dim = pressure_mesh->Dimension();

    TPZMaterial *mat = pressure_mesh->FindMaterial(skeletonMatId);
    if (!mat) DebugStop();
    int nstate = mat->NStateVariables();

    TPZGeoElSide node_gelside(node_celside.Reference());
    TPZGeoEl *gel = node_gelside.Element();
    int side = node_gelside.Side();

    // celstack will contain all zero dimensional sides connected to the side
    TPZStack<TPZCompElSide> celstack;
    int onlyinterpolated = 1;
    int removeduplicates = 0;

    node_gelside.ConnectedCompElementList(celstack, onlyinterpolated, removeduplicates);

    // Impose solution of the first node connect in the celstack we can find

    celstack.Push(node_celside);

    // This map stores the connects, the weight associated with the element
    // and the solution of that connect. The weight of Dirichlet condition is
    // higher and will be used later to impose the value of the BC in the
    // connects when needed
    std::map<int64_t, std::pair<REAL, TPZVec<STATE>>> connects;
    TPZBlock &block =  pressure_mesh->Block();
    TPZFMatrix<STATE> &solMatrix = pressure_mesh->Solution();
    for (int elc = 0; elc < celstack.size(); elc++) {
        TPZCompElSide celside = celstack[elc];
        TPZGeoElSide gelside = celside.Reference();
        if (gelside.Element()->Dimension() != dim - 1) {
            continue;
        }

        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(celside.Element());
        if (!intel) DebugStop();
        int64_t index = intel->Index();
        REAL weight = fPressureweights[index];

        if (IsZero(weight)) {
            TPZMaterial *mat = intel->Material();
            TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
            if (!bc) DebugStop();

            continue;
        }

        int64_t conindex = intel->ConnectIndex(celside.Side());
        if (connects.find(conindex) != connects.end()) {
            DebugStop();
        }

        TPZConnect &c = intel->Connect(celside.Side());
        int64_t seqnum = c.SequenceNumber();
        if (c.NState() != nstate || c.NShape() != 1) DebugStop();
        TPZManVector<REAL,3> pt(0), x(3);
        TPZManVector<STATE, 3> sol(nstate, 0.);
        for (int istate = 0; istate < nstate; istate++) {
            sol[istate] = solMatrix.at(block.at(seqnum, 0, istate, 0));
        }
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "conindex " << conindex << " weight " << weight << " state " << sol;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        connects.insert({conindex, {weight, sol}});
    }

    TPZManVector<STATE, 3> averageSol(nstate, 0);
    REAL sum_weight = 0.;
    for (auto it = connects.begin(); it != connects.end(); it++) {
        REAL weight = it->second.first;
        sum_weight += weight;
        for (int istate = 0; istate < nstate; istate++) {
            averageSol[istate] += it->second.second[istate] * weight;
        }
    }

    for (int istate = 0; istate < nstate; istate++) {
        averageSol[istate] /= sum_weight;
    }

#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Comparing connect indexes ";
        for (auto it : connects) {
            sout << it.first << " ";
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    for (auto it = connects.begin(); it != connects.end(); it++) {
        int64_t conindex = it->first;
        TPZConnect &c = pressure_mesh->ConnectVec()[conindex];
        int64_t seqnum = c.SequenceNumber();

        if (c.NState() != nstate || c.NShape() != 1) DebugStop();
        for (int istate = 0; istate < nstate; istate++) {
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "value before " << pressure_mesh->Block()(seqnum, 0, istate, 0) <<
                    " value after " << averageSol[istate] << " diff "
                     << pressure_mesh->Block()(seqnum, 0, istate, 0) - averageSol[istate];
                //                        res2.Print("Residual",sout);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            solMatrix.at(block.at(seqnum, 0, istate, 0))= averageSol[istate];
        }
    }
}

// compute the average of an element iel in the pressure mesh looking at its neighbours
void TPZHybridH1CreateH1Reconstruction::ComputeAverage(TPZCompMesh *pressuremesh, int64_t iel)
{
    TPZGeoMesh *gmesh = fPressureMesh->Reference();
    int dim = gmesh->Dimension();
    gmesh->ResetReference();
    fPressureMesh->LoadReferences();
    TPZCompEl *cel = fPressureMesh->Element(iel);
    if (!cel) DebugStop();
    TPZGeoEl *gel = cel->Reference();
    if (!gel) DebugStop();
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
    if (!intel) DebugStop();
    int nc = cel->NConnects();
    int order = cel->Connect(nc - 1).Order();

#ifdef ERRORESTIMATION_DEBUG
    for (int ic = 0; ic < nc; ic++) {
        if (cel->Connect(ic).HasDependency()) DebugStop();
    }
#endif

    int target_dim = gel->Dimension();
    if (target_dim == dim - 1 && gel->MaterialId() != fHybridH1EE->fPressureSkeletonMatId) {
        DebugStop();
    }

    // We denote large the skeleton element passed to this method, which is whether the skeleton between two elements
    // of the same level or the larger skeleton in the region of a hanging node.
    TPZGeoElSide largeSkeletonSide(gel, gel->NSides() - 1);

    // Stack of TPZCompElSides containing the volumetric neighbours from which we obtain the pressure solution to
    // be used in the average computation.
    TPZStack<TPZCompElSide> volumeNeighSides;
    largeSkeletonSide.EqualLevelCompElementList(volumeNeighSides, 1, 0);
    if (volumeNeighSides.size() < 1) DebugStop();

    // If we are on a hanging side, we need to store the smaller skeleton sides and their respective volumetric
    // neighbours
    TPZStack<TPZGeoElSide> smallerSkelSides;
    bool noHangingSide = true;
    if (volumeNeighSides.size() == 1) {
        noHangingSide = false;
        TPZGeoElSidePartition partition(largeSkeletonSide);
        if (!partition.HigherLevelNeighbours(smallerSkelSides, fHybridH1EE->fPressureSkeletonMatId)) {
            DebugStop();
        }
        if (smallerSkelSides.size() == 1) DebugStop();
        for (int iskel = 0; iskel < smallerSkelSides.size(); iskel++) {
            smallerSkelSides[iskel].EqualLevelCompElementList(volumeNeighSides, 1, 0);
        }
    }

#ifdef ERRORESTIMATION_DEBUG
    for (int ineigh = 0; ineigh < volumeNeighSides.size(); ineigh++) {
        if (volumeNeighSides[ineigh].Element()->Dimension() != dim) DebugStop();
    }
#endif

    // Calculate average weights for left/right elements
    REAL sum_weights = 0.;
    TPZGeoEl *leftVolumeGel = volumeNeighSides[0].Element()->Reference();
    int matId = leftVolumeGel->MaterialId();
    sum_weights += fMatid_weights[matId];
    TPZGeoEl *rightVolumeGel = volumeNeighSides[1].Element()->Reference();
    matId = rightVolumeGel->MaterialId();
    sum_weights += fMatid_weights[matId];
    fPressureweights[cel->Index()] = sum_weights / 2;
    if (!noHangingSide) {
        for (int i = 0; i < smallerSkelSides.size(); i++) {
            TPZGeoElSide smallSkelGelside = smallerSkelSides[i];
            TPZCompEl* smallSkelCel = smallerSkelSides[i].Element()->Reference();
            fPressureweights[smallSkelCel->Index()] = sum_weights / 2;
        }
    }

    // TODO refactor this part so the transformations aren't stored but rather calculated inside the main loop at line 1480

    // Fills vectors of TPZTransform with the transformation from the skeleton to each left/right volumetric neighbours
    // and from the the small (right) to the large (left) skel if we are on a hanging side.
    int64_t nSkelsToIntegrate = volumeNeighSides.size() - 1;
    TPZManVector<TPZTransform<REAL>, 5> leftSkelToVolumeTrans(nSkelsToIntegrate);
    TPZManVector<TPZTransform<REAL>, 5> rightSkelToVolumeTrans(nSkelsToIntegrate);
    TPZManVector<TPZTransform<REAL>, 5> smallSkelToLargeSkelTrans(nSkelsToIntegrate);

    // Transformation between the large skeleton to the target_dim dimensional side of the volume element
    TPZTransform<REAL> largeSkelToVolumeTrans =
        largeSkeletonSide.NeighbourSideTransform(volumeNeighSides[0].Reference());
    TPZGeoEl *volumeNeigh = volumeNeighSides[0].Element()->Reference();
    // Transformation of the target_dim dimensional  side of the volume element to its higher dimension side
    TPZTransform<REAL> tmp = volumeNeigh->SideToSideTransform(volumeNeighSides[0].Side(), volumeNeigh->NSides() - 1);
    largeSkelToVolumeTrans = tmp.Multiply(largeSkelToVolumeTrans);

    if (noHangingSide) {
        // In this case the largeSkeletonSide is the same for both left/right regions.
        // We denote left the side at volumeNeighSides[0] and the remaining right.
        leftSkelToVolumeTrans[0] = largeSkelToVolumeTrans;

        rightSkelToVolumeTrans[0] = largeSkeletonSide.NeighbourSideTransform(volumeNeighSides[1].Reference());
        volumeNeigh = volumeNeighSides[1].Element()->Reference();
        tmp = volumeNeigh->SideToSideTransform(volumeNeighSides[1].Side(), volumeNeigh->NSides() - 1);
        rightSkelToVolumeTrans[0] = tmp.Multiply(rightSkelToVolumeTrans[0]);
    } else {
        // If we are on a hanging side we need to build left/right transformations for every small skeleton
        for (int iskel = 0; iskel < smallerSkelSides.size(); iskel++) {
            rightSkelToVolumeTrans[iskel] =
                smallerSkelSides[iskel].NeighbourSideTransform(volumeNeighSides[iskel + 1].Reference());
            volumeNeigh = volumeNeighSides[iskel + 1].Element()->Reference();
            tmp = volumeNeigh->SideToSideTransform(volumeNeighSides[iskel + 1].Side(), volumeNeigh->NSides() - 1);
            rightSkelToVolumeTrans[iskel] = tmp.Multiply(rightSkelToVolumeTrans[iskel]);

            smallSkelToLargeSkelTrans[iskel] = TPZTransform<REAL>(target_dim);
            smallerSkelSides[iskel].SideTransform3(largeSkeletonSide, smallSkelToLargeSkelTrans[iskel]);
            leftSkelToVolumeTrans[iskel] = largeSkelToVolumeTrans.Multiply(smallSkelToLargeSkelTrans[iskel]);
        }
    }

    int nshape = intel->NShapeF();
    TPZFNMatrix<20, REAL> L2Mat(nshape, nshape, 0.), L2Rhs(nshape, 1, 0.);
    TPZFNMatrix<220, REAL> phi(nshape, 1, 0.), dphi(dim, nshape);
    for (int iskel = 0; iskel < nSkelsToIntegrate; iskel++) {
        TPZGeoElSide integrationGeoElSide;
        if (noHangingSide) {
            integrationGeoElSide = largeSkeletonSide;
        } else {
            integrationGeoElSide = smallerSkelSides[iskel];
        }
        std::unique_ptr<TPZIntPoints> intpoints(gel->CreateSideIntegrationRule(integrationGeoElSide.Side(), 2 * order));

        REAL left_weight = fMatid_weights[leftVolumeGel->MaterialId()] / sum_weights;
        rightVolumeGel = volumeNeighSides[iskel + 1].Element()->Reference();
        REAL right_weight = fMatid_weights[rightVolumeGel->MaterialId()] / sum_weights;

        if(!noHangingSide){
            AdditionalAverageWeights(largeSkeletonSide.Element(), integrationGeoElSide.Element(), left_weight, right_weight,sum_weights);
        }

        int64_t nintpoints = intpoints->NPoints();
        for (int64_t ip = 0; ip < nintpoints; ip++) {
            TPZManVector<REAL, 3> pt_left_skel(target_dim, 0.), pt_right_skel(target_dim, 0.);
            TPZManVector<REAL, 3> pt_left_vol(target_dim + 1, 0.), pt_right_vol(target_dim + 1, 0.);

            REAL weight;
            intpoints->Point(ip, pt_right_skel, weight);
            if (noHangingSide) {
                pt_left_skel = pt_right_skel;
            } else {
                smallSkelToLargeSkelTrans[iskel].Apply(pt_right_skel, pt_left_skel);
            }

            // Get shape at integration point
            intel->Shape(pt_left_skel, phi, dphi);

            // Get solution from left/right sides
            leftSkelToVolumeTrans[iskel].Apply(pt_right_skel, pt_left_vol);
            rightSkelToVolumeTrans[iskel].Apply(pt_right_skel, pt_right_vol);
            TPZVec<STATE> left_sol, right_sol;
            volumeNeighSides[0].Element()->Solution(pt_left_vol, 0, left_sol);
            volumeNeighSides[iskel + 1].Element()->Solution(pt_right_vol, 0, right_sol);

            STATE average_sol = left_weight * left_sol[0] + right_weight * right_sol[0];

            TPZFNMatrix<9, REAL> jac(dim, dim), jacinv(dim, dim), axes(dim, 3);
            REAL detjac;
            integrationGeoElSide.Jacobian(pt_right_skel, jac, axes, detjac, jacinv);

#ifdef LOG4CXX
            if(logger->isDebugEnabled()) {
                std::stringstream ss;
                phi.Print(ss);
                ss << detjac << '\n';
                LOGPZ_DEBUG(logger, ss.str())
            }
#endif
            for (int ishape = 0; ishape < nshape; ishape++) {
                L2Rhs(ishape, 0) += weight * phi(ishape, 0) * detjac * average_sol;
            }
            for (int ishape = 0; ishape < nshape; ishape++) {
                for (int jshape = 0; jshape < nshape; jshape++) {
                    L2Mat(ishape, jshape) += weight * detjac * phi(ishape, 0) * phi(jshape, 0);
                }
            }
        }
    }

#ifdef LOG4CXX
    if(logger->isDebugEnabled()) {
        std::stringstream ss;
        L2Rhs.Print("Rhs =", ss, EMathematicaInput);
        L2Mat.Print("Stiffness =", ss, EMathematicaInput);
        LOGPZ_DEBUG(logger, ss.str())
    }
#endif

    L2Mat.SolveDirect(L2Rhs, ECholesky);
    // Stores solution in the computational mesh

    TPZFMatrix<STATE> &sol = fPressureMesh->Solution();
    int count = 0;
    for (int ic = 0; ic < nc; ic++) {
        TPZConnect &c = cel->Connect(ic);
        int64_t seqnum = c.SequenceNumber();
        int64_t pos = fPressureMesh->Block().Position(seqnum);
        int ndof = c.NShape() * c.NState();
        for (int idf = 0; idf < ndof; idf++) {
            sol(pos + idf, 0) = L2Rhs(count++);
        }
    }
}

void TPZHybridH1CreateH1Reconstruction::AdditionalAverageWeights(TPZGeoEl *large, TPZGeoEl *small, REAL &large_weight, REAL &small_weight, REAL sum){
    REAL diff_weight = large_weight+small_weight - 1.;
    if(!IsZero(diff_weight)) DebugStop();
    switch(fAverageMode){
    case 0:
        return;
    case 1:
        large_weight = 0.;
        small_weight = 1.;
        break;
    case 2: // w ~ 1/h^p
        REAL large_size, small_size;
        int large_order,small_order;

        //Compute h
        large_size = large->CharacteristicSize();
        small_size = small->CharacteristicSize();

        //Compute k
        large_order = large->Reference()->Connect(large->NSides()-1).Order();
        small_order = small->Reference()->Connect(small->NSides()-1).Order();

        //Compute 1/h^(k)
        REAL large_precision = 1.;
        REAL small_precision = 1.;
        for(int pindex = 0; pindex < large_order; pindex++){
            large_precision *= large_size;
        }
        for(int pindex = 0; pindex < small_order; pindex++){
            small_precision *= small_size;
        }
        large_precision = 1./large_precision;
        small_precision = 1./small_precision;

        REAL new_large_weight = large_weight*large_precision;
        REAL new_small_weight = small_weight*small_precision;
        REAL new_sum = new_large_weight+new_small_weight;

        //std::cout << "large0: " << large_weight << "; large_h: " << large_size << "; large_p: " <<  large_order << "; large_precision: " << large_precision << "; large1: " << new_large_weight/new_sum <<"\n";
        //std::cout << "small0: " << small_weight << "; small_h: " << small_size << "; small_p: " <<  small_order << "; small_precision: " << small_precision << "; small1: " << new_small_weight/new_sum <<"\n";

        large_weight = new_large_weight/new_sum;
        small_weight = new_small_weight/new_sum;

        if(!IsZero(large_weight + small_weight - 1.)) DebugStop();
    }
}

/// compute the pressure weights and material weights
//  fills in the data structure pressureweights and matid_weights
void TPZHybridH1CreateH1Reconstruction::ComputePressureWeights() {

    int dim = fPressureMesh->Dimension();
    int64_t nel = fPressureMesh->NElements();
    fPressureweights.Resize(nel, 0);
    fMatid_weights.clear();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fPressureMesh->Element(el);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        TPZMaterial *mat = fOriginal->FindMaterial(matid);
        // TODO: Case HybridSquared: insert if( ... || matid == fLagrangeInterface)
        if (matid == fHybridH1EE->fPressureSkeletonMatId ) {
            continue;
        }
        if (!mat)  DebugStop();

        TPZBndCondT<STATE> *bcmat = dynamic_cast<TPZBndCondT<STATE> *>(mat);

        if (bcmat) {
            if (bcmat->Type() == 0) {
                fPressureweights[el] = 1.e12;
                fMatid_weights[matid] = 1.e12;
                continue;
            }
            else if (bcmat->Type() == 4){
                fPressureweights[el] = bcmat->Val1()(0,0);//valor de Km
                fMatid_weights[matid] = bcmat->Val1()(0,0);

            }
            // Neumann
            else {
                fPressureweights[el] = 0.;
                fMatid_weights[matid] = 0.;
                continue;
            }
        }

        else{
            TPZMatLaplacianHybrid *matlaplacian = dynamic_cast<TPZMatLaplacianHybrid *>(mat);
            if (!matlaplacian) DebugStop();

            REAL perm;
            perm = matlaplacian->GetPermeability({0,0,0});
            if (IsZero(perm)) DebugStop();
            fPressureweights[el] = perm;
            fMatid_weights[matid] = perm;
        }
    }
}

// Update the mesh solution at dirichlet boundary locations with exact values
void TPZHybridH1CreateH1Reconstruction::ComputeBoundaryL2Projection(TPZCompMesh *pressuremesh, int target_dim){
    {
        std::ofstream out("PressureBeforeL2Projection.txt");
        fPressureMesh->Print(out);
    }

    TPZGeoMesh *gmesh = fPressureMesh->Reference();
    gmesh->ResetReference();
    int64_t nel = fPressureMesh->NElements();

    TPZAdmChunkVector<TPZCompEl *> &elementvec = fPressureMesh->ElementVec();

    TPZElementMatrixT<STATE> ekbc, efbc;
    for (int iel = 0; iel < nel; iel++) {
        TPZCompEl *cel = elementvec[iel];
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();

        int matid = gel->MaterialId();
        TPZMaterial *mat = fPressureMesh->FindMaterial(matid);

        TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
        if (!bc || (bc->Type() != 0)) continue;

        cel->CalcStiff(ekbc, efbc);

        ekbc.fMat.SolveDirect(efbc.fMat, ELU);

#ifdef ERRORESTIMATION_DEBUG43
        if (bc->HasForcingFunctionBC()) {
            TPZManVector<STATE> res(3),resex(3);
            int dim=2;
            TPZFNMatrix<9, STATE> gradu(dim, 1);
            auto bcstate = dynamic_cast<TPZBndCondT<STATE>*>(bc);
            if(!bcstate)
                DebugStop();

            TPZMaterialDataT<STATE> data;
            auto intel = dynamic_cast<TPZInterpolationSpace*> (cel);
            intel->InitMaterialData(data);

            int dimbc = intel->Dimension();
            TPZManVector<REAL,3> intpoint(dimbc,0.);
            intpoint[0] =-3;
            for(int inode =0; inode <2;inode++){
                intpoint[0] +=2;

                intel->ComputeRequiredData(data, intpoint);
                data.phi.Print(std::cout);
                efbc.fMat.Print(std::cout);
                double sol =0;
                for(int i=0; i < data.phi.Rows(); i++){
                    sol += data.phi[i]*efbc.fMat(i);
                }
                std::cout << "x[0]:\t"<<data.x[0]<<"\nx[1]:\t"<<data.x[1]<<"\n";
                std::cout << "solution:\t"<< sol <<"\n";

                bcstate->ForcingFunctionBC()(data.x, res, gradu);
                fHybridH1EE->fProblemConfig.exact->ExactSolution()(data.x, resex, gradu);
                std::cout << "bc value: " << res[0] <<"\n";
                std::cout << "exact value: " << resex[0] <<"\n";
            }
        }
#endif


        int count = 0;
        int nc = cel->NConnects();
        TPZFMatrix<STATE> &mesh_sol = fPressureMesh->Solution();
        for (int ic = 0; ic < nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            int64_t seqnum = c.SequenceNumber();
            int64_t pos = fPressureMesh->Block().Position(seqnum);
            int ndof = c.NShape() * c.NState();
            for (int idf = 0; idf < ndof; idf++) {
                mesh_sol(pos + idf, 0) = efbc.fMat(count++);
            }
        }
    }

    {
        std::ofstream out("PressureAfterL2Projection.txt");
        fPressureMesh->Print(out);
    }
}

void TPZHybridH1CreateH1Reconstruction::CopySolutionFromSkeleton() {

    fPressureMesh->Reference()->ResetReference();
    fPressureMesh->LoadReferences();
    int dim = fPressureMesh->Dimension();
    int64_t nel = fPressureMesh->NElements();

    TPZBlock &block = fPressureMesh->Block();
    TPZFMatrix<STATE> &sol = fPressureMesh->Solution();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fPressureMesh->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() != dim) continue;
        int nsides = gel->NSides();

        for (int is = 0; is < nsides; is++) {
            //
            TPZGeoElSide gelside(gel, is);
            TPZConnect &c = intel->Connect(is);
            int64_t c_seqnum = c.SequenceNumber();
            int c_blocksize = c.NShape() * c.NState();
            //TPZGeoElSide gelside(gel,is);
            TPZStack<TPZCompElSide> celstack;
            gelside.EqualLevelCompElementList(celstack, 1, 0);
            //std::cout << "celstack size: "<< celstack.size() <<"\n";
            int nst = celstack.NElements();
            for (int ist = 0; ist < nst; ist++) {
                TPZCompElSide cneigh = celstack[ist];
                TPZGeoElSide gneigh = cneigh.Reference();
                //std::cout << "Gel/side: " << intel->Index() << "/" << gelside.Side() << "\nNeighbour Gel/side: " <<cneigh.Element()->Index() << "/" << cneigh.Side() << "\n\n";
                if ((gneigh.Element()->MaterialId() == fHybridH1EE->fPressureSkeletonMatId)||
                    IsDirichletCondition(gneigh)) {
                    TPZInterpolatedElement *intelneigh = dynamic_cast<TPZInterpolatedElement *>(cneigh.Element());
                    if (!intelneigh) DebugStop();
                    TPZConnect &con_neigh = intelneigh->Connect(cneigh.Side());
                    int64_t con_seqnum = con_neigh.SequenceNumber();
                    int con_size = con_neigh.NState() * con_neigh.NShape();
                    if (con_size != c_blocksize) DebugStop();
                    for (int ibl = 0; ibl < con_size; ibl++) {
                        sol.at(block.at(c_seqnum,0,ibl,0)) = sol.at(block.at(con_seqnum,0,ibl,0));
                    }
                    break;
                }
                // all elements must have at least one neighbour of type skeleton--> esta premissa nao vale para reconstrucao Hdiv-H1
                if (ist == nst - 1) {
#ifdef LOG4CXX
                    std::stringstream sout;
                    sout << "Connect " << is << " from element el " << el
                         << " was not updated \n";
                    LOGPZ_DEBUG(logger, sout.str())
#endif
                }
            }
        }
    }
}

/// returns true if the material associated with the element is a boundary condition
/// and if the boundary condition is dirichlet type
bool TPZHybridH1CreateH1Reconstruction::IsDirichletCondition(TPZGeoElSide gelside) {
    TPZGeoEl *gel = gelside.Element();
    int matid = gel->MaterialId();
    TPZMaterial *mat = fMultiphysicsH1reconstructionMesh->FindMaterial(matid);
    TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
    if (!bc) return false;
    int typ = bc->Type();
    //  std::cout<<"type "<< typ<<"\n";
    if ((typ == 0)||(typ == 4)) return true;//for copy the robin boundary too
    return false;
}

void TPZHybridH1CreateH1Reconstruction::BoundaryPressureProjection(TPZCompMesh *pressuremesh, int target_dim){
    // TODO remove unused target_dim variable
    //    std::ofstream out("PressureProjBefore.txt");
    //    pressuremesh->Print(out);

    TPZGeoMesh *gmesh = fMultiphysicsH1reconstructionMesh->Reference();
    std::map<int,TPZMaterial *> matvec = fMultiphysicsH1reconstructionMesh->MaterialVec();
    for (auto matit=matvec.begin(); matit != matvec.end(); matit++) {
        TPZHybridH1ErrorEstimateMaterial *castmat = dynamic_cast<TPZHybridH1ErrorEstimateMaterial *>(matit->second);
        if(castmat)
        {
            TPZPressureProjection *pressproj = new TPZPressureProjection(*castmat);
            fMultiphysicsH1reconstructionMesh->MaterialVec()[matit->first] = pressproj;
        }
    }
    for (auto matit=matvec.begin(); matit != matvec.end(); matit++) {
        TPZBndCond *bndcond = dynamic_cast<TPZBndCond *>(matit->second);
        if(bndcond)
        {
            TPZMaterial *refmat = bndcond->Material();
            int matid = refmat->Id();
            bndcond->SetMaterial(fMultiphysicsH1reconstructionMesh->FindMaterial(matid));
        }
    }

    gmesh->ResetReference();


    //    {
    //        TPZCompMesh *pressmesh = fMultiphysicsH1reconstructionMesh->MeshVector()[1];
    //        std::ofstream out("pressure.txt");
    //        pressmesh->Print(out);
    //    }

    TPZElementMatrixT<STATE> ekbc,efbc;

    int64_t nel = fMultiphysicsH1reconstructionMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel_orig = fMultiphysicsH1reconstructionMesh->Element(el);
        if(!cel_orig) continue;
        TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>(cel_orig);
        if(!cond) continue;
        TPZCompEl *condref = cond->ReferenceCompEl();
        TPZElementGroup *celgr = dynamic_cast<TPZElementGroup *>(condref);
        if(!celgr) continue;//DebugStop();
        TPZVec<TPZCompEl *> elgrST = celgr->GetElGroup();
        int nelst = elgrST.size();

        for (int elst = 0; elst<nelst; elst++)
        {
            TPZCompEl *cel = elgrST[elst];
            TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement*>(cel);
            if(mphys == 0) DebugStop();
            TPZMaterial *mat = cel->Material();
            TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
            if(!bnd) continue;
            cel->CalcStiff(ekbc, efbc);

            if(ekbc.HasDependency()) DebugStop();
            ekbc.fMat.SolveDirect(efbc.fMat, ECholesky);
            //   efbc.fMat.Print("Solution",std::cout);


            TPZCompEl *celpressure = mphys->Element(1);
            int nc = celpressure->NConnects();
            int count = 0;
            TPZFMatrix<STATE> &mesh_sol = fPressureMesh->Solution();
            for (int ic = 0; ic < nc; ic++) {
                TPZConnect &c = celpressure->Connect(ic);
                int64_t seqnum = c.SequenceNumber();
                int64_t pos = fPressureMesh->Block().Position(seqnum);
                int ndof = c.NShape() * c.NState();
                for (int idf = 0; idf < ndof; idf++) {
                    mesh_sol(pos + idf, 0) = efbc.fMat(count++);
                    //                    std::cout<<" connect seqnum "<< seqnum<<" sol "<<pressuremesh->Solution()(pos + idf, 0)<<"\n";
                }
            }


        }
    }

    for (auto matit=matvec.begin(); matit != matvec.end(); matit++) {
        TPZBndCond *bndcond = dynamic_cast<TPZBndCond *>(matit->second);
        if(bndcond)
        {
            TPZMaterial *refmat = bndcond->Material();
            int matid = refmat->Id();
            bndcond->SetMaterial(matvec[matid]);
        }
    }

    {
        std::ofstream out("PressureProjectionAfter.txt");
        fPressureMesh->Print(out);
    }

}



void TPZHybridH1CreateH1Reconstruction::NewComputeBoundaryL2Projection(
    TPZCompMesh *pressuremesh, int target_dim) {

    //{
    //    std::ofstream out("PressureBeforeL2Projection.txt");
    //    pressuremesh->Print(out);
    //}

    if (target_dim == 2) {
        std::cout << "Not implemented for 2D interface.\n";
        DebugStop();
    }

    TPZGeoMesh *gmesh = fPressureMesh->Reference();
    gmesh->ResetReference();
    fOriginal->LoadReferences();

    int64_t nel = fPressureMesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZCompEl *cel = fPressureMesh->ElementVec()[iel];
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        int dim = gel->Dimension();

        int matid = gel->MaterialId();
        TPZMaterial *mat = fPressureMesh->FindMaterial(matid);

        TPZBndCondT<STATE> *bc = dynamic_cast<TPZBndCondT<STATE> *>(mat);
        if (!bc) continue;

        if(bc->Type()==4 && IsZero(bc->Val1()(0,0))) continue;

        int nc = cel->NConnects();
        int order = cel->Connect(nc - 1).Order();
        TPZGeoElSide gelside(gel, gel->NSides() - 1);
        TPZAutoPointer<TPZIntPoints> intp =
            gel->CreateSideIntegrationRule(gel->NSides() - 1, 2 * order);

        TPZInterpolatedElement *intel =
            dynamic_cast<TPZInterpolatedElement *>(cel);
        int nshape = intel->NShapeF();

        TPZFNMatrix<20, REAL> ekbc(nshape, nshape, 0.), efbc(nshape, 1, 0.);
        TPZFNMatrix<220, REAL> phi(nshape, 1, 0.), dshape(dim, nshape);

        int64_t npoints = intp->NPoints();
        for (int64_t ip = 0; ip < npoints; ip++) {
            TPZManVector<REAL, 3> pt(dim, 0.);
            REAL weight;
            intp->Point(ip, pt, weight);
            intel->Shape(pt, phi, dshape);
            REAL u_D = 0.;
            REAL g = 0.;

            if (bc->HasForcingFunctionBC()) {
                TPZManVector<STATE> result(3);
                TPZFNMatrix<9, STATE> gradu(dim, 1);
                TPZManVector<REAL, 3> x;
                gel->X(pt, x);
                bc->ForcingFunctionBC()(x, result, gradu);
                u_D = result[0];
                //tem que calcular o g aqui..que é Kgradu.n

                TPZFNMatrix<9,REAL> PermTensor;
                TPZFNMatrix<9,REAL> InvPermTensor;



            }

            int bcType = bc->Type();
            switch (bcType) {
            case 0: {
                for (int ishape = 0; ishape < nshape; ishape++) {
                    efbc(ishape, 0) += weight * phi(ishape, 0) * u_D;
                    for (int jshape = 0; jshape < nshape; jshape++) {
                        ekbc(ishape, jshape) += weight * phi(ishape, 0) * phi(jshape, 0);
                    }
                }
                break;
            }
            case 4: {
                TPZCompEl *origCel = gel->Reference();
                TPZMultiphysicsElement *multiphysicsCel = dynamic_cast<TPZMultiphysicsElement *>(origCel);
                if (!multiphysicsCel) {
                    DebugStop();
                }
                TPZCompEl *fluxCel = multiphysicsCel->Element(0);
                TPZManVector<REAL> sol;
                fluxCel->Solution(pt, 0, sol);

                REAL InvKm = 1./ bc->Val1()(0, 0);
                REAL g = bc->Val1()(1, 0);//tem que ser Kgradu
                for (int ishape = 0; ishape < nshape; ishape++) {
                    efbc(ishape, 0) += weight * (InvKm * (sol[0] + g) + u_D) * phi(ishape, 0);
                    for (int jshape = 0; jshape < nshape; jshape++) {
                        ekbc(ishape, jshape) += weight * phi(ishape, 0) * phi(jshape, 0);
                    }
                }
                break;
            }
            default:
                std::cout << "Invalid BC type.\n";
                DebugStop();
            }
        }

        //ekbc.Print(std::cout);
        //efbc.Print(std::cout);

        ekbc.SolveDirect(efbc, ELU);
        //        efbc.Print(std::cout << "Solution ");

        int count = 0;
        TPZFMatrix<STATE> &mesh_sol = fPressureMesh->Solution();
        for (int ic = 0; ic < nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            int64_t seqnum = c.SequenceNumber();
            int64_t pos = fPressureMesh->Block().Position(seqnum);
            int ndof = c.NShape() * c.NState();
            for (int idf = 0; idf < ndof; idf++) {
                mesh_sol(pos + idf, 0) = efbc(count++);
            }
        }
    }

    {
        std::ofstream out("PressureAfterL2Projection_2.txt");
        fPressureMesh->Print(out);
    }
}

void TPZHybridH1CreateH1Reconstruction::VerifyAverage(int target_dim) {

    TPZCompMesh *pressure_mesh = fMultiphysicsH1reconstructionMesh->MeshVector()[1];
    TPZGeoMesh *gmesh = pressure_mesh->Reference();

    //    std::ofstream out("PressureToAverage.txt");
    //    pressure_mesh->Print(out);

    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    int64_t nel = pressure_mesh->NElements();
    // load the pressure elements of dimension target_dim+1
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() != target_dim + 1 && gel->MaterialId() != fHybridH1EE->fPressureSkeletonMatId) continue;
        gel->SetReference(cel);
    }

    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        if (!cel) {
            continue;
        }

        TPZGeoEl *gel = cel->Reference();

        if (!gel || gel->Dimension() != target_dim) {
            continue;
        }

        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);

        // Skip calculation if the element is a boundary condition
        int matid = gel->MaterialId();
        TPZMaterial *mat = pressure_mesh->FindMaterial(matid);
        // TODO change this. Look for matIDs in bcMatIds instead. Only cast in debug mode for further checks
#ifdef ERRORESTIMATION_DEBUG
        TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
        if (bc) continue;

#endif
        // Skip calculation if the element is a small skeleton
        bool largeSideExists = false;
        if (cel->Connect(0).HasDependency()) largeSideExists = true;

#ifdef ERRORESTIMATION_DEBUG
        int nsides = gel->NSides();
        TPZGeoElSide side(gel, nsides - 1);
        TPZGeoElSideAncestors ancestors(side);
        TPZGeoElSide largerNeigh = ancestors.HasLarger(fPressureSkeletonMatId);
        if (largeSideExists && !largerNeigh) DebugStop();
#endif
        if (largeSideExists) continue;

        // Start veryfication

        TPZGeoMesh *gmesh = pressure_mesh->Reference();
        int dim = gmesh->Dimension();
        gmesh->ResetReference();
        pressure_mesh->LoadReferences();
        int nc = cel->NConnects();
        int order = cel->Connect(nc - 1).Order();

#ifdef ERRORESTIMATION_DEBUG
        for (int ic = 0; ic < nc; ic++) {
            if (cel->Connect(ic).HasDependency()) DebugStop();
        }
#endif

        std::cout << "Computing average for compel " << el << ", gel: " << cel->Reference()->Index() << "\n";

        int target_dim = gel->Dimension();
        if (target_dim == dim - 1 && gel->MaterialId() != fHybridH1EE->fPressureSkeletonMatId) {
            DebugStop();
        }

        // We denote large the skeleton element passed to this method, which is whether the skeleton between two elements
        // of the same level or the larger skeleton in the region of a hanging node.
        TPZGeoElSide largeSkeletonSide(gel, gel->NSides() - 1);

        // Stack of TPZCompElSides containing the volumetric neighbours from which we obtain the pressure solution to
        // be used in the average computation.
        TPZCompElSide SkelSide = largeSkeletonSide.Reference();
        TPZStack<TPZCompElSide> volumeNeighSides;
        largeSkeletonSide.EqualLevelCompElementList(volumeNeighSides, 1, 0);
        if (volumeNeighSides.size() < 1) DebugStop();

        // If we are on a hanging side, we need to store the smaller skeleton sides and their respective volumetric
        // neighbours
        TPZStack<TPZGeoElSide> smallerSkelSides;
        bool noHangingSide = true;
        if (volumeNeighSides.size() == 1) {
            noHangingSide = false;
            TPZGeoElSidePartition partition(largeSkeletonSide);
            partition.HigherLevelNeighbours(smallerSkelSides, fHybridH1EE->fPressureSkeletonMatId);
            if (smallerSkelSides.size() == 1) DebugStop();
            for (int iskel = 0; iskel < smallerSkelSides.size(); iskel++) {
                smallerSkelSides[iskel].EqualLevelCompElementList(volumeNeighSides, 1, 0);
            }
        }

#ifdef ERRORESTIMATION_DEBUG
        for (int ineigh = 0; ineigh < volumeNeighSides.size(); ineigh++) {
            if (volumeNeighSides[ineigh].Element()->Dimension() != dim) DebugStop();
        }
#endif

        // Calculate average weights for left/right elements
        REAL sum_weights = 0.;
        TPZGeoEl *leftVolumeGel = volumeNeighSides[0].Element()->Reference();
        int matId = leftVolumeGel->MaterialId();
        sum_weights += fMatid_weights[matId];
        TPZGeoEl *rightVolumeGel = volumeNeighSides[1].Element()->Reference();
        matId = rightVolumeGel->MaterialId();
        sum_weights += fMatid_weights[matId];
        fPressureweights[cel->Index()] = sum_weights / 2;
        if (!noHangingSide) {
            for (int i = 0; i < smallerSkelSides.size(); i++) {
                TPZGeoElSide smallSkelGelside = smallerSkelSides[i];
                TPZCompEl *smallSkelCel = smallerSkelSides[i].Element()->Reference();
                fPressureweights[smallSkelCel->Index()] = sum_weights / 2;
            }
        }

        // TODO refactor this part so the transformations aren't stored but rather calculated inside the main loop at line 1480

        // Fills vectors of TPZTransform with the transformation from the skeleton to each left/right volumetric neighbours
        // and from the the small (right) to the large (left) skel if we are on a hanging side.
        int64_t nSkelsToIntegrate = volumeNeighSides.size() - 1;
        TPZManVector<TPZTransform<REAL>, 5> leftSkelToVolumeTrans(nSkelsToIntegrate);
        TPZManVector<TPZTransform<REAL>, 5> rightSkelToVolumeTrans(nSkelsToIntegrate);
        TPZManVector<TPZTransform<REAL>, 5> smallSkelToLargeSkelTrans(nSkelsToIntegrate);

        // Transformation between the large skeleton to the target_dim dimensional side of the volume element
        TPZTransform<REAL> largeSkelToVolumeTrans =
            largeSkeletonSide.NeighbourSideTransform(volumeNeighSides[0].Reference());
        TPZGeoEl *volumeNeigh = volumeNeighSides[0].Element()->Reference();
        // Transformation of the target_dim dimensional  side of the volume element to its higher dimension side
        TPZTransform<REAL> tmp = volumeNeigh->SideToSideTransform(volumeNeighSides[0].Side(),
                                                                  volumeNeigh->NSides() - 1);
        largeSkelToVolumeTrans = tmp.Multiply(largeSkelToVolumeTrans);

        if (noHangingSide) {
            // In this case the largeSkeletonSide is the same for both left/right regions.
            // We denote left the side at volumeNeighSides[0] and the remaining right.
            leftSkelToVolumeTrans[0] = largeSkelToVolumeTrans;

            rightSkelToVolumeTrans[0] = largeSkeletonSide.NeighbourSideTransform(volumeNeighSides[1].Reference());
            volumeNeigh = volumeNeighSides[1].Element()->Reference();
            tmp = volumeNeigh->SideToSideTransform(volumeNeighSides[1].Side(), volumeNeigh->NSides() - 1);
            rightSkelToVolumeTrans[0] = tmp.Multiply(rightSkelToVolumeTrans[0]);
        } else {
            // If we are on a hanging side we need to build left/right transformations for every small skeleton
            for (int iskel = 0; iskel < smallerSkelSides.size(); iskel++) {
                rightSkelToVolumeTrans[iskel] =
                    smallerSkelSides[iskel].NeighbourSideTransform(volumeNeighSides[iskel + 1].Reference());
                volumeNeigh = volumeNeighSides[iskel + 1].Element()->Reference();
                tmp = volumeNeigh->SideToSideTransform(volumeNeighSides[iskel + 1].Side(), volumeNeigh->NSides() - 1);
                rightSkelToVolumeTrans[iskel] = tmp.Multiply(rightSkelToVolumeTrans[iskel]);

                smallSkelToLargeSkelTrans[iskel] = TPZTransform<REAL>(target_dim);
                smallerSkelSides[iskel].SideTransform3(largeSkeletonSide, smallSkelToLargeSkelTrans[iskel]);
                leftSkelToVolumeTrans[iskel] = largeSkelToVolumeTrans.Multiply(smallSkelToLargeSkelTrans[iskel]);
            }
        }

        int nshape = intel->NShapeF();
        TPZFNMatrix<20, REAL> L2Mat(nshape, nshape, 0.), L2Rhs(nshape, 1, 0.);
        TPZFNMatrix<220, REAL> phi(nshape, 1, 0.), dphi(dim, nshape);
        for (int iskel = 0; iskel < nSkelsToIntegrate; iskel++) {
            TPZGeoElSide integrationGeoElSide;
            if (noHangingSide) {
                integrationGeoElSide = largeSkeletonSide;
            } else {
                integrationGeoElSide = smallerSkelSides[iskel];
            }
            std::unique_ptr<TPZIntPoints> intpoints(
                gel->CreateSideIntegrationRule(integrationGeoElSide.Side(), 2 * order));

            REAL left_weight = fMatid_weights[leftVolumeGel->MaterialId()] / sum_weights;
            rightVolumeGel = volumeNeighSides[iskel + 1].Element()->Reference();
            REAL right_weight = fMatid_weights[rightVolumeGel->MaterialId()] / sum_weights;

            if (!noHangingSide) {
                AdditionalAverageWeights(largeSkeletonSide.Element(), integrationGeoElSide.Element(), left_weight,
                                         right_weight, sum_weights);
            }

            int64_t nintpoints = intpoints->NPoints();
            //std::cout <<"########\n\n skelId = " << gel->Index() << "; leftVolId: " <<  leftVolumeGel->Index() <<"; rightVolId: " <<  rightVolumeGel->Index() << "\n";
            for (int64_t ip = 0; ip < nintpoints; ip++) {
                TPZManVector<REAL, 3> pt_left_skel(target_dim, 0.), pt_right_skel(target_dim, 0.);
                TPZManVector<REAL, 3> pt_left_vol(target_dim + 1, 0.), pt_right_vol(target_dim + 1, 0.);
                TPZManVector<REAL, 3> x_left(3, 0.);

                REAL weight;
                intpoints->Point(ip, pt_right_skel, weight);
                if (noHangingSide) {
                    pt_left_skel = pt_right_skel;
                } else {
                    smallSkelToLargeSkelTrans[iskel].Apply(pt_right_skel, pt_left_skel);
                }

                // Get shape at integration point
                intel->Shape(pt_left_skel, phi, dphi);

                // Get solution from left/right sides
                leftSkelToVolumeTrans[iskel].Apply(pt_right_skel, pt_left_vol);
                rightSkelToVolumeTrans[iskel].Apply(pt_right_skel, pt_right_vol);
                TPZVec<STATE> left_sol, right_sol,skel_sol;
                volumeNeighSides[0].Element()->Solution(pt_left_vol, 0, left_sol);
                volumeNeighSides[iskel + 1].Element()->Solution(pt_right_vol, 0, right_sol);

                STATE average_sol = left_weight * left_sol[0] + right_weight * right_sol[0];

                TPZFNMatrix<9, REAL> jac(dim, dim), jacinv(dim, dim), axes(dim, 3);
                REAL detjac;
                integrationGeoElSide.Jacobian(pt_right_skel, jac, axes, detjac, jacinv);

                SkelSide.Element()->Solution(pt_left_skel,0,skel_sol);

                leftVolumeGel->X(pt_left_vol,x_left);

                //std::cout <<"X: [" << x_left[0] << "," <<x_left[1] << "," <<x_left[2] << "]: skelSol: " << skel_sol[0] << "; leftSol/weight: " << left_sol[0] << "/" << left_weight << "; rightSol/weight: " << right_sol[0] << "/" << right_weight <<"\n";

                if(!IsZero(skel_sol[0]- (left_sol[0]+right_sol[0])/2)) {
                    DebugStop();
                }
            }
            //std::cout <<"\n";
        }
    }
}


/// compute the average pressures of the hybridized form of the H(div) mesh
void TPZHybridH1CreateH1Reconstruction::ComputeAverageFacePressures() {
    DebugStop();
    /*TPZCompMesh *pressure = MeshVector()[1];
    TPZCompMesh *pressure_mesh = fMultiphysicsH1reconstructionMesh->MeshVector()[1];
    int fInterfaceMatid = fHybridizer.fLagrangeInterface;
    TPZGeoMesh *gmesh = pressure->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    pressure->LoadReferences();
    int64_t nel = pressure_mesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        if (!cel || !cel->Reference() || cel->Reference()->Dimension() != dim - 1) {
            continue;
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        TPZGeoEl *gel = cel->Reference();
        if (gel->MaterialId() != fInterfaceMatid) {
            continue;
        }
        if (!intel || gel->Dimension() != dim - 1) {
            DebugStop();
        }
        int nc = cel->NConnects();
        int order = cel->Connect(nc - 1).Order();
        TPZGeoElSide gelside(gel, gel->NSides() - 1);
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 1, 0);
        TPZManVector<TPZTransform<REAL>, 2> tr(2);
        tr[0] = gelside.NeighbourSideTransform(celstack[0].Reference());
        {
            TPZGeoEl *right = celstack[0].Element()->Reference();
            TPZTransform<REAL> tmp = right->SideToSideTransform(celstack[0].Side(), right->NSides() - 1);
            tr[0] = tmp.Multiply(tr[0]);
        }
        if (celstack.size() == 1) {
            TPZCompElSide lowlevel = gelside.LowerLevelCompElementList2(1);
            if (!lowlevel) {
                DebugStop();
            }
            celstack.Push(lowlevel);
            tr[1] = TPZTransform<REAL>(gelside.Dimension());
            gelside.SideTransform3(lowlevel.Reference(), tr[1]);
        } else if (celstack.size() == 2) {
            tr[1] = gelside.NeighbourSideTransform(celstack[1].Reference());
        } else {
            DebugStop();
        }
        {
            TPZGeoEl *right = celstack[1].Element()->Reference();
            TPZTransform<REAL> tmp = right->SideToSideTransform(celstack[1].Side(), right->NSides() - 1);
            tr[1] = tmp.Multiply(tr[1]);
        }

        std::unique_ptr<TPZIntPoints> intp(gel->CreateSideIntegrationRule(gel->NSides() - 1, 2 * order));
        int nshape = intel->NShapeF();
        TPZFNMatrix<20, REAL> L2Mat(nshape, nshape, 0.), L2Rhs(nshape, 1, 0.);
        TPZFNMatrix<220, REAL> phi(nshape, 1, 0.), dshape(dim, nshape);
        int64_t npoints = intp->NPoints();
        for (int64_t ip = 0; ip < npoints; ip++) {
            TPZManVector<REAL, 3> pt(dim - 1, 0.), pt1(dim, 0.), pt2(dim, 0.), sol1(1), sol2(1);
            REAL weight;
            intp->Point(ip, pt, weight);
            intel->Shape(pt, phi, dshape);
            tr[0].Apply(pt, pt1);
            tr[1].Apply(pt, pt2);
            celstack[0].Element()->Solution(pt1, 0, sol1);//solucao a esquerda
            celstack[1].Element()->Solution(pt2, 0, sol2);//solucao a direita
                                                          //           std::cout << "Values " << sol1 << " " << sol2 << std::endl;
                                                          //projecao L2 da media das soluceos no espaco Lh, do esqueleto da malha
            for (int ishape = 0; ishape < nshape; ishape++) {
                L2Rhs(ishape, 0) += weight * phi(ishape, 0) * (sol1[0] + sol2[0]) / 2.;
                for (int jshape = 0; jshape < nshape; jshape++) {
                    L2Mat(ishape, jshape) += weight * phi(ishape, 0) * phi(jshape, 0);
                }
            }
        }
        L2Mat.SolveDirect(L2Rhs, ECholesky);
        //apos este passo temos uma pressao que é continua ao longo das interfaces dos elementos, nos esqueletos. Falta suavizar nos vértices
        // L2Rhs.Print("Average pressure");
        int count = 0;
        for (int ic = 0; ic < nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            int64_t seqnum = c.SequenceNumber();
            int64_t pos = pressure_mesh->Block().Position(seqnum);
            int ndof = c.NShape() * c.NState();
            for (int idf = 0; idf < ndof; idf++) {
                pressure_mesh->Solution()(pos + idf, 0) = L2Rhs(count++);
            }
        }
    }
    TPZManVector<TPZCompMesh *, 2> meshvec(2);
    meshvec[0] = fMultiphysicsH1reconstructionMesh->MeshVector()[0];
    meshvec[1] = fMultiphysicsH1reconstructionMesh->MeshVector()[1];
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, fMultiphysicsH1reconstructionMesh);
*/}