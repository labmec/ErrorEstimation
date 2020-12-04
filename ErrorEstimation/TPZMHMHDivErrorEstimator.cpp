
#include "TPZMHMHDivErrorEstimator.h"
#include "TPZMHMHDivErrorEstimatorMaterial.h"
#include "TPZNullMaterial.h"
#include "TPZVTKGeoMesh.h"
#include "TPZVecL2.h"
#include "pzbndcond.h"
#include "pzintel.h"
#include "pzsubcmesh.h"
#include <Material/TPZHDivErrorEstimateMaterial.h>
#include <Mesh/pzmultiphysicscompel.h>
#include <Mesh/TPZCompMeshTools.h>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("HDivErrorEstimator"));
#endif

// a method for generating the hybridized multiphysics post processing mesh
void TPZMHMHDivErrorEstimator::CreatePostProcessingMesh()
{

    // initialize the post processing mesh
    fPostProcMesh.SetReference(fOriginal->Reference());

    fOriginal->CopyMaterials(fPostProcMesh);
    // Switch the material from mixed to TPZMHMHDivErrorEstimationMaterial
    SwitchMaterialObjects();

    TPZManVector<TPZCompMesh *> meshvec(4);
    TPZManVector<int,4> active(4,0);
    active[1] = 1;

    meshvec[0] = 0;
    meshvec[1] = CreatePressureMesh();
    meshvec[2] = fOriginal->MeshVector()[0];
    meshvec[3] = fOriginal->MeshVector()[1];

    if (fPostProcesswithHDiv) {
        meshvec[0] = CreateFluxMesh();
        active[0] = 1;
    }

    if (!fPostProcesswithHDiv) {
        CreateSkeletonElements(meshvec[1]);
        CreateSkeletonApproximationSpace(meshvec[1]);
    }

    //enriquecer no MHM tbem?
    {
        IncreasePressureSideOrders(meshvec[1]);//malha da pressao
        if(fPostProcesswithHDiv)
        {
            IncreaseSideOrders(meshvec[0]);//malha do fluxo
        }
    }

    {
        std::ofstream out("PressureMeshAfterIncreaseSideOrders.txt");
        meshvec[1]->Print(out);
    }

    RemoveMaterialObjects(fPostProcMesh.MaterialVec());
    fPostProcMesh.BuildMultiphysicsSpace(active, meshvec);
    bool groupelements = false;

    //for reconstruction on H1 hybridezed mesh is not applied
    //just for Hdiv reconstruction the hybridizer will be use
    if(fPostProcesswithHDiv){
        fHybridizer.HybridizeGivenMesh(fPostProcMesh,groupelements);
    }

    {
        std::ofstream out3("MalhaComPressureSkeletonMat.txt");
        fPostProcMesh.Print(out3);
    }
#ifdef PZDEBUG2
    {
        std::ofstream out1("fluxbeforeSub.txt");
        meshvec[0]->Print(out1);
        std::ofstream out2("pressurebeforeSub.txt");
        meshvec[1]->Print(out2);
        std::ofstream out3("mphyspost_before.txt");
        fPostProcMesh.Print(out3);
    }
#endif

    if (fPostProcesswithHDiv) {
       // fPressureSkeletonMatId = fHybridizer.fLagrangeInterface;
    }

    SubStructurePostProcessingMesh();

    {
        std::ofstream out("CompletePressureMeshMHM.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(meshvec[1], out, {}, false, false);
    }

    ComputePressureWeights();

#ifdef PZDEBUG2
    {
        std::ofstream out1("fluxbePostSub.txt");
        meshvec[0]->Print(out1);
    }
#endif
}

// a method for transferring the multiphysics elements in submeshes
void TPZMHMHDivErrorEstimator::SubStructurePostProcessingMesh()
{
    TPZGeoMesh *gmesh = fOriginal->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    fOriginal->LoadReferences();
    int64_t numgel = gmesh->NElements();
    TPZVec<TPZSubCompMesh *> ReferredMesh(numgel,0);
    // associate with each geometric element a subcmesh object
    // we do this because the original mesh and post processing mesh share the same geometric mesh
    // if an element of the original mesh is in a given subcmesh then the corresponding element in the
    // post processing mesh will be put in the corresponding subcmesh
    for(int64_t el=0; el<numgel; el++)
    {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) continue;
        TPZCompEl *cel = gel->Reference();
        if(!cel) continue;
        TPZCompMesh *mesh = cel->Mesh();
        TPZSubCompMesh *ref = dynamic_cast<TPZSubCompMesh *>(mesh);
        ReferredMesh[el] = ref;
    }
#ifdef PZDEBUG2
    {
        std::ofstream file("GmeshSub.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
        std::ofstream out("gmesh_sub.txt");
        gmesh->Print(out);
        std::ofstream out2("original.txt");
        fOriginal->Print(out2);
        
    }

#endif

    // create the sub comp meshes and create a mapping data structure

    int64_t nel = fPostProcMesh.NElements();
    TPZVec<TPZSubCompMesh *> ElementMesh(nel,0);
    std::map<TPZSubCompMesh *,TPZSubCompMesh *> submeshmap;
    // associate the submesh with the volumetric elements
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fPostProcMesh.Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        // we are only building the data structure for elements with mesh dimension
        // this excludes interface elements, wrappers, lagrange multipliers etc

        TPZSubCompMesh *submesh = ReferredMesh[gel->Index()];
        if(!submesh) continue;

        auto iter = submeshmap.find(submesh);
        if (iter == submeshmap.end()) {
            int64_t index;
            TPZSubCompMesh *subcmesh = new TPZSubCompMesh(fPostProcMesh,index);
            submeshmap[submesh] = subcmesh;
            ElementMesh[el] = subcmesh;
        }
        else
        {
            ElementMesh[el] = iter->second;
        }
    }
    // create a data structure associating each element with a group
    // the value of elementgroup is the index of the computational element that
    // will nucleate the group
    // this only works if the post processing mesh is H(div) hybridized
    
    if(fPostProcesswithHDiv) {
    
        TPZVec<int64_t> elementgroup;
        fHybridizer.AssociateElements(&fPostProcMesh, elementgroup);
        // transfer the elements in the submesh indicated by elementgroup
        // associate the submesh with the volumetric elements
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = fPostProcMesh.Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if(!gel) DebugStop();
            if(elementgroup[el] == -1)
            {
                // this element is not assocated with a volumetric element
                continue;
            }
            int64_t volume_element = elementgroup[el];
            if(ElementMesh[volume_element] == 0) DebugStop();
            TPZSubCompMesh *submesh = ElementMesh[volume_element];
            submesh->TransferElement(&fPostProcMesh, el);
        }
        fPostProcMesh.ComputeNodElCon();
        for(auto iter : submeshmap)
        {
            iter.second->MakeAllInternal();
        }
    #ifdef PZDEBUG2
        {
            int64_t nel = fPostProcMesh.NElements();
            for(int64_t el = 0; el<nel; el++)
            {
                TPZCompEl *cel = fPostProcMesh.Element(el);
                TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
                if(sub)
                {
                    std::stringstream sout;
                    sout << "postproc_part_submesh_" << el << ".vtk";
                    std::ofstream file(sout.str());
                    TPZVTKGeoMesh::PrintCMeshVTK(sub, file,true);
                }
            }

        }
#endif

    } else {
        std::set<int64_t> connectlist;
        ComputeConnectsNextToSkeleton(connectlist);

        {
            std::ofstream out("MalhaTesteBeforeTransfer.txt");
            fPostProcMesh.Print(out);
        }

        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl* cel = fPostProcMesh.Element(el);
            if (!cel) continue;
            TPZGeoEl* gel = cel->Reference();
            if (!gel) DebugStop();
            TPZSubCompMesh* submesh = ElementMesh[el];
            if (!submesh)continue;

            submesh->TransferElement(&fPostProcMesh, el);
            cel = fPostProcMesh.Element(el);
            if (cel) DebugStop();
        }

        fPostProcMesh.ComputeNodElCon();

        // {
        //     std::cout << "connectlist ";
        //     for (auto it : connectlist) {
        //         std::cout << it << " ";
        //     }
        //     std::cout << '\n';
        // }

        for (auto it : connectlist) {
            fPostProcMesh.ConnectVec()[it].IncrementElConnected();
        }

        for (auto iter : submeshmap) {
            iter.second->ExpandSolution();
        }

//        {
//            std::ofstream out("MalhaTeste1.txt");
//            fPostProcMesh.Print(out);
//        }

        for (auto iter : submeshmap) {
            iter.second->MakeAllInternal();
        }

    }

    fPostProcMesh.ComputeNodElCon();
    fPostProcMesh.CleanUpUnconnectedNodes();
    {
        std::ofstream out("MalhaTeste2.txt");
        fPostProcMesh.Print(out);
    }

    // set an analysis type for the submeshes
    {
        int64_t nel = fPostProcMesh.NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = fPostProcMesh.Element(el);
            TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
            if (sub) {
//                {
//                    std::ofstream out2("Sub.txt");
//                    sub->Print(out2);
//                }

                if (fPostProcesswithHDiv) {
 //                   fHybridizer.GroupandCondenseElements(sub);
                    sub->CleanUpUnconnectedNodes();
                }
                int numthreads = 0;
                int preconditioned = 0;
                TPZAutoPointer<TPZGuiInterface> guiInterface;

                sub->SetAnalysisSkyline(numthreads, preconditioned, guiInterface);
            }
        }

#ifdef PZDEBUG2
        {
            int64_t nel = fPostProcMesh.NElements();
            for(int64_t el = 0; el<nel; el++)
            {
                TPZCompEl *cel = fPostProcMesh.Element(el);
                TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
                if(sub)
                {
                    std::stringstream sout;
                    sout << "postproc_submesh_" << el << ".vtk";
                    std::ofstream file(sout.str());
                    TPZVTKGeoMesh::PrintCMeshVTK(sub, file,true);
                }
            }
            fPostProcMesh.LoadReferences();
            std::ofstream out("postproc_cmesh_substruct.vtk");
            TPZVTKGeoMesh::PrintCMeshVTK(&fPostProcMesh, out,true);
        }
#endif
    }

}


// a method for generating the HDiv mesh
TPZCompMesh *TPZMHMHDivErrorEstimator::CreateFluxMesh()
{
    
    TPZCompMesh *OrigFlux = fOriginal->MeshVector()[0];
    TPZGeoMesh *gmesh = OrigFlux->Reference();
    gmesh->ResetReference();
    TPZCompMesh *fluxmesh = new TPZCompMesh(OrigFlux->Reference());
    fluxmesh->SetAllCreateFunctionsHDiv();
    fluxmesh->SetDefaultOrder(OrigFlux->GetDefaultOrder());
    OrigFlux->CopyMaterials(*fluxmesh);
    RemoveMaterialObjects(fluxmesh->MaterialVec());
    int dim = gmesh->Dimension();
    int64_t nel = OrigFlux->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = OrigFlux->Element(el);
        if(!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        TPZGeoEl *gel = cel->Reference();
        bool isbcmat = fMHM->fMaterialBCIds.find(gel->MaterialId()) != fMHM->fMaterialBCIds.end();
        // if the element is of lower dimension and is not a boundary
        // don't create a flux element
        if(gel->Dimension() != dim && !isbcmat) continue;
        TPZMaterial *mat = fluxmesh->FindMaterial(gel->MaterialId());
        if(!mat) DebugStop();
        int64_t index;
        TPZCompEl *create_cel = fluxmesh->CreateCompEl(gel, index);
        TPZInterpolatedElement *hdiv_cel = dynamic_cast<TPZInterpolatedElement *>(create_cel);
        /// equate the order of the connects with the order of the original mesh
        for (int is = gel->NCornerNodes(); is<gel->NSides(); is++) {
            int nside_connects = intel->NSideConnects(is);
            // if the number of side connects is zero do nothing
            if(nside_connects == 0) continue;
            // get the last connect of the side (in the case of hdiv there is always a single connect)
            TPZConnect &corig = intel->SideConnect(nside_connects-1, is);
            int loc_index = intel->SideConnectLocId(nside_connects-1, is);
            
            TPZConnect &cnew = hdiv_cel->Connect(loc_index);
            if(cnew.Order() != corig.Order());
            hdiv_cel->SetSideOrder(is, corig.Order());
        }
    }
    fluxmesh->ExpandSolution();
    return fluxmesh;
}

// method for creating a discontinuous pressure mesh
TPZCompMesh *TPZMHMHDivErrorEstimator::CreatePressureMesh() {
    if (fPostProcesswithHDiv) {
        return CreateDiscontinuousPressureMesh();
    } else {
        return CreateInternallyContinuousPressureMesh();
    }
}

// a method for creating the pressure mesh
TPZCompMesh *TPZMHMHDivErrorEstimator::CreateDiscontinuousPressureMesh()
{
    TPZCompMesh *OrigPressure = fOriginal->MeshVector()[1];
    TPZGeoMesh *gmesh = OrigPressure->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    TPZCompMesh *pressmesh = new TPZCompMesh(gmesh);
    OrigPressure->CopyMaterials(*pressmesh);
    RemoveMaterialObjects(pressmesh->MaterialVec());
    pressmesh->SetDefaultOrder(OrigPressure->GetDefaultOrder());
    pressmesh->SetAllCreateFunctionsContinuous();
    pressmesh->ApproxSpace().CreateDisconnectedElements(true);
    int64_t nel = OrigPressure->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = OrigPressure->Element(el);
        if(!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if(!intel) DebugStop();
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() != dim)
        {
            std::cout << __PRETTY_FUNCTION__ << " dimension of pressure element unexpected\n";
            continue;
        }
        int64_t index;
        TPZCompEl *celnew = pressmesh->CreateCompEl(gel, index);
        TPZInterpolatedElement *newcel = dynamic_cast<TPZInterpolatedElement *>(celnew);
        int nc = cel->NConnects();
        int order = cel->Connect(nc-1).Order();
        newcel->PRefine(order);
        gel->ResetReference();
    }
    pressmesh->ExpandSolution();
    int64_t nc = pressmesh->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        pressmesh->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    return pressmesh;
}

TPZCompMesh *TPZMHMHDivErrorEstimator::CreateInternallyContinuousPressureMesh() {
    TPZCompMesh *original_pressure = fOriginal->MeshVector()[1];
    TPZGeoMesh *gmesh = original_pressure->Reference();
    gmesh->ResetReference();
    original_pressure->LoadReferences();

    // We need to fill a tuple with the information of the MHM domain that each geometric element belongs to
    // and the corresponding computational element in the original pressure mesh.
    // The MHM domain info allows the creation of continuous space inside a MHM domain, but not globally.
    // The original pressure comp. element is used to retrieve the original approximation order.
    int64_t nel = gmesh->NElements();
    auto geoToMHM = fMHM->GetGeoToMHMDomain();
    TPZManVector<std::tuple<int64_t, int64_t, TPZCompEl*>> MHMOfEachGeoEl(nel);
    for (int i = 0; i < nel; i++) {
        TPZGeoEl * gel = gmesh->Element(i);
        if (!gel) {
            MHMOfEachGeoEl[i] = {-1, i, nullptr};
            continue;
        }
        TPZCompEl * orig_cel = gel->Reference();
        if (!orig_cel) {
            MHMOfEachGeoEl[i] = {-1, i, nullptr};
            continue;
        }
        MHMOfEachGeoEl[i] = std::make_tuple(geoToMHM[i], i, orig_cel);
    }

    int dim = gmesh->Dimension();
    std::map<int64_t, int64_t> bcToMHM;
    int64_t nElem = gmesh->NElements();
    for (int64_t el = 0; el < nElem; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel) continue;
        if (gel->Dimension() != dim) continue;

        for (int iside = gel->NCornerNodes(); iside < gel->NSides() - 1; iside++) {
            TPZGeoElSide gelside(gel, iside);
            TPZGeoElSide neighbour;
            neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                if (neighbour.Element()->Dimension() != dim) {
                    TPZMaterial *neighMat = fPostProcMesh.FindMaterial(neighbour.Element()->MaterialId());
                    if (neighMat) {
                        TPZBndCond *bc = dynamic_cast<TPZBndCond *>(neighMat);
                        if (bc) {
                            int64_t bcId = neighbour.Element()->Index();
                            MHMOfEachGeoEl[bcId] = {geoToMHM[el], bcId, nullptr};
                        }
                    }
                }
                neighbour = neighbour.Neighbour();
            }
        }
    }

    std::sort(&MHMOfEachGeoEl[0], &MHMOfEachGeoEl[nel - 1] + 1);

    // Create pressure mesh
    TPZCompMesh *reconstruction_pressure = new TPZCompMesh(gmesh);

    // Copies volume materials
    original_pressure->CopyMaterials(*reconstruction_pressure);
    RemoveMaterialObjects(reconstruction_pressure->MaterialVec());

    // Copies BC materials
    std::set<int> bcMatIDs;
    for (auto it : fOriginal->MaterialVec()) {
        TPZMaterial *mat = it.second;
        TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
        if (bc) {
            int bcID = bc->Material()->Id();
            TPZMaterial *pressure_mat = original_pressure->FindMaterial(bcID);
            TPZMaterial *bc_mat = pressure_mat->CreateBC(pressure_mat, mat->Id(), bc->Type(), bc->Val1(), bc->Val2());
            if (fExact) {
                bc_mat->SetForcingFunction(fExact->Exact());
            }
            reconstruction_pressure->InsertMaterialObject(bc_mat);
        }
    }

    reconstruction_pressure->SetDefaultOrder(original_pressure->GetDefaultOrder());
    reconstruction_pressure->SetAllCreateFunctionsContinuous();
    reconstruction_pressure->ApproxSpace().CreateDisconnectedElements(false);
    gmesh->ResetReference();

    // Creates elements in pressure mesh
    int64_t previousMHMDomain = -1;
    int64_t firstElemInMHMDomain = -1;
    for (int i = 0; i < MHMOfEachGeoEl.size(); i++) {
        int64_t MHMDomain = std::get<0>(MHMOfEachGeoEl[i]);
        int64_t elIndex = std::get<1>(MHMOfEachGeoEl[i]);

        if (MHMDomain == -1) continue;

        if (MHMDomain != previousMHMDomain) {
            if (previousMHMDomain != -1) {
                for (int j = firstElemInMHMDomain; j < i; j++) {
                    gmesh->Element(std::get<1>(MHMOfEachGeoEl[j]))->ResetReference();
                }
            }
            firstElemInMHMDomain = i;
            previousMHMDomain = MHMDomain;
        }

        // Create the pressure element
        TPZGeoEl *gel = gmesh->Element(elIndex);
        if (!gel || gel->HasSubElement()) continue;

        int64_t index;
        TPZCompEl *new_cel = reconstruction_pressure->CreateCompEl(gel, index);
        TPZInterpolatedElement *new_intel = dynamic_cast<TPZInterpolatedElement *>(new_cel);

        TPZCompEl * orig_cel = std::get<2>(MHMOfEachGeoEl[i]);
        if (orig_cel) {
            int nc = gel->Reference()->NConnects();
            int order = gel->Reference()->Connect(nc - 1).Order();
            new_intel->PRefine(order);
        }
        else {
            // TODO: review these choices
            // There are no BC elements in original pressure, so I'm not sure what to set as default order in this case.
            // What I'm doing right now is to set as the same order of the neighbour. I think PZ does this automatically
            // but I'm not sure. (Gustavo, 9/11/2020)
            TPZGeoElSide gelside(gel);
            TPZStack<TPZCompElSide> celstack;
            int onlyinterpolated = 1;
            int removeduplicates = 0;

            gelside.EqualLevelCompElementList(celstack, onlyinterpolated, removeduplicates);
            // A BC element should have only one neighbour from its highest dimension side
            if (celstack.size() != 1) DebugStop();

            int neigh_side_id = celstack[0].Side();
            int order = celstack[0].Element()->Connect(neigh_side_id).Order();
            new_intel->PRefine(order);
        }
    }

    // Resets references of last MHM domain
    for (int j = firstElemInMHMDomain; j < MHMOfEachGeoEl.size(); j++) {
        gmesh->Element(std::get<1>(MHMOfEachGeoEl[j]))->ResetReference();
    }

    return reconstruction_pressure;
}

// remove the materials that are not listed in MHM
void TPZMHMHDivErrorEstimator::RemoveMaterialObjects(std::map<int,TPZMaterial *> &matvec)
{
    bool changed = true;
    while(changed)
    {
        changed = false;
        for(auto iter : matvec)
        {
            int matid = iter.first;
            if(fMHM->fMaterialIds.find(matid) == fMHM->fMaterialIds.end() &&
               fMHM->fMaterialBCIds.find(matid) == fMHM->fMaterialBCIds.end())
            {
                std::cout << __PRETTY_FUNCTION__  << ": Removing material " << matid << '\n';
                TPZMaterial *mat = iter.second;
                delete mat;
                matvec.erase(matid);
                changed = true;
                break;
            }
        }
    }

}

// transfer embedded elements to submeshes
// the substructuring method will only transfer the H(div) and surrounding elements to the submesh
// it does not detect that the boundary pressure elements belong to the submesh. This is done in the
// following method
void TPZMHMHDivErrorEstimator::TransferEmbeddedElements()
{
    int64_t nel = fPostProcMesh.NElements();
    int64_t ncon = fPostProcMesh.NConnects();
    TPZVec<int> numelcon(ncon,0);
    TPZVec<TPZSubCompMesh *> connect_submesh(ncon,0);
    // numelcon : number of subcmesh sharing a connect
    // connect_submesh : pointer to the submesh if there is only one subdomain that owns the connect
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fPostProcMesh.Element(el);
        TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
        if(!sub) continue;
        int nconnect = sub->NConnects();
        for(int ic=0; ic<nconnect; ic++)
        {
            int64_t conindex = cel->ConnectIndex(ic);
            numelcon[conindex]++;
            if(numelcon[conindex] == 1)
            {
                connect_submesh[conindex] = sub;
            }
            else
            {
                connect_submesh[conindex] = 0;
            }
        }
    }
    // transfer the elements whose connects belong to a single subdomain
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fPostProcMesh.Element(el);
        if(!cel) continue;
        TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
        if(sub) continue;
        std::set<TPZSubCompMesh *> submeshes;
        int nconnect = cel->NConnects();
        for(int ic=0; ic<nconnect; ic++)
        {
            int64_t conindex = cel->ConnectIndex(ic);
            submeshes.insert(connect_submesh[conindex]);
        }
        if(submeshes.size() == 1)
        {
            TPZSubCompMesh *sub = *submeshes.begin();
            if(sub)
            {
                sub->TransferElement(&fPostProcMesh, el);
            }
        }
    }
    fPostProcMesh.ComputeNodElCon();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fPostProcMesh.Element(el);
        TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
        if(!sub) continue;
        sub->MakeAllInternal();
    }
}

// a method for computing the pressures between subdomains as average pressures
/// compute the average pressures of across edges of the H(div) mesh
void TPZMHMHDivErrorEstimator::ComputeAveragePressures(int target_dim)
{
    // load the pressure elements of the finite element approximation
    TPZCompMesh *OrigPressure = fOriginal->MeshVector()[1];
    TPZGeoMesh *gmesh = OrigPressure->Reference();
    gmesh->ResetReference();
    int dim = fPostProcMesh.Dimension();
    
    TPZCompMesh *postpressuremesh =  fPostProcMesh.MeshVector()[1];
    TPZCompMesh *loadmesh = OrigPressure;
    if(target_dim < dim-1) {
        loadmesh = postpressuremesh ;
    }
    int64_t nel = loadmesh->NElements();
    // load all elements of dimension target_dim+1
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = loadmesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() == target_dim+1) {
            cel->LoadElementReference();
        }
    }
    // compute the averages one element at a time
    nel = postpressuremesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = postpressuremesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        if(gel->Dimension() == target_dim)
        {
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            if(!intel) DebugStop();
            
            int matId = gel->MaterialId();
            if(matId != fPressureSkeletonMatId) continue;
            int64_t index = intel->Index();
            ComputeAverage(postpressuremesh,index);
        }
    }
    
}
/// compute the average pressure over corners
/// set the cornernode values equal to the averages
void TPZMHMHDivErrorEstimator::ComputeNodalAverages()
{
    // load the one dimensional interface elements
    // load the pressure elements of the finite element approximation
    TPZCompMesh *pressuremesh = fPostProcMesh.MeshVector()[1];
    TPZGeoMesh *gmesh = pressuremesh->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    int64_t nel = pressuremesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = pressuremesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        if (gel->Dimension() == dim-1) {
            cel->LoadElementReference();
        }
    }

    nel = pressuremesh->NElements();
    // compute the averages
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = pressuremesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel || !gel->Reference()) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if(!intel) DebugStop();
        if (gel->Dimension() == dim-1) {
            int ncorner = gel->NCornerNodes();
            for (int side = 0; side<ncorner; side++) {
                TPZCompElSide celside(cel,side);
                ComputeNodalAverage(celside);
            }
        }
    }
}

//void TPZMHMHDivErrorEstimator::SwitchMaterialObjects() {
//    // switch the material
//    for (auto matid : fPostProcMesh.MaterialVec()) {
//        TPZMixedPoisson *mixpoisson = dynamic_cast<TPZMixedPoisson *>(matid.second);
//        if (mixpoisson) {
//
//            TPZHDivErrorEstimateMaterial *newmat = new TPZHDivErrorEstimateMaterial(*mixpoisson);
//            newmat->fNeumannLocalProblem = false;
//
//            if (fExact) {
//                newmat->SetForcingFunction(fExact->Exact());
//                newmat->SetForcingFunction(fExact->ForcingFunction());
//            }
//
//            for (auto bcmat : fPostProcMesh.MaterialVec()) {
//                TPZBndCond *bc = dynamic_cast<TPZBndCond *>(bcmat.second);
//                if (bc) {
//                    bc->SetMaterial(newmat);
//                }
//            }
//            fPostProcMesh.MaterialVec()[newmat->Id()] = newmat;
//            delete mixpoisson;
//        }
//    }
//}

void TPZMHMHDivErrorEstimator::CreateSkeletonElements(TPZCompMesh * pressure_mesh) {

    TPZCompMesh* cmesh = fOriginal;
    TPZGeoMesh* gmesh = fOriginal->Reference();
    gmesh->ResetReference();
    cmesh->LoadReferences();

#ifdef PZDEBUG
    {
        std::ofstream fileVTK("GeoMeshBeforePressureSkeleton.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, fileVTK);
        std::ofstream fileTXT("GeoMeshBeforePressureSkeleton.txt");
        gmesh->Print(fileTXT);
    }
#endif

    // Assigns a material ID that has not been used yet
    int maxMatId = std::numeric_limits<int>::min();
    const int nel = gmesh->NElements();

    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->Element(iel);
        if(gel) maxMatId = std::max(maxMatId, gel->MaterialId());
    }

    if (maxMatId == std::numeric_limits<int>::min()) maxMatId = 0;

    fPressureSkeletonMatId = maxMatId + 1;
    int dim = gmesh->Dimension();

    for (int iel = 0; iel < nel; iel++) {

        TPZGeoEl *gel = gmesh->Element(iel);
        if (!gel) continue;
        TPZCompEl* cel = gel->Reference();

        if (!cel) continue;
        if (gel->Dimension() != dim) continue;

        // Iterates through the sides of the element
        int nsides = gel->NSides();
        for (int iside = 0; iside < nsides; iside++) {
            TPZGeoElSide gelside(gel, iside);

            // Filters boundary sides
            if (gelside.Dimension() != dim - 1) continue;

            // Gets compel sides of equal and lower (if existing) level linked to the gelside
            TPZStack<TPZCompElSide> celstack;
            gelside.EqualLevelCompElementList3(celstack, 1, 0);

            TPZCompElSide large = gelside.LowerLevelCompElementList2(1);
            if (large) celstack.Push(large);

            int nstack = celstack.size();
            if (nstack == 0) continue;

            for (int ist = 0; ist < nstack; ist++) {
                TPZCompElSide cneighbour = celstack[ist];
                if (!cneighbour) continue;

                TPZGeoElSide neighbour = cneighbour.Reference();

                // Filters neighbour sides that belong to volume elements
                if (neighbour.Element()->Dimension() != dim) continue;

                if (cneighbour.Element()->Mesh() != gel->Reference()->Mesh()) {

                    // This lambda checks if a skeleton has already been created over gelside
                    auto hasSkeletonNeighbour = [&] () -> bool {
                        neighbour = gelside.Neighbour();
                        while (neighbour != gelside) {
                            int neighbourMatId = neighbour.Element()->MaterialId();
                            if (neighbourMatId == fPressureSkeletonMatId) {
                                return true;
                            }
                            neighbour = neighbour.Neighbour();
                        }
                        return false;
                    };

                    if (!hasSkeletonNeighbour()) {
                        TPZGeoElBC(gelside, fPressureSkeletonMatId);
                        break;
                    }
                }
            }
        }
    }

#ifdef PZDEBUG
    {
        std::ofstream fileVTK("GeoMeshAfterPressureSkeleton.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, fileVTK);
        std::ofstream fileTXT("GeoMeshAfterPressureSkeleton.txt");
        gmesh->Print(fileTXT);
    }
#endif
}

void TPZMHMHDivErrorEstimator::CreateSkeletonApproximationSpace(TPZCompMesh *pressure_mesh) {

    TPZGeoMesh *gmesh = pressure_mesh->Reference();
    int dim = gmesh->Dimension();

    // Create skeleton elements in pressure mesh
    TPZNullMaterial *skeletonMat = new TPZNullMaterial(fPressureSkeletonMatId);
    skeletonMat->SetDimension(dim - 1);
    pressure_mesh->InsertMaterialObject(skeletonMat);

    set<int> matIdSkeleton = { fPressureSkeletonMatId };
    gmesh->ResetReference();

    pressure_mesh->ApproxSpace().CreateDisconnectedElements(true);
    pressure_mesh->AutoBuild(matIdSkeleton);
    pressure_mesh->ExpandSolution();
}

void TPZMHMHDivErrorEstimator::CopySolutionFromSkeleton() {

    // Pressure skeleton
    TPZCompMesh *pressuremesh = PressureMesh();
//    {
//        std::ofstream out("MeshBeforeCopySkeletonMeshVector1.txt");
//        pressuremesh->Print(out);
//    }

    PlotState("PressureBeforeCopyskeleton.vtk", pressuremesh->Dimension(), &fPostProcMesh);

    pressuremesh->Reference()->ResetReference();

    if(pressuremesh->Reference()!= fPostProcMesh.Reference()) DebugStop();
    
    int dim = pressuremesh->Dimension();
    int64_t nel = pressuremesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl* cel = pressuremesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!cel) continue;
        if(!intel) DebugStop();
        // load just d dimensional elements
        if (cel->Dimension() != dim) continue;
        cel->LoadElementReference();
    }

    {
        std::ofstream file("GmeshCopySkelton.txt");
        pressuremesh->Reference()->Print(file);
    }

    nel = pressuremesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl* cel = pressuremesh->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement* intel = dynamic_cast<TPZInterpolatedElement*>(cel);
        if (!intel) DebugStop();
        TPZGeoEl* gel = cel->Reference();
        // filters just (d-1) dimensional elements
        if (gel->Dimension() == dim) continue;

        int nsides = gel->NSides();
        for (int is = 0; is < nsides; is++) {
            TPZGeoElSide gelside(gel, is);
            
            int matgelSide = gelside.Element()->MaterialId();
            
            //std::cout<<"MatIdgelSide  "<<matgelSide<<"\n";
            TPZConnect &c = intel->Connect(is);
            int64_t c_gelSide_seqnum  = c.SequenceNumber();
            int c_blocksize = c.NShape() * c.NState();
            TPZStack<TPZCompElSide> celstack;

            gelside.EqualLevelCompElementList(celstack, 1, 0);

            int nst = celstack.NElements();
            if(nst==0) DebugStop();
            for (int ist = 0; ist < nst; ist++) {
                TPZCompElSide cneigh = celstack[ist];
                TPZGeoElSide gneigh = cneigh.Reference();

                TPZInterpolatedElement *intelneigh = dynamic_cast<TPZInterpolatedElement *>(cneigh.Element());
                if (!intelneigh) DebugStop();
                TPZConnect &con_neigh = intelneigh->Connect(cneigh.Side());
                int64_t c_neigh_seqnum = con_neigh.SequenceNumber();
                int con_size = con_neigh.NState() * con_neigh.NShape();
                if (con_size != c_blocksize) DebugStop();
                for (int ibl = 0; ibl < con_size; ibl++) {
                    //std::cout<<"valor da pressao connect neigh (d-dimensional) "<<c_neigh_seqnum<<" = "<<pressuremesh->Block()(c_neigh_seqnum, 0, ibl, 0)<<"\n";
                    //std::cout<<"valor da pressao connect "<<c_gelSide_seqnum<<" = "<<pressuremesh->Block()(c_gelSide_seqnum, 0, ibl, 0)<<"\n";
                    pressuremesh->Block()(c_neigh_seqnum, 0, ibl, 0) = pressuremesh->Block()(c_gelSide_seqnum, 0, ibl, 0);
                }
            }
        }
    }
//    {
//        std::ofstream out("MultiphysicsAfterCopySkeleton.txt");
//        fPostProcMesh.Print(out);
//        std::string file("PressureAfterCopyskeleton.vtk");
//        PlotState(file, 2, &fPostProcMesh);
//    }

    std::set<int64_t> connectList;
    ComputeConnectsNextToSkeleton(connectList);
    
}

void TPZMHMHDivErrorEstimator::VerifySolutionConsistency(TPZCompMesh* cmesh) {
//    {
//        std::ofstream outvtk("MeshToVerifyConsistency.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(), outvtk);
//        std::ofstream outtxt("MeshToVerifyConsistency.txt");
//        cmesh->Print(outtxt);
//    }

    TPZGeoMesh* gmesh = fOriginal->Reference();
    gmesh->ResetReference();
    cmesh->LoadReferences();

    int dim = gmesh->Dimension();

    int64_t nel = cmesh->NElements();
    // Iterates through all elements of the mesh
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if (!cel) continue;

        // Filters elements of highest dimension (2 or 3)
        TPZGeoEl *gel = cel->Reference();

        if (gel->Dimension() != dim) continue;

        // Iterates through the sides of the element
        int nsides = gel->NSides();
        for (int iside = 0; iside < nsides; iside++) {
            TPZGeoElSide gelside(gel, iside);

            // Filters sides of lower dimension
            if (gelside.Dimension() != dim - 1) continue;

            // Gets compel sides of equal and lower (if existing) level linked to the gelside
            TPZStack<TPZCompElSide> celstack;
            gelside.EqualLevelCompElementList(celstack, 1, 0);

            TPZCompElSide large = gelside.LowerLevelCompElementList2(1);
            if (large) celstack.Push(large);

            if (celstack.size() == 0) continue;

            int intOrder = 2;

            TPZIntPoints *intRule = gelside.CreateIntegrationRule(intOrder);

            // Iterates through the comp sides connected to the reference gelside
            int nstack = celstack.size();
            for (int ist = 0; ist < nstack; ist++) {
                TPZCompElSide cneighbour = celstack[ist];
                if (!cneighbour) continue;

                TPZGeoElSide neighbour = cneighbour.Reference();

                // Filters comp sides in elements of highest dimension (2 or 3)
                if (neighbour.Element()->Dimension() != dim) continue;

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

                    // Maps pt0 and pt1 to volume and gets solution on this points
                    TPZTransform<REAL> sideToVolume(dim, dim);
                    sideToVolume = gelside.Element()->SideToSideTransform(iside, nsides - 1);

                    TPZManVector<REAL> pt0_vol(dim, 0);
                    sideToVolume.Apply(pt0, pt0_vol);
                    TPZManVector<STATE> sol0(1);
                    cel->Solution(pt0_vol, 0, sol0);

                    TPZTransform<REAL> neighSideToVolume(dim, dim);
                    neighSideToVolume = neighbour.Element()->SideToSideTransform(cneighbour.Side(), neighbour.Element()->NSides() - 1);

                    TPZManVector<REAL> pt1_vol(dim, 0);
                    neighSideToVolume.Apply(pt1, pt1_vol);
                    TPZManVector<STATE> sol1(1);
                    cneighbour.Element()->Solution(pt1_vol, 0, sol1);

#ifdef LOG4CXX
                    if (logger->isDebugEnabled()) {
                        std::stringstream sout;
                        sout << "\nSide Element =  " << gelside.Element()->Index() << "\n";
                        sout << "Neighbour Element =  " << neighbour.Element()->Index() << "\n";
                        sout << "Side solution =  " << sol0[0] << "\n";
                        sout << "Neigh solution = " << sol1[0] << "\n";
                        sout << "Diff = " << sol1[0] - sol0[0] << "\n";
                        sout << "Side coord:  [" << x0[0] << ", " << x0[1] << ", " << x0[2] << "]\n";
                        sout << "Neigh coord: [" << x1[0] << ", " << x1[1] << ", " << x1[2] << "]\n";

                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif

                    // Checks pressure value on these nodes
                    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cneighbour.Element());
                    if (!intel) DebugStop();
                }
            }
            delete intRule;
        }
    }
}

void TPZMHMHDivErrorEstimator::ComputeConnectsNextToSkeleton(std::set<int64_t>& connectList) {

    TPZCompMesh *pressure_reconstruction = fPostProcMesh.MeshVector()[1];
    TPZGeoMesh * gmesh = pressure_reconstruction->Reference();
    int dim = gmesh->Dimension();
    gmesh->ResetReference();
    pressure_reconstruction->LoadReferences();

    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl* gel = gmesh->Element(el);
        if (!gel) continue;
        if (gel->Dimension() == dim) continue;
        if (gel->MaterialId() != this->fPressureSkeletonMatId) continue;
        TPZCompEl* cel = gel->Reference();
        if (!cel) continue;

        int nsides = gel->NSides();
        for (int iside = 0; iside < nsides; iside++) {
            TPZGeoElSide gelside(gel, iside);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                if (neighbour.Element()->Dimension() == dim) {
                    TPZCompEl * neigh = neighbour.Element()->Reference();
                    if (neigh) {
                        int sideId = neighbour.Side();
                        connectList.insert(neigh->ConnectIndex(sideId));
                    }
                }
                neighbour = neighbour.Neighbour();
            }
        }
    }
}
