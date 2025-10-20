#include "TPZPostProcessErrorSBFem.h"
#include "TPZH1ErrorHybridH1EstimateMaterial.h"
#include "TPZNullMaterialCS.h"
#include "TPZSBFemElementGroupPostProcess.h"
#include "TPZRenumbering.h"
#include "TPZMultiPhysicsMeshWindow.h"

TPZPostProcessErrorSBFem::TPZPostProcessErrorSBFem(TPZCompMesh *origin, TPZCompMesh *postprocess, TPZBuildSBFemHybrid &buildSBFemHybrid) : TPZPostProcessError(origin), fBuildSBFemHybrid(buildSBFemHybrid), fInterfaceMaterials() {
    this->fMatWrap = buildSBFemHybrid.GetSkeletonMatid();
    // Initialize members if needed
    fMeshVector[Epressure] = postprocess;
    InitializeInterfaceMaterialObjects();
    IdentifyElementGroups();
    CreatePartitionofUnityMesh();
    if(0)
    {
        std::ofstream out("partition_of_unity.txt");
        fMeshVector[Epatch]->Print(out);
    }
    AddLevel0SkeletonElements(fMeshVector[Epatch]);
    if(0)
    {
        std::ofstream out("partition_of_unity.txt");
        fMeshVector[Epatch]->Print(out);
    }
    RestrainCentralConnectofPatches();
    if(0)
    {
        std::ofstream out("partition_of_unity.txt");
        fMeshVector[Epatch]->Print(out);
    }

    
}

// Destructor
TPZPostProcessErrorSBFem::~TPZPostProcessErrorSBFem() {
    // Clean up resources if needed
}


void TPZPostProcessErrorSBFem::BuildPatchStructures2() {
    // Implementation for building patch structures specific to SBFem
    // This could involve grouping patches based on SBFem characteristics
    // For now, we will call the base class implementation
    
    // a set of all geometric cornernodes in the mesh
    std::set<int64_t> cornernodes;
    std::set<int> materialids;
    auto matidtranslation = fBuildSBFemHybrid.GetMatIdTranslation();
    for(auto it : matidtranslation)
    {
        materialids.insert(it.first);
    }

    TPZGeoMesh *gmesh = fMeshVector[Eorigin]->Reference();
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel) {
            continue;
        }
        if(gel->HasSubElement())
        {
            continue;
        }
        if(materialids.find(gel->MaterialId()) == materialids.end())
        {
            continue;
        }
        int ncorners = gel->NCornerNodes();
        for(int in=0; in<ncorners; in++)
        {
            cornernodes.insert(gel->NodeIndex(in));
        }
    }
//    std::cout << "Corner nodes size " << cornernodes.size() << std::endl;
    #include <iterator>

//    std::copy(cornernodes.begin(), cornernodes.end(),
//          std::ostream_iterator<int64_t>(std::cout, " "));
//    std::cout << std::endl;
    TPZCompMesh *partitionunity = fMeshVector[Epatch];
    TPZVec<int64_t> connecttogeonode(partitionunity->NConnects(),-1);
    {
        int64_t nel = partitionunity->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = partitionunity->Element(el);
            if(!cel) DebugStop();
            TPZGeoEl *gel = cel->Reference();
            int ncorner = gel->NCornerNodes();
            for(int ic=0; ic<ncorner; ic++) {
                int64_t geonode = gel->NodeIndex(ic);
                int64_t cindex = cel->ConnectIndex(ic);
                connecttogeonode[cindex] = geonode;
            }
        }
    }
    gmesh->ResetReference();
    // the number of patches is equal to the number of independent connects
    // which is also equal to the number of equations of the mesh
    // the patch mesh is already configured with connect dependencies
    // this includes the center node of each patch
    int64_t npatches = partitionunity->NEquations();
    TPZBlock &block = partitionunity->Block();
    partitionunity->Block().Resequence();
    fElementPatches.Resize(npatches);
    // we have the corner nodes. Find the the full set of element groups linked to the node
    // 1 : find the elements connected to each connect of the partition of unity mesh
    // the graph will only include connects with nonzero block sizes.
    TPZVec<int64_t> nodtoelgraph,nodtoelgraphindex;
    ComputeConnectToElementGraph(nodtoelgraph, nodtoelgraphindex);
    // 2 : add all the elements connected to each node to the patches
    int64_t nodegraphindexsize = nodtoelgraphindex.size();
    if(nodegraphindexsize != partitionunity->NConnects()+1) {
        DebugStop();
    }
    int64_t patchcount = 0;
    for(int64_t nod = 0; nod < nodegraphindexsize-1; nod++) {
        int64_t firstind = nodtoelgraphindex[nod];
        int64_t lastind = nodtoelgraphindex[nod+1];
        // there is no element connected
        if(lastind == firstind) continue;
        TPZConnect &c = partitionunity->ConnectVec()[nod];
        if(c.NShape() == 0) DebugStop();
        int64_t seqnum = c.SequenceNumber();
        int64_t patchindex = block.Position(seqnum);
        // invert the connect dependency to the list of connects contributing to a connect
        // this will include the root connect (with value 1)
        // if the connect has no dependency deptree will have a single entry (itself)
        // deptree will only point to connects without dependency, meaning connects with associated patches
        // the equation associated with the connect is equal to the index of the patch
        // dependent connects have equation larger than the equations of the independent connects
        std::map<int64_t,REAL> deptree;
        ComputeConnectMultiplicators(*partitionunity, nod, 1., deptree);
        if(0) {
            std::cout << "Connect index " << nod << " patch number " << patchindex << " deptree ";
            for (const auto& [key, value] : deptree) {
                std::cout << key << ":" << value << ' ';
            }
            std::cout << std::endl;
        }
        // add the contribution of all connects to the different patches
        for(auto it : deptree) {
            int64_t cindex = it.first;
            TPZConnect &cdep = partitionunity->ConnectVec()[cindex];
            int64_t seqnum = partitionunity->ConnectVec()[cindex].SequenceNumber();
            int64_t patchindex_loc = block.Position(seqnum);
            if(patchindex_loc >= npatches) {
                DebugStop();
            }
            TPZGeoPatch &patch = fElementPatches[patchindex_loc];
            patch.AddConnect(nod, it.second, connecttogeonode[nod]);
//            int64_t vecsize = patch.ConnectIndexes().size();
//            patch.ConnectIndexes().Resize(vecsize+1,nod);
//            patch.HatFunctionValues().Resize(vecsize+1, it.second);
        }
        // for connects without dependency, initialize the datastructure of geometric elements
        // contained in the domain of the patch
        if(!c.HasDependency()) {
            fElementPatches[patchindex].ElementIndexes().Resize(lastind-firstind);
            for(int64_t elind = firstind; elind < lastind; elind++) {
                int64_t element = nodtoelgraph[elind];
                TPZCompEl *cel = partitionunity->Element(element);
                TPZGeoEl *gel = cel->Reference();
                fElementPatches[patchindex].ElementIndexes()[elind-firstind] = gel->Index();
            }
        }
    }
//    PrintPatchInformation();
    {
        std::ofstream out("GeoMesh.txt");
        gmesh->Print(out);
    }
    ExpandGeoPatchMeshes();
//    PrintPatchInformation();
}
#include "TPZNullMaterial.h"


/// add coarse scale skeleton elements so that the intermediate nodes will be restrained
void TPZPostProcessErrorSBFem::AddLevel0SkeletonElements(TPZCompMesh *partition_unitiy_mesh) {
    int skeleton_matid = fBuildSBFemHybrid.GetSkeletonMatid();
    std::set<int> allmatids = fBuildSBFemHybrid.GetBoundaryMatIds();
    TPZGeoMesh *gmesh = partition_unitiy_mesh->Reference();
    // identify the relevant boundary material elements
    std::set<int64_t,int> boundary_groupindex;
    {
        int64_t nel = gmesh->NElements();
        // the groupindex of a boundary father element is
        // the groupindex of its sons
        // the groupindex of the neighbours
    }
    allmatids.insert(skeleton_matid);
    for(auto it : allmatids) {
        int dim = partition_unitiy_mesh->Dimension();
        int nstate = 1;
        TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial<STATE> (it, dim, nstate);
        partition_unitiy_mesh->InsertMaterialObject(nullmat);
    }
    const auto &translate = fBuildSBFemHybrid.GetMatIdTranslation();
    std::set<int> origmatid;
    for(auto &it : translate) origmatid.insert(it.first);
    gmesh->ResetReference();
    partition_unitiy_mesh->LoadReferences();
    int64_t nel = gmesh->NElements();
    for(int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) continue;
        int gelmatid = gel->MaterialId();
//        std::cout << "gel matid " << gelmatid << std::endl;
        if(allmatids.find(gelmatid) == allmatids.end()) continue;
        TPZGeoElSide gelside(gel);
        TPZGeoElSide volume_el = gelside.HasNeighbour(origmatid);
        if(!volume_el) continue;
        if(volume_el.HasSubElement()) continue;
//        if(!gel->HasSubElement()) continue;
//        std::cout << "generating comp element from geo index " << gel->Index() << " matid " << gelmatid << std::endl;
        partition_unitiy_mesh->ApproxSpace().CreateCompEl(gel, *partition_unitiy_mesh);
    }
    partition_unitiy_mesh->ExpandSolution();
    partition_unitiy_mesh->CleanUpUnconnectedNodes();
}

/// restrain the central node of each patch
void TPZPostProcessErrorSBFem::RestrainCentralConnectofPatches() {
    TPZCompMesh *patch = fMeshVector[Epatch];
    TPZGeoMesh *gmesh = patch->Reference();
    if(gmesh->Reference() != patch) DebugStop();
    int64_t ngroups = fElementGroupsHybrid.size();
    for(auto igr : fElementGroupsHybrid) {
        std::set<int64_t> connectids;
        auto celvec = igr->GetElGroup();
        int64_t centralconnectindex = -1;
        {
            TPZGeoEl *gel = celvec[0]->Reference();
            // identify the element of the partition of unity mesh
            TPZCompEl *cel = gel->Reference();
            if(!cel) DebugStop();
            int64_t nc = gel->NCornerNodes();
            centralconnectindex = cel->ConnectIndex(nc-1);
        }
        for(auto cel : celvec) {
            std::set<int64_t> loclist;
            TPZGeoEl *gel = cel->Reference();
            TPZCompEl *celpatch = gel->Reference();
            celpatch->BuildConnectList(loclist);
            for(auto it : loclist) {
                if(it == centralconnectindex) continue;
                TPZConnect &c = patch->ConnectVec()[it];
                if(c.NShape() == 0 || c.HasDependency()) continue;
                connectids.insert(it);
            }
        }
        if(0) {
            std::cout << "Central connect " << centralconnectindex << " surrounding con ";
            for(auto it : connectids) std::cout << it << " ";
            std::cout << std::endl;
        }
        REAL fact = 1./connectids.size();
        TPZConnect &ccentral = patch->ConnectVec()[centralconnectindex];
        TPZFNMatrix<1,REAL> rest(1, 1, fact);
        for(auto cindex : connectids) {
            TPZConnect &c = patch->ConnectVec()[cindex];
            ccentral.AddDependency(centralconnectindex, cindex, rest, 0, 0, 1, 1);
        }
    }
    patch->CleanUpUnconnectedNodes();
}

/// Expand the geometric patch data to include the skeleton elements, the interface elements and flux elements
/// These will be used to create a TPZMultiphysicsWindow object
void TPZPostProcessErrorSBFem::ExpandGeoPatchMeshes() {
    TPZCompMesh *partition_unity_mesh = fMeshVector[Epatch];
    TPZGeoMesh *gmesh = partition_unity_mesh->Reference();
    int dim = gmesh->Dimension();
    int skelmatid = fBuildSBFemHybrid.GetSkeletonMatid();
    auto intfacematidpair = fBuildSBFemHybrid.GetInterfaceMaterialIds();
    std::set<int> intfacematid = {intfacematidpair.first,intfacematidpair.second};
    int fluxmatid = fBuildSBFemHybrid.GetFluxMaterialId();
    std::set<int> bcmatids = fBuildSBFemHybrid.GetBoundaryMatIds();
    bcmatids.insert(fluxmatid);
    // loop over the patches
    for(auto &patch : fElementPatches) {
        auto & elind = patch.ElementIndexes();
        std::set<int64_t> elindices;
        elindices.insert(&elind[0],(&elind[0])+elind.size());
        // there should not be duplicate element indexes
        if(elindices.size() != elind.size()) DebugStop();
        // loop over the elements contained in the patches
        for(auto el : elind) {
            TPZGeoEl *gel = gmesh->Element(el);
            // as we are dealing with SBFem Volume elements, only along the first side of dim-1 there
            // will be fluxes
            int firstside = gel->FirstSide(dim-1);
            TPZGeoElSide gelside(gel,firstside);
            // verify if the dim-1 elements have to be included
            bool include = false;
            // loop over the cornernodes
            for(int ic = 0; ic <gelside.NSideNodes(); ic++) {
                int64_t nodeindex = gelside.SideNodeIndex(ic);
                if(patch.isInternalNode(nodeindex)) {
                    include = true;
                    break;
                }
            }
            // the dim-1 side contains and internal node
//            include = true;
            if(include) {
                
                TPZGeoElSide neighskel = gelside.Neighbour();
                if(neighskel.Element()->MaterialId() != skelmatid) DebugStop();
                elindices.insert(neighskel.Element()->Index());
                TPZGeoElSide neighintface = neighskel.Neighbour();
                int neighmatid = neighintface.Element()->MaterialId();
                if(intfacematid.find(neighmatid) == intfacematid.end()) DebugStop();
                elindices.insert(neighintface.Element()->Index());
                TPZGeoElSide neighflux = neighintface.HasNeighbour(bcmatids);
                if(!neighflux) DebugStop();
                elindices.insert(neighflux.Element()->Index());
            }
        }
        elind.Resize(elindices.size());
        auto it = elindices.begin();
        for(int i = 0; i<elindices.size(); i++, it++) {
            elind[i] = *it;
        }
    }
}



void TPZPostProcessErrorSBFem::CreateMultiphysicsMesh() {
    TPZGeoMesh *gmesh = fMeshVector[Eorigin]->Reference();
    CreateAveragePressureMesh();
    fMeshVector[Eflux] = new TPZCompMesh(gmesh);
    HideSBFemVolume();
    if(fMeshVector[Emulti]) delete fMeshVector[Emulti];
    TPZMultiphysicsCompMesh *mphysics;
    mphysics = new TPZMultiphysicsCompMesh(gmesh);
    fMeshVector[Emulti] = mphysics;
    InsertPostProcessingMaterials(mphysics);
    TPZManVector<TPZCompMesh *> meshvec(5,0);
    meshvec[0] = fMeshVector[Eflux];
    meshvec[1] = fMeshVector[Epressure];
    meshvec[2] = fMeshVector[Epatch];
    meshvec[3] = fMeshVector[Eorigin];
    meshvec[4] = fMeshVector[Epressureaverage];
    //Fazendo auto build
    TPZManVector<int> active = {0,1,0,0,1};
    mphysics->SetMeshVectorAndActiveSpaces(meshvec, active);
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->AutoBuild();
    ExposeSBFemVolume();
    
    InsertInterfaceMaterial(mphysics);
    AddInterfaceElements(mphysics);
}

// create a multiphysics window mesh
void TPZPostProcessErrorSBFem::CreateMultiphysicsWindowMesh(){
    TPZGeoMesh *gmesh = fMeshVector[Eorigin]->Reference();
    CreateAveragePressureMesh();
    fMeshVector[Eflux] = new TPZCompMesh(gmesh);
    HideSBFemVolume();
    TPZMultiPhysicsMeshWindow *mphysics = new TPZMultiPhysicsMeshWindow(gmesh);
    fMeshVector[Emulti] = mphysics;
    InsertPostProcessingMaterials(mphysics);
    TPZManVector<TPZCompMesh *> meshvec(5,0);
    meshvec[0] = fMeshVector[Eflux];
    meshvec[1] = fMeshVector[Epressure];
    meshvec[2] = fMeshVector[Epatch];
    meshvec[3] = fMeshVector[Eorigin];
    meshvec[4] = fMeshVector[Epressureaverage];
    //Fazendo auto build
    TPZManVector<int> active = {0,1,0,0,1};
    mphysics->SetMeshVectorAndActiveSpaces(meshvec, active);

}

#include "TPZLagrangeMultiplierCS.h"
/// insert post processing materials in the multiphysics mesh
void TPZPostProcessErrorSBFem::InsertPostProcessingMaterials(TPZMultiphysicsCompMesh *mphys) {
    TPZCompMesh *origin = fMeshVector[Eorigin];
    int meshdim = origin->Dimension();
    auto &matvec = origin->MaterialVec();
    for(auto matit : matvec) {
        TPZMaterial *mat = matit.second;
        TPZDarcyFlow *darcy = dynamic_cast<TPZDarcyFlow *>(mat);
        if(darcy) {
            auto newmat = new TPZH1ErrorHybridH1EstimateMaterial(*darcy);
            mphys->InsertMaterialObject(newmat);
        }
    }
    for(auto matit : matvec) {
        TPZMaterial *mat = matit.second;
        TPZBndCondT<STATE> *bc = dynamic_cast<TPZBndCondT<STATE> *>(mat);
        if(bc) {
            int bcid = bc->Id();
            auto matref = bc->Material();
            int matid = matref->Id();
            TPZMaterialT<STATE> *newmat = dynamic_cast<TPZMaterialT<STATE> *>(mphys->FindMaterial(matid));
            auto *newbc = newmat->CreateBC(newmat, bcid, bc->Type(), bc->Val1(), bc->Val2());
            int order = bc->ForcingFunctionBCPOrder();
            auto forcefunc = bc->ForcingFunctionBC();
            newbc->SetForcingFunctionBC(forcefunc, order);
            mphys->InsertMaterialObject(newbc);
        }
    }
    int skelmatid = fBuildSBFemHybrid.GetSkeletonMatid();
    int nstate = 1;
//    TPZNullMaterialCS(int matid, int dimension, int nstate)
    auto skelmat = new TPZNullMaterialCS<STATE>(skelmatid,meshdim-1,nstate);
    int fluxmatid = fBuildSBFemHybrid.GetFluxMaterialId();
    auto fluxmat = new TPZNullMaterialCS<>(fluxmatid,meshdim-1,nstate);
    mphys->InsertMaterialObject(skelmat);
    mphys->InsertMaterialObject(fluxmat);
}


void TPZPostProcessErrorSBFem::CreateAveragePressureMesh() {
    TPZPostProcessError::CreateAveragePressureMesh();
    TPZCompMesh *meshpress = fMeshVector[Epressureaverage];
    int meshdim = meshpress->Reference()->Dimension();
    int64_t nel = meshpress->NElements();
    for(int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = meshpress->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int dim = gel->Dimension();
        if(dim != meshdim) {
            delete cel;
        } else if (cel->NConnects() == 1) {
            cel->SetConnectIndex(0, 0);
        }
    }
    meshpress->ComputeNodElCon();
    meshpress->CleanUpUnconnectedNodes();
}


/// Identify the element groups in the H1 and Hybrid H1 mesh
void TPZPostProcessErrorSBFem::IdentifyElementGroups() {
    {
        TPZCompMesh *cmesh = fMeshVector[Eorigin];
        int64_t nel = cmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
            if(elgr) fElementGroupsH1.Push(elgr);
        }
    }
    {
        TPZCompMesh *cmesh = fMeshVector[Epressure];
        int64_t nel = cmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
            if(elgr) fElementGroupsHybrid.Push(elgr);
        }
    }
}

/// Expose the SBFemVolume elements and hide the element groups
/// this method is necessary to build the multiphysics mesh
void TPZPostProcessErrorSBFem::HideSBFemVolume() {
    for(auto elgr : fElementGroupsH1) {
        auto &elvec = elgr->GetElGroup();
        for(auto elsbfem : elvec) {
            int64_t elindex = elsbfem->Index();
            if(elsbfem->Mesh() != fMeshVector[Eorigin]) DebugStop();
            fMeshVector[Eorigin]->ElementVec()[elindex] = elsbfem;
        }
        int64_t grindex = elgr->Index();
        if(elgr->Mesh() != fMeshVector[Eorigin]) DebugStop();
        fMeshVector[Eorigin]->ElementVec()[grindex] = 0;
    }
    for(auto elgr : fElementGroupsHybrid) {
        auto &elvec = elgr->GetElGroup();
        for(auto elsbfem : elvec) {
            int64_t elindex = elsbfem->Index();
            if(elsbfem->Mesh() != fMeshVector[Epressure]) DebugStop();
            fMeshVector[Epressure]->ElementVec()[elindex] = elsbfem;
        }
        int64_t grindex = elgr->Index();
        if(elgr->Mesh() != fMeshVector[Epressure]) DebugStop();
        fMeshVector[Epressure]->ElementVec()[grindex] = 0;
    }
}

void TPZPostProcessErrorSBFem::ExposeSBFemVolume() {
    for(auto elgr : fElementGroupsH1) {
        auto &elvec = elgr->GetElGroup();
        for(auto elsbfem : elvec) {
            int64_t elindex = elsbfem->Index();
            if(elsbfem->Mesh() != fMeshVector[Eorigin]) DebugStop();
            fMeshVector[Eorigin]->ElementVec()[elindex] = 0;
        }
        int64_t grindex = elgr->Index();
        if(elgr->Mesh() != fMeshVector[Eorigin]) DebugStop();
        fMeshVector[Eorigin]->ElementVec()[grindex] = elgr;
    }
    for(auto elgr : fElementGroupsHybrid) {
        auto &elvec = elgr->GetElGroup();
        for(auto elsbfem : elvec) {
            int64_t elindex = elsbfem->Index();
            if(elsbfem->Mesh() != fMeshVector[Epressure]) DebugStop();
            fMeshVector[Epressure]->ElementVec()[elindex] = 0;
        }
        int64_t grindex = elgr->Index();
        if(elgr->Mesh() != fMeshVector[Epressure]) DebugStop();
        fMeshVector[Epressure]->ElementVec()[grindex] = elgr;
    }
}

/// put the multiphysics elements corresponding to TPZSBFemVolume into postprocess groups
void TPZPostProcessErrorSBFem::GroupSBFemMultiphysics(TPZMultiphysicsCompMesh *mphys, bool hasboundary) {
    std::map<TPZSBFemElementGroup *, std::set<TPZCompEl *> > AtomicToMPhysis;
    int64_t nel = mphys->NElements();
    for(int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = mphys->Element(el);
        if(!cel) continue;
        TPZMultiphysicsElement *mcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if(!mcel) continue;
        TPZCompEl *press = mcel->Element(1);
        TPZSBFemVolume *sbvol = dynamic_cast<TPZSBFemVolume *>(press);
        if(sbvol) {
            TPZSBFemElementGroup *elgr = sbvol->ElementGroup();
            AtomicToMPhysis[elgr].insert(cel);
        }
    }
    // we create a group in the multiphysics mesh that will compute the siffness of SBFem and the
    // right hand side of hybrid H1 reconstruction
    for(auto &it : AtomicToMPhysis) {
        std::set<TPZCompEl *> &grp = it.second;
        TPZSBFemElementGroup *origgroup = it.first;
        TPZSBFemElementGroupPostProcess *grppp = new TPZSBFemElementGroupPostProcess(*mphys,origgroup);
        grppp->SetHasBoundary(hasboundary);
        for(auto sbvol : grp) {
            grppp->AddElement(sbvol);
        }
    }
}

/// @brief identify connect multiplicator values
/// if the connect has no dependency, the multiplicator value is 1
/// other wise it will be a list of connects and scalar values
///  the sum of the multiplicator values should always be one
void TPZPostProcessErrorSBFem::ComputeConnectMultiplicators(TPZCompMesh &cmesh, int64_t cindex, REAL depval, std::map<int64_t,REAL> &mult) {
    TPZConnect &c = cmesh.ConnectVec()[cindex];
    if(c.HasDependency()) {
        TPZConnect::TPZDepend<STATE> *dep = dynamic_cast<TPZConnect::TPZDepend<STATE> *>(c.FirstDepend());
        std::map<int64_t,REAL> depmult;
        while(dep) {
            STATE depvalloc = dep->fDepMatrix(0,0);
            ComputeConnectMultiplicators(cmesh, dep->fDepConnectIndex, depval*depvalloc, depmult);
            dep = dynamic_cast<TPZConnect::TPZDepend<STATE> *>(dep->fNext);
        }
        REAL sum = 0.;
        for(auto it : depmult) {
            mult[it.first] += it.second;
            sum += it.second;
        }
        if(!IsZero(sum - depval)) {
            
        }
    } else {
        mult[cindex] = depval;
    }
}

/// Compute the element graph as a function of connect indices for the partition of unity mesh
void TPZPostProcessErrorSBFem::ComputeElementGraph(TPZStack<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex) {
    
    TPZCompMesh &cmesh = *fMeshVector[Epatch];
    int dim = cmesh.Dimension();
    int64_t nel = cmesh.NElements();
    elgraphindex.Resize(nel+1);
    elgraphindex[0] = 0;
    for (int el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh.Element(el);
        if(!cel) DebugStop();
        if(cel->Reference()->Dimension() != dim) {
            elgraphindex[el+1] = elgraph.size();
            continue;
        }
        std::set<int64_t> connects;
        cel->BuildConnectList(connects);
        for(auto it : connects) {
            TPZConnect &c = cmesh.ConnectVec()[it];
            if(c.NShape() == 0) continue;
            elgraph.Push(it);
        }
        if(0) {
            std::cout << "el " << el << " cindexes ";
            for(int64_t i = elgraphindex[el]; i< elgraph.size(); i++) std::cout << elgraph[i] << " ";
            std::cout << std::endl;
        }
        elgraphindex[el+1] = elgraph.size();
    }
}

/// compute Connect to element graph
void TPZPostProcessErrorSBFem::ComputeConnectToElementGraph(TPZVec<int64_t> &connectgraph, TPZVec<int64_t> &connectgraphindex) {
    
    TPZStack<int64_t> elgraph;
    TPZVec<int64_t> elgraphindex;
    ComputeElementGraph(elgraph, elgraphindex);
    TPZCompMesh *cmesh = fMeshVector[Epatch];
    int64_t nconnects = cmesh->NConnects();
    TPZRenumbering renum(elgraphindex.NElements() -1 ,nconnects);
    renum.SetElementGraph(elgraph,elgraphindex);
    renum.NodeToElGraph(elgraph,elgraphindex,connectgraph, connectgraphindex);

}

// print the relevant information of the patches
void TPZPostProcessErrorSBFem::PrintPatchInformation(std::ostream &out) {
    out << "Number of patches " << fElementPatches.size() << std::endl;
    int64_t count = 0;
    for(auto &it : fElementPatches) {
        out << "Patch index " << count << std::endl;
        it.Print(out);
        count++;
    }
}

// plot the relevant patch information to a plot file
void TPZPostProcessErrorSBFem::PlotPatches(const std::string rootname) {
    int64_t count = 0;
    TPZMultiphysicsCompMesh *multi = dynamic_cast<TPZMultiphysicsCompMesh *>(fMeshVector[Emulti]);
    for(auto &it : fElementPatches) {
        it.PlotPatch(rootname, count, *multi);
        count++;
    }
    HideInterfaceMaterial(multi);
}

#include "pzfstrmatrix.h"
#include "TPZVTKGenerator.h"


void TPZPostProcessErrorSBFem::ComputeElementErrors(TPZVec<STATE> &errors) {
    /// @brief Compute the difference between two meshes
    /// @param cmesh1 reference to the first computational mesh
    /// @param cmesh2 reference to the second computational mesh
//    void CompareMeshes(TPZCompMesh &cmesh1, TPZCompMesh &cmesh2, TPZVec<REAL> &errors)
//    {
    if(!fMeshVector[Emulti]) DebugStop();
    TPZCompMesh &cmeshmulti = *fMeshVector[Emulti];
    std::set<int> matpostprocess = fBuildSBFemHybrid.GetMaterialIds();
    cmeshmulti.EvaluateError(false, errors, matpostprocess);


    int matid = *matpostprocess.begin();
    TPZH1ErrorHybridH1EstimateMaterial *posmat = dynamic_cast<TPZH1ErrorHybridH1EstimateMaterial *>(cmeshmulti.FindMaterial(matid));
    if(!posmat) DebugStop();
    if(errors.size() != posmat->NEvalErrors()) {
        errors.Resize(posmat->NEvalErrors());
    }
    int64_t NErrors = errors.size();
    errors.Fill(0.);
    int64_t nel = cmeshmulti.NElements();
    cmeshmulti.ElementSolution().Redim(nel, NErrors);
    TPZFMatrix<REAL> &elementsolution = cmeshmulti.ElementSolution();
    int H1Pos = TPZH1ErrorHybridH1EstimateMaterial::EH1;
    int EstPos = TPZH1ErrorHybridH1EstimateMaterial::EEstimate;
    int HybridPos = TPZH1ErrorHybridH1EstimateMaterial::EHybrid;

    for(int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmeshmulti.Element(el);
        if(!cel) continue;
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if(mfcel) {
            TPZGeoEl *gel = mfcel->Reference();
            if(gel) {
                int matid = gel->MaterialId();
                if(matpostprocess.find(matid) != matpostprocess.end()) {
                    TPZManVector<REAL> elerrors(NErrors,0.);
                    mfcel->EvaluateError(elerrors, false);
                    elementsolution(el,0) = elerrors[H1Pos];
                    elementsolution(el,1) = elerrors[HybridPos];
                    elementsolution(el,2) = elerrors[EstPos];
                    // storing the effectivity index
                    if(elerrors[H1Pos] > 1.e-8) elementsolution(el,3) = elementsolution(el,2)/elementsolution(el,0);
                    for (int i = 0; i < NErrors; i++) {
                        errors[i] += elerrors[i]*elerrors[i];
                    }
                }
            }
        }
    }
    for (int i = 0; i < NErrors; i++) {
        errors[i] = sqrt(errors[i]);
    }
    // storing the effectivity index in the 4th position
    errors[3] = errors[EstPos]/errors[H1Pos];

}

// compute the estimated H1 seminorm errors
void TPZPostProcessErrorSBFem::ReconstructHybridH1() {
    if(fMeshVector[Emulti]) {
        HideInterfaceMaterial(fMeshVector[Emulti]);
        delete fMeshVector[Emulti];
//        fInterfaceMaterials.Resize(0);
        fMeshVector[Emulti] = 0;
    }

//    std::cout << "Creating multiphysics mesh\n";
    
    CreateMultiphysicsWindowMesh();
    TPZMultiPhysicsMeshWindow *meshw = dynamic_cast<TPZMultiPhysicsMeshWindow *>(fMeshVector[Emulti]);

//    std::cout << "mesh created\n";

    TPZCompMesh *hybridsbfem = fMeshVector[Epressure];
    TPZFMatrix<STATE> &hybridsol = hybridsbfem->Solution();
    int64_t hybrid_neq = hybridsbfem->Solution().Rows();
    TPZFMatrix<STATE> accumulateSolution(hybrid_neq,1,0.);
    hybridsbfem->ComputeNodElCon();
    TPZGeoMesh *gmesh = meshw->Reference();
    if(!meshw) DebugStop();
    std::set<int> boundarymatids = fBuildSBFemHybrid.GetBoundaryMatIds();
    std::set<int> eliminat;
    for(auto it : boundarymatids) {
        TPZMaterial *mat = meshw->FindMaterial(it);
        TPZBndCondT<STATE> *bnd = dynamic_cast<TPZBndCondT<STATE> *>(mat);
        if(!bnd) DebugStop();
        if(bnd->Type() == 1) eliminat.insert(mat->Id());
    }
    for(auto it : eliminat) boundarymatids.erase(it);
    TPZCompMesh *patchmesh = fMeshVector[Epatch];
    {
        TPZFMatrix<STATE> &sol = patchmesh->Solution();
        sol.Zero();
    }
    int fluxmatid = fBuildSBFemHybrid.GetFluxMaterialId();
    int64_t npatches = fElementPatches.size();
//    std::cout << "Post processing\n";
    extern std::complex<STATE> integrateF;
    integrateF = 0.;
//    std::cout << __PRETTY_FUNCTION__ << "********************************** Remove this\n";
//    npatches = 1;
    for(int64_t ip = 0; ip < npatches; ip++) {
        if(ip%20 == 0) std::cout << "*";
//        std::cout << "Processing patch " << ip << std::endl;
        meshw->BuildMultiphysicsSpace(fElementPatches[ip].ElementIndexes());
        fElementPatches[ip].LoadPatchVelues(patchmesh);
        bool hasboundary = fElementPatches[ip].HasBoundary(*gmesh, boundarymatids);
        // eliminate the average pressure connect
        if(hasboundary) {
            std::map<int64_t,int64_t> connectmap;
            meshw->BuildConnectMap(4, connectmap);
            int64_t mphysconnectindex = connectmap[0];
            TPZConnect &c = meshw->ConnectVec()[mphysconnectindex];
            c.SetNShape(0);
            int64_t seqnum = c.SequenceNumber();
            meshw->Block().Set(seqnum, 0);
            meshw->ExpandSolution();
        }
        InsertInterfaceMaterial(meshw);
        AddInterfaceElements(meshw);
        if(0)
        {
            std::ofstream out("MultiphysicsWindow.txt");
            meshw->Print(out);
        }
        GroupSBFemMultiphysics(meshw,hasboundary);
        if(0)
        {
            std::ofstream out("MultiphysicsWindow.txt");
            meshw->Print(out);
        }
        std::set<int64_t> exclude_eq;
        {
            int64_t nel = meshw->NElements();
            for(int64_t el = 0; el<nel; el++) {
                TPZCompEl *cel = meshw->Element(el);
                if(!cel) continue;
                TPZGeoEl *gel = cel->Reference();
                if(!gel) continue;
                if(gel->MaterialId() == fluxmatid) {
                    if(cel->NConnects() != 1) DebugStop();
                    TPZMultiphysicsElement *mphysel = dynamic_cast<TPZMultiphysicsElement *>(cel);
                    TPZCompEl *fluxel = mphysel->Element(1);
                    if(fluxel->NConnects() != 1) DebugStop();
                    int nelconmphys = cel->Connect(0).NElConnected();
                    int nelconhybrid = fluxel->Connect(0).NElConnected();
                    if(nelconmphys != nelconhybrid) {
//                        std::cout << "gel " << gel->Index() << " is patch boundary " << nelconmphys
//                        << " " << nelconhybrid << std::endl;
                        TPZConnect &c = mphysel->Connect(0);
                        int64_t seqnum = c.SequenceNumber();
                        int64_t pos = meshw->Block().Position(seqnum);
                        int blockdim = meshw->Block().Size(seqnum);
                        for(int64_t i = pos; i<pos+blockdim; i++) exclude_eq.insert(i);
                    }
                }
            }
        }
        int64_t nexclude = exclude_eq.size();
        int64_t neq = meshw->NEquations() - nexclude;
//        std::cout << "Patch created now assembling\n";
        TPZFStructMatrix<STATE> strmat(meshw);
        strmat.EquationFilter().ExcludeEquations(exclude_eq);
        TPZFMatrix<STATE> stiff(neq,neq,0), rhslarge(neq+nexclude,1,0.), rhs(neq,1,0.);
        
        extern std::complex<STATE> integrateF;
        integrateF = 0.;
        
        strmat.Assemble(stiff, rhslarge);
        
//        std::cout << "Integrated residual for hat function " << integrateF << std::endl;
        integrateF = 0.;
        strmat.EquationFilter().Gather(rhslarge, rhs);
//        std::cout << "Assemble finished, now inverting the system of " << stiff.Rows() << " equations\n";
        stiff.SolveDirect(rhs, ELDLt);
//        std::cout <<  "Inversion finished - computing the error\n";
        strmat.EquationFilter().Scatter(rhs, rhslarge);
        meshw->LoadSolution(rhslarge);
        // eliminate the average pressure connect
        if(hasboundary) {
            std::map<int64_t,int64_t> connectmap;
            meshw->BuildConnectMap(4, connectmap);
            int64_t mphysconnectindex = connectmap[0];
            TPZConnect &c = meshw->ConnectVec()[mphysconnectindex];
            c.SetNShape(1);
            int64_t seqnum = c.SequenceNumber();
            meshw->Block().Set(seqnum, 1);
            meshw->ExpandSolution();
        }
        hybridsol.Zero();
        meshw->TransferMultiphysicsSolution();
        for(int64_t eq = 0; eq < hybrid_neq; eq++) {
            accumulateSolution(eq,0) += hybridsol(eq,0);
        }
        if(1)
        {
            std::string plotname("PatchReconstruct");
            TPZStack<std::string> fields;
            fields.Push("Partition");
            fields.Push("DistFlux");
            fields.Push("Pressure");
            fields.Push("SolH1Hat");
            TPZVTKGenerator vtk(meshw, fields, plotname, 3);
            int step = (int)ip;
            vtk.SetStep(step);
            vtk.Do();

        }
        fElementPatches[ip].ZeroPatchValues(patchmesh);
        meshw->CleanElementsConnects();
        HideInterfaceMaterial(meshw);
//        std::cout << "Patch finished processing\n";
    }
    std::cout << std::endl;
//    accumulateSolution.Print(std::cout);
    hybridsol = accumulateSolution;
    ExposeSBFemVolume();
    hybridsbfem->LoadSolution(hybridsol);
    {
        std::string plotname("Reconstructed");
        TPZStack<std::string> fields;
        fields.Push("Pressure");
        TPZVTKGenerator vtk(hybridsbfem, fields, plotname, 3);
//        int step = (int) npatches;
        vtk.SetStep(0);
        vtk.Do();
    }
}

#include "TPZMultiphysicsInterfaceEl.h"

/// Add the Interface elements to the multiphysics mesh
void TPZPostProcessErrorSBFem::AddInterfaceElements(TPZMultiphysicsCompMesh *mfmesh) {
    mfmesh->LoadReferences();
    int64_t nel = mfmesh->NElements();
    int skeletonmatid = this->fBuildSBFemHybrid.GetSkeletonMatid();
    std::set<int> fluxmatids = this->fBuildSBFemHybrid.GetBoundaryMatIds();
    fluxmatids.insert(fBuildSBFemHybrid.GetFluxMaterialId());
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = mfmesh->Element(el);
        if(!cel) DebugStop();
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        if(gel->MaterialId() == skeletonmatid) {
            TPZGeoElSide skelside(gel);
            TPZCompElSide skelsideC = skelside.Reference();
            TPZGeoElSide intfaceside = skelside.Neighbour();
            
            TPZGeoElSide fluxside = intfaceside.HasNeighbour(fluxmatids);
            TPZCompElSide fluxsideC = fluxside.Reference();
            new TPZMultiphysicsInterfaceElement(*mfmesh,intfaceside.Element(),skelsideC,fluxsideC);
            
        }
    }
}

/// create and initialize the interface material objects
void TPZPostProcessErrorSBFem::InitializeInterfaceMaterialObjects() {
    fInterfaceMaterials.Resize(2, 0);
    int meshdim = fMeshVector[Eorigin]->Dimension();
    auto interfacematids = fBuildSBFemHybrid.GetInterfaceMaterialIds();
    TPZLagrangeMultiplierCS<STATE> *lagrange = new TPZLagrangeMultiplierCS<STATE>(interfacematids.first, meshdim-1, 1);
    lagrange->SetLinear(true);
    fInterfaceMaterials[0] = lagrange;
    lagrange = new TPZLagrangeMultiplierCS<STATE>(interfacematids.second, meshdim-1, 1);
    lagrange->SetMultiplier(-1.);
    lagrange->SetLinear(true);
    fInterfaceMaterials[1] = lagrange;
//    std::cout << "Create Interface matids " << fInterfaceMaterials[0]->Id() << " " << fInterfaceMaterials[1]->Id() << std::endl;
}

/// HideInterfaceMaterial objects
void TPZPostProcessErrorSBFem::HideInterfaceMaterial(TPZCompMesh *cmesh) {
    for(auto it : fInterfaceMaterials) {
        int matid = it->Id();
        if(!cmesh->FindMaterial(matid)) {
            std::cout << "Material id not found in mesh matid = " << matid << std::endl;
            DebugStop();
        }
        cmesh->MaterialVec().erase(matid);
    }
//    std::cout << "Hide Interface matids " << fInterfaceMaterials[0]->Id() << " " << fInterfaceMaterials[1]->Id() << std::endl;

}

/// Insert the interface material objects
void TPZPostProcessErrorSBFem::InsertInterfaceMaterial(TPZCompMesh *cmesh) {
//    std::cout << "insert Interface matids " << fInterfaceMaterials[0]->Id() << " " << fInterfaceMaterials[1]->Id() << std::endl;
    for(auto it : fInterfaceMaterials) {
        cmesh->InsertMaterialObject(it);
    }
}

