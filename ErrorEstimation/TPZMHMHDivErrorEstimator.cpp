
#include "TPZMHMHDivErrorEstimator.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#include "TPZVTKGeoMesh.h"
#include "pzintel.h"
#include "pzsubcmesh.h"
#include <Material/TPZLagrangeMultiplierCS.h>
#include <Mesh/TPZCompMeshTools.h>
#include <Mesh/TPZGeoElSideAncestors.h>
#include <Mesh/pzmultiphysicscompel.h>
#include <Mesh/TPZMultiphysicsInterfaceEl.h>
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <Elasticity/TPZMixedElasticityND.h>
#include "Projection/TPZL2Projection.h"
#include "TPZHDivErrorEstimateElasticityMaterial.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("HDivErrorEstimator"));
#endif

// a method for generating the hybridized multiphysics post processing mesh
template <typename MixedMaterial>
void TPZMHMHDivErrorEstimator<MixedMaterial>::CreatePostProcessingMesh()
{
    // initialize the post processing mesh
    this->fPostProcMesh.SetReference(this->fOriginal->Reference());

    this->fOriginal->CopyMaterials(this->fPostProcMesh);
    // Switch the material from mixed to TPZMHMHDivErrorEstimationMaterial
    this->SwitchMaterialObjects();

    TPZManVector<TPZCompMesh *> meshvec(4);
    TPZManVector<int,4> active(4,0);
    active[1] = 1;

    meshvec[0] = 0;
    meshvec[1] = CreatePrimalMesh();
    meshvec[2] = this->fOriginal->MeshVector()[0];
    meshvec[3] = this->fOriginal->MeshVector()[1];

    if (this->fPostProcesswithHDiv) {
        meshvec[0] = CreateHDivMesh();
        active[0] = 1;
    }

    CreateSkeletonElements(meshvec[1]);
    CreateSkeletonApproximationSpace(meshvec[1]);

    // If we reconstruct in H(div) we need to create an additional skeleton for the multiphysics interfaces
    if (this->fPostProcesswithHDiv) {
        CreateFluxSkeletonElements(meshvec[0]);
        // TODO move this to InsertMaterials method
        TPZNullMaterialCS<> *skeletonMat = new TPZNullMaterialCS<>(this->fPrimalSkeletonMatId);
        skeletonMat->SetDimension(this->GMesh()->Dimension() - 1);
        this->fPostProcMesh.InsertMaterialObject(skeletonMat);

        TPZNullMaterialCS<> *wrapMat = new TPZNullMaterialCS<>(this->fHDivWrapMatId);
        skeletonMat->SetDimension(this->GMesh()->Dimension() - 1);
        this->fPostProcMesh.InsertMaterialObject(wrapMat);
    }

    //enriquecer no MHM tbem?
    this->IncreasePrimalSideOrders(meshvec[1]);//malha da pressao
    if(this->fPostProcesswithHDiv) {
        this->IncreaseSideOrders(meshvec[0]);//malha do fluxo
    }

    this->RemoveMaterialObjects(this->fPostProcMesh.MaterialVec());
    this->fPostProcMesh.ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;
    this->fPostProcMesh.BuildMultiphysicsSpace(active, meshvec);

    if(this->fPostProcesswithHDiv) {
        // Create multiphysics interface computational elements
        CreateMultiphysicsInterfaces();

        TPZCompMesh* fluxMesh = this->fPostProcMesh.MeshVector()[0];
        fluxMesh->ComputeNodElCon();
        for (int i = 0; i < fluxMesh->NElements(); i++) {
            TPZCompEl *cel = fluxMesh->Element(i);
            if (!cel) continue;
            TPZInterpolatedElement * intel = dynamic_cast<TPZInterpolatedElement*>(cel);
            if (!intel) continue;

            TPZGeoEl * gel = cel->Reference();
            if (!gel) continue;

            if (gel->MaterialId() != fHDivWrapMatId) continue;
            int gelid = gel->Index();

            for (int iside = 0; iside < gel->NSides(); iside++) {

                TPZGeoElSide side(gel, iside);

                for (int ic = 0; ic < intel->NSideConnects(iside);ic++) {
                    TPZConnect &c = intel->SideConnect(ic, iside);
                    if (c.NElConnected() != 2) DebugStop();
                }
            }
        }
    }

    SubStructurePostProcessingMesh();
    this->ComputePrimalWeights();

    {
        std::ofstream file("GmeshAfterCreatePostProcessing.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(this->GMesh(), file);
    }
}

// a method for transferring the multiphysics elements in submeshes
template <typename MixedMaterial>
void TPZMHMHDivErrorEstimator<MixedMaterial>::SubStructurePostProcessingMesh()
{
    TPZGeoMesh *gmesh = this->GMesh();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    this->fOriginal->LoadReferences();
    int64_t ngel = gmesh->NElements();

    // associate with each geometric element a subcmesh object
    // we do this because the original mesh and post processing mesh share the same geometric mesh
    // if an element of the original mesh is in a given subcmesh then the corresponding element in the
    // post processing mesh will be put in the corresponding subcmesh
    TPZVec<TPZSubCompMesh *> orig_submesh_per_gel(ngel,nullptr);
    for (int64_t iel = 0; iel < ngel; iel++) {
        TPZGeoEl *gel = gmesh->Element(iel);
        if(!gel) continue;
        TPZCompEl *cel = gel->Reference();
        if(!cel) continue;
        TPZCompMesh *mesh = cel->Mesh();
        TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *>(mesh);
        if (submesh) {
            orig_submesh_per_gel[iel] = submesh;
        }
    }

    // create new submeshes and store a pointer of the submeshes to which the compels in 'fPostProcMesh' belongs.
    // the map 'submeshmap' is create to store the original submesh and its correspondent in 'fPostProcMesh'.
    std::map<TPZSubCompMesh *,TPZSubCompMesh *> submeshmap;
    int64_t nel = this->fPostProcMesh.NElements();
    TPZVec<TPZSubCompMesh *> new_submesh_per_cel(nel,nullptr);

    // this block creates subcmeshs in fPostProcMesh without association to any elements
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = this->fPostProcMesh.Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();

        if (gel->Dimension() != dim) continue;
        // we are only building the data structure for elements with mesh dimension
        // this excludes interface elements, wrappers, lagrange multipliers etc
        TPZSubCompMesh *submesh_orig = orig_submesh_per_gel[gel->Index()];
        if(!submesh_orig) continue;
        auto iter = submeshmap.find(submesh_orig);
        if (iter == submeshmap.end()) {
            TPZSubCompMesh *new_submesh = new TPZSubCompMesh(this->fPostProcMesh);
            submeshmap[submesh_orig] = new_submesh;
            new_submesh_per_cel[el] = new_submesh;
        } else {
            new_submesh_per_cel[el] = iter->second;
        }
    }

    this->fPostProcMesh.ComputeNodElCon();
    // TODO  refactor here
    std::cout << "Transferring volumetric elements...\n";
    if(1) {

        TPZVec<TPZSubCompMesh *> connectToSubcmesh(this->fPostProcMesh.NConnects(), nullptr);
        // transfer the elements in the submesh indicated by elementgroup
        // associate the submesh with the volumetric elements
        for (int64_t iel = 0; iel < nel; iel++) {
            TPZCompEl *cel = this->fPostProcMesh.Element(iel);
            if(!cel) continue;

            TPZGeoEl *gel = cel->Reference();
            if(!gel) DebugStop();
            if(new_submesh_per_cel[iel] == nullptr) {
                //std::cout << "Not transfering compel: " << iel << ", dim" << gel->Dimension()
                //          << ", matid: " << gel->MaterialId() << '\n';
                continue;
            }
            TPZSubCompMesh *submesh = new_submesh_per_cel[iel];

            TPZStack<int64_t> connectlist;
            cel->BuildConnectList(connectlist);
            auto nc = connectlist.size();
            for (int ic = 0; ic < nc; ic++) {
                auto conindex = connectlist[ic];
                connectToSubcmesh[conindex] = submesh;
            }
            
            //std::cout << "Transfer compel: " << iel << ", dim" << gel->Dimension() << ", matid: " << gel->MaterialId() <<  ' ';
            //for(int i=0; i< nc; i++) std::cout << connectlist[i] << " ";
            //std::cout  << '\n';
            submesh->TransferElement(&(this->fPostProcMesh), iel);
        }

        std::cout << "Transferring 'dim -1' elements...\n";
        for (int64_t iel = 0; iel < nel; iel++) {
            TPZCompEl *cel = this->fPostProcMesh.Element(iel);
            if(!cel) continue;

            TPZGeoEl *gel = cel->Reference();
            if(!gel) {
                continue;
            }

            std::set<TPZSubCompMesh*> celDomain;
            int nc = cel->NConnects();
            for (int ic = 0; ic < nc; ic++) {
                auto conindex = cel->ConnectIndex(ic);
                if (connectToSubcmesh[conindex] != nullptr) {
                    celDomain.insert(connectToSubcmesh[conindex]);
                }
            }

#ifdef ERRORESTIMATION_DEBUG
            if (celDomain.size() > 1) DebugStop();
#endif

            if (celDomain.size() == 1) {
                TPZSubCompMesh *submesh = *celDomain.begin();
                //std::cout << "Transfer element: " << iel << ", gel " << gel->Index()
                //          << ", matId: " << gel->MaterialId();
                //if (nc) std::cout << " connect index " << cel->ConnectIndex(0);

                //std::cout << '\n';
                submesh->TransferElement(&(this->fPostProcMesh), iel);
                // log: elementos e matId
            }
        }
        std::cout << "Finished transfering dim-1 elems\n";
        std::cout << "Printing gmeshsub\n";

        std::ofstream file("GmeshSub.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);

        std::cout << "Finished  printing gmesh\n";

        this->fPostProcMesh.ComputeNodElCon();
        if (!this->fPostProcesswithHDiv) {
            // TODO increment el connected of connects next to the skeleton
            std::set<int64_t> connectlist;
            ComputeConnectsNextToSkeleton(connectlist);
            for (auto it : connectlist) {
                this->fPostProcMesh.ConnectVec()[it].IncrementElConnected();
            }
        }

        for(auto iter : submeshmap)
        {
            iter.second->MakeAllInternal();

            if (this->fPostProcesswithHDiv) {
                iter.second->ComputeNodElCon();
                bool keeplagrangian = true;
                bool keepmatrix = false;
                TPZCompMeshTools::CreatedCondensedElements(iter.second, keeplagrangian, keepmatrix);
            }
        }

    } else {
        std::set<int64_t> connectlist;
        ComputeConnectsNextToSkeleton(connectlist);

        {
            std::ofstream out("MalhaTesteBeforeTransfer.txt");
            this->fPostProcMesh.Print(out);
        }

        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl* cel = this->fPostProcMesh.Element(el);
            if (!cel) continue;
            TPZGeoEl* gel = cel->Reference();
            if (!gel) DebugStop();
            TPZSubCompMesh* submesh = new_submesh_per_cel[el];
            if (!submesh)continue;

            submesh->TransferElement(&(this->fPostProcMesh), el);
            cel = this->fPostProcMesh.Element(el);
            if (cel) DebugStop();
        }

        this->fPostProcMesh.ComputeNodElCon();

        for (auto it : connectlist) {
            this->fPostProcMesh.ConnectVec()[it].IncrementElConnected();
        }

        for (auto iter : submeshmap) {
            iter.second->ExpandSolution();
        }

        for (auto iter : submeshmap) {
            iter.second->MakeAllInternal();
        }

    }

    this->fPostProcMesh.ComputeNodElCon();
    this->fPostProcMesh.CleanUpUnconnectedNodes();

    // set an analysis type for the submeshes
    {
        int64_t nel = this->fPostProcMesh.NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = this->fPostProcMesh.Element(el);
            TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
            if (sub) {

                if (this->fPostProcesswithHDiv) {
                    sub->CleanUpUnconnectedNodes();
                }
                int numthreads = 0;
                int preconditioned = 0;

                sub->SetAnalysisSkyline(numthreads, preconditioned);
            }
        }
    }
    std::cout << "Finished substructuring post proc mesh\n";
}

template <typename MixedMaterial>
TPZCompMesh *TPZMHMHDivErrorEstimator<MixedMaterial>::CreateHDivMesh()
{
    TPZCompMesh *OrigFlux = this->fOriginal->MeshVector()[0];
    TPZGeoMesh *gmesh = OrigFlux->Reference();
    gmesh->ResetReference();
    TPZCompMesh *fluxmesh = OrigFlux->Clone();
    RemoveMaterialObjects(fluxmesh->MaterialVec());
    int dim = gmesh->Dimension();
    // creating a copy of all flux elements except for flux elements for the skeleton
    // this will generate a consistent H(div) space over the complete domain
    // the fluxes should be discontinuous between subdomains?
    int64_t nel = fluxmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fluxmesh->Element(el);
        if(!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        TPZGeoEl *gel = cel->Reference();
        bool isbcmat = fMHM->fMaterialBCIds.find(gel->MaterialId()) != fMHM->fMaterialBCIds.end();
        // if the element is of lower dimension and is not a boundary
        // don't create a flux element
        TPZMaterial *mat = fluxmesh->FindMaterial(gel->MaterialId());
        if(!mat)
        {
            // this deletes the skeleton elements
            delete cel;
            continue;
        }
        if(gel->Dimension() != dim && !isbcmat) DebugStop();
        if(gel->Dimension() != dim) continue;
        /// equate the order of the connects with the order of the original mesh
        for (int is = gel->NCornerNodes(); is<gel->NSides(); is++) {
            // we are only interested in sides of dimension dim-1
            if(gel->SideDimension(is) != dim-1) continue;
            TPZGeoElSide gelside(gel,is);
            TPZGeoElSideAncestors gelsideancestor(gelside);
            // if the element is not neighbour of the skeleton mesh, continue
            if(!gelside.HasNeighbour(fMHM->fSkeletonMatId) && !gelsideancestor.HasLarger(fMHM->fSkeletonMatId)) continue;
            
            int nside_connects = intel->NSideConnects(is);
            // if the number of side connects is zero do nothing
            if(nside_connects != 1) DebugStop();
            // get the last connect of the side (in the case of hdiv there is always a single connect)
            TPZConnect &corig = intel->SideConnect(nside_connects-1, is);
            if(!corig.HasDependency()) DebugStop();
            corig.RemoveDepend();
        }
    }
    fluxmesh->ExpandSolution();
    return fluxmesh;
}

template <typename MixedMaterial>
TPZCompMesh *TPZMHMHDivErrorEstimator<MixedMaterial>::CreatePrimalMesh() {
    if (this->fPostProcesswithHDiv) {
        return CreateDiscontinuousPressureMesh();
    } else {
        return CreateInternallyContinuousPressureMesh();
    }
}

template <typename MixedMaterial>
TPZCompMesh *TPZMHMHDivErrorEstimator<MixedMaterial>::CreateDiscontinuousPressureMesh()
{
    TPZCompMesh *OrigPressure = this->fOriginal->MeshVector()[1];
    TPZGeoMesh *gmesh = OrigPressure->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    TPZCompMesh *pressmesh = OrigPressure->Clone();
    RemoveMaterialObjects(pressmesh->MaterialVec());
    int64_t nel = pressmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = pressmesh->Element(el);
        if(!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if(!intel) DebugStop();
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        if(!pressmesh->FindMaterial(matid))
        {
            delete cel;
        }
    }
    pressmesh->ExpandSolution();
    return pressmesh;
}

template <typename MixedMaterial>
TPZCompMesh *TPZMHMHDivErrorEstimator<MixedMaterial>::CreateInternallyContinuousPressureMesh() {
    TPZCompMesh *original_pressure = this->fOriginal->MeshVector()[1];
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
                    TPZMaterial *neighMat = this->fPostProcMesh.FindMaterial(neighbour.Element()->MaterialId());
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
    // original_pressure->CopyMaterials(*reconstruction_pressure);
    // RemoveMaterialObjects(reconstruction_pressure->MaterialVec());

    TPZL2Projection<STATE>* l2p = new TPZL2Projection<STATE>(1,2,2);
    reconstruction_pressure->InsertMaterialObject(l2p);

    // Copies BC materials
    std::set<int> bcMatIDs;
    for (auto it : this->fOriginal->MaterialVec()) {
        TPZMaterial *mat = it.second;
        TPZBndCondT<STATE> *bc = dynamic_cast<TPZBndCondT<STATE> *>(mat);
        if (bc) {
            int bcID = bc->Material()->Id();
            TPZMaterial *pressure_abs = original_pressure->FindMaterial(bcID);
            TPZMaterialT<STATE> *pressure_mat = dynamic_cast<TPZMaterialT<STATE>*>(pressure_abs);
            TPZBndCondT<STATE> *bc_mat = l2p->CreateBC(l2p, mat->Id(), bc->Type(), bc->Val1(), bc->Val2());
            if (this->fExact) {
                bc_mat->SetForcingFunctionBC(this->fExact->ExactSolution(),4);
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

        TPZCompEl *new_cel = reconstruction_pressure->CreateCompEl(gel);
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
            if (celstack.size() != 1) {
                continue;
                //DebugStop();
            }

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
template <typename MixedMaterial>
void TPZMHMHDivErrorEstimator<MixedMaterial>::RemoveMaterialObjects(std::map<int,TPZMaterial *> &matvec)
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
                if (!this->fPostProcesswithHDiv || (this->fPostProcesswithHDiv
                    && matid != fMultiPhysicsInterfaceMatId
                    && matid != fHDivWrapMatId
                    && matid != this->fPrimalSkeletonMatId)) {
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
}

// transfer embedded elements to submeshes
// the substructuring method will only transfer the H(div) and surrounding elements to the submesh
// it does not detect that the boundary pressure elements belong to the submesh. This is done in the
// following method
template <typename MixedMaterial>
void TPZMHMHDivErrorEstimator<MixedMaterial>::TransferEmbeddedElements()
{
    int64_t nel = this->fPostProcMesh.NElements();
    int64_t ncon = this->fPostProcMesh.NConnects();
    TPZVec<int> numelcon(ncon,0);
    TPZVec<TPZSubCompMesh *> connect_submesh(ncon,0);
    // numelcon : number of subcmesh sharing a connect
    // connect_submesh : pointer to the submesh if there is only one subdomain that owns the connect
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = this->fPostProcMesh.Element(el);
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
        TPZCompEl *cel = this->fPostProcMesh.Element(el);
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
                sub->TransferElement(&(this->fPostProcMesh), el);
            }
        }
    }
    this->fPostProcMesh.ComputeNodElCon();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = this->fPostProcMesh.Element(el);
        TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
        if(!sub) continue;
        sub->MakeAllInternal();
    }
}

// a method for computing the pressures/displacements between subdomains as average pressures/displacements
/// compute the average pressures/displacements of across edges of the H(div) mesh
template <typename MixedMaterial>
void TPZMHMHDivErrorEstimator<MixedMaterial>::ComputeAveragePrimal(int target_dim)
{
    // load the pressure/displacement elements of the finite element approximation
    TPZCompMesh *OrigPressure = this->fOriginal->MeshVector()[1];
    TPZGeoMesh *gmesh = OrigPressure->Reference();
    gmesh->ResetReference();
    int dim = this->fPostProcMesh.Dimension();
    
    TPZCompMesh *postprimalmesh =  this->fPostProcMesh.MeshVector()[1];
    TPZCompMesh *loadmesh = OrigPressure;
    if(target_dim < dim-1) {
        loadmesh = postprimalmesh ;
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
    nel = postprimalmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = postprimalmesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        if(gel->Dimension() == target_dim)
        {
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            if(!intel) DebugStop();
            
            int matId = gel->MaterialId();
            if(matId != this->fPrimalSkeletonMatId) continue;
            int64_t index = intel->Index();
            this->ComputeAverage(postprimalmesh,index);
        }
    }
    
}
/// compute the average pressure over corners
/// set the cornernode values equal to the averages
template <typename MixedMaterial>
void TPZMHMHDivErrorEstimator<MixedMaterial>::ComputeNodalAverages()
{
    // load the one dimensional interface elements
    // load the pressure elements of the finite element approximation
    TPZCompMesh *pressuremesh = this->fPostProcMesh.MeshVector()[1];
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
                this->ComputeNodalAverage(celside);
            }
        }
    }
}

// TODO we dont need to pass pressure_mesh as an argument here, if I divide the method in the father class into
//  creation of gels and compels I can remove it
template <typename MixedMaterial>
void TPZMHMHDivErrorEstimator<MixedMaterial>::CreateSkeletonElements(TPZCompMesh * pressure_mesh) {

    TPZCompMesh *cmesh = this->fOriginal;
    TPZGeoMesh *gmesh = this->fOriginal->Reference();
    gmesh->ResetReference();
    cmesh->LoadReferences();

#ifdef ERRORESTIMATION_DEBUG
    {
        std::ofstream fileVTK("GeoMeshBeforePressureSkeleton.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, fileVTK);
    }
#endif

    if (this->fPrimalSkeletonMatId == 0) {
        this->fPrimalSkeletonMatId = this->FindFreeMatId(this->GMesh());
        std::cout << "Created new pressure skeleton material of index " << this->fPrimalSkeletonMatId << '\n';
    }

    const TPZManVector<int64_t> geoToMHM = fMHM->GetGeoToMHMDomain();

    const int nel = gmesh->NElements();
    int dim = gmesh->Dimension();
    for (int iel = 0; iel < nel; iel++) {

        TPZGeoEl *gel = gmesh->Element(iel);
        if (!gel) continue;
        TPZCompEl *cel = gel->Reference();

        if (!cel) continue;
        if (gel->Dimension() != dim) continue;

        // Iterates through the sides of the element
        int nsides = gel->NSides();
        for (int iside = 0; iside < nsides; iside++) {
            TPZGeoElSide gelside(gel, iside);

            // Filters boundary sides
            if (gelside.Dimension() != dim - 1) continue;

            for (TPZGeoElSide neighbour = gelside.Neighbour(); neighbour != gelside; neighbour++) {
                TPZGeoEl *neigh_gel = neighbour.Element();
                if (neigh_gel->Dimension() != dim) continue;

                int64_t gel_index = gel->Index();
                int64_t neigh_gel_index = neighbour.Element()->Index();
                if (geoToMHM[gel_index] != geoToMHM[neigh_gel_index]) {
                    if (!gelside.HasNeighbour(this->fPrimalSkeletonMatId)) {
                        TPZGeoElBC gbc(gelside, this->fPrimalSkeletonMatId);
                        break;
                    }
                }
            }
        }
    }

#ifdef ERRORESTIMATION_DEBUG
    {
        std::ofstream fileVTK("GeoMeshAfterPressureSkeleton.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, fileVTK);
    }
#endif
}

template <typename MixedMaterial>
void TPZMHMHDivErrorEstimator<MixedMaterial>::CreateSkeletonApproximationSpace(TPZCompMesh *pressure_mesh) {

    TPZGeoMesh *gmesh = pressure_mesh->Reference();
    int dim = gmesh->Dimension();

    // Create skeleton elements in pressure mesh
    TPZNullMaterial<> *skeletonMat = new TPZNullMaterial<>(this->fPrimalSkeletonMatId);
    skeletonMat->SetDimension(dim - 1);
    if (fMHM->ProblemType() == TPZMHMeshControl::MProblemType::EElasticity2D) skeletonMat->SetNStateVariables(dim);
    // if(fConfig.problemtype == ProblemConfig::TProbType::EElasticity) 
    
    pressure_mesh->InsertMaterialObject(skeletonMat);

    std::set<int> matIdSkeleton = { this->fPrimalSkeletonMatId };
    gmesh->ResetReference();

    pressure_mesh->ApproxSpace().CreateDisconnectedElements(true);
    pressure_mesh->AutoBuild(matIdSkeleton);
    pressure_mesh->ExpandSolution();
}

template <typename MixedMaterial>
void TPZMHMHDivErrorEstimator<MixedMaterial>::CopySolutionFromSkeleton() {

    TPZCompMesh *pressuremesh = this->PrimalMesh();

    pressuremesh->Reference()->ResetReference();

    if(pressuremesh->Reference()!= this->fPostProcMesh.Reference()) DebugStop();
    
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


    TPZBlock &block =  pressuremesh->Block();
    TPZFMatrix<STATE> &sol = pressuremesh->Solution();
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
                    sol.at(block.at(c_neigh_seqnum, 0, ibl, 0)) = sol.at(block.at(c_gelSide_seqnum, 0, ibl, 0));
                }
            }
        }
    }
}

template <typename MixedMaterial>
void TPZMHMHDivErrorEstimator<MixedMaterial>::VerifySolutionConsistency(TPZCompMesh* cmesh) {
//    {
//        std::ofstream outvtk("MeshToVerifyConsistency.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(), outvtk);
//        std::ofstream outtxt("MeshToVerifyConsistency.txt");
//        cmesh->Print(outtxt);
//    }

    TPZGeoMesh* gmesh = this->fOriginal->Reference();
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
                    TPZManVector<STATE> sol0(dim);
                    cel->Solution(pt0_vol, 1, sol0);

                    TPZTransform<REAL> neighSideToVolume(dim, dim);
                    neighSideToVolume = neighbour.Element()->SideToSideTransform(cneighbour.Side(), neighbour.Element()->NSides() - 1);

                    TPZManVector<REAL> pt1_vol(dim, 0);
                    neighSideToVolume.Apply(pt1, pt1_vol);
                    TPZManVector<STATE> sol1(dim);
                    cneighbour.Element()->Solution(pt1_vol, 1, sol1);

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
                    
{
    std::cout << "\nSide Element =  " << gelside.Element()->Index() << "\n";
    std::cout << "Neighbour Element =  " << neighbour.Element()->Index() << "\n";
    std::cout << "Side solutionX =  " << sol0[0] << " Side solutionY =  " << sol0[1] <<"\n";
    std::cout << "Neigh solutionX = " << sol1[0] << " Neigh solutionY = " << sol1[1] <<"\n";
    std::cout << "DiffX = " << sol1[0] - sol0[0] << " DiffY = " << sol1[1] - sol0[1] <<"\n";
    std::cout << "Side coord:  [" << x0[0] << ", " << x0[1] << ", " << x0[2] << "]\n";
    std::cout << "Neigh coord: [" << x1[0] << ", " << x1[1] << ", " << x1[2] << "]\n";
}
                    

                    // Checks pressure value on these nodes
                    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cneighbour.Element());
                    if (!intel) DebugStop();
                }
            }
            delete intRule;
        }
    }
}

template <typename MixedMaterial>
void TPZMHMHDivErrorEstimator<MixedMaterial>::ComputeConnectsNextToSkeleton(std::set<int64_t>& connectList) {

    TPZCompMesh *pressure_reconstruction = this->fPostProcMesh.MeshVector()[1];
    TPZGeoMesh * gmesh = pressure_reconstruction->Reference();
    int dim = gmesh->Dimension();
    gmesh->ResetReference();
    pressure_reconstruction->LoadReferences();

    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl* gel = gmesh->Element(el);
        if (!gel) continue;
        if (gel->Dimension() == dim) continue;
        if (gel->MaterialId() != this->fPrimalSkeletonMatId) continue;
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

template <typename MixedMaterial>
void TPZMHMHDivErrorEstimator<MixedMaterial>::CreateFluxSkeletonElements(TPZCompMesh *flux_mesh) {

    std::map<TPZGeoEl *, TPZCompEl *> geltocel;
    {
        int64_t nel = flux_mesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = flux_mesh->Element(el);
            if(cel)
            {
                TPZGeoEl *gel = cel->Reference();
                geltocel[gel] = cel;
            }
        }
    }
    TPZGeoMesh *gmesh = this->GMesh();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();

    if (this->fHDivWrapMatId == 0) {
        this->fHDivWrapMatId = this->FindFreeMatId(this->GMesh());
    }

    auto wrapMat = new TPZNullMaterial(this->fHDivWrapMatId);
    wrapMat->SetDimension(dim - 1);
    wrapMat->SetNStateVariables(1);
    flux_mesh->InsertMaterialObject(wrapMat);
    std::cout << "Created new HDivWrap material of index " << this->fHDivWrapMatId << '\n';

    // Iterates over geometric elements
    int64_t nel = gmesh->NElements();
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->Element(iel);
        if (!gel) continue;

        // Continue if element is not of pressure skeleton mat id
        if (gel->MaterialId() != this->fPrimalSkeletonMatId) continue;

        // Get side of largest dimension of the skeleton
        TPZGeoElSide skelSide = TPZGeoElSide(gel);

        int count = 0;
        for (TPZGeoElSide neighbour = skelSide.Neighbour(); neighbour != skelSide; neighbour++) {
            TPZGeoEl *neighGel = neighbour.Element();
            if (!neighGel) continue;
            if (neighGel->Dimension() != dim) continue;
            TPZCompEl *neighCel = geltocel[neighGel];
            if (!neighCel) DebugStop();
            TPZInterpolatedElement *neighIntel = dynamic_cast<TPZInterpolatedElement*>(neighCel);
            if (!neighIntel) DebugStop();

            // TODO sometimes the order here is 1, I think it should always be 3 for the case I'm running. Need to
            //  think if there are better ways to obtain the side connect order
            TPZConnect &c = neighIntel->SideConnect(0, neighbour.Side());
            int order = c.Order();
            int64_t cindex = neighIntel->SideConnectIndex(0, neighbour.Side());

            neighIntel->SetPreferredOrder(order);
            neighGel->SetReference(neighCel);

            TPZGeoElBC gbc(neighbour, this->fHDivWrapMatId);
            TPZCompEl * wrap = flux_mesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *flux_mesh);
            TPZInterpolatedElement *wrapintel = dynamic_cast<TPZInterpolatedElement *>(wrap);
            neighIntel->SetSideOrient(neighbour.Side(), 1);
            wrapintel->SetSideOrient(gbc.CreatedElement()->NSides()-1,1);
            int64_t newc_index = wrap->ConnectIndex(0);
            if(newc_index != cindex) DebugStop();
            neighGel->ResetReference();
            gbc.CreatedElement()->ResetReference();
            count++;
        }
        if (count != 2) DebugStop();
    }

}

template <typename MixedMaterial>
void TPZMHMHDivErrorEstimator<MixedMaterial>::CreateMultiphysicsInterfaces() {

    this->fPostProcMesh.LoadReferences();
    TPZGeoMesh *gmesh = this->GMesh();
    int dim = gmesh->Dimension();

    if (this->fMultiPhysicsInterfaceMatId == 0) {
        this->fMultiPhysicsInterfaceMatId = this->FindFreeMatId(this->GMesh());
    }

    TPZLagrangeMultiplierCS<> *interfaceMat = new TPZLagrangeMultiplierCS<>(fMultiPhysicsInterfaceMatId, dim - 1, 1);
    this->fPostProcMesh.InsertMaterialObject(interfaceMat);
    std::cout << "Created interface material of index " << fMultiPhysicsInterfaceMatId << '\n';

    // Iterates over geometric elements
    int64_t nel = gmesh->NElements();
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->Element(iel);

        if (gel->MaterialId() != this->fPrimalSkeletonMatId) continue;

        // Get side of largest dimension of the skeleton
        TPZGeoElSide skelSide = TPZGeoElSide(gel);
        TPZCompEl *skelCel = gel->Reference();
        if (!skelCel) DebugStop();

        int count = 0;
        for (TPZGeoElSide neigh = skelSide.Neighbour(); neigh != skelSide; neigh++) {
            TPZCompElSide skelCelSide = skelSide.Reference();
            if (!skelCelSide) DebugStop();

            TPZGeoEl *neighGel = neigh.Element();
            if (neighGel->MaterialId() != this->fHDivWrapMatId) continue;

            TPZGeoElBC gbc(skelSide, this->fMultiPhysicsInterfaceMatId);
            count++;

            TPZCompElSide neighSide = neigh.Reference();
            if (!neighSide) DebugStop();

            auto *interface = new TPZMultiphysicsInterfaceElement(this->fPostProcMesh, gbc.CreatedElement(), skelCelSide, neighSide);
        }
        if (count != 2) DebugStop();
    }
    {
        std::ofstream file("GmeshAfterInterfaces.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(this->GMesh(), file);
    }
}

///// compute the pressure weights and material weights
//// fills in the data structure fPressureweights and fMatid_weights
//void TPZMHMHDivErrorEstimator<TPZMixedDarcyFlow>::ComputePrimalWeights() {
//    std::cout << "Computing pressure weights\n";
//    TPZCompMesh *primalMesh = fPostProcMesh.MeshVector()[1];
//    const int dim = primalMesh->Dimension();
//    const int64_t nel = primalMesh->NElements();
//    fPrimalWeights.Resize(nel, 0);
//    fMatid_weights.clear();
//    for (int64_t el = 0; el < nel; el++) {
//        TPZCompEl *cel = primalMesh->Element(el);
//        if (!cel) continue;
//        TPZGeoEl *gel = cel->Reference();
//        int matid = gel->MaterialId();
//        TPZMaterial *mat = this->fOriginal->FindMaterial(matid);
//        if (matid == fPrimalSkeletonMatId || matid == fHybridizer.fLagrangeInterface) {
//            fPrimalWeights[el] = 0.;
//            fMatid_weights[matid] = 0.;
//            continue;
//        }
//        if (!mat) DebugStop();
//
//        TPZBndCondT<STATE> *bcmat = dynamic_cast<TPZBndCondT<STATE> *>(mat);
//        if (bcmat) {
//            switch(bcmat->Type()) {
//                case 0: // Dirichlet BC
//                    fPrimalWeights[el] = 1.e12;
//                    fMatid_weights[matid] = 1.e12;
//                    break;
//                case 1: // Neumann BC
//                    fPrimalWeights[el] = 0;
//                    fMatid_weights[matid] = 0;
//                    break;
//                case 4: // Robin BC, weight = Km
//                    fPrimalWeights[el] = bcmat->Val1()(0, 0);
//                    fMatid_weights[matid] = bcmat->Val1()(0, 0);
//                    break;
//                default:
//                    DebugStop();
//            }
//        } else {
//            MixedMaterial *mixpoisson = dynamic_cast<MixedMaterial *>(mat);
//            if (!mixpoisson) DebugStop();
//
//            REAL perm;
//            TPZVec<REAL> xi(gel->Dimension(), 0.);
//            gel->CenterPoint(gel->NSides() - 1, xi);
//            TPZVec<REAL> x(3, 0.);
//            gel->X(xi, x);
//            perm = mixpoisson->GetPermeability(x);
//
//            if (IsZero(perm)) DebugStop();
//            this->fPrimalWeights[el] = perm;
//            fMatid_weights[matid] = perm;
//        }
//    }
//    std::cout << "Finished computing pressure weights\n";
//}

template class TPZMHMHDivErrorEstimator<TPZMixedDarcyFlow>;
template class TPZMHMHDivErrorEstimator<TPZMixedElasticityND>;
