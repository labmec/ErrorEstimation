
#include "TPZMHMHDivErrorEstimator.h"
#include "pzbndcond.h"
#include "pzintel.h"
#include "pzsubcmesh.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMHMHDivErrorEstimatorMaterial.h"
#include "TPZNullMaterial.h"

// a method for generating the hybridized multiphysics post processing mesh
void TPZMHMHDivErrorEstimator::CreatePostProcessingMesh()
{
    
    // initialize the post processing mesh
    fPostProcMesh.SetReference(fOriginal->Reference());

    fOriginal->CopyMaterials(fPostProcMesh);
        // switch the material from mixed to TPZMHMHDivErrorEstimationMaterial...
    SwitchMaterialObjects();
    

    TPZManVector<TPZCompMesh *> meshvec(4);
    meshvec[0] = CreateFluxMesh();
    meshvec[1] = CreatePressureMesh();
    meshvec[2] = fOriginal->MeshVector()[0];
    meshvec[3] = fOriginal->MeshVector()[1];
    
    TPZManVector<int,4> active(4,0);
    active[0] = 1;
    active[1] = 1;
 
    RemoveMaterialObjects(fPostProcMesh.MaterialVec());
    fPostProcMesh.BuildMultiphysicsSpace(active, meshvec);
    bool groupelements = false;
#ifdef PZDEBUG2
    {
        std::ofstream out1("fluxpostNH.txt");
        meshvec[0]->Print(out1);
        std::ofstream out2("pressurepostNH.txt");
        meshvec[1]->Print(out2);
        std::ofstream out3("mphyspost_beforeNH.txt");
        fPostProcMesh.Print(out3);
    }
#endif
    
    fHybridizer.HybridizeGivenMesh(fPostProcMesh,groupelements);

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

    SubStructurePostProcessingMesh();
    
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
    fOriginal->Reference()->ResetReference();
    TPZGeoMesh *gmesh = fOriginal->Reference();
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
#ifdef PZDEBUG
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
        if(gel->Dimension() != dim) continue;
        TPZSubCompMesh *submesh = ReferredMesh[gel->Index()];
        if(!submesh) DebugStop();
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
        fPostProcMesh.LoadReferences();
        std::ofstream out("postproc_part_cmesh_substruct.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(&fPostProcMesh, out,true);
    }
#endif
    // transfer the elements this procedure didn't recognize
    // this will transfer the wrappers, interface elements and lagrange multipliers
    // transfer the elements whose connects belong to a unique submesh
    TransferEmbeddedElements();
    fPostProcMesh.ComputeNodElCon();
    fPostProcMesh.CleanUpUnconnectedNodes();
    // set an analysis type for the submeshes
    {
        int64_t nel = fPostProcMesh.NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = fPostProcMesh.Element(el);
            TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
            if(sub)
            {
                fHybridizer.GroupandCondenseElements(sub);
                sub->CleanUpUnconnectedNodes();
                int numthreads = 0;
                int preconditioned = 0;
                TPZAutoPointer<TPZGuiInterface> guiInterface;
                sub->SetAnalysisSkyline(numthreads, preconditioned, guiInterface);
            }
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



// a method for generating the HDiv mesh
TPZCompMesh *TPZMHMHDivErrorEstimator::CreateFluxMesh()
{
    if(fPostProcesswithHDiv == false)
    {
        return 0;
    }
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

// method fro creating a discontinuous pressure mesh
TPZCompMesh *TPZMHMHDivErrorEstimator::CreatePressureMesh()
{
    if(fPostProcesswithHDiv)
    {
        return CreateDiscontinuousPressureMesh();
    }
    else
    {
        return CreateContinousPressureMesh();
    }
}

TPZCompMesh *TPZMHMHDivErrorEstimator::CreateContinousPressureMesh()
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
    }
    pressmesh->ExpandSolution();
    return pressmesh;

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
    if(target_dim < dim-1) OrigPressure = fPostProcMesh.MeshVector()[1];
    int64_t nel = OrigPressure->NElements();
    // load all elements of dimension target_dim+1
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = OrigPressure->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() == target_dim+1) {
            cel->LoadElementReference();
        }
    }
    // compute the averages one element at a time
    nel = fPostProcMesh.NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fPostProcMesh.Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        if(gel->Dimension() == target_dim)
        {
            TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if(!mphys) DebugStop();
            TPZCompEl *pressel = mphys->Element(1);
            int64_t index = pressel->Index();
            ComputeAverage(pressel->Mesh(),index);
        }
    }
    
}
/// compute the average pressure over corners
/// set the cornernode values equal to the averages
void TPZMHMHDivErrorEstimator::ComputeNodalAverages()
{
    // load the one dimensional interface elements
    // load the pressure elements of the finite element approximation
    TPZCompMesh *OrigPressure = fPostProcMesh.MeshVector()[1];
    TPZGeoMesh *gmesh = OrigPressure->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    int64_t nel = fPostProcMesh.NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fPostProcMesh.Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        if (gel->Dimension() == dim-1) {
            TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if(!mphys) DebugStop();
            TPZCompEl *pressel = mphys->Element(1);
            pressel->LoadElementReference();
        }
    }

    nel = OrigPressure->NElements();
    // compute the averages
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = OrigPressure->Element(el);
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


void TPZMHMHDivErrorEstimator::SwitchMaterialObjects()
{
    // switch the material
    for(auto matid : fPostProcMesh.MaterialVec())
    {
        TPZMixedPoisson *mixpoisson = dynamic_cast<TPZMixedPoisson *> (matid.second);
        if(mixpoisson)
        {
            
            TPZMHMHDivErrorEstimateMaterial *newmat = new TPZMHMHDivErrorEstimateMaterial(*mixpoisson);
            
            if(fExact)
            {
                newmat->SetForcingFunctionExact(fExact->Exact());
                newmat->SetForcingFunction(fExact->ForcingFunction());
                
            }
            
            for (auto bcmat : fPostProcMesh.MaterialVec()) {
                TPZBndCond *bc = dynamic_cast<TPZBndCond *>(bcmat.second);
                if (bc) {
                    bc->SetMaterial(newmat);
                }
            }
            fPostProcMesh.MaterialVec()[newmat->Id()] = newmat;
            delete mixpoisson;
        }
    }
    
}


void TPZMHMHDivErrorEstimator::CreatePressureSkeleton() {

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

#ifdef PZDEBUG2
    {
        std::ofstream fileVTK("GeoMeshAfterPressureSkeleton.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, fileVTK);
        std::ofstream fileTXT("GeoMeshAfterPressureSkeleton.txt");
        gmesh->Print(fileTXT);
    }
#endif
}

TPZCompMesh *TPZMHMHDivErrorEstimator::CreateSkeletonApproximationSpace() {
    TPZCompMesh *cmesh = new TPZCompMesh(fOriginal->Reference());

    TPZNullMaterial *skeletonMat = new TPZNullMaterial(fSkeletonMatId);
    cmesh->InsertMaterialObject(skeletonMat);

    cmesh->SetDefaultOrder(fOriginal->GetDefaultOrder());

    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    cmesh->AutoBuild();

    return cmesh;
}