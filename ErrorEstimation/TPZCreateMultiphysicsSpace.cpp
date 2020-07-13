//
//  TPZCreateMultiphysicsSpace.cpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 13/07/19.
//

#include "TPZCreateMultiphysicsSpace.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoelbc.h"
#include "pzintel.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZNullMaterial.h"
#include "TPZLagrangeMultiplier.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZHybridizeHDiv.h"

TPZCreateMultiphysicsSpace::TPZCreateMultiphysicsSpace(TPZGeoMesh *gmesh, MSpaceType spacetype) :
            fSpaceType(spacetype), fGeoMesh(gmesh) {
    fDimension = gmesh->Dimension();
}

/// copy constructor
TPZCreateMultiphysicsSpace::TConfigH1Hybrid::TConfigH1Hybrid(const TConfigH1Hybrid &copy)
{
    
}

/// copy operator
TPZCreateMultiphysicsSpace::TConfigH1Hybrid &TPZCreateMultiphysicsSpace::TConfigH1Hybrid::operator=(const TConfigH1Hybrid &copy)
{
    return *this;
}


/// copy constructor
TPZCreateMultiphysicsSpace::TPZCreateMultiphysicsSpace(const TPZCreateMultiphysicsSpace &copy)
{
    
}

/// = operator
TPZCreateMultiphysicsSpace & TPZCreateMultiphysicsSpace::operator=(const TPZCreateMultiphysicsSpace &copy)
{
    return *this;
}

/// Indicate to create Hybridized H1 meshes
void TPZCreateMultiphysicsSpace::SetH1Hybridized(const TConfigH1Hybrid &config)
{
    fSpaceType = EH1Hybrid;
    fH1Hybrid = config;
}

/// create meshes and elements for all geometric elements
void TPZCreateMultiphysicsSpace::CreateAtomicMeshes(TPZVec<TPZCompMesh *> &meshvec, int pressureOrder, int lagrangeorder)
{
    AddGeometricWrapElements();
    SetPOrder(pressureOrder);
    SetLagrangeOrder(lagrangeorder);
    TPZCompMesh *pressure = CreatePressureMesh();
//    CreateLagrangeGeometricElements(pressure);
    TPZCompMesh *fluxmesh = CreateBoundaryFluxMesh();
    TPZCompMesh *gspace = new TPZCompMesh(fGeoMesh);
    {
        InsertNullSpaceMaterialIds(gspace);
        gspace->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
        gspace->SetDefaultOrder(0);//sao espacos de pressao media 
        gspace->AutoBuild();
        int64_t nconnects = gspace->NConnects();
        for (int ic = 0; ic<nconnects; ic++) {
            gspace->ConnectVec()[ic].SetLagrangeMultiplier(2);
        }
    }
    TPZCompMesh *average = new TPZCompMesh(fGeoMesh);
    {
        InsertNullSpaceMaterialIds(average);
        average->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
        average->SetDefaultOrder(0);
        average->AutoBuild();
        int64_t nconnects = average->NConnects();
        for (int ic = 0; ic<nconnects; ic++) {
            average->ConnectVec()[ic].SetLagrangeMultiplier(5);
        }
    }

    meshvec.Resize(4,0);
    meshvec[0] = pressure;
    meshvec[1] = fluxmesh;
    meshvec[2] = gspace;
    meshvec[3] = average;
}


/// if there a neighbouring element with matid == lagrangematid -> return true
bool TPZCreateMultiphysicsSpace::ShouldCreateFluxElement(TPZGeoElSide &gelside, int lagrangematid)
{
    TPZGeoElSide neighbour(gelside.Neighbour());
    while(neighbour != gelside)
    {
        int matid = neighbour.Element()->MaterialId();
        if(matid == lagrangematid) return false;
        if(fBCMaterialIds.find(matid) != fBCMaterialIds.end()) return false;
        
        neighbour = neighbour.Neighbour();
    }
    return true;
}
/// create the geometric elements for the lagrange multipliers
// these elements will go with the largest H1 element
void TPZCreateMultiphysicsSpace::CreateLagrangeGeometricElements(TPZCompMesh *pressure)
{
    // this method shouldn t be called anymore. All geometric elements are created in
    // AddGeometricWrapElements
    DebugStop();
    TPZGeoMesh *gmesh = pressure->Reference();
    gmesh->ResetReference();
    pressure->LoadReferences();
    int64_t nel = pressure->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = pressure->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        if(gel->Dimension() != fDimension)
        {
            continue;
        }
        int matid = gel->MaterialId();
        if(fMaterialIds.find(matid) == fMaterialIds.end()) DebugStop();
        int nsides = gel->NSides();
        int ncorner = gel->NCornerNodes();
        for (int side = ncorner; side < nsides-1; side++) {
            if(gel->SideDimension(side) != fDimension-1) continue;
            TPZGeoElSide gelside(gel,side);
            TPZStack<TPZCompElSide> celstack_eq, celstack_sm;
            TPZCompElSide large;
            gelside.EqualLevelCompElementList(celstack_eq, 1, 0);
            gelside.HigherLevelCompElementList2(celstack_sm, 1, 0);
            large = gelside.LowerLevelCompElementList2(1);
            int ndif = 0;
            if(celstack_eq.size()) ndif++;
            if(celstack_sm.size()) ndif++;
            if(large) ndif++;
//            if(ndif > 1) DebugStop();
            // do not create a geometric flux element if there is a larger element connected
            if(large) continue;
            if(celstack_eq.size())
            {
                // figure out if there is a lagrange material element neighbour
                if(ShouldCreateFluxElement(gelside, fH1Hybrid.fFluxMatId) == true)
                {
                    TPZGeoElBC(gelside,fH1Hybrid.fFluxMatId);
                }
            }
            if(celstack_sm.size())
            {
                // figure out if there is a lagrange material element neighbour
                if(ShouldCreateFluxElement(gelside, fH1Hybrid.fFluxMatId) == true)
                {
                    TPZGeoElBC(gelside,fH1Hybrid.fFluxMatId);
                }
            }
            // if large exists we do not create a lagrange element
        }
    }
}

/// if there a neighbouring element with matid == lagrangematid -> return true
static TPZGeoElSide HasBCNeighbour(const TPZGeoElSide &gelside, const std::set<int> &matbc)
{
    TPZGeoElSide neighbour(gelside.Neighbour());
    while(neighbour != gelside)
    {
        int matid = neighbour.Element()->MaterialId();
        if(matbc.find(matid) != matbc.end()) return neighbour;
        neighbour = neighbour.Neighbour();
    }
    return TPZGeoElSide();
}

/// create the pressure boundary elements if the boundary is not hybridized
void TPZCreateMultiphysicsSpace::CreatePressureBoundaryElements(TPZCompMesh *pressure)
{
    if (fSpaceType != EH1Hybrid && fSpaceType != EH1HybridSquared) {
        DebugStop();
    }
#ifdef PZDEBUG
    std::map<int,int> numcreated;
#endif
    // create elements connected to the pressure elements, the neighbour of each pressure
    // element should be either matwrap or boundary condition
    pressure->Reference()->ResetReference();
    pressure->ApproxSpace().CreateDisconnectedElements(false);
    int64_t nel = pressure->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure->Element(el);
        if(!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        int porder = intel->GetPreferredOrder();
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        int matid = gel->MaterialId();
        // loop only over volumetric elements
        if(fMaterialIds.find(matid) == fMaterialIds.end()) continue;
        if (gel->Dimension() != fDimension) {
            DebugStop();
        }
        int nsides = gel->NSides();
        int ncorner = gel->NCornerNodes();
        for (int side = ncorner; side < nsides; side++) {
            if(gel->SideDimension(side) != fDimension-1) continue;
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            TPZGeoElSide bcneighbour = HasBCNeighbour(gelside, fBCMaterialIds);
            // the boundary condition element should be the first neighbour
            if(fH1Hybrid.fHybridizeBCLevel == 0 && bcneighbour && bcneighbour != neighbour) DebugStop();
            if(!bcneighbour && neighbour.Element()->MaterialId() != fH1Hybrid.fMatWrapId) DebugStop();
            pressure->SetDefaultOrder(porder);
            int dir = gel->NormalOrientation(side);
            int wrapmat = fH1Hybrid.fMatWrapId;
            // if neighbour exists then create the wrap conditionally
            {
                // load the element reference so that the created element will share the connects
                cel->LoadElementReference();
                int64_t index;
                TPZCompEl *bc_cel = pressure->ApproxSpace().CreateCompEl(neighbour.Element(), *pressure, index);
                // reset the references so that future elements will not share connects
#ifdef PZDEBUG
                numcreated[neighbour.Element()->MaterialId()]++;
#endif
                gel->ResetReference();
                bc_cel->Reference()->ResetReference();
            }
        }
    }
    if(fSpaceType == EH1HybridSquared && fH1Hybrid.fHybridizeBCLevel == 2)
    {
        // be sure to verify on first execution
        // create boundary elements in the pressure space
        int64_t nel = fGeoMesh->NElements();
        pressure->SetDefaultOrder(fDefaultLagrangeOrder);
        for(int64_t el = 0; el<nel; el++)
        {
            TPZGeoEl *gel = fGeoMesh->Element(el);
            int matid = gel->MaterialId();
            if(fBCMaterialIds.find(matid) == fBCMaterialIds.end()) continue;
            int64_t index;
            TPZCompEl *cel = pressure->ApproxSpace().CreateCompEl(gel, *pressure, index);
#ifdef PZDEBUG
            numcreated[matid]++;
#endif
            int nc = cel->NConnects();
            for(int ic=0; ic<nc; ic++)
            {
                cel->Connect(ic).SetLagrangeMultiplier(6);
            }
            gel->ResetReference();
        }
        pressure->ExpandSolution();
    }
#ifdef PZDEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    for(auto it : numcreated) std::cout << "Material ID " << it.first << " number of elements created " << it.second << std::endl;
        
#endif
}

/// insert the pressure material ids
void TPZCreateMultiphysicsSpace::InsertPressureMaterialIds(TPZCompMesh *pressure)
{
    for (auto matid:fMaterialIds) {
        TPZNullMaterial *nullmat = new TPZNullMaterial(matid);
        nullmat->SetDimension(fDimension);
        nullmat->SetNStateVariables(1);
        pressure->InsertMaterialObject(nullmat);
    }
    if((fH1Hybrid.fHybridizeBCLevel == 2) || (fH1Hybrid.fHybridizeBCLevel == 0))
    {
        for (auto matid:fBCMaterialIds) {
            TPZNullMaterial *nullmat = new TPZNullMaterial(matid);
            nullmat->SetDimension(fDimension);
            nullmat->SetNStateVariables(1);
            pressure->InsertMaterialObject(nullmat);
        }
    }
    {
        TPZNullMaterial *nullmat = new TPZNullMaterial(fH1Hybrid.fMatWrapId);
        nullmat->SetDimension(fDimension);
        nullmat->SetNStateVariables(1);
        pressure->InsertMaterialObject(nullmat);
    }
    if(fSpaceType == EH1HybridSquared)
    {
        TPZNullMaterial *nullmat = new TPZNullMaterial(fH1Hybrid.fInterfacePressure);
        nullmat->SetDimension(fDimension);
        nullmat->SetNStateVariables(1);
        pressure->InsertMaterialObject(nullmat);

    }
}

/// insert flux material ids
void TPZCreateMultiphysicsSpace::InsertFluxMaterialIds(TPZCompMesh *fluxmesh)
{
    if (fSpaceType == EH1Hybrid || fSpaceType == EH1HybridSquared) {
        int matid = fH1Hybrid.fFluxMatId;
        TPZNullMaterial *nullmat = new TPZNullMaterial(matid);
        nullmat->SetDimension(fDimension-1);
        nullmat->SetNStateVariables(1);
        fluxmesh->InsertMaterialObject(nullmat);
    }
    else {
        DebugStop();
    }
    if(fH1Hybrid.fHybridizeBCLevel == 1)
    {
        for (auto matid:fBCMaterialIds) {
            TPZNullMaterial *nullmat = new TPZNullMaterial(matid);
            nullmat->SetDimension(fDimension);
            nullmat->SetNStateVariables(1);
            fluxmesh->InsertMaterialObject(nullmat);
        }
    }

}

/// insert materialids for the null space
void TPZCreateMultiphysicsSpace::InsertNullSpaceMaterialIds(TPZCompMesh *nullspace)
{
    for (auto matid:fMaterialIds) {
        TPZNullMaterial *nullmat = new TPZNullMaterial(matid);
        nullmat->SetDimension(fDimension);
        nullmat->SetNStateVariables(1);
        nullspace->InsertMaterialObject(nullmat);
    }
}

/// Create the pressure mesh
TPZCompMesh *TPZCreateMultiphysicsSpace::CreatePressureMesh()
{
    
    // create the pressure mesh
    TPZCompMesh *pressure = new TPZCompMesh(fGeoMesh);
    InsertPressureMaterialIds(pressure);
    pressure->ApproxSpace().SetAllCreateFunctionsContinuous();
    pressure->ApproxSpace().CreateDisconnectedElements(true);
    
    pressure->SetDefaultOrder(fDefaultPOrder);
    pressure->AutoBuild(fMaterialIds);
    int64_t nelem = pressure->NElements();
    for (int64_t el = 0; el<nelem; el++) {
        TPZCompEl *cel = pressure->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int nconnects = cel->NConnects();
        cel->Connect(0).SetLagrangeMultiplier(3);
        for (int ic=1; ic<nconnects; ic++) {
            cel->Connect(ic).SetLagrangeMultiplier(1);
        }
    }
#ifdef PZDEBUG
    std::cout << "Number of volumetric pressure elements created " << nelem << std::endl;
#endif
    // se nao condensar tem que mudar o nivel de lagrange multiplier de um connect
    if(fSpaceType == EH1HybridSquared)
    {
        std::set<int> matids;
        if(fH1Hybrid.fHybridizeBCLevel == 2) matids = fBCMaterialIds;
        matids.insert(fH1Hybrid.fInterfacePressure);
        pressure->SetDefaultOrder(fDefaultLagrangeOrder);
        pressure->AutoBuild(matids);
    }
    int64_t nelem_big = pressure->NElements();
    for (int64_t el = nelem; el<nelem_big; el++) {
        TPZCompEl *cel = pressure->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int nconnects = cel->NConnects();
        for (int ic=0; ic<nconnects; ic++) {
            cel->Connect(ic).SetLagrangeMultiplier(6);
        }
    }
#ifdef PZDEBUG
    std::cout << "Number of lower dimensional pressure elements created " << nelem_big-nelem << std::endl;
#endif
    CreatePressureBoundaryElements(pressure);
    return pressure;
}

/// Create the flux mesh for 1D-elements
TPZCompMesh *TPZCreateMultiphysicsSpace::CreateBoundaryFluxMesh()
{
    TPZCompMesh *fluxmesh = new TPZCompMesh(fGeoMesh);
    InsertFluxMaterialIds(fluxmesh);
    fluxmesh->ApproxSpace().SetAllCreateFunctionsHDiv(fDimension);
    fluxmesh->ApproxSpace().CreateDisconnectedElements(true);
    fluxmesh->SetDefaultOrder(fDefaultLagrangeOrder);
    fluxmesh->AutoBuild();
    int64_t nconnects = fluxmesh->NConnects();
    for (int ic=0; ic<nconnects; ic++) {
        fluxmesh->ConnectVec()[ic].SetLagrangeMultiplier(4);
    }
    return fluxmesh;
}

/// add interface elements to the multiphysics space
void TPZCreateMultiphysicsSpace::AddInterfaceElements(TPZMultiphysicsCompMesh *mphys)
{
#ifdef PZDEBUG
    std::map<int,int> numcreated;
#endif
    TPZGeoMesh *gmesh = mphys->Reference();
    gmesh->ResetReference();
    mphys->LoadReferences();
    int64_t nel = mphys->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = mphys->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
//        if(matid != fH1Hybrid.fMatWrapId.first && matid != fH1Hybrid.fMatWrapId.second) continue;
        if(matid == fH1Hybrid.fMatWrapId)
        {
            TPZCompEl *fluxel = FindFluxElement(cel);
            TPZGeoEl *fluxgel = fluxel->Reference();
            TPZGeoElSide gelside(gel);
            TPZGeoElSide neighbour = gelside.Neighbour();
            int neighmat = neighbour.Element()->MaterialId();
            if(neighmat != fH1Hybrid.fLagrangeMatid.first && neighmat != fH1Hybrid.fLagrangeMatid.second)
            {
                DebugStop();
            }
            // determine if the interface should be positive or negative...
            int interfacematid = neighmat;
            int64_t index;
            TPZCompElSide celwrap(cel,gel->NSides()-1);
            TPZGeoElSide fluxgelside(fluxgel);
            TPZCompElSide fluxside = fluxgelside.Reference();
            std::cout << "Creating interface from wrap element " << gel->Index() << " using neighbour " << neighbour.Element()->Index() <<
             " and flux element " << fluxgel->Index() << std::endl;
            if(neighbour.Element()->Reference()) DebugStop();
            new TPZMultiphysicsInterfaceElement(*mphys,neighbour.Element(),index,celwrap,fluxside);
#ifdef PZDEBUG
            numcreated[neighmat]++;
#endif
        }
        if(fSpaceType == EH1HybridSquared && matid == fH1Hybrid.fFluxMatId)
        {
            TPZGeoElSide gelsideflux(gel);
            TPZGeoElSide neighbour = gelsideflux.Neighbour();
            // we only handle the flux element neighbour to second lagrange multiplier
            if(neighbour.Element()->MaterialId() != fH1Hybrid.fSecondLagrangeMatid) continue;
            TPZGeoElSide firstlagrange = neighbour;
            TPZGeoElSide pressureinterface = firstlagrange.Neighbour();
            if(pressureinterface.Element()->MaterialId() != fH1Hybrid.fInterfacePressure)
            {
                int pressmatid = pressureinterface.Element()->MaterialId();
                if(fBCMaterialIds.find(pressmatid) == fBCMaterialIds.end()) DebugStop();
            }
            else {
                TPZGeoElSide secondlagrange = pressureinterface.Neighbour();
                if(secondlagrange.Element()->MaterialId() != fH1Hybrid.fSecondLagrangeMatid) DebugStop();
                // now we have to find the second flux element
                TPZGeoElSide fluxcandidate = secondlagrange.HasNeighbour(fH1Hybrid.fFluxMatId);
                // if the fluxelement found is the first flux element
                if(fluxcandidate == gelsideflux) {
                    // we have to find a larger (lower level) flux element
                    fluxcandidate = gelsideflux.HasLowerLevelNeighbour(fH1Hybrid.fFluxMatId);
                    if(!fluxcandidate) DebugStop();
                }
                {
                    TPZCompElSide celflux = fluxcandidate.Reference();
                    TPZCompElSide pressure = pressureinterface.Reference();
                    if(!celflux || !pressure) DebugStop();
                    int64_t index;
                    if(secondlagrange.Element()->Reference()) DebugStop();
                    new TPZMultiphysicsInterfaceElement(*mphys,secondlagrange.Element(),index,celflux,pressure);
#ifdef PZDEBUG
                    numcreated[secondlagrange.Element()->MaterialId()]++;
#endif
                }
            }
            {
                TPZCompElSide celflux = gelsideflux.Reference();
                TPZCompElSide pressure = pressureinterface.Reference();
                if(!celflux || !pressure) DebugStop();
                int64_t index;
#ifdef PZDEBUG
                numcreated[firstlagrange.Element()->MaterialId()]++;
#endif
                if(firstlagrange.Element()->Reference()) DebugStop();
                new TPZMultiphysicsInterfaceElement(*mphys,firstlagrange.Element(),index,celflux,pressure);
            }
        }
    }
    std::cout << __PRETTY_FUNCTION__ << "Number of computational interface elements created by material id\n";
    for(auto it : numcreated) std::cout << "Material id " << it.first << " number of elements created " << it.second << std::endl;
}

/// group and condense the elements
void TPZCreateMultiphysicsSpace::GroupandCondenseElements(TPZMultiphysicsCompMesh *cmesh)
{
    /// same procedure as hybridize hdiv
    int64_t nel = cmesh->NElements();
    TPZVec<int64_t> groupnumber(nel,-1);
    /// compute a groupnumber associated with each element
    AssociateElements(cmesh, groupnumber);
    std::map<int64_t, TPZElementGroup *> groupmap;
    //    std::cout << "Groups of connects " << groupindex << std::endl;
    for (int64_t el = 0; el<nel; el++) {
        int64_t groupnum = groupnumber[el];
        if(groupnum == -1) continue;
        auto iter = groupmap.find(groupnum);
        if (groupmap.find(groupnum) == groupmap.end()) {
            int64_t index;
            TPZElementGroup *elgr = new TPZElementGroup(*cmesh,index);
            groupmap[groupnum] = elgr;
            elgr->AddElement(cmesh->Element(el));
        }
        else
        {
            iter->second->AddElement(cmesh->Element(el));
        }
        //        std::cout << std::endl;
    }
    cmesh->ComputeNodElCon();
    int64_t nconnects = cmesh->NConnects();
    for (int64_t ic = 0; ic<nconnects; ic++) {
        TPZConnect &c = cmesh->ConnectVec()[ic];
        if(c.LagrangeMultiplier() == 5) c.IncrementElConnected();
    }
    nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *> (cel);
        if (elgr) {
            TPZCondensedCompEl *cond = new TPZCondensedCompEl(elgr);
            cond->SetKeepMatrix(false);
        }
    }
}

/// Find the neighbouring flux element
TPZCompEl *TPZCreateMultiphysicsSpace::FindFluxElement(TPZCompEl *wrapelement)
{
    TPZGeoEl *gel = wrapelement->Reference();
    int nsides = gel->NSides();
    TPZGeoElSide gelside(gel,nsides-1);
    TPZStack<TPZCompElSide> celstack;
    gelside.EqualLevelCompElementList(celstack, 0, 0);
    int nelstack = celstack.size();
    for (int st = 0; st<nelstack; st++) {
        TPZCompElSide celside = celstack[st];
        TPZCompEl *cel = celside.Element();
        TPZGeoEl *gelneigh = cel->Reference();
        int matid = gelneigh->MaterialId();
        if (matid == fH1Hybrid.fFluxMatId) {
            return cel;
        }
        if (fH1Hybrid.fHybridizeBCLevel == 1 && fBCMaterialIds.find(matid) != fBCMaterialIds.end()) {
            return cel;
        }
    }
    TPZCompElSide cellarge = gelside.LowerLevelCompElementList2(0);
    if(!cellarge) DebugStop();
    {
        TPZStack<TPZCompElSide> celstack;
        TPZGeoElSide gellarge = cellarge.Reference();
        if(gellarge.Element()->MaterialId() == fH1Hybrid.fFluxMatId)
        {
            return cellarge.Element();
        }
        gellarge.EqualLevelCompElementList(celstack, 0, 0);
        int nelst = celstack.size();
        for (int ist = 0; ist < nelst; ist++) {
            TPZCompElSide celside = celstack[ist];
            TPZGeoElSide gelside = celside.Reference();
            TPZGeoEl *candidate = gelside.Element();
            if(candidate->MaterialId() == fH1Hybrid.fFluxMatId)
            {
                return celside.Element();
            }
        }
    }
    return NULL;
}

/// Compute Periferal Material ids
// the material ids will be computed from a number whose modulus by base is zero
void TPZCreateMultiphysicsSpace::ComputePeriferalMaterialIds(int base)
{
    if(base < 2) base = 2;
    int max_matid = 0;
    int64_t nel = fGeoMesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if(!gel) continue;
        max_matid = std::max(max_matid,gel->MaterialId());
    }
    int remain = max_matid % base;
    int matid_base = max_matid-remain + base;
    fH1Hybrid.fFluxMatId = matid_base;
//        fH1Hybrid.fMatWrapId.first = matid_base+base;
//        fH1Hybrid.fMatWrapId.second = matid_base+base+1;
    fH1Hybrid.fMatWrapId = matid_base+base;
    fH1Hybrid.fLagrangeMatid.first = matid_base + 2*base;
    fH1Hybrid.fLagrangeMatid.second = matid_base + 2*base+1;
    fH1Hybrid.fSecondLagrangeMatid = matid_base + 3*base;
    fH1Hybrid.fInterfacePressure = matid_base + 4*base;
    
}

static void InsertNullMaterial(int matid, int dim, int nstate, TPZCompMesh *cmesh)
{
    TPZMaterial *mat = cmesh->FindMaterial(matid);
    if(mat) return;
    TPZNullMaterial *nullmat = new TPZNullMaterial(matid);
    nullmat->SetDimension(dim);
    nullmat->SetNStateVariables(nstate);
    cmesh->InsertMaterialObject(nullmat);
    
}
/// Insert the periferal material objects (for wrapmatid, fluxmatid and lagrange matid
void TPZCreateMultiphysicsSpace::InsertPeriferalMaterialObjects(TPZMultiphysicsCompMesh *mphys)
{
    if(fSpaceType == EH1Hybrid)
    {
        InsertNullMaterial(fH1Hybrid.fFluxMatId, fDimension-1, 1, mphys);
        InsertNullMaterial(fH1Hybrid.fMatWrapId, fDimension-1, 1, mphys);
    }
    else if (fSpaceType == EH1HybridSquared) {
        InsertNullMaterial(fH1Hybrid.fMatWrapId, fDimension-1, 1, mphys);
        InsertNullMaterial(fH1Hybrid.fFluxMatId, fDimension-1, 1, mphys);
        InsertNullMaterial(fH1Hybrid.fInterfacePressure, fDimension-1, 1, mphys);
    }
    else {
        DebugStop();
    }

}

void TPZCreateMultiphysicsSpace::InsertLagranceMaterialObjects(TPZMultiphysicsCompMesh *mphys)
{
    if(fSpaceType == EH1Hybrid)
    {
        TPZLagrangeMultiplier *lag1 = new TPZLagrangeMultiplier(fH1Hybrid.fLagrangeMatid.first, fDimension-1, 1);
        mphys->InsertMaterialObject(lag1);
        TPZLagrangeMultiplier *lag2 = new TPZLagrangeMultiplier(fH1Hybrid.fLagrangeMatid.second, fDimension-1, 1);
        lag2->SetMultiplier(-1.);
        mphys->InsertMaterialObject(lag2);
    }
    else if (fSpaceType == EH1HybridSquared) {
        TPZLagrangeMultiplier *lag1 = new TPZLagrangeMultiplier(fH1Hybrid.fLagrangeMatid.first, fDimension-1, 1);
        mphys->InsertMaterialObject(lag1);
        TPZLagrangeMultiplier *lag2 = new TPZLagrangeMultiplier(fH1Hybrid.fSecondLagrangeMatid, fDimension-1, 1);
        lag2->SetMultiplier(-1.);
        mphys->InsertMaterialObject(lag2);
    }
    else {
        DebugStop();
    }

}

/// Create geometric elements needed for the computational elements
void TPZCreateMultiphysicsSpace::AddGeometricWrapElements()
{
#ifdef PZDEBUG
    std::map<int,int> numcreated;
#endif
    int64_t nel = fGeoMesh->NElements();
    int dim = fGeoMesh->Dimension();
    for(int64_t el = 0; el<nel; el++)
    {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if(!gel || gel->HasSubElement() || gel->Dimension() != dim) continue;
        int nsides = gel->NSides();
        int side = 0;
        for(int d=0; d<dim-1; d++) side += gel->NSides(d);
        // loop over the sides of dimension dim-1
        for(; side < nsides-1; side++)
        {
            TPZGeoElSide gelside(gel,side);
            // we want to create side elements of type
            // first fMatWrapId
            TPZGeoElSide neighbour = gelside.Neighbour();
            int neighmat = neighbour.Element()->MaterialId();
#ifdef PZDEBUG
            {
                if(neighmat == fH1Hybrid.fMatWrapId)
                {
                    std::cout << __PRETTY_FUNCTION__ << " should be called only once!\n";
                    DebugStop();
                }
            }
#endif
            // if the neighbour is a boundary condition and no hybridization is applied
            // do not create the wrap layers
            bool HasBCNeighbour = (fBCMaterialIds.find(neighmat) != fBCMaterialIds.end());
            if(fH1Hybrid.fHybridizeBCLevel == 0 && HasBCNeighbour)
            {
                // no interface will be created between the element and a flux space
                continue;
            }
            // first fH1Hybrid.fMatWrapId
            TPZGeoElBC(gelside, fH1Hybrid.fMatWrapId);
            neighbour = gelside.Neighbour();
#ifdef PZDEBUG
            numcreated[fH1Hybrid.fMatWrapId]++;
            if(neighbour.Element()->MaterialId() != fH1Hybrid.fMatWrapId)
            {
                DebugStop();
            }
#endif
            if(fSpaceType == EH1Hybrid)
            {
                // then, depending on orientation fH1Hybrid.fLagrangeMatid.first or second
                int dir = gel->NormalOrientation(side);
                if(dir == 1)
                {
                    TPZGeoElBC(neighbour,fH1Hybrid.fLagrangeMatid.first);
#ifdef PZDEBUG
                    numcreated[fH1Hybrid.fLagrangeMatid.first]++;
#endif
                } else {
                    TPZGeoElBC(neighbour,fH1Hybrid.fLagrangeMatid.second);
#ifdef PZDEBUG
                    numcreated[fH1Hybrid.fLagrangeMatid.second]++;
#endif
                }
            } else if(fSpaceType == EH1HybridSquared) {
                TPZGeoElBC(neighbour,fH1Hybrid.fLagrangeMatid.first);
#ifdef PZDEBUG
                numcreated[fH1Hybrid.fLagrangeMatid.first]++;
#endif
                neighbour = neighbour.Neighbour();
                // if a second hybridization is applied
                // add a flux element fH1Hybrid.fFluxMatId
                // add a flux element
                if(!HasBCNeighbour || fH1Hybrid.fHybridizeBCLevel == 2)
                {
                    TPZGeoElBC(neighbour,fH1Hybrid.fFluxMatId);
#ifdef PZDEBUG
                    numcreated[fH1Hybrid.fFluxMatId]++;
#endif
                }
            } else {
                // only two cases handled so far
                DebugStop();
            }
        }
    }
    nel = fGeoMesh->NElements();
    for(int64_t el = 0; el<nel; el++)
    {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if(!gel || gel->HasSubElement() || gel->Dimension() != dim-1) continue;
        int matid = gel->MaterialId();
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        auto lagrange = fH1Hybrid.fLagrangeMatid;
        bool islagrange = (matid == lagrange.first || matid == lagrange.second);
        bool isflux = matid == fH1Hybrid.fFluxMatId;
        if(fSpaceType == EH1Hybrid && islagrange)
        {
            TPZGeoElSidePartition partition(gelside);
            if(partition.HasHigherLevelNeighbour(lagrange.first) ||
               partition.HasHigherLevelNeighbour(lagrange.second))
            {
                // create the flux element
                TPZGeoElBC(gelside,fH1Hybrid.fFluxMatId);
#ifdef PZDEBUG
                numcreated[fH1Hybrid.fFluxMatId]++;
#endif
            }
            else if(!gelside.HasNeighbour(fH1Hybrid.fFluxMatId))
            {
                // create the flux element
                TPZGeoElBC(gelside,fH1Hybrid.fFluxMatId);
#ifdef PZDEBUG
                numcreated[fH1Hybrid.fFluxMatId]++;
#endif
            }
        }
        if(fSpaceType == EH1HybridSquared && isflux)
        {
            // we need to include three layers of geometric elements
            // if there is fH1Hybrid.fInterfacePressure element as neighbour dont create
            if(gelside.HasNeighbour(fH1Hybrid.fInterfacePressure)) continue;
            // if there is a higher level connected element, dont create
            TPZGeoElSidePartition partition(gelside);
            if(partition.HasHigherLevelNeighbour(fH1Hybrid.fFluxMatId)) continue;
            // now we can create the elements
            TPZGeoElSide neighbour = gelside.Neighbour();
            int neighmat = neighbour.Element()->MaterialId();
            bool HasBCNeighbour = (fBCMaterialIds.find(neighmat) != fBCMaterialIds.end());

            // lagrange element
            TPZGeoElBC gbc1(gelside, fH1Hybrid.fSecondLagrangeMatid);
            // if the flux is neighbour, the lagrange multiplier will be between the flux and the
            // boundary condition
#ifdef PZDEBUG
            numcreated[fH1Hybrid.fSecondLagrangeMatid]++;
#endif
            if(HasBCNeighbour) continue;
            neighbour = gelside.Neighbour();
            // pressure element
            TPZGeoElBC gbc2(neighbour,fH1Hybrid.fInterfacePressure);
            neighbour = neighbour.Neighbour();
            // lagrange element
            TPZGeoElBC gbc3(neighbour,fH1Hybrid.fSecondLagrangeMatid);
#ifdef PZDEBUG
            numcreated[fH1Hybrid.fInterfacePressure]++;
            numcreated[fH1Hybrid.fSecondLagrangeMatid]++;
#endif
        }
    }
#ifdef PZDEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    for(auto it : numcreated)
    {
        std::cout << "Material id " << it.first << " number of elements created " << it.second << std::endl;
    }
#endif
}

/// Associate elements with a volumetric element
// elementgroup[el] = index of the element with which the element should be grouped
// this method only gives effective result for hybridized hdiv meshes
void TPZCreateMultiphysicsSpace::AssociateElements(TPZCompMesh *cmesh, TPZVec<int64_t> &elementgroup)
{
    int64_t nel = cmesh->NElements();
    elementgroup.Resize(nel, -1);
    elementgroup.Fill(-1);
    int64_t nconnects = cmesh->NConnects();
    TPZVec<int64_t> groupindex(nconnects, -1);
    int dim = cmesh->Dimension();
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference() || cel->Reference()->Dimension() != dim) {
            continue;
        }
        elementgroup[cel->Index()] = cel->Index();
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        for (auto cindex : connectlist) {
#ifdef PZDEBUG
            if (groupindex[cindex] != -1) {
                DebugStop();
            }
#endif
            groupindex[cindex] = cel->Index();
        }
    }
    std::cout << "Groups of connects " << groupindex << std::endl;
    int numloops = 1;
    if(fSpaceType == EH1HybridSquared) numloops = 2;
    // this loop will associate a first layer of interface elements to the group
    // this loop will associate the wrap elements with the group
    // if HybridSquared the connects of interface elements with matid fLagrangeMatId will be added to the group
    // in the second pass :
        // incorporate the flux elements in the group
        // incorporate the interface elements to the pressure lagrange DOFs in the group
    for (int iloop = 0; iloop < numloops; iloop++) for (TPZCompEl *cel : cmesh->ElementVec())
    {
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        int matid = cel->Reference()->MaterialId();
        int64_t celindex = cel->Index();
        std::cout << "Analysing element " << celindex << " matid " << matid;
        std::cout << " connect list " << connectlist << std::endl;
        TPZVec<int> connectgroup(connectlist.size());
        for(int i=0; i<connectlist.size(); i++) connectgroup[i] = groupindex[connectlist[i]];
        std::cout << "groupindexes " << connectgroup << std::endl;
        int64_t groupfound = -1;
        for (auto cindex : connectlist) {
            if (groupindex[cindex] != -1) {
                elementgroup[celindex] = groupindex[cindex];
                if(groupfound != -1 && groupfound != groupindex[cindex])
                {
                    DebugStop();
                }
//                if(groupfound == -1)
//                {
//                    std::cout << " added to " << groupindex[cindex];
//                }
                groupfound = groupindex[cindex];
            }
        }
        if(fSpaceType == EH1HybridSquared && matid == fH1Hybrid.fLagrangeMatid.first)
        {
            std::cout << "Changing connect indexes group for element " << celindex;
            for(auto cindex : connectlist)
            {
                std::cout << " cindex " << cindex << " from " << groupindex[cindex] << " to " << groupfound << std::endl;
                groupindex[cindex] = groupfound;
            }
        }

//        std::cout << std::endl;
    }
}

