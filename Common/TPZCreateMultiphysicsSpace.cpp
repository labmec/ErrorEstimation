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
void TPZCreateMultiphysicsSpace::CreateAtomicMeshes(TPZVec<TPZCompMesh *> &meshvec)
{
    TPZCompMesh *pressure = CreatePressureMesh();
    CreateLagrangeGeometricElements(pressure);
    TPZCompMesh *fluxmesh = CreateFluxMesh();
    TPZCompMesh *gspace = new TPZCompMesh(fGeoMesh);
    {
        InsertNullSpaceMaterialIds(gspace);
        gspace->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
        gspace->SetDefaultOrder(0);
        gspace->AutoBuild();
        int64_t nconnects = gspace->NConnects();
        for (int ic = 0; ic<nconnects; ic++) {
            gspace->ConnectVec()[ic].SetLagrangeMultiplier(1);
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
            average->ConnectVec()[ic].SetLagrangeMultiplier(4);
        }
}

    if(fSpaceType == EH1Hybrid)
    {
        meshvec.Resize(4,0);
        meshvec[0] = pressure;
        meshvec[1] = fluxmesh;
        meshvec[2] = gspace;
        meshvec[3] = average;
    }
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
    if (fSpaceType != EH1Hybrid) {
        DebugStop();
    }
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
        if (gel->Dimension() != fDimension) {
            DebugStop();
        }
        int nsides = gel->NSides();
        int ncorner = gel->NCornerNodes();
        for (int side = ncorner; side < nsides; side++) {
            if(gel->SideDimension(side) != fDimension-1) continue;
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = HasBCNeighbour(gelside, fBCMaterialIds);
            pressure->SetDefaultOrder(porder);
            TPZGeoEl * gelcreate = 0;
            int dir = gel->NormalOrientation(side);
            int wrapmat = 0;
            if(dir == 1)
            {
                wrapmat = fH1Hybrid.fMatWrapId.first;
            }
            else if(dir == -1)
            {
                wrapmat = fH1Hybrid.fMatWrapId.second;
            }
            else
            {
                DebugStop();
            }
            // if neighbour exists then create the wrap conditionally
            if(neighbour && fH1Hybrid.fHybridizeBC == false)
            {
                gelcreate = neighbour.Element();
            }
            else
            {
                // create a geometric element neighbouring to the element
                TPZGeoElBC gbc(gelside,wrapmat);
                gelcreate = gbc.CreatedElement();
            }
            if(gelcreate)
            {
                // load the element reference so that the created element will share the connects
                cel->LoadElementReference();
                int64_t index;
                TPZCompEl *bc_cel = pressure->ApproxSpace().CreateCompEl(gelcreate, *pressure, index);
                // reset the references so that future elements will not share connects
                gel->ResetReference();
                bc_cel->Reference()->ResetReference();
            }
        }
    }
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
    if(fSpaceType == EH1Hybrid && fH1Hybrid.fHybridizeBC == false)
    {
        for (auto matid:fBCMaterialIds) {
            TPZNullMaterial *nullmat = new TPZNullMaterial(matid);
            nullmat->SetDimension(fDimension);
            nullmat->SetNStateVariables(1);
            pressure->InsertMaterialObject(nullmat);
        }
    }
    {
        TPZNullMaterial *nullmat = new TPZNullMaterial(fH1Hybrid.fMatWrapId.first);
        nullmat->SetDimension(fDimension);
        nullmat->SetNStateVariables(1);
        pressure->InsertMaterialObject(nullmat);
    }
    {
        TPZNullMaterial *nullmat = new TPZNullMaterial(fH1Hybrid.fMatWrapId.second);
        nullmat->SetDimension(fDimension);
        nullmat->SetNStateVariables(1);
        pressure->InsertMaterialObject(nullmat);
    }
}

/// insert flux material ids
void TPZCreateMultiphysicsSpace::InsertFluxMaterialIds(TPZCompMesh *fluxmesh)
{
    if (fSpaceType == EH1Hybrid) {
        int matid = fH1Hybrid.fFluxMatId;
        TPZNullMaterial *nullmat = new TPZNullMaterial(matid);
        nullmat->SetDimension(fDimension-1);
        nullmat->SetNStateVariables(1);
        fluxmesh->InsertMaterialObject(nullmat);
    }
    if(fSpaceType == EH1Hybrid && fH1Hybrid.fHybridizeBC == true)
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
    CreatePressureBoundaryElements(pressure);
    int64_t nelem = pressure->NElements();
    if(fDimension != 2)
    {
        std::cout << "fix me\n";
        DebugStop();
    }
    for (int64_t el = 0; el<nelem; el++) {
        TPZCompEl *cel = pressure->Element(el);
        int nconnects = cel->NConnects();
        for (int ic=2; ic<nconnects; ic++) {
            cel->Connect(ic).SetLagrangeMultiplier(2);
        }
    }
    return pressure;
}

/// Create the flux mesh
TPZCompMesh *TPZCreateMultiphysicsSpace::CreateFluxMesh()
{
    TPZCompMesh *fluxmesh = new TPZCompMesh(fGeoMesh);
    InsertFluxMaterialIds(fluxmesh);
    fluxmesh->ApproxSpace().SetAllCreateFunctionsHDiv(fDimension);
    fluxmesh->SetDefaultOrder(fDefaultLagrangeOrder);
    fluxmesh->AutoBuild();
    int64_t nconnects = fluxmesh->NConnects();
    for (int ic=0; ic<nconnects; ic++) {
        fluxmesh->ConnectVec()[ic].SetLagrangeMultiplier(3);
    }
    return fluxmesh;
}

/// add interface elements to the multiphysics space
void TPZCreateMultiphysicsSpace::AddInterfaceElements(TPZMultiphysicsCompMesh *mphys)
{
    TPZGeoMesh *gmesh = mphys->Reference();
    gmesh->ResetReference();
    mphys->LoadReferences();
    int64_t nel = mphys->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = mphys->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        if(matid != fH1Hybrid.fMatWrapId.first && matid != fH1Hybrid.fMatWrapId.second) continue;
        TPZCompEl *fluxel = FindFluxElement(cel);
        TPZGeoEl *fluxgel = fluxel->Reference();
        // determine if the interface should be positive or negative...
        int interfacematid = 0;
        if (matid == fH1Hybrid.fMatWrapId.first) {
            interfacematid = fH1Hybrid.fLagrangeMatid.first;
        }
        else
        {
            interfacematid = fH1Hybrid.fLagrangeMatid.second;
        }
        int64_t index;
        TPZCompElSide celwrap(cel,gel->NSides()-1);
        TPZCompElSide fluxside(fluxel,fluxgel->NSides()-1);
        TPZGeoElBC gbc(gel,gel->NSides()-1,interfacematid);
        new TPZMultiphysicsInterfaceElement(*mphys,gbc.CreatedElement(),index,celwrap,fluxside);
    }
}

/// group and condense the elements
void TPZCreateMultiphysicsSpace::GroupandCondenseElements(TPZMultiphysicsCompMesh *cmesh)
{
    /// same procedure as hybridize hdiv
    int64_t nel = cmesh->NElements();
    TPZVec<int64_t> groupnumber(nel,-1);
    /// compute a groupnumber associated with each element
    TPZHybridizeHDiv::AssociateElements(cmesh, groupnumber);
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
        if(c.LagrangeMultiplier() == 4) c.IncrementElConnected();
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
        if (fH1Hybrid.fHybridizeBC == true && fBCMaterialIds.find(matid) != fBCMaterialIds.end()) {
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
    if (fSpaceType == EH1Hybrid) {
        fH1Hybrid.fFluxMatId = matid_base;
        fH1Hybrid.fMatWrapId.first = matid_base+base;
        fH1Hybrid.fMatWrapId.second = matid_base+base+1;
        fH1Hybrid.fLagrangeMatid.first = matid_base + 2*base;
        fH1Hybrid.fLagrangeMatid.second = matid_base + 2*base+1;
    }
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
        InsertNullMaterial(fH1Hybrid.fMatWrapId.first, fDimension-1, 1, mphys);
        InsertNullMaterial(fH1Hybrid.fMatWrapId.second, fDimension-1, 1, mphys);
    }
    TPZLagrangeMultiplier *lag1 = new TPZLagrangeMultiplier(fH1Hybrid.fLagrangeMatid.first, fDimension-1, 1);
    mphys->InsertMaterialObject(lag1);
    TPZLagrangeMultiplier *lag2 = new TPZLagrangeMultiplier(fH1Hybrid.fLagrangeMatid.second, fDimension-1, 1);
    lag2->SetMultiplier(-1.);
    mphys->InsertMaterialObject(lag2);

}

