//
//  TPZCreateMultiPhysicsSpace.cpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 13/07/19.
//

#include "TPZCreateMultiPhysicsSpace.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoelbc.h"
#include "pzintel.h"
#include "TPZNullMaterial.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZHybridizeHDiv.h"

TPZCreateMultiPhysicsSpace::TPZCreateMultiPhysicsSpace(TPZGeoMesh *gmesh) : fGeoMesh(gmesh) {
    fDimension = gmesh->Dimension();
}

/// copy constructor
TPZCreateMultiPhysicsSpace::TConfigH1Hybrid::TConfigH1Hybrid(const TConfigH1Hybrid &copy)
{
    
}

/// copy operator
TPZCreateMultiPhysicsSpace::TConfigH1Hybrid &TPZCreateMultiPhysicsSpace::TConfigH1Hybrid::operator=(const TConfigH1Hybrid &copy)
{
    return *this;
}


/// copy constructor
TPZCreateMultiPhysicsSpace::TPZCreateMultiPhysicsSpace(const TPZCreateMultiPhysicsSpace &copy)
{
    
}

/// = operator
TPZCreateMultiPhysicsSpace & TPZCreateMultiPhysicsSpace::operator=(const TPZCreateMultiPhysicsSpace &copy)
{
    return *this;
}

/// Indicate to create Hybridized H1 meshes
void TPZCreateMultiPhysicsSpace::SetH1Hybridized(const TConfigH1Hybrid &config)
{
    fSpaceType = EH1Hybrid;
    fH1Hybrid = config;
}

/// create meshes and elements for all geometric elements
void TPZCreateMultiPhysicsSpace::CreateAtomicMeshes(TPZVec<TPZCompMesh *> &meshvec)
{
    TPZCompMesh *pressure = CreatePressureMesh();
    TPZCompMesh *fluxmesh = CreateFluxMesh();
    TPZCompMesh *gspace = new TPZCompMesh(fGeoMesh);
    InsertNullSpaceMaterialIds(gspace);
    gspace->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    gspace->SetDefaultOrder(0);
    gspace->AutoBuild();
    TPZCompMesh *average = new TPZCompMesh(fGeoMesh);
    InsertNullSpaceMaterialIds(average);
    average->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    average->SetDefaultOrder(0);
    average->AutoBuild();

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
static bool HasLagrangeNeighbour(TPZGeoElSide gelside, int lagrangematid)
{
    TPZGeoElSide neighbour(gelside.Neighbour());
    while(neighbour != gelside)
    {
        if(neighbour.Element()->MaterialId() == lagrangematid) return true;
        neighbour = neighbour.Neighbour();
    }
    return false;
}
/// create the geometric elements for the lagrange multipliers
// these elements will go with the largest H1 element
void TPZCreateMultiPhysicsSpace::CreateLagrangeGeometricElements(TPZCompMesh *H1mesh)
{
    TPZGeoMesh *gmesh = H1mesh->Reference();
    gmesh->ResetReference();
    H1mesh->LoadReferences();
    int64_t nel = H1mesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = H1mesh->Element(el);
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
            if(ndif > 1) DebugStop();
            if(celstack_eq.size())
            {
                // figure out if there is a lagrange material element neighbour
                if(HasLagrangeNeighbour(gelside, fH1Hybrid.fFluxMatId) == false)
                {
                    TPZGeoElBC(gelside,fH1Hybrid.fFluxMatId);
                }
            }
            if(celstack_sm.size())
            {
                TPZGeoElBC(gelside,fH1Hybrid.fFluxMatId);
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
void TPZCreateMultiPhysicsSpace::CreatePressureBoundaryElements(TPZCompMesh *pressure)
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
void TPZCreateMultiPhysicsSpace::InsertPressureMaterialIds(TPZCompMesh *pressure)
{
    for (auto matid:fMaterialIds) {
        TPZNullMaterial *nullmat = new TPZNullMaterial(matid);
        nullmat->SetDimension(fDimension);
        nullmat->SetNStateVariables(1);
        pressure->InsertMaterialObject(nullmat);
    }
    if(fSpaceType == EH1Hybrid && fH1Hybrid.fHybridizeBC == false)
    {
        for (auto matid:fMaterialIds) {
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
void TPZCreateMultiPhysicsSpace::InsertFluxMaterialIds(TPZCompMesh *fluxmesh)
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
        for (auto matid:fMaterialIds) {
            TPZNullMaterial *nullmat = new TPZNullMaterial(matid);
            nullmat->SetDimension(fDimension);
            nullmat->SetNStateVariables(1);
            fluxmesh->InsertMaterialObject(nullmat);
        }
    }

}

/// insert materialids for the null space
void TPZCreateMultiPhysicsSpace::InsertNullSpaceMaterialIds(TPZCompMesh *nullspace)
{
    for (auto matid:fMaterialIds) {
        TPZNullMaterial *nullmat = new TPZNullMaterial(matid);
        nullmat->SetDimension(fDimension);
        nullmat->SetNStateVariables(1);
        nullspace->InsertMaterialObject(nullmat);
    }
}

/// Create the pressure mesh
TPZCompMesh *TPZCreateMultiPhysicsSpace::CreatePressureMesh()
{
    // create the pressure mesh
    TPZCompMesh *pressure = new TPZCompMesh(fGeoMesh);
    InsertPressureMaterialIds(pressure);
    pressure->ApproxSpace().SetAllCreateFunctionsContinuous();
    pressure->ApproxSpace().CreateDisconnectedElements(true);
    pressure->SetDefaultOrder(fDefaultPOrder);
    pressure->AutoBuild(fMaterialIds);
    CreatePressureBoundaryElements(pressure);
    return pressure;
}

/// Create the flux mesh
TPZCompMesh *TPZCreateMultiPhysicsSpace::CreateFluxMesh()
{
    TPZCompMesh *fluxmesh = new TPZCompMesh(fGeoMesh);
    InsertFluxMaterialIds(fluxmesh);
    fluxmesh->ApproxSpace().SetAllCreateFunctionsHDiv(fDimension);
    fluxmesh->SetDefaultOrder(fDefaultLagrangeOrder);
    fluxmesh->AutoBuild();
    return fluxmesh;
}

/// add interface elements to the multiphysics space
void TPZCreateMultiPhysicsSpace::AddInterfaceElements(TPZMultiphysicsCompMesh *mphys)
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
void TPZCreateMultiPhysicsSpace::GroupandCondenseElements(TPZMultiphysicsCompMesh *mphys)
{
    /// same procedure as hybridize hdiv
    TPZHybridizeHDiv::GroupandCondenseElements(mphys);
}

/// Find the neighbouring flux element
TPZCompEl *TPZCreateMultiPhysicsSpace::FindFluxElement(TPZCompEl *wrapelement)
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
        if (gelneigh->MaterialId() == fH1Hybrid.fFluxMatId) {
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
