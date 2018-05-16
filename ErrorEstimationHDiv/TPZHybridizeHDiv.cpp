//
//  TPZHybridizeHDiv.cpp
//  ErrorEstimate
//
//  Created by Philippe Devloo on 16/05/18.
//

#include "TPZHybridizeHDiv.h"
#include "pzconnect.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZMaterial.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "pzintel.h"
#include "pzgeoelbc.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"


/// split the connect between two neighbouring elements
void TPZHybridizeHDiv::SplitConnects(const TPZCompElSide &left, const TPZCompElSide &right, TPZVec<TPZCompMesh *> &meshvec)
{
    TPZCompMesh *fluxmesh = meshvec[0];
    TPZCompMesh *pressuremesh = meshvec[1];
    TPZGeoElSide gleft(left.Reference());
    TPZGeoElSide gright(right.Reference());
    TPZInterpolatedElement *intelleft = dynamic_cast<TPZInterpolatedElement *>(left.Element());
    TPZInterpolatedElement *intelright = dynamic_cast<TPZInterpolatedElement *>(right.Element());
    intelleft->SetSideOrient(left.Side(), 1);
    intelright->SetSideOrient(right.Side(), 1);
    TPZConnect &cleft = intelleft->SideConnect(0, left.Side());
    if (cleft.HasDependency()) {
        cleft.RemoveDepend();
    }
    int64_t index = fluxmesh->ConnectVec().AllocateNewElement();
    fluxmesh->ConnectVec()[index] = cleft;
    cleft.DecrementElConnected();
    fluxmesh->ConnectVec()[index].ResetElConnected();
    fluxmesh->ConnectVec()[index].IncrementElConnected();
    
    int rightlocindex = intelright->SideConnectLocId(0, right.Side());
    intelright->SetConnectIndex(rightlocindex, index);
    int sideorder = cleft.Order();
    // create HDivBound on the sides of the elements
    {
        intelright->Reference()->ResetReference();
        intelleft->LoadElementReference();
        TPZGeoElBC gbc(gleft,HDivWrapMatid);
        int64_t index;
        TPZCompEl *wrap = fluxmesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *fluxmesh, index);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(wrap);
        int wrapside = gbc.CreatedElement()->NSides()-1;
        intel->SetSideOrient(wrapside, 1);
        intelleft->Reference()->ResetReference();
        wrap->Reference()->ResetReference();
    }
    {
        intelleft->Reference()->ResetReference();
        intelright->LoadElementReference();
        TPZGeoElBC gbc(gleft,HDivWrapMatid);
        int64_t index;
        TPZCompEl *wrap = fluxmesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *fluxmesh, index);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(wrap);
        int wrapside = gbc.CreatedElement()->NSides()-1;
        intel->SetSideOrient(wrapside, 1);
        intelright->Reference()->ResetReference();
        wrap->Reference()->ResetReference();
    }
    {
        TPZGeoElBC gbc(gleft,LagrangeInterface);
        int64_t index;
        TPZCompEl *press = pressuremesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *pressuremesh, index);
        press->Reference()->ResetReference();
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(press);
        intel->SetPreferredOrder(sideorder);
    }
    intelleft->LoadElementReference();
    intelright->LoadElementReference();
}

/// find an element which shares the connect and has the same dimension
// else return 0
static TPZCompElSide RightElement(TPZInterpolatedElement *intel, int side)
{
    bool isrestrained = false;
    {
        TPZConnect &c = intel->SideConnect(0, side);
        if (c.HasDependency()) {
            isrestrained = true;
        }
    }
    TPZGeoEl *gel = intel->Reference();
    TPZGeoElSide gelside(gel,side);
    
    if (isrestrained == true)
    {
        /// if the side is restrained we will hybridize between the element and the larger element
        TPZCompElSide celside = gelside.LowerLevelCompElementList2(1);
        if(!celside) DebugStop();
        TPZGeoEl *neigh = celside.Element()->Reference();
        /// we assume that a larger element should not be a boundary element
        if (neigh->Dimension() != gel->Dimension()) {
            DebugStop();
        }
        return celside;
    }
    else
    {
        // if the connect is not restrained
        // - there should be only one neighbour
        //   if there is more than one neighbour the element side has already been hybridized
        // - the neighbour should be of the same dimension
        //   if the neighbour is of lower dimension it is a boundary element
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 1, 0);
        if (celstack.size() == 1) {
            TPZGeoEl *neigh = celstack[0].Element()->Reference();
            if (neigh->Dimension() == gel->Dimension()) {
                return celstack[0];
            }
        }
    }
    return TPZCompElSide();
}


/// split the connects between flux elements and create a dim-1 pressure element
void TPZHybridizeHDiv::HybridizeInternalSides(TPZVec<TPZCompMesh *> &meshvec)
{
    TPZCompMesh *fluxmesh = meshvec[0];
    TPZGeoMesh *gmesh = fluxmesh->Reference();
    int dim = gmesh->Dimension();
    gmesh->ResetReference();
    fluxmesh->LoadReferences();
    int64_t nel = fluxmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fluxmesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel || intel->Reference()->Dimension() != dim) {
            continue;
        }
        // loop over the side of dimension dim-1
        TPZGeoEl *gel = intel->Reference();
        for (int side = gel->NCornerNodes(); side < gel->NSides()-1; side++) {
            if (gel->SideDimension(side) != dim-1) {
                continue;
            }
            TPZCompElSide celside(intel,side);
            TPZCompElSide neighcomp = RightElement(intel, side);
            if (neighcomp) {
                SplitConnects(celside, neighcomp, meshvec);
            }
        }
    }
}

void TPZHybridizeHDiv::CreateInterfaceElements(TPZCompMesh *cmesh, TPZVec<TPZCompMesh *> &meshvec)
{
    TPZCompMesh *pressuremesh = meshvec[1];
    int dim = pressuremesh->Dimension();
    TPZVec<TPZCompEl *> celpressure(pressuremesh->NElements(),0);
    for (int64_t el=0; el<pressuremesh->NElements(); el++) {
        TPZCompEl *cel = pressuremesh->Element(el);
        if (cel && cel->Reference() && cel->Reference()->Dimension() == dim-1) {
            celpressure[el] = cel;
        }
    }
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    for(auto cel:celpressure)
    {
        if(!cel) continue;
        TPZStack<TPZCompElSide> celstack;
        TPZGeoEl *gel = cel->Reference();
        TPZCompEl *mphysics = gel->Reference();
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        gelside.EqualLevelCompElementList(celstack, 1, 0);
        for (auto &celstackside:celstack) {
            if (celstackside.Reference().Element()->Dimension() == dim-1) {
                TPZCompElSide celside(mphysics,gel->NSides()-1);
                TPZGeoElBC gbc(gelside,InterfaceMatid);
                int64_t index;
                TPZMultiphysicsInterfaceElement *intface = new TPZMultiphysicsInterfaceElement(*cmesh,gbc.CreatedElement(),index,celside,celstackside);
            }
        }
    }
}

TPZCompMesh *TPZHybridizeHDiv::CreateMultiphysicsMesh(const TPZVec<TPZMaterial *> &matvec, TPZVec<TPZCompMesh *> &meshvector)
{
    TPZCompMesh *cmesh = new TPZCompMesh(meshvector[0]->Reference());
    for (int imat = 0; imat<matvec.size(); imat++) {
        cmesh->MaterialVec()[matvec[imat]->Id()] = matvec[imat];
    }
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    cmesh->AutoBuild();
    
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh);
    cmesh->LoadReferences();
    CreateInterfaceElements(cmesh, meshvector);
    GroupElements(cmesh);
    return cmesh;
}

/// group and condense the elements
void TPZHybridizeHDiv::GroupElements(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    int64_t nconnects = cmesh->NConnects();
    TPZVec<TPZElementGroup *> groupindex(nconnects,0);
    int dim = cmesh->Dimension();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel || !cel->Reference() || cel->Reference()->Dimension() != dim) {
            continue;
        }
        int64_t index;
        TPZElementGroup *elgr = new TPZElementGroup(*cmesh,index);
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        for (auto cindex:connectlist)
        {
#ifdef PZDEBUG
            if (groupindex[cindex] != 0) {
                DebugStop();
            }
#endif
            groupindex[cindex] = elgr;
        }
        elgr->AddElement(cel);
    }
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        for (auto cindex:connectlist) {
            if (groupindex[cindex] != 0) {
                groupindex[cindex]->AddElement(cel);
                break;
            }
        }
    }
    cmesh->ComputeNodElCon();
    nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cmesh->Element(el));
        if (elgr) {
            TPZCondensedCompEl *cond = new TPZCondensedCompEl(elgr);
            cond->SetKeepMatrix(false);
        }
    }
}

