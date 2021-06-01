
//
//  TPZHDivErrorEstimator.cpp
//  ErrorEstimateHDiv
//
//  Created by Philippe Devloo on 10/06/18.
//

#include <Material/TPZHDivErrorEstimateMaterial.h>
#include "TPZElementMatrixT.h"
#include "TPZHDivErrorEstimator.h"
#include "TPZGeoElSideAncestors.h"
#include "TPZGeoElSidePartition.h"
#include "TPZInterfaceEl.h"
#include "TPZMixedHdivErrorEstimate.h"
#include "TPZNullMaterial.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZAnalysis.h"
#include "TPZBndCond.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"
#include "pzelmat.h"
#include "pzintel.h"
//#include "pzmat1dlin.h"
#include "pzsubcmesh.h"
#include "pzstepsolver.h"
#include "TPZVTKGeoMesh.h"
#include "pzmultiphysicscompel.h"
#include "TPZHDivErrorEstimateMaterial.h"
#include "TPZCompMeshTools.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("HDivErrorEstimator"));
#endif

static void VerifyConnectConsistency(TPZCompMesh *cmesh) {
    cmesh->ComputeNodElCon();

    for (int64_t icon = 0; icon < cmesh->NConnects(); icon++) {
        TPZConnect &con = cmesh->ConnectVec()[icon];
        if (con.NElConnected() == 0 && con.HasDependency()) {
            con.RemoveDepend();
        }
    }
}

TPZHDivErrorEstimator::~TPZHDivErrorEstimator() {
    TPZVec<TPZCompMesh *> meshvec = fPostProcMesh.MeshVector();
    
    fPostProcMesh.Reference()->ResetReference();
    for (int i = 0; i < 2; i++) {
        if (!meshvec[i]) continue;
        TPZCompMesh *cmesh = meshvec[i];
        cmesh->SetReference(nullptr);
        for (int64_t ic = 0; ic < meshvec[i]->NConnects(); ic++) {
            TPZConnect &c = meshvec[i]->ConnectVec()[ic];
            c.RemoveDepend();
        }
        delete meshvec[i];
    }
    for (int64_t ic = 0; ic < fPostProcMesh.NConnects(); ic++) {
        TPZConnect &c = fPostProcMesh.ConnectVec()[ic];
        c.RemoveDepend();
    }
}

/// compute the element errors comparing the reconstructed solution based on average pressures
/// with the original solution
void TPZHDivErrorEstimator::ComputeErrors(TPZVec<REAL>&errorVec, TPZVec<REAL>& elementErrors, std::string& vtkPath) {
    TPZLinearAnalysis an(&fPostProcMesh, false);
    
    if (fExact) {
        an.SetExact(fExact->ExactSolution());
    }

    int64_t nErrorCols = 6;
    errorVec.resize(nErrorCols);
    for (int64_t i = 0; i < nErrorCols; i++) {
        errorVec[i] = 0;
    }
    
    int64_t nelem = fPostProcMesh.NElements();
    fPostProcMesh.LoadSolution(fPostProcMesh.Solution());
    fPostProcMesh.ExpandSolution();
    fPostProcMesh.ElementSolution().Redim(nelem, 5);
    for(int64_t el = 0; el<nelem; el++)
    {
        TPZCompEl *cel = fPostProcMesh.Element(el);
        TPZSubCompMesh *subc = dynamic_cast<TPZSubCompMesh *>(cel);
        if(subc)
        {
            int64_t nelsub = subc->NElements();
            subc->ElementSolution().Redim(nelsub, 5);
        }
    }
    
    // Calculate error and store in element solution
    bool store_error = true;
    an.PostProcessError(errorVec, store_error);
    
    TPZCompMeshTools::UnCondensedElements(&fPostProcMesh);
    TPZCompMeshTools::UnGroupElements(&fPostProcMesh);

    if (fExact) ComputeEffectivityIndices();

    if (!vtkPath.empty()) {
        PostProcessing(an, vtkPath);
    }

    elementErrors.resize(fPostProcMesh.Reference()->NElements());
    for (REAL & elementerror : elementErrors) {
        elementerror = 0;
    }

    TPZFMatrix<STATE> &elsol = fPostProcMesh.ElementSolution();
    for (int64_t i = 0; i < nelem; i++) {
        TPZCompEl *cel = fPostProcMesh.Element(i);
        if (!cel) continue;
        TPZGeoEl* gel = fPostProcMesh.Element(i)->Reference();
        if (!gel) continue;
        elementErrors[gel->Index()] = elsol(i, 3);
    }
    
}

void TPZHDivErrorEstimator::PostProcessing(TPZAnalysis &an, std::string &out) {

    TPZMaterial *mat = fPostProcMesh.FindMaterial(1);
    int varindex = -1;
    if (mat) varindex = mat->VariableIndex("PressureFem");
    if (varindex != -1) {
        TPZStack<std::string> scalnames, vecnames;
        if (fExact) {
            scalnames.Push("PressureExact");
            scalnames.Push("PressureErrorExact");
            scalnames.Push("EnergyErrorExact");
            scalnames.Push("PressureEffectivityIndex");
            scalnames.Push("EnergyEffectivityIndex");
            vecnames.Push("FluxExact");
        }
        scalnames.Push("PressureFem");
        scalnames.Push("PressureReconstructed");
        scalnames.Push("PressureErrorEstimate");
        scalnames.Push("EnergyErrorEstimate");
        vecnames.Push("FluxFem");
        vecnames.Push("FluxReconstructed");
        scalnames.Push("POrder");
        
        int dim = fPostProcMesh.Reference()->Dimension();

        an.DefineGraphMesh(dim, scalnames, vecnames, out);
        an.PostProcess(0, dim);
    }
    else {
        std::cout << __PRETTY_FUNCTION__ << "\nPost Processing variable not found!\n";
    }
}

// a method for generating the HDiv mesh
TPZCompMesh *TPZHDivErrorEstimator::CreateFluxMesh()
{
    return fOriginal->MeshVector()[0]->Clone();
}

// a method for creating the pressure mesh
TPZCompMesh *TPZHDivErrorEstimator::CreatePressureMesh() {
    if (fPostProcesswithHDiv) {
        TPZCompMesh *pressure = fOriginal->MeshVector()[1]->Clone();
        return pressure;
    }
    
    // For H1 reconstruction, we need to build BC materials
    else {
        TPZCompMesh *mult = fOriginal;
        TPZCompMesh *pressureMesh = fOriginal->MeshVector()[1]->Clone();
        TPZGeoMesh *gmesh = pressureMesh->Reference();
        gmesh->ResetReference();
        int dim = gmesh->Dimension();
        
        // Delete compels of dimension dim - 1
        for (int64_t el = 0; el < pressureMesh->NElements(); el++) {
            TPZCompEl *cel = pressureMesh->Element(el);
            if (!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if (!gel) DebugStop();
            if (gel->Dimension() == dim - 1) {
                delete cel;
            }
        }
        
        pressureMesh->ComputeNodElCon();
        pressureMesh->CleanUpUnconnectedNodes();
        pressureMesh->LoadReferences();
        pressureMesh->ApproxSpace().SetAllCreateFunctionsContinuous();
        pressureMesh->ApproxSpace().CreateDisconnectedElements(true);

        // Insert BC materials in pressure reconstruction mesh
        std::set<int> bcMatIDs = GetBCMatIDs(&fPostProcMesh);
        for (auto bcID : bcMatIDs) {
            TPZMaterial *mat = mult->FindMaterial(bcID);
            TPZBndCondT<STATE> *bc = dynamic_cast<TPZBndCondT<STATE> *>(mat);
            if (!bc) DebugStop();

            int volumetricMatId = bc->Material()->Id();
            TPZMatBase<STATE> *pressuremat = dynamic_cast<TPZMatBase<STATE>*>(pressureMesh->FindMaterial(volumetricMatId));
            if (!pressuremat) DebugStop();

            TPZBndCondT<STATE> *newbc = pressuremat->CreateBC(pressuremat, bc->Id(), bc->Type(), bc->Val1(), bc->Val2());
            if (bc->HasForcingFunctionBC()) {
                newbc->SetForcingFunctionBC(bc->ForcingFunctionBC());
            }
            pressureMesh->InsertMaterialObject(newbc);
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
        
        // Create BC condition elements. We reset references at each step to allow for a discontinuous mesh globally
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
            pressureMesh->SetDefaultOrder(order);

            int64_t id;
            TPZCompEl *bcCel = pressureMesh->CreateCompEl(bcGeoEl, id);
            // Reset references so that future elements will not share connects with these two elements
            bcCel->Reference()->ResetReference();
            neighCel->Reference()->ResetReference();
        }
        
#ifdef PZDEBUG
        {
            std::ofstream outTXT("PostProcPressureMesh.txt");
            std::ofstream outVTK("PostProcPressureMesh.vtk");
            pressureMesh->Print(outTXT);
            TPZVTKGeoMesh::PrintCMeshVTK(pressureMesh, outVTK);
        }
#endif

        CreateSkeletonElements(pressureMesh);
        
        return pressureMesh;
    }
}

/// create the post processed multiphysics mesh (which is necessarily
/// hybridized)
void TPZHDivErrorEstimator::CreatePostProcessingMesh() {

#ifdef PZDEBUG
    {
        std::ofstream out("OriginalFlux.txt");
        fOriginal->MeshVector()[0]->Print(out);
        std::ofstream out2("OriginalPotential.txt");
        fOriginal->MeshVector()[1]->Print(out2);
        std::ofstream out3("OriginalMeshHybrid.txt");
        fPostProcMesh.Print(out3);
    }
#endif
    
    // initialize the post processing mesh
    fPostProcMesh.SetReference(fOriginal->Reference());
    int dim = fOriginal->Dimension();
    fOriginal->CopyMaterials(fPostProcMesh);
    
    {
        fPostProcMesh.DeleteMaterial(fHybridizer.fHDivWrapMatid);
        fPostProcMesh.DeleteMaterial(fHybridizer.fInterfaceMatid.first);
        fPostProcMesh.DeleteMaterial(fHybridizer.fInterfaceMatid.second);
        fPostProcMesh.DeleteMaterial(fHybridizer.fLagrangeInterface);
    }
    
    // switch the material from mixed to TPZMixedHdivErrorEstimate...
    SwitchMaterialObjects();
    
    TPZManVector<TPZCompMesh *> meshvec(4, 0);
    meshvec[0] = 0;
    meshvec[2] = fOriginal->MeshVector()[0];//flux
    meshvec[3] = fOriginal->MeshVector()[1];//potential
    // 
    meshvec[1] = CreatePressureMesh();//potential reconstructed

    if(fPostProcesswithHDiv)
    {
        //flux reconstructed just using Hdiv reconstruction
        meshvec[0] = CreateFluxMesh();
    }

    // TODO (Gustavo) I commented these lines since the materials are not being identified correctly and
    //      we don't need when we set fHybridizer externally.
    //if (!fOriginalIsHybridized) {
    //    fHybridizer.ComputePeriferalMaterialIds(meshvec);
    //    fHybridizer.ComputeNState(meshvec);
    //    fHybridizer.HybridizeInternalSides(meshvec);
    //} else {
    //    IdentifyPeripheralMaterialIds();
    //}
    
    if (fPostProcesswithHDiv) {
        IncreaseSideOrders(meshvec[0]);//malha do fluxo
    }
    
    if (dim == 3) {
        CreateEdgeSkeletonMesh(meshvec[1]);
    }
#ifdef PZDEBUG2
    {
        if(fPostProcesswithHDiv)
        {
            std::ofstream out("EnrichedFluxBorder.txt");
            meshvec[0]->Print(out);
        }
        std::ofstream out2("EnrichedPressure.txt");
        meshvec[1]->Print(out2);
    }
#endif
    
    
    TPZManVector<int> active(4, 0);
    active[1] = 1;
    
    if(fPostProcesswithHDiv)
    {
        // the flux mesh is active only if we postprocess with an H(div) approximation
        active[0] = 1;
    }
    
    fPostProcMesh.BuildMultiphysicsSpace(active, meshvec);
    {
        std::ofstream out("multiphysicsCreated.txt");
        fPostProcMesh.Print(out);
    }

    VerifyConnectConsistency(&fPostProcMesh); // TODO

    if(fPostProcesswithHDiv) {
        // construction of the multiphysics mesh
        //cria elementos de interface
        fHybridizer.CreateInterfaceElements(&fPostProcMesh);
        fHybridizer.GroupandCondenseElements(&fPostProcMesh);
        fPostProcMesh.CleanUpUnconnectedNodes();
    }
    else {
        PrepareElementsForH1Reconstruction();
    }

    ComputePressureWeights();

#ifdef PZDEBUG2
    {
        std::ofstream out("multiphysicsgrouped.txt");
        fPostProcMesh.Print(out);
        //            std::ofstream outvtk("multiphysics.vtk");
        //            TPZVTKGeoMesh::PrintCMeshVTK(cmesh_Hybrid,outvtk);
        std::ofstream outgvtk("postprocessgmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fPostProcMesh.Reference(), outgvtk);
    }
#endif
}

/// computing the element stifnesses will "automatically" compute the condensed form of the matrices
void TPZHDivErrorEstimator::ComputeElementStiffnesses() {
    std::cout << "Solving local Dirichlet problem " << std::endl;

    for (auto cel:fPostProcMesh.ElementVec()) {
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
#ifdef PZDEBUG
        if(subcmesh && condense)
        {
            DebugStop();
        }
#endif
    }
}

/// increase the side orders of the post processing flux mesh
void TPZHDivErrorEstimator::IncreaseSideOrders(TPZCompMesh *mesh) {
    int64_t nel = mesh->NElements();
    int dim = mesh->Dimension();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = mesh->Element(el);
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != dim) {
            continue;
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        int nc = cel->NConnects();
        int order = cel->Connect(nc - 1).Order();
        int nsides = gel->NSides();
        int ncorner = gel->NCornerNodes();
        intel->SetPreferredOrder(order);
        for (int side = ncorner; side < nsides - 1; side++) {
            if (intel->NSideConnects(side)) {
                TPZConnect &c = intel->SideConnect(0, side);
                if (c.Order() != order) {
                    intel->SetSideOrder(side, order);
                }
            }
        }
    }
    mesh->InitializeBlock();
}

void TPZHDivErrorEstimator::IncreasePressureSideOrders(TPZCompMesh *cmesh) {

    TPZGeoMesh *gmesh = cmesh->Reference();
    
    gmesh->ResetReference();
    cmesh->LoadReferences();
    
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
        TPZMaterial *mat = cel->Material();

        TPZGeoElSide gelside(gel, gel->NSides() - 1);
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 1, 0);

        int nneigh = celstack.NElements();
        if (nneigh == 1) {
            TPZCompElSide celside = gelside.LowerLevelCompElementList2(1);
            if (!celside) continue;//DebugStop();/// para nao incremenentar ordem na condicao de contorno
            celstack.Push(celside);
            nneigh++;
        } else if (nneigh != 2) {
            DebugStop();
        }
        
        int maxOrder = 0;

        for (int ineigh = 0; ineigh < nneigh; ineigh++) {
            TPZInterpolatedElement *intelS = dynamic_cast<TPZInterpolatedElement *>(celstack[ineigh].Element());
            int orderEl = intelS->GetPreferredOrder();
            maxOrder = (orderEl > maxOrder) ? orderEl : maxOrder;
        }

        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        int nsides = gel->NSides();
        int ncorner = gel->NCornerNodes();
        
        for (int side = ncorner; side < nsides; side++) {
            if (intel->NSideConnects(side)) {
                intel->SetSideOrder(side, maxOrder);
            }
        }
    }
    cmesh->InitializeBlock();
}

/// searches for a neighbour whose element has the proper dimension and materialid
static TPZGeoElSide HasNeighbour(const TPZGeoElSide &gelside, int matid) {
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

/// create dim-2 skeleton mesh based on the dim-1 faces
// will do nothing if the dimension of the mesh == 2
void TPZHDivErrorEstimator::CreateEdgeSkeletonMesh(TPZCompMesh *pressuremesh) {
    
    if (pressuremesh->MaterialVec().find(fPressureSkeletonMatId) != pressuremesh->MaterialVec().end()) {
        DebugStop();
    }
    TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial<STATE>(fPressureSkeletonMatId);
    pressuremesh->InsertMaterialObject(nullmat);
    int dim = fPostProcMesh.Dimension();
    int64_t nel = pressuremesh->NElements();
    std::map<int64_t, int> gelpressures;
    // create the geometrical elements
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressuremesh->Element(el);
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
            TPZGeoElSide hasneigh = HasNeighbour(gelside, fPressureSkeletonMatId);
            if (!hasneigh) {
                TPZGeoElBC gbc(gelside, fPressureSkeletonMatId);
                TPZGeoEl *createdelement = gbc.CreatedElement();
                hasneigh = TPZGeoElSide(createdelement, createdelement->NSides() - 1);
                gelpressures[createdelement->Index()] = polynomialorder;
            } else {
                int64_t gelindex = hasneigh.Element()->Index();
#ifdef PZDEBUG
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
    pressuremesh->Reference()->ResetReference();
    TPZGeoMesh *gmesh = pressuremesh->Reference();
    nel = gmesh->NElements();
    for (auto indexpair : gelpressures) {
        int64_t index = indexpair.first;
        int polynomialorder = indexpair.second;
        TPZGeoEl *gel = gmesh->Element(index);
        if (!gel) DebugStop();
        TPZCompEl *cel = nullptr;
        int64_t celindex = -1;
        pressuremesh->SetDefaultOrder(polynomialorder);
        cel = pressuremesh->ApproxSpace().CreateCompEl(gel, *pressuremesh, celindex);
#ifdef PZDEBUG
        {
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            if (!intel) DebugStop();
            int porder = intel->GetPreferredOrder();
            if (porder != polynomialorder) DebugStop();
        }
#endif
        gel->ResetReference();
    }
    AdjustNeighbourPolynomialOrders(pressuremesh);
    pressuremesh->ExpandSolution();
    RestrainSmallEdges(pressuremesh);
}

/// restrain the edge elements that have larger elements as neighbours
void TPZHDivErrorEstimator::RestrainSmallEdges(TPZCompMesh *pressuremesh) {
    //    TPZCompMesh *pressuremesh = fPostProcMesh.MeshVector()[1];
    TPZGeoMesh *gmesh = pressuremesh->Reference();
    gmesh->ResetReference();
    int dim = fPostProcMesh.Dimension();
    int64_t nel = pressuremesh->NElements();
    // load the face and edge elements
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressuremesh->Element(el);
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
        TPZCompEl *cel = pressuremesh->Element(el);
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
void TPZHDivErrorEstimator::AdjustNeighbourPolynomialOrders(TPZCompMesh *pressureHybrid) {
    //    TPZCompMesh *pressureHybrid = fPostProcMesh.MeshVector()[1];
    TPZGeoMesh *gmesh = pressureHybrid->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    // load the elements of lower dimension than dim
    int64_t nel = pressureHybrid->NElements();
    std::map<std::pair<int64_t, int>, int> polynomialorders;
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressureHybrid->Element(el);
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
            TPZCompEl *cel = pressureHybrid->Element(el);
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
        TPZCompEl *cel = pressureHybrid->Element(index);
        if (!cel) DebugStop();
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        intel->SetSideOrder(side, porder);
    }
}

/// return a pointer to the pressure mesh
TPZCompMesh *TPZHDivErrorEstimator::PressureMesh() {
    return fPostProcMesh.MeshVector()[1];
}

/// compute the average pressures of across edges of the H(div) mesh
void TPZHDivErrorEstimator::ComputeAveragePressures(int target_dim) {
    
    TPZCompMesh *pressure_mesh = PressureMesh();
    TPZGeoMesh *gmesh = pressure_mesh->Reference();
    int dim = gmesh->Dimension();
    gmesh->ResetReference();

    int64_t nel = pressure_mesh->NElements();
    // load the pressure elements of dim target_dim + 1 and of skeleton mat IDs
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() != target_dim + 1 && gel->MaterialId() != fPressureSkeletonMatId) continue;
        gel->SetReference(cel);
    }
    
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        TPZGeoEl *gel = cel->Reference();

        if (!cel || !cel->Reference() || cel->Reference()->Dimension() != target_dim) {
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

#ifdef PZDEBUG
        int nsides = gel->NSides();
        TPZGeoElSide side(gel, nsides - 1);
        TPZGeoElSideAncestors ancestors(side);
        TPZGeoElSide largerNeigh = ancestors.HasLarger(fPressureSkeletonMatId);
        if (largeSideExists && !largerNeigh) DebugStop();
#endif
        if (largeSideExists) continue;
        ComputeAverage(pressure_mesh, el);
    }
    
    // Loads solution into the connects of the smaller skeletons
    pressure_mesh->LoadSolution(pressure_mesh->Solution());

    // apply the restraints to the edge connects
    if (target_dim == dim - 2) {
        pressure_mesh->LoadSolution(pressure_mesh->Solution());
        TransferEdgeSolution();
    }
    
}

//compute de L2 projection of Dirichlet boundary condition for Hdi-H1 reconstruction
void TPZHDivErrorEstimator::ComputeBoundaryL2Projection(int target_dim){
    std::cout << "Computing boundary L2 projection\n";
    TPZCompMesh* pressuremesh = PressureMesh();
    {
        std::ofstream out("PressureBeforeL2Projection.txt");
        pressuremesh->Print(out);
    }
    if (target_dim == 2) {
        std::cout << "Not implemented for 2D interface" << std::endl;
        DebugStop();
    }
    
    TPZGeoMesh *gmesh = pressuremesh->Reference();
    gmesh->ResetReference();
    int64_t nel = pressuremesh->NElements();
    
    TPZAdmChunkVector<TPZCompEl *> &elementvec = pressuremesh->ElementVec();
    
    TPZElementMatrixT<STATE> ekbc, efbc;
    TPZFMatrix<STATE> &mesh_sol = pressuremesh->Solution();
    for (int iel = 0; iel < nel; iel++) {
        TPZCompEl *cel = elementvec[iel];
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        
        int matid = gel->MaterialId();
        TPZMaterial *mat = pressuremesh->FindMaterial(matid);
        
        TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
        if (!bc || (bc->Type() != 0)) continue;
        
        cel->CalcStiff(ekbc, efbc);
        ekbc.fMat.SolveDirect(efbc.fMat, ELU);
        int count = 0;
        int nc = cel->NConnects();
        for (int ic = 0; ic < nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            int64_t seqnum = c.SequenceNumber();
            int64_t pos = pressuremesh->Block().Position(seqnum);
            int ndof = c.NShape() * c.NState();
            for (int idf = 0; idf < ndof; idf++) {
                mesh_sol(pos + idf, 0) = efbc.fMat(count++);
            }
        }
    }
    
    {
        std::ofstream out("PressureAfterL2Projection.txt");
        pressuremesh->Print(out);
    }
    std::cout << "Finished computing boundary L2 projection\n";
}

// compute the average of an element iel in the pressure mesh looking at its neighbours
void TPZHDivErrorEstimator::ComputeAverage(TPZCompMesh *pressuremesh, int64_t iel) {
    TPZGeoMesh *gmesh = pressuremesh->Reference();
    int dim = gmesh->Dimension();
    TPZCompEl *cel = pressuremesh->Element(iel);
    if (!cel) DebugStop();
    TPZGeoEl *gel = cel->Reference();
    if (!gel) DebugStop();
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
    if (!intel) DebugStop();
    int nc = cel->NConnects();
    int order = cel->Connect(nc - 1).Order();

#ifdef PZDEBUG
    for (int ic = 0; ic < nc; ic++) {
        if (cel->Connect(ic).HasDependency()) DebugStop();
    }
#endif
    
    //std::cout << "Computing average for compel " << iel << ", gel: " << cel->Reference()->Index() << "\n";

    TPZFMatrix<STATE> &mesh_sol = pressuremesh->Solution();

    int target_dim = gel->Dimension();
    if (target_dim == dim - 1 && gel->MaterialId() != fPressureSkeletonMatId) {
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
        partition.HigherLevelNeighbours(smallerSkelSides, fPressureSkeletonMatId);
        if (smallerSkelSides.size() == 1) DebugStop();
        for (int iskel = 0; iskel < smallerSkelSides.size(); iskel++) {
            smallerSkelSides[iskel].EqualLevelCompElementList(volumeNeighSides, 1, 0);
        }
    }

#ifdef PZDEBUG
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
    
    L2Mat.SolveDirect(L2Rhs, ECholesky);
    // Stores solution in the computational mesh
    int count = 0;
    for (int ic = 0; ic < nc; ic++) {
        TPZConnect &c = cel->Connect(ic);
        int64_t seqnum = c.SequenceNumber();
        int64_t pos = pressuremesh->Block().Position(seqnum);
        int ndof = c.NShape() * c.NState();
        for (int idf = 0; idf < ndof; idf++) {
            mesh_sol(pos + idf, 0) = L2Rhs(count++);
        }
    }
}


/// transfer the solution of the edge functions to the face functions
void TPZHDivErrorEstimator::TransferEdgeSolution() {
    // copy the solution associated with one-d edge connect to the corresponding side connect of the face mesh
    TPZCompMesh *pressureHybrid = PressureMesh();
    TPZGeoMesh *gmesh = pressureHybrid->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    if (dim != 3) {
        std::cout << __PRETTY_FUNCTION__ << " should not be called for mesh dimension " << dim << std::endl;
        return;
    }
    int lagrangematid = fHybridizer.fLagrangeInterface;
    TPZMaterial *mat = pressureHybrid->FindMaterial(lagrangematid);
    if (!mat) DebugStop();
    int nstate = mat->NStateVariables();
    int64_t nel = pressureHybrid->NElements();
    // load the pressure elements of dimension 1 and 2
    pressureHybrid->Reference()->ResetReference();

    TPZBlock &block =  pressureHybrid->Block();
    TPZFMatrix<STATE> &sol = pressureHybrid->Solution();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressureHybrid->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() > 2) continue;
        gel->SetReference(cel);
    }
    
    // loop over the edge elements
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressureHybrid->Element(el);
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
            int nshape_edge = pressureHybrid->Block().Size(edge_seqnum);
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
                int nshape_neigh = pressureHybrid->Block().Size(neighblock);
                if (nshape_edge != nshape_neigh) DebugStop();
                for (int i = 0; i < nshape_neigh; i++) {
                    sol.at(block.at(neighblock, 0, i, 0)) = sol.at(block.at(edge_seqnum, 0, i, 0));
                }
            }
        }
    }
}

/// set the cornernode values equal to the averages
void TPZHDivErrorEstimator::ComputeNodalAverages() {
    TPZCompMesh *pressure_mesh = PressureMesh();
    TPZGeoMesh *gmesh = pressure_mesh->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
#ifdef PZDEBUG
    TPZMaterial *mat = pressure_mesh->FindMaterial(fPressureSkeletonMatId);
    if (!mat) DebugStop();
#endif
    int64_t nel = pressure_mesh->NElements();
    // load the pressure elements of dimension dim-1
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() != dim - 1) continue;
        gel->SetReference(cel);
    }

    // We create a stack of TPZCompElSide to store skeleton nodes that are adjacent to hanging nodes, but their elements
    // aren't part of the set of small skeletons. In these node connects, we don't calculate an average but, instead,
    // impose the solution of the restrained neighbours.
    TPZStack<TPZCompElSide> nodesToImposeSolution;
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        
        if (gel->Dimension() != dim - 1 || gel->MaterialId() != fPressureSkeletonMatId) {
            // MatId of element not important for nodal average computation
            continue;
        }
        
        // Compute nodal average for nodal sides of skeleton element
        int ncorners = gel->NCornerNodes();
        for (int side = 0; side < ncorners; side++) {
            TPZGeoElSide gelside(gel, side);
            TPZCompElSide celside(intel,side);

            if (IsAdjacentToHangingNode(celside)) {
                nodesToImposeSolution.Push(celside);
                continue;
            }

            int64_t conindex = intel->ConnectIndex(side);
            TPZConnect &c = pressure_mesh->ConnectVec()[conindex];
            if (!c.HasDependency()) {
                ComputeNodalAverage(celside);
            }
        }
    }

    pressure_mesh->LoadSolution(pressure_mesh->Solution());

    TPZBlock &block =  pressure_mesh->Block();
    TPZFMatrix<STATE> &sol = pressure_mesh->Solution();

    // Impose solution on nodes adjacent to hanging nodes
    for (int64_t i = 0; i < nodesToImposeSolution.size(); i++) {
        TPZCompElSide node_celside = nodesToImposeSolution[i];
        TPZGeoElSide node_gelside(node_celside.Reference());

        // celstack will contain all zero dimensional sides connected to the side
        TPZStack<TPZCompElSide> celstack;
        int onlyinterpolated = 1;
        int removeduplicates = 0;

        node_gelside.ConnectedCompElementList(celstack, onlyinterpolated, removeduplicates);

        for (int elc = 0; elc < celstack.size(); elc++) {
            TPZCompElSide neigh_celside = celstack[elc];
            if (neigh_celside.Reference().Dimension() != 0) continue;

            // Get solution of the neighbour
            TPZInterpolatedElement *neigh_intel = dynamic_cast<TPZInterpolatedElement *>(neigh_celside.Element());
            if (!neigh_intel) DebugStop();

            int64_t neigh_conindex = neigh_intel->ConnectIndex(neigh_celside.Side());
            TPZConnect &neigh_c = pressure_mesh->ConnectVec()[neigh_conindex];

            int64_t neigh_seqnum = neigh_c.SequenceNumber();
            int nstate = 1;
            if (neigh_c.NState() != nstate || neigh_c.NShape() != 1) DebugStop();
            TPZManVector<STATE, 3> neigh_sol(nstate, 0.);
            for (int istate = 0; istate < nstate; istate++) {
                neigh_sol[istate] = sol.at(block.at(neigh_seqnum, 0, istate, 0));
            }

            // Set solution to given connect
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(node_celside.Element());
            if (!intel) DebugStop();

            int side = node_gelside.Side();
            int64_t conindex = intel->ConnectIndex(side);
            TPZConnect &c = pressure_mesh->ConnectVec()[conindex];

            int64_t seqnum = c.SequenceNumber();
            if (c.NState() != nstate || c.NShape() != 1) DebugStop();
            for (int istate = 0; istate < nstate; istate++) {
                sol.at(block.at(seqnum, 0, istate, 0)) = neigh_sol[istate];
            }
            break;
        }
        pressure_mesh->LoadSolution(pressure_mesh->Solution());
    }
}

/// compute the nodal average of all elements that share a point
void TPZHDivErrorEstimator::ComputeNodalAverage(TPZCompElSide &node_celside)
{
    int skeletonMatId = fPressureSkeletonMatId;
    TPZCompMesh *pressure_mesh = PressureMesh();
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
    celstack.Push(node_celside);

    TPZBlock &block =  pressure_mesh->Block();
    TPZFMatrix<STATE> &solMatrix = pressure_mesh->Solution();
    
    // This map stores the connects, the weight associated with the element
    // and the solution of that connect. The weight of Dirichlet condition is
    // higher and will be used later to impose the value of the BC in the
    // connects when needed
    std::map<int64_t, std::pair<REAL, TPZVec<STATE>>> connects;
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
  //          std::cout << conindex << std::endl;
 //           DebugStop();
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
            solMatrix.at(block.at(seqnum, 0, istate, 0)) = averageSol[istate];
        }
    }
    
}

/// compute the effectivity indices of the pressure error and flux error and store in the element solution
void TPZHDivErrorEstimator::ComputeEffectivityIndices(TPZSubCompMesh *subcmesh) {
    TPZFMatrix<STATE> &elsol = subcmesh->ElementSolution();
    int64_t nrows = elsol.Rows();
    int64_t ncols = 5; // subcmesh->ElementSolution().Cols();

    if (subcmesh->ElementSolution().Cols() != 7) {
        // TODO I made some changes to be able to run the code again.
        //  Sometimes subcmesh->ElementSolution().Cols() equals 5  sometimes it's already 7.
        //  I'm not sure if this behaviour is expected.
        subcmesh->ElementSolution().Resize(nrows, 7);
    }

    int64_t nel = subcmesh->NElements();
    TPZManVector<REAL, 5> errors(5, 0.);
    for (int64_t el = 0; el < nel; el++) {
        for (int i = 0; i < ncols; i++) {
            errors[i] += elsol(el, i) * elsol(el, i);
        }
    }
    for (int i = 0; i < ncols; i++) {
        errors[i] = sqrt(errors[i]);
    }

    for (int64_t el = 0; el < nel; el++) {
        for (int i = 0; i < ncols; i++) {
            elsol(el, i) = errors[i];
        }
    }

    for (int64_t el = 0; el < nrows; el++) {
        TPZCompEl *cel = subcmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        REAL hk = gel->CharacteristicSize();

        for (int i = 0; i < 3; i += 2) {

            REAL tol = 1.e-10;
            REAL ErrorEstimate = errors[i + 1];
            REAL ErrorExact = errors[i];

            REAL oscillatoryterm = 0;
            if (i == 2) {
                oscillatoryterm = elsol(el, i + 2);
                oscillatoryterm *= (hk / M_PI);
            }

            if (abs(ErrorEstimate) < tol) {
                elsol(el, ncols + i / 2) = 1.;

            } else {
                REAL EfIndex = (ErrorEstimate + oscillatoryterm) / ErrorExact;
                elsol(el, ncols + i / 2) = EfIndex;
            }
        }
    }
}

/// compute the effectivity indices of the pressure error and flux error and store in the element solution
void TPZHDivErrorEstimator::ComputeEffectivityIndices() {
    /**The  ElementSolution() is a matrix with 4 cols,
     col 0: pressure exact error
     col 1: pressure estimate error
     col 2: flux exact error
     col 3: flux estimate error
     Is increased 2 cols on ElementSolution() to store the efectivity index for pressure and flux
     **/

    TPZCompMesh *cmesh = &fPostProcMesh;
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    TPZFMatrix<STATE> &elsol = cmesh->ElementSolution();
    int64_t nrows = elsol.Rows();
    int64_t ncols = elsol.Cols();

    // std::ostream &out;
    //    cmesh->ElementSolution().Print("ElSolution",std::cout);

    TPZFMatrix<REAL> dataIeff(nrows, 1);
    dataIeff.Zero();
    TPZFMatrix<REAL> InnerEstimated(nrows, 1);
    InnerEstimated.Zero();
    TPZFMatrix<REAL> InnerExact(nrows, 1);
    InnerExact.Zero();
    TPZFMatrix<REAL> BoundEstimated(nrows, 1);
    BoundEstimated.Zero();
    TPZFMatrix<REAL> BoundExact(nrows, 1);
    BoundExact.Zero();

    int dim = cmesh->Dimension();
    elsol.Resize(nrows, ncols + 2);

    std::set<int> bcMatIDs = GetBCMatIDs(&fPostProcMesh);

    for (int64_t el = 0; el < nrows; el++) {

        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) continue;
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (subcmesh) {
            ComputeEffectivityIndices(subcmesh);
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) continue;
        if (gel->Dimension() != dim) continue;
        int nsides = gel->NSides();
        for (int is = 0; is < nsides; is++) {
            if (gel->SideDimension(is) != dim - 1) continue;
            TPZStack<TPZCompElSide> equal;
            TPZGeoElSide gelside(gel, is);
            gelside.EqualLevelCompElementList(equal, 0, 0);
            //            if(equal.size() != 1){
            //                std::cout<<"Number of neighbour "<<equal.size()<<"\n";
            //                DebugStop();
            //            }
            TPZGeoElSide neighbour;
            TPZCompElSide selected;
            for (int i = 0; i < equal.size(); i++) {
                TPZGeoEl *gequal = equal[i].Element()->Reference();
                int eldim = gequal->Dimension();
                if (eldim != dim - 1) continue;
                int elmatid = gequal->MaterialId();
                if (bcMatIDs.find(elmatid) != bcMatIDs.end()) {
                    neighbour = equal[i].Reference();
                    selected = equal[i];
                    break;
                }
            }

            if (!neighbour) continue;
            if (neighbour.Element()->Dimension() != dim - 1) DebugStop();
            int64_t neighindex = selected.Element()->Index();
            for (int i = 0; i < 3; i += 2) {

                // std::cout << "linha = " << el << " col = " << 4 + i / 2 << " neinEl " << neighindex << std::endl;

                if (neighindex > nrows) {
                    std::cout << " neighindex= " << neighindex << " nrows " << nrows << "\n";
                    DebugStop();
                }

                REAL NeighbourErrorEstimate = elsol(neighindex, i + 1);
                REAL NeighbourErrorExact = elsol(neighindex, i);
                REAL ErrorEstimate = elsol(el, i + 1);
                REAL ErrorExact = elsol(el, i);

                InnerEstimated(el, 0) = ErrorEstimate;
                InnerExact(el, 0) = ErrorExact;
                BoundExact(el, 0) = NeighbourErrorExact;
                BoundEstimated(el, 0) = NeighbourErrorEstimate;

#ifdef LOG4CXX
                if (logger->isDebugEnabled()) {
                    std::stringstream sout;
                    sout << "El " << el << " dim " << dim << " ErrorEstimate " << ErrorEstimate << " ErrorExact "
                         << ErrorExact << "\n";
                    sout << "neighbour " << neighindex << " dim " << neighbour.Element()->Dimension()
                         << " NeighbourErrorEstimate " << NeighbourErrorEstimate << " NeighbourErrorExact "
                         << NeighbourErrorExact << "\n";

                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif

                REAL sumErrorExact = sqrt(NeighbourErrorExact * NeighbourErrorExact + ErrorExact * ErrorExact);
                REAL sumErrorEstimate =
                    sqrt(NeighbourErrorEstimate * NeighbourErrorEstimate + ErrorEstimate * ErrorEstimate);
                elsol(neighindex, i + 1) = 0.;
                elsol(neighindex, i) = 0.;
                elsol(el, i) = sumErrorExact;
                elsol(el, i + 1) = sumErrorEstimate;
            }
        }
    }
    for (int64_t el = 0; el < nrows; el++) {

        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) continue;
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (subcmesh) {
            ComputeEffectivityIndices(subcmesh);
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) continue;
        REAL hk = gel->CharacteristicSize();

        for (int i = 0; i < 3; i += 2) {

            //  std::cout<<"linha = "<<el<< "col = "<<4 + i / 2<<std::endl;

            REAL tol = 1.e-10;
            REAL ErrorEstimate = elsol(el, i + 1);
            REAL ErrorExact = elsol(el, i);

#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "El " << el << " dim " << gel->Dimension() << " ErrorEstimate " << ErrorEstimate
                     << " ErrorExact " << ErrorExact << "\n";
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif

            TPZGeoEl *gel = cel->Reference();

            REAL hk = gel->CharacteristicSize();

            REAL oscilatorytherm = 0;
            if (i == 2) {
                oscilatorytherm = elsol(el, i + 2);
                oscilatorytherm *= (hk / M_PI);
            }

            if (abs(ErrorEstimate) < tol) {
                elsol(el, ncols + i / 2) = 1.;
                dataIeff(el, 0) = 1.;
            } else {
                REAL EfIndex = (ErrorEstimate + oscilatorytherm) / ErrorExact;
                dataIeff(el, 0) = EfIndex;
                elsol(el, ncols + i / 2) = EfIndex;
            }
        }
    }

    {
        std::ofstream out("IeffPerElement.nb");
        dataIeff.Print("Ieff = ", out, EMathematicaInput);
        std::ofstream out1("InnerEstimated.nb");
        InnerEstimated.Print("InnerEstimated = ", out1, EMathematicaInput);
        std::ofstream out2("InnerExact.nb");
        InnerExact.Print("InnerExact = ", out2, EMathematicaInput);
        std::ofstream out3("BoundExact.nb");
        BoundExact.Print("BoundExact = ", out3, EMathematicaInput);
        std::ofstream out4("BoundEstimated.nb");
        BoundEstimated.Print("BoundEstimated = ", out4, EMathematicaInput);
    }
}

/// returns true if the material associated with the element is a boundary condition
/// and if the boundary condition is dirichlet type
bool TPZHDivErrorEstimator::IsDirichletCondition(TPZGeoElSide gelside) {
    TPZGeoEl *gel = gelside.Element();
    int matid = gel->MaterialId();
    TPZMaterial *mat = fPostProcMesh.FindMaterial(matid);
    TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
    if (!bc) return false;
    int typ = bc->Type();
    //  std::cout<<"type "<< typ<<"\n";
    if ((typ == 0)||(typ == 4)) return true;//for copy the robin boundary too
    return false;
}

void TPZHDivErrorEstimator::PotentialReconstruction() {
    
    if (fPostProcMesh.MeshVector().size()) {
        DebugStop();
    }

#ifdef PZDEBUG
    // Create directories to store debugging files
    std::string command;
    command = "mkdir -p ReconstructionSteps";
    system(command.c_str());
    command = "mkdir -p DebuggingTransfer";
    system(command.c_str());
#endif

    CreatePostProcessingMesh();

#ifdef PZDEBUG
    {
        REAL presnorm = Norm(fPostProcMesh.MeshVector()[1]->Solution());
        if(std::isnan(presnorm)) DebugStop();
    }
#endif

    // L2 projection for Dirichlet and Robin boundary condition for H1 reconstruction
    if (!fPostProcesswithHDiv) {
        int target_dim = 1;
        // TODO ver se fica igual para dimensao maior
        ComputeBoundaryL2Projection(target_dim);
        //BoundaryPressureProjection(pressuremesh, target_dim);
    }

#ifdef PZDEBUG
    {
        PlotPressureSkeleton("ReconstructionSteps/SkelBoundaryProjection");
    }
#endif

    // Calculates average pressure on interface edges and vertices
    int dim = fPostProcMesh.Dimension();
    ComputeAveragePressures(dim - 1);
    // in three dimensions make the one-d polynomials compatible
    if (dim == 3) {
        ComputeAveragePressures(1);
    }

#ifdef PZDEBUG
    {
        PlotPressureSkeleton("ReconstructionSteps/SkelInterfaceAverage");
    }
#endif

    ComputeNodalAverages();

#ifdef PZDEBUG
    {
        PlotPressureSkeleton("ReconstructionSteps/SkelNodalAverage");
    }
#endif

    //PlotState("ReconstructionSteps/VolumePressureBeforeCopyFromSkel", 2, fPostProcMesh.MeshVector()[1]);
    //PlotState("ReconstructionSteps/VolumeMFPressureBeforeCopyFromSkel", 2, &fPostProcMesh, false);
    CopySolutionFromSkeleton();
    //PlotState("ReconstructionSteps/VolumePressureAfterCopyFromSkel", 2, fPostProcMesh.MeshVector()[1]);
    //PlotState("ReconstructionSteps/VolumeMFPressureAfterCopyFromSkel", 2, &fPostProcMesh, false);

    // transfer the continuous pressures to the multiphysics space
    TPZManVector<TPZCompMesh *, 2> meshvec(2);
    // fPostProcMesh[0] is the H(div) mesh when post processing with H(div)
    // fPostProcMesh[1] is the L2 mesh
    meshvec[0] = fPostProcMesh.MeshVector()[0];
    meshvec[1] = fPostProcMesh.MeshVector()[1];

    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, &fPostProcMesh);

    std::ofstream out("ReconstructionSteps/MFMeshBeforeManualTransfer.txt");
    fPostProcMesh.Print(out);

#ifdef PZDEBUG
    {
        PlotState("ReconstructionSteps/VolumeMFPressureAfterTransferFromMeshes", 2, &fPostProcMesh, false);
        PlotState("ReconstructionSteps/VolumePressureAfterTransferFromMeshes", 2, fPostProcMesh.MeshVector()[1]);
        std::ofstream out("DebuggingTransfer/PressureAfterTransferFromMeshes.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(fPostProcMesh.MeshVector()[1], out);
        std::ofstream outMultiphysics ("DebuggingTransfer/MultiphysicsAfterTransferFromMeshes.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(&fPostProcMesh, outMultiphysics);
    }
#endif

    std::ofstream outafter("ReconstructionSteps/MFMeshAfterManualTransfer.txt");
    fPostProcMesh.Print(outafter);
    ComputeElementStiffnesses();

    fPostProcMesh.LoadSolution(fPostProcMesh.Solution());
    //PlotState("ReconstructionSteps/VolumeMFPressureAfterLoadSolution", 2, &fPostProcMesh, false);
    //PlotState("ReconstructionSteps/VolumePressureAfterLoadSolution", 2, fPostProcMesh.MeshVector()[1]);


    {
        std::ofstream out("DebuggingTransfer/PressureBeforeTransferFromMult.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(fPostProcMesh.MeshVector()[1], out);
        std::ofstream outMultiphysics("DebuggingTransfer/MultiphysicsBeforeTransferFromMult.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(&fPostProcMesh, outMultiphysics);
    }
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, &fPostProcMesh);
    {
        //PlotState("ReconstructionSteps/VolumeMFPressureAfterTransferFromMult", 2, &fPostProcMesh, false);
        //PlotState("ReconstructionSteps/VolumePressureAfterTransferFromMult", 2, fPostProcMesh.MeshVector()[1]);
        std::ofstream out("DebuggingTransfer/PressureAfterTransferFromMult.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(fPostProcMesh.MeshVector()[1], out);
        std::ofstream outMultiphysics("DebuggingTransfer/MultiphysicsAfterTransferFromMult.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(&fPostProcMesh, outMultiphysics);
    }


#ifdef PZDEBUG
    VerifySolutionConsistency(PressureMesh());
#endif

    PlotPressureSkeleton("ReconstructionSteps/FinalSkeletonPressure");

    if (fPostProcesswithHDiv) {
        PlotInterfaceFluxes("ReconstructedInterfaceFluxes", true);
    }
}

void TPZHDivErrorEstimator::PlotPressureSkeleton(const std::string &filename, bool reconstructed) {
    
    TPZCompMesh *pressure = nullptr;

    TPZStack<std::string> scalnames, vecnames;

    if (!reconstructed) {
        pressure = fOriginal->MeshVector()[1];
        scalnames.Push("State");
    } else {
        pressure = PressureMesh();
        scalnames.Push("State");
    }

    TPZLinearAnalysis an(pressure, false);

    {
        int dim = pressure->Reference()->Dimension() - 1;
        std::string plotname;
        {
            std::stringstream out;
            out << filename << ".vtk";
            plotname = out.str();
        }
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(5, dim);
    }
}

void TPZHDivErrorEstimator::PlotInterfaceFluxes(const std::string &filename, bool reconstructed) {
    TPZCompMesh *flux_mesh = nullptr;
    if (reconstructed) {
        if (!fPostProcesswithHDiv) DebugStop();
        flux_mesh = fPostProcMesh.MeshVector()[0];
    } else {
        flux_mesh = fPostProcMesh.MeshVector()[2];
    }

    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("State");

    TPZLinearAnalysis an(flux_mesh, false);

    {
        int dim = flux_mesh->Reference()->Dimension() - 1;
        std::string plotname;
        {
            std::stringstream out;
            out << filename << ".vtk";
            plotname = out.str();
        }
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(1, dim);
    }
}

/// switch material object from mixed poisson to TPZMixedHdivErrorEstimate
void TPZHDivErrorEstimator::SwitchMaterialObjects() {

    for (auto mat : fPostProcMesh.MaterialVec()) {
        TPZMixedDarcyFlow *mixpoisson = dynamic_cast<TPZMixedDarcyFlow *>(mat.second);
        if (mixpoisson) {
            TPZMixedDarcyFlow *newmat;
            if (fPostProcesswithHDiv) {
                newmat = new TPZMixedHDivErrorEstimate<TPZMixedDarcyFlow>(*mixpoisson);
            } else {
                newmat = new TPZHDivErrorEstimateMaterial(*mixpoisson);
            }

            if (mixpoisson->HasForcingFunction()) {
                newmat->SetForcingFunction(mixpoisson->ForcingFunction());
            }
            if (mixpoisson->HasExactSol()) {
                newmat->SetExactSol(mixpoisson->ExactSol(), mixpoisson->PolynomialOrderExact());
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

void TPZHDivErrorEstimator::VerifySolutionConsistency(TPZCompMesh *cmesh) {
    {
        std::ofstream outvtk("MeshToVerifyConsistency.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(), outvtk);
    }
    
    int64_t nel = cmesh->NElements();
    int dim = cmesh->Reference()->Dimension();
    
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
            if (gelside.Dimension() == dim) continue;
            
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
                        if (!IsZero(sol1[0] - sol0[0])) {

                            std::stringstream sout;
                            sout << "\nSide solution =  " << sol0[0] << "\n";
                            sout << "Neigh solution = " << sol1[0] << "\n";
                            sout << "Diff = " << sol1[0] - sol0[0] << "\n";
                            sout << "Side coord:  [" << x0[0] << ", " << x0[1] << ", " << x0[2] << "]\n";
                            sout << "Neigh coord: [" << x1[0] << ", " << x1[1] << ", " << x1[2] << "]\n";
                            std::cout << sout.str(); // TODO remove
                            LOGPZ_DEBUG(logger, sout.str())
                        }
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

void TPZHDivErrorEstimator::PrepareElementsForH1Reconstruction() {
    
    // This vector stores the connects from elements which have a neighbour of
    // an internal boundary material. We don't want to condense these connects,
    // so we are later incrementing the number of elements connected to them.
    // Then we compute the stiffness matrix and load the solution of the
    // internal degrees of freedom.
    TPZManVector<int64_t> connectsToIncrement(fPostProcMesh.NConnects(), -1);
    fPostProcMesh.ComputeNodElCon();
    
    TPZCompMesh *pressureMesh = fPostProcMesh.MeshVector()[1];
    pressureMesh->LoadReferences();
    
    for (int64_t el = 0; el < pressureMesh->NElements(); el++) {
        TPZCompEl *cel = pressureMesh->Element(el);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (!gel) continue;
        
        if (gel->MaterialId() != fPressureSkeletonMatId) continue;
        
        TPZGeoElSide skelSide(gel, gel->NSides() - 1);
        TPZStack<TPZCompElSide> compNeighSides;
        skelSide.EqualLevelCompElementList(compNeighSides, 1, 0);
        if (compNeighSides.size() == 1) {
            TPZCompElSide large = skelSide.LowerLevelCompElementList2(true);
            if (large) {
                std::cout << "Gel: " << gel->Index() << " Side: " << skelSide.Side() << '\n';
                //compNeighSides.Push(large);
            }
        }
        
        for (int i = 0; i < compNeighSides.size(); i++) {
            TPZCompEl *neighCel = compNeighSides[i].Element();
            TPZInterpolatedElement *neighIntEl = dynamic_cast<TPZInterpolatedElement *>(neighCel);
            if (!neighIntEl) DebugStop();
            
            int sideNum = compNeighSides[i].Side();
            int nCon = neighIntEl->NSideConnects(sideNum);
            
            for (int iCon = 0; iCon < nCon; iCon++) {
                TPZConnect &con = neighIntEl->SideConnect(iCon, sideNum);
                int64_t seqNum = con.SequenceNumber();
                int64_t conindex = neighIntEl->SideConnectIndex(iCon, sideNum);
                connectsToIncrement[conindex] = 1;
            }
        }
    }
    
    fPostProcMesh.LoadReferences();
    int64_t nel = fPostProcMesh.NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fPostProcMesh.Element(el);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (!gel) DebugStop();
        
        if (gel->Dimension() != fPostProcMesh.Dimension()) continue;
        TPZMultiphysicsElement *mcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mcel) DebugStop();
        
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
            int64_t index;
            TPZElementGroup *elGroup = new TPZElementGroup(fPostProcMesh, index);
            for (const auto &it : elementsToGroup) {
                elGroup->AddElement(fPostProcMesh.Element(it));
            }
        }
    }
    
    // Increments NElConnected of connects that should not be condensed
    fPostProcMesh.ComputeNodElCon();
    for (int64_t i = 0; i < fPostProcMesh.NConnects(); i++) {
        if (connectsToIncrement[i] == 1) {
            fPostProcMesh.ConnectVec()[i].IncrementElConnected();
        }
    }
    
#ifdef PZDEBUG2
    {
        std::ofstream txtPostProcMesh("PostProcMeshAfterIncrementingConnects.txt");
        fPostProcMesh.Print(txtPostProcMesh);
        std::ofstream txtPressureMesh("PressureMeshAfterIncrementingConnects.txt");
        pressureMesh->Print(txtPressureMesh);
    }
#endif
    
    for (int64_t el = 0; el < fPostProcMesh.NElements(); el++) {
        TPZCompEl *cel = fPostProcMesh.Element(el);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            TPZElementGroup *group = dynamic_cast<TPZElementGroup *>(cel);
            if (!group) DebugStop();
        }
        if (gel && gel->Dimension() != fPostProcMesh.Dimension()) continue;
        TPZCondensedCompEl *condense = new TPZCondensedCompEl(cel, false);
    }
    
    // @TODO what is the meaning of this? phil
    for (auto matit : fPostProcMesh.MaterialVec()) {
        TPZMaterial *mat = matit.second;
        TPZHDivErrorEstimateMaterial *errormat = dynamic_cast<TPZHDivErrorEstimateMaterial *>(mat);
        if (errormat) {
            errormat->fNeumannLocalProblem = false;
        }
    }
    
    fPostProcMesh.CleanUpUnconnectedNodes();
    
#ifdef PZDEBUG
    std::ofstream outTXT("PostProcessMeshAfterPreparingElements.txt");
    fPostProcMesh.Print(outTXT);
#endif
}

void TPZHDivErrorEstimator::CopySolutionFromSkeleton() {
    
    TPZCompMesh *pressuremesh = PressureMesh();
    pressuremesh->Reference()->ResetReference();
    pressuremesh->LoadReferences();
    int dim = pressuremesh->Dimension();
    int64_t nel = pressuremesh->NElements();

    TPZBlock &block =  pressuremesh->Block();
    TPZFMatrix<STATE> &sol = pressuremesh->Solution();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressuremesh->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() != dim) continue;
        
        int nsides = gel->NSides();
        for (int is = 0; is < nsides; is++) {
            TPZGeoElSide gelside(gel, is);
            TPZConnect &c = intel->Connect(is);
            int64_t c_seqnum = c.SequenceNumber();
            int c_blocksize = c.NShape() * c.NState();
            
            TPZStack<TPZCompElSide> celstack;
            gelside.EqualLevelCompElementList(celstack, 1, 0);
            int nst = celstack.NElements();
            for (int ist = 0; ist < nst; ist++) {
                TPZCompElSide cneigh = celstack[ist];
                TPZGeoElSide gneigh = cneigh.Reference();
                if (gneigh.Element()->MaterialId() == this->fPressureSkeletonMatId || IsDirichletCondition(gneigh)) {
                    TPZInterpolatedElement *intelneigh = dynamic_cast<TPZInterpolatedElement *>(cneigh.Element());
                    if (!intelneigh) DebugStop();
                    TPZConnect &con_neigh = intelneigh->Connect(cneigh.Side());
                    int64_t con_seqnum = con_neigh.SequenceNumber();
                    int con_size = con_neigh.NState() * con_neigh.NShape();
                    if (con_size != c_blocksize) DebugStop();
                    for (int ibl = 0; ibl < con_size; ibl++) {

                        sol.at(block.at(c_seqnum, 0, ibl, 0)) = sol.at(block.at(con_seqnum, 0, ibl, 0));
                    }
                    break;
                }
#ifdef LOG4CXX
                // all elements must have at least one neighbour of type skeleton--> esta premissa nao vale para reconstrucao Hdiv-H1
                if (ist == nst - 1) {
                    std::stringstream sout;
                    sout << "Connect " << is << " from element el " << el
                    << " was not updated \n";
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
            }
        }
    }
}


/// compute the pressure weights and material weights
// fills in the data structure fPressureweights and fMatid_weights
void TPZHDivErrorEstimator::ComputePressureWeights() {
    std::cout << "Computing pressure weights\n";
    TPZCompMesh *pressuremesh = fPostProcMesh.MeshVector()[1];
    int dim = pressuremesh->Dimension();
    int64_t nel = pressuremesh->NElements();
    fPressureweights.Resize(nel, 0);
    fMatid_weights.clear();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressuremesh->Element(el);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        TPZMaterial *mat = this->fOriginal->FindMaterial(matid);
        if (matid == fPressureSkeletonMatId || matid == fHybridizer.fLagrangeInterface) {
            fPressureweights[el] = 0.;
            fMatid_weights[matid] = 0.;
            continue;
        }
        if (!mat) DebugStop();

        TPZBndCond *bcmat = dynamic_cast<TPZBndCond *>(mat);
        if (bcmat) {
            switch(bcmat->Type()) {
                case 0: // Dirichlet BC
                    fPressureweights[el] = 1.e12;
                    fMatid_weights[matid] = 1.e12;
                    break;
                case 1: // Neumann BC
                    fPressureweights[el] = 0;
                    fMatid_weights[matid] = 0;
                    break;
                case 4: // Robin BC, weight = Km
                    fPressureweights[el] = bcmat->Val1()(0, 0);
                    fMatid_weights[matid] = bcmat->Val1()(0, 0);
                    break;
                default:
                    DebugStop();
            }
        } else {
            auto *mixpoisson = dynamic_cast<TPZMixedDarcyFlow*>(mat);
            if (!mixpoisson) DebugStop();

            REAL perm;
            TPZVec<REAL> xi(gel->Dimension(), 0.);
            gel->CenterPoint(gel->NSides() - 1, xi);
            TPZVec<REAL> x(3, 0.);
            gel->X(xi, x);

            TPZFNMatrix<1, STATE> K(1, 1, 0.), InvK(1, 1, 0.);
            mixpoisson->GetPermeabilities(x, K, InvK);
            perm = K(0, 0);
            if (IsZero(perm)) DebugStop();
            this->fPressureweights[el] = perm;
            fMatid_weights[matid] = perm;
        }
    }
    std::cout << "Finished computing pressure weights\n";
}

void TPZHDivErrorEstimator::PlotState(const std::string& filename, int targetDim, TPZCompMesh* cmesh, bool atomic) {
    
    std::ofstream outTXT("PressuretoStateGraph.txt");
    cmesh->Print(outTXT);
    
    {
        TPZAnalysis an(cmesh, false);
        TPZStack<std::string> scalnames, vecnames;
        if (atomic) {
            scalnames.Push("State");
        } else {
            scalnames.Push("PressureReconstructed");
        }

        std::string plotname;
        {
            std::stringstream out;
            out << filename << ".vtk";
            plotname = out.str();
        }
        an.DefineGraphMesh(targetDim, scalnames, vecnames, plotname);
        an.PostProcess(4, targetDim);
    }
}

// TODO separate geometric element creation into new method and improve documentation
void TPZHDivErrorEstimator::CreateSkeletonElements(TPZCompMesh *pressure_mesh) {

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

    // Assigns a material ID that has not been used yet // TODO move this to IdentifyPeripheralMatIds
    const int nel = gmesh->NElements();

    if (fPressureSkeletonMatId == 0) {
        fPressureSkeletonMatId = FindFreeMatId(gmesh);
    }

    std::set<int> bcMatIDs = GetBCMatIDs(&fPostProcMesh);
    bcMatIDs.insert(fPressureSkeletonMatId);

    int dim = gmesh->Dimension();

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

            // Create Geometric element if there is no boundary neighbour and no skeleton neighbour were created
            TPZGeoElSide neighSide = gelside.HasNeighbour(bcMatIDs);
            if(!neighSide.Exists())
            {
                TPZGeoElBC gbc(gelside, fPressureSkeletonMatId);
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

    // Create skeleton elements in pressure mesh
    TPZNullMaterial *skeletonMat = new TPZNullMaterial(fPressureSkeletonMatId);
    skeletonMat->SetDimension(dim - 1);
    pressure_mesh->InsertMaterialObject(skeletonMat);

    std::set<int> matIdSkeleton = { fPressureSkeletonMatId };
    gmesh->ResetReference();

    pressure_mesh->ApproxSpace().CreateDisconnectedElements(true);
    pressure_mesh->AutoBuild(matIdSkeleton);
    pressure_mesh->ExpandSolution();

    // increase the order of the dim-1 elements to the maximum of both neighbouring elements
    IncreasePressureSideOrders(pressure_mesh);//malha da pressao

    // Restrain the spaces of smaller to larger skeleton elements on hanging nodes
    RestrainSkeletonSides(pressure_mesh);

}

void TPZHDivErrorEstimator::RestrainSkeletonSides(TPZCompMesh *pressure_mesh) {

    TPZGeoMesh *gmesh = pressure_mesh->Reference();
    gmesh->ResetReference();
    pressure_mesh->LoadReferences();

#ifdef PZDEBUG
    {
        std::ofstream out("MeshBeforeRestrainSkeleton.txt");
        pressure_mesh->Print(out);
    }
#endif

    int64_t nel = gmesh->NElements();
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->Element(iel);
        TPZCompEl *cel = gel->Reference();
        if (!cel) continue;
        if (gel->MaterialId() != fPressureSkeletonMatId) continue;

        // If the element is a small skeleton, restrain its highest dimension side and then its subsides
        int nsides = gel->NSides();
        TPZGeoElSide small(gel, nsides - 1);
        TPZGeoElSideAncestors ancestors(small);
        TPZGeoElSide largerNeigh = ancestors.HasLarger(fPressureSkeletonMatId);
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
        std::cout << "Restriction @ [" << xcenter << "]:"
                  << "  Small El: " << small.Element()->Index() << ", Side: " << small.Side()
                  << "  Large El: " << largerNeigh.Element()->Index() << ", Side: " << largerNeigh.Side() << "\n";
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
            std::cout << "SubRestriction @ [" << xcenter << "]:"
                      << "  Small El: " << small.Element()->Index() << ", Side: " << subsmall.Side()
                      << "  Large El: " << largerNeigh.Element()->Index() << ", Side: " << subLargerNeigh.Side() << "\n";
            smallIntel->RestrainSide(subsmall.Side(), largeIntel, subLargerNeigh.Side());
        }
    }
    
    pressure_mesh->CleanUpUnconnectedNodes();

#ifdef PZDEBUG
    {
        std::ofstream out("MeshAfterRestrainSkeleton.txt");
        pressure_mesh->Print(out);
    }
#endif
}

// If all nodal (0-dimensional) connects linked to this node have dependency, it lies* on a hanging side and
// instead of calculating the average we should impose the solution of the other connects in it.
// * Although it seems to work for now, we are not sure if this logic will work for every case/ref. pattern,
// specially in 3D. So I'm leaving this disclaimer. (Gustavo Batistela, 14/10/2020)
// TODO I just realized this logic has flaws. Imagine the case where you two squares, a left and a right one.
//  Assume we refine only the right square, but instead of diving it uniformly, we divide (by 2) only the edge that
//  touches the left square, resulting in 3 triangles. We now have two edges which aren't restrained.
//  Thus, I think that we shouldn't check if 'all' connected nodes are restrained, but if 'n' connects are restained,
//  where 'n' is the number of elements by which the larger (father) edge side is divided.
bool TPZHDivErrorEstimator::IsAdjacentToHangingNode(const TPZCompElSide &celside) {

    TPZCompMesh *pressure_mesh = PressureMesh();
    int dim = pressure_mesh->Dimension();
    TPZGeoElSide gelside(celside.Reference());

    // celstack will contain all zero dimensional sides connected to the side
    TPZStack<TPZCompElSide> celstack;
    int onlyinterpolated = 1;
    int removeduplicates = 0;

    gelside.ConnectedCompElementList(celstack, onlyinterpolated, removeduplicates);

    bool allConnectedNodesAreRestrained = true;
    for (int elc = 0; elc < celstack.size(); elc++) {
        TPZCompElSide neigh_celside = celstack[elc];
        TPZGeoElSide neigh_gelside = neigh_celside.Reference();
        if (neigh_gelside.Element()->Dimension() != dim - 1) {
            continue;
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(neigh_celside.Element());
        if (!intel) DebugStop();

        int64_t conindex = intel->ConnectIndex(neigh_celside.Side());
        TPZConnect &c = pressure_mesh->ConnectVec()[conindex];
        if (neigh_gelside.Dimension() == 0 && !c.HasDependency()) allConnectedNodesAreRestrained = false;
    }

    return allConnectedNodesAreRestrained;
}

// Finds a material ID that has not been used yet
int TPZHDivErrorEstimator::FindFreeMatId(TPZGeoMesh *gmesh) {

    int maxMatId = std::numeric_limits<int>::min();
    const int nel = gmesh->NElements();

    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->Element(iel);
        if(gel) maxMatId = std::max(maxMatId, gel->MaterialId());
    }

    if (maxMatId == std::numeric_limits<int>::min()) maxMatId = 0;

    return maxMatId + 1;
}

std::set<int> TPZHDivErrorEstimator::GetBCMatIDs(const TPZCompMesh* cmesh) {
    std::set<int> bc_mat_ids;
    const auto mat_vec = cmesh->MaterialVec();
    for (const auto mat : mat_vec) {
        if (dynamic_cast<TPZBndCond*>(mat.second)) {
            bc_mat_ids.insert(mat.first);
        }
    }
    return bc_mat_ids;
}