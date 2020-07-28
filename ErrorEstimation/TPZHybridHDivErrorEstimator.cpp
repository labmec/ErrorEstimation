
//
//  TPZHybridHDivErrorEstimator.cpp
//  ErrorEstimateHDiv
//
//  Created by Philippe Devloo on 10/06/18.
//

#include "TPZHybridHDivErrorEstimator.h"
#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include "pzcompel.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"
#include "pzintel.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "mixedpoisson.h"
#include "TPZMixedHdivErrorEstimate.h"
#include "pzbuildmultiphysicsmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzanalysis.h"
#include "pzmat1dlin.h"
#include "TPZVecL2.h"
#include "TPZNullMaterial.h"
#include "pzelementgroup.h"
#include "TPZInterfaceEl.h"

#include "TPZParFrontStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzskylstrmatrix.h"

#include "TPZMultiphysicsCompMesh.h"
#include "pzmultiphysicscompel.h"

#include "TPZVTKGeoMesh.h"

#include "TPZHDivErrorEstimateMaterial.h"

#include "TPZCompMeshTools.h"
#include "TPZPressureProjection.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("HDivErrorEstimator"));
#endif

TPZHybridHDivErrorEstimator::~TPZHybridHDivErrorEstimator() {
    TPZVec<TPZCompMesh *> meshvec = fPostProcMesh.MeshVector();
    for (int i = 0; i < 2; i++) {
        delete meshvec[i];
    }
}

/// compute the element errors comparing the reconstructed solution based on average pressures
/// with the original solution
void TPZHybridHDivErrorEstimator::ComputeErrors(TPZVec<REAL>& errorVec, TPZVec<REAL> &elementerrors, bool store = true) {
    TPZAnalysis an(&fPostProcMesh, false);
    
    if (fExact) {
        an.SetExact(fExact->ExactSolution());
    }
    
#ifdef PZDEBUG2
    {
        std::ofstream out("PressureRecMeshComputeError.txt");
        fPostProcMesh.MeshVector()[1]->Print(out);
        std::ofstream out2("PressureMeshComputeError.txt");
        fPostProcMesh.MeshVector()[3]->Print(out2);
        //        std::ofstream out3("FluxRecMeshComputeError.txt");
        //        fPostProcMesh.MeshVector()[0]->Print(out3);
        std::ofstream out4("FluxMeshComputeError.txt");
        fPostProcMesh.MeshVector()[2]->Print(out);
        
    }
#endif
    
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
    
#ifdef PZDEBUG1
    {
        std::ofstream out("MeshToComputeError2.txt");
        fPostProcMesh.Print(out);
        
    }
#endif
    
    an.PostProcessError(errorVec, store);//calculo do erro com sol exata e aprox e armazena no elementsolution
    
    std::cout << "Computed errors " << errorVec << std::endl;
    
    TPZCompMeshTools::UnCondensedElements(&fPostProcMesh);
    TPZCompMeshTools::UnGroupElements(&fPostProcMesh);
    
    //Erro global
    ofstream myfile;
    myfile.open("ErrorsReconstruction.txt", ios::app);
    myfile << "\n\n Estimator errors for Problem " << fProblemConfig.problemname;
    myfile << "\n-------------------------------------------------- \n";
    myfile << "Ndiv = " << fProblemConfig.ndivisions << "AdaptativStep "<<fProblemConfig.adaptivityStep<<" Order k= " << fProblemConfig.porder << " Order n= "<< fProblemConfig.hdivmais<< " K_R = "<<fProblemConfig.Km<<"\n";
    myfile << "DOF Total = " << fPostProcMesh.NEquations() << "\n";
    myfile << "Global estimator = " << errorVec[3] << "\n";
    myfile << "Global exact error = " << errorVec[2] << "\n";
    myfile <<"|uex-ufem|= "<< errorVec[0] << "\n";
    myfile <<"|ufem-urec| = "<< errorVec[1] << "\n";
    myfile <<"Residual ErrorL2= "<< errorVec[4] << "\n";
    myfile <<"Global Index = "<< sqrt(errorVec[3]+ errorVec[4])/ errorVec[2];
    
    myfile.close();
    
    ComputeEffectivityIndices();
    
    //    GlobalEffectivityIndex();
    
    PostProcessing(an);
    
    elementerrors.resize(fPostProcMesh.Reference()->NElements());
    for (REAL & elementerror : elementerrors) {
        elementerror = 0;
    }
    for (int64_t i = 0; i < nelem; i++) {
        TPZCompEl *cel = fPostProcMesh.Element(i);
        if (!cel) continue;
        TPZGeoEl* gel = fPostProcMesh.Element(i)->Reference();
        if (!gel) continue;
        elementerrors[gel->Index()] = fPostProcMesh.ElementSolution()(i, 3);
    }
    
}

void TPZHybridHDivErrorEstimator::GlobalEffectivityIndex(){
    
    //    error[0] - error computed with exact pressure
    //    error[1] - error computed with reconstructed pressure
    //    error[2] - energy error computed with exact solution
    //    error[3] - energy error computed with reconstructed solution
    //    error[4] - residual data error
    
    int dim = fPostProcMesh.Dimension();
    int64_t nelem = fPostProcMesh.NElements();
    TPZManVector<REAL, 10> globalerrors(5, 0.);
    REAL Ieff_global = 0.;
    
    for (int64_t el = 0; el < nelem; el++) {
        
        TPZCompEl *cel = fPostProcMesh.ElementVec()[el];
        TPZGeoEl *gel = cel->Reference();
        REAL hk = gel->CharacteristicSize();
        if (cel->Reference()->Dimension() != dim) continue;
        TPZManVector<REAL, 10> elerror(10, 0.);
        elerror.Fill(0.);
        cel->EvaluateError(fExact->ExactSolution(), elerror, false);
        int nerr = elerror.size();
        
        for (int i = 0; i < nerr; i++) {
            
            globalerrors[i] += elerror[i] * elerror[i];
        }
        
        REAL coef1 = (hk / M_PI) * (hk / M_PI);
        
        globalerrors[nerr - 1] += coef1 * elerror[nerr - 1] * elerror[nerr - 1];
    }
    
    if (sqrt(globalerrors[2]) < 1.e-10) {
        Ieff_global = 1.;
    } else {
        std::cout << "residual " << globalerrors[4] << " estimated "
        << globalerrors[3] << " exact " << globalerrors[2] << "\n";
        
        Ieff_global =
        sqrt(globalerrors[4] + globalerrors[3]) / sqrt(globalerrors[2]);
    }
    
    std::cout << "Ieff_global " << Ieff_global << std::endl;
    
    ofstream myfile;
    myfile.open("ErrorEstimationResults.txt", ios::app);
    myfile << "\n\n Estimator errors for Problem "
    << fProblemConfig.problemname;
    myfile << "\n-------------------------------------------------- \n";
    myfile << "AdaptativStep " << fProblemConfig.adaptivityStep
    << " Order k= " << fProblemConfig.porder
    << " Order n= " << fProblemConfig.hdivmais << "\n";
    myfile << "DOF Total = " << fPostProcMesh.NEquations() << "\n";
    myfile << "Global exact error = " << sqrt(globalerrors[2]) << "\n";
    myfile << "Global estimator = " << sqrt(globalerrors[3]) << "\n";
    myfile << "Global residual error = " << sqrt(globalerrors[4]) << "\n";
    myfile << "Ieff_global = " << Ieff_global << "\n";
    myfile.close();
}

void TPZHybridHDivErrorEstimator::PostProcessing(TPZAnalysis &an) {
    
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
        
        std::stringstream out;
        out << fProblemConfig.dir_name << "/" << fProblemConfig.problemname
        << "_POrder" << fProblemConfig.porder << "_HdivMais_"
        << fProblemConfig.hdivmais;
        if (fProblemConfig.ndivisions != -1) {
            out << "_Ndiv_" << fProblemConfig.ndivisions;
        }
        if (fProblemConfig.adaptivityStep != -1) {
            out << "_AdaptivityStep_" << fProblemConfig.adaptivityStep;
        }
        out << ".vtk";
        
        an.DefineGraphMesh(dim, scalnames, vecnames, out.str());
        an.PostProcess(0, dim);
    }
    else {
        std::cout << __PRETTY_FUNCTION__ << "\nPost Processing variable not found!\n";
    }
}

// a method for generating the HDiv mesh
TPZCompMesh *TPZHybridHDivErrorEstimator::CreateFluxMesh()
{
    return fOriginal->MeshVector()[0]->Clone();
}
// a method for creating the pressure mesh
TPZCompMesh *TPZHybridHDivErrorEstimator::CreatePressureMesh()
{
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
        pressureMesh->LoadReferences();
        
        pressureMesh->ApproxSpace().SetAllCreateFunctionsContinuous();
        pressureMesh->ApproxSpace().CreateDisconnectedElements(true);
        
        std::set<int> bcMatIDs;
        for (auto it : mult->MaterialVec()) {
            TPZMaterial *mat = it.second;
            TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
            if (bc) {
                int pressureMatId = bc->Material()->Id();
                TPZMaterial *pressuremat =
                pressureMesh->FindMaterial(pressureMatId);
                TPZMaterial *newbc = pressuremat->CreateBC(
                                                           pressuremat, bc->Id(), bc->Type(), bc->Val1(), bc->Val2());
                if (bc->HasForcingFunction()) {
                    newbc->SetForcingFunction(bc->ForcingFunction());
                }
                pressureMesh->InsertMaterialObject(newbc);
                bcMatIDs.insert(bc->Id());
            }
        }
        
        // Create computational elements for BCs in pressureMesh mesh
        std::map<int64_t, TPZCompEl*> bcGeoElToNeighCompEl;
        for (int64_t el = 0; el < gmesh->NElements(); el++) {
            // This map stores the ID of the BC geometric element and a pointer
            // to the computational element of its volumetric neighbour
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
                
                bcGeoElToNeighCompEl.insert({gel->Index(), cel});
            }
        }
        gmesh->ResetReference();
        for (const auto &it : bcGeoElToNeighCompEl) {
            int64_t bcElID = it.first;
            TPZCompEl *neighCel = it.second;
            TPZGeoEl *bcGeoEl = gmesh->Element(bcElID);
            neighCel->LoadElementReference();
            int64_t id;
            TPZCompEl *bcCel = pressureMesh->CreateCompEl(bcGeoEl, id);
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
        
        return pressureMesh;
    }
}

/// create the post processed multiphysics mesh (which is necessarily
/// hybridized)
void TPZHybridHDivErrorEstimator::CreatePostProcessingMesh() {
    if(!fOriginalIsHybridized && !fPostProcesswithHDiv)
    {
        // we can not post process with H1 if the original mesh is not hybridized
        DebugStop();
    }
    
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
    
    TPZManVector<TPZCompMesh *> mesh_vectors(4, 0);
    mesh_vectors[0] = 0;
    mesh_vectors[2] = fOriginal->MeshVector()[0];//flux
    mesh_vectors[3] = fOriginal->MeshVector()[1];//potential
    mesh_vectors[1] = CreatePressureMesh();//potential reconstructed
    
    
    if(fPostProcesswithHDiv)
    {
        //flux reconstructed just using Hdiv reconstruction
        mesh_vectors[0] = CreateFluxMesh();
    }
    
    if (!fOriginalIsHybridized) {
        fHybridizer.ComputePeriferalMaterialIds(mesh_vectors);
        fHybridizer.ComputeNState(mesh_vectors);
        fHybridizer.HybridizeInternalSides(mesh_vectors);
        
        
    } else {
        IdentifyPeripheralMaterialIds();
        
    }
    
    // increase the order of the dim-1 elements to the maximum of both neighbouring elements
    IncreasePressureSideOrders(mesh_vectors[1]);//malha da pressao
    if(fPostProcesswithHDiv)
    {
        IncreaseSideOrders(mesh_vectors[0]);//malha do fluxo
    }
    
    if (dim == 3) {
        CreateEdgeSkeletonMesh(mesh_vectors[1]);
    }
#ifdef PZDEBUG2
    {
        if(fPostProcesswithHDiv)
        {
            std::ofstream out("EnrichedFluxBorder.txt");
            mesh_vectors[0]->Print(out);
        }
        std::ofstream out2("EnrichedPressure.txt");
        mesh_vectors[1]->Print(out2);
    }
#endif
    
    
    TPZManVector<int> active(4, 0);
    active[1] = 1;
    
    if(fPostProcesswithHDiv)
    {
        // the flux mesh is active only if we postprocess with an H(div) approximation
        active[0] = 1;
    }
    
    fPostProcMesh.BuildMultiphysicsSpace(active, mesh_vectors);
    {
        std::ofstream out("multiphysicsCreated.txt");
        fPostProcMesh.Print(out);
    }
    
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
void TPZHybridHDivErrorEstimator::ComputeElementStiffnesses() {
    std::cout << "Solving local Dirichlet problem " << std::endl;
#ifdef PZDEBUG2
    
    {
        std::ofstream out("MeshToComputeStiff.txt");
        fPostProcMesh.Print(out);
    }
#endif
    
    
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
    
    //    {
    //        std::ofstream out("MeshPostComputeStiff.txt");
    //        fPostProcMesh.Print(out);
    //    }
    
    
}


/// increase the side orders of the post processing flux mesh
void TPZHybridHDivErrorEstimator::IncreaseSideOrders(TPZCompMesh *mesh) {
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
        //        intel->Print();
    }
    mesh->InitializeBlock();
    
}

void TPZHybridHDivErrorEstimator::IncreasePressureSideOrders(TPZCompMesh *cmesh) {
    
    
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
        TPZMaterial *mat=cel->Material();
        
        //   std::cout<<"material "<<mat->Id()<<std::endl;
        TPZGeoElSide gelside(gel, gel->NSides() - 1);
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 1, 0);
        
        
        int nneigh = celstack.NElements();
        if (nneigh == 1) {
            TPZCompElSide celside = gelside.LowerLevelCompElementList2(1);
            if (!celside) continue;//DebugStop();/// para nao incremenentar ordem na condicao de contorno
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
void TPZHybridHDivErrorEstimator::CreateEdgeSkeletonMesh(TPZCompMesh *pressuremesh) {
    
    if (pressuremesh->MaterialVec().find(fPressureSkeletonMatId) != pressuremesh->MaterialVec().end()) {
        DebugStop();
    }
    TPZNullMaterial *nullmat = new TPZNullMaterial(fPressureSkeletonMatId);
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
                    polorder = max(polorder, polynomialorder);
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
        TPZCompEl *cel = 0;
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
void TPZHybridHDivErrorEstimator::RestrainSmallEdges(TPZCompMesh *pressuremesh) {
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
void TPZHybridHDivErrorEstimator::AdjustNeighbourPolynomialOrders(TPZCompMesh *pressureHybrid) {
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
TPZCompMesh *TPZHybridHDivErrorEstimator::PressureMesh() {
    return fPostProcMesh.MeshVector()[1];
}

/// compute the average pressures of the hybridized form of the H(div) mesh
void TPZHybridHDivErrorEstimator::ComputeAverageFacePressures() {
    DebugStop();
    TPZCompMesh *pressure = fOriginal->MeshVector()[1];
    TPZCompMesh *pressureHybrid = fPostProcMesh.MeshVector()[1];
    int fInterfaceMatid = fHybridizer.fLagrangeInterface;
    TPZGeoMesh *gmesh = pressure->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    pressure->LoadReferences();
    int64_t nel = pressureHybrid->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressureHybrid->Element(el);
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
            int64_t pos = pressureHybrid->Block().Position(seqnum);
            int ndof = c.NShape() * c.NState();
            for (int idf = 0; idf < ndof; idf++) {
                pressureHybrid->Solution()(pos + idf, 0) = L2Rhs(count++);
            }
        }
    }
    TPZManVector<TPZCompMesh *, 2> meshvec(2);
    meshvec[0] = fPostProcMesh.MeshVector()[0];
    meshvec[1] = fPostProcMesh.MeshVector()[1];
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, &fPostProcMesh);
}

/// compute the average pressures of across edges of the H(div) mesh
void TPZHybridHDivErrorEstimator::ComputeAveragePressures(int target_dim) {
    
    TPZCompMesh *pressureHybrid = PressureMesh();
    
    //    std::ofstream out("PressureToAverage.txt");
    //    pressureHybrid->Print(out);
    
    
    TPZGeoMesh *gmesh = pressureHybrid->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    int64_t nel = pressureHybrid->NElements();
    // load the pressure elements of dimension target_dim+1
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressureHybrid->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() != target_dim + 1) continue;
        gel->SetReference(cel);
    }
    
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressureHybrid->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        TPZMaterial *mat = pressureHybrid->FindMaterial(matid);
        TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
        //nao calcular a media para condicao de contorno
        if (bc) {
            // std::cout<<" No average processor for bc element"<<std::endl;
            continue ;
        }
        if (!cel || !cel->Reference() || cel->Reference()->Dimension() != target_dim) {
            continue;
        }
        ComputeAverage(pressureHybrid, el);
    }
    // apply the restraints to the edge connects
    if (target_dim == dim - 2) {
        pressureHybrid->LoadSolution(pressureHybrid->Solution());
        TransferEdgeSolution();
    }
    
}
//compute de L2 projection of Dirichlet boundary condition for Hdi-H1 reconstruction
void TPZHybridHDivErrorEstimator::ComputeBoundaryL2Projection(TPZCompMesh *pressuremesh, int target_dim){
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
    
    TPZElementMatrix ekbc, efbc;
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
        //efbc.fMat.Print("Solution",std::cout);
        int count = 0;
        int nc = cel->NConnects();
        for (int ic = 0; ic < nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            int64_t seqnum = c.SequenceNumber();
            int64_t pos = pressuremesh->Block().Position(seqnum);
            int ndof = c.NShape() * c.NState();
            for (int idf = 0; idf < ndof; idf++) {
                pressuremesh->Solution()(pos + idf, 0) = efbc.fMat(count++);
            }
        }
    }
    
    {
        std::ofstream out("PressureAfterL2Projection.txt");
        pressuremesh->Print(out);
    }
}

void TPZHybridHDivErrorEstimator::BoundaryPressureProjection(TPZCompMesh *pressuremesh, int target_dim){
    // TODO remove unused target_dim variable
    //    std::ofstream out("PressureProjBefore.txt");
    //    pressuremesh->Print(out);
    
    TPZGeoMesh *gmesh = fPostProcMesh.Reference();
    std::map<int,TPZMaterial *> matvec = fPostProcMesh.MaterialVec();
    for (auto matit=matvec.begin(); matit != matvec.end(); matit++) {
        TPZHDivErrorEstimateMaterial *castmat = dynamic_cast<TPZHDivErrorEstimateMaterial *>(matit->second);
        if(castmat)
        {
            TPZPressureProjection *pressproj = new TPZPressureProjection(*castmat);
            fPostProcMesh.MaterialVec()[matit->first] = pressproj;
        }
    }
    for (auto matit=matvec.begin(); matit != matvec.end(); matit++) {
        TPZBndCond *bndcond = dynamic_cast<TPZBndCond *>(matit->second);
        if(bndcond)
        {
            TPZMaterial *refmat = bndcond->Material();
            int matid = refmat->Id();
            bndcond->SetMaterial(fPostProcMesh.FindMaterial(matid));
        }
    }
    
    gmesh->ResetReference();
    

    //    {
    //        TPZCompMesh *pressmesh = fPostProcMesh.MeshVector()[1];
    //        std::ofstream out("pressure.txt");
    //        pressmesh->Print(out);
    //    }
    
    TPZElementMatrix ekbc,efbc;

    int64_t nel = fPostProcMesh.NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel_orig = fPostProcMesh.Element(el);
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
            // ekbc.fMat.Print(std::cout);
            // efbc.fMat.Print(std::cout);
            
            
            
            
            
            if(ekbc.HasDependency()) DebugStop();
            ekbc.fMat.SolveDirect(efbc.fMat, ECholesky);
            //   efbc.fMat.Print("Solution",std::cout);
            
            
            TPZCompEl *celpressure = mphys->Element(1);
            int nc = celpressure->NConnects();
            int count = 0;
            for (int ic = 0; ic < nc; ic++) {
                TPZConnect &c = celpressure->Connect(ic);
                int64_t seqnum = c.SequenceNumber();
                int64_t pos = pressuremesh->Block().Position(seqnum);
                int ndof = c.NShape() * c.NState();
                for (int idf = 0; idf < ndof; idf++) {
                    pressuremesh->Solution()(pos + idf, 0) = efbc.fMat(count++);
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
    for (auto matit=matvec.begin(); matit != matvec.end(); matit++) {
        TPZHDivErrorEstimateMaterial *castmat = dynamic_cast<TPZHDivErrorEstimateMaterial *>(matit->second);
        if(castmat)
        {
            int matid = matit->first;
            delete fPostProcMesh.MaterialVec()[matit->first];
            fPostProcMesh.MaterialVec()[matid] = matvec[matid];
        }
    }
    
    
    {
        std::ofstream out("PressureProjectionAfter.txt");
        pressuremesh->Print(out);
    }
    
}



void TPZHybridHDivErrorEstimator::NewComputeBoundaryL2Projection(
                                                                 TPZCompMesh *pressuremesh, int target_dim) {
    
    //{
    //    std::ofstream out("PressureBeforeL2Projection.txt");
    //    pressuremesh->Print(out);
    //}
    
    if (target_dim == 2) {
        std::cout << "Not implemented for 2D interface.\n";
        DebugStop();
    }
    
    TPZGeoMesh *gmesh = pressuremesh->Reference();
    gmesh->ResetReference();
    fOriginal->LoadReferences();
    
    int64_t nel = pressuremesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZCompEl *cel = pressuremesh->ElementVec()[iel];
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        int dim = gel->Dimension();
        
        int matid = gel->MaterialId();
        TPZMaterial *mat = pressuremesh->FindMaterial(matid);
        
        TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
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
            
            if (bc->HasForcingFunction()) {
                TPZManVector<STATE> result(3);
                TPZFNMatrix<9, STATE> gradu(dim, 1);
                TPZManVector<REAL, 3> x;
                gel->X(pt, x);
                bc->ForcingFunction()->Execute(x, result, gradu);
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
        for (int ic = 0; ic < nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            int64_t seqnum = c.SequenceNumber();
            int64_t pos = pressuremesh->Block().Position(seqnum);
            int ndof = c.NShape() * c.NState();
            for (int idf = 0; idf < ndof; idf++) {
                pressuremesh->Solution()(pos + idf, 0) = efbc(count++);
            }
        }
    }
    
    {
        std::ofstream out("PressureAfterL2Projection_2.txt");
        pressuremesh->Print(out);
    }
}

// compute the average of an element iel in the pressure mesh looking at its neighbours
void TPZHybridHDivErrorEstimator::ComputeAverage(TPZCompMesh *pressuremesh, int64_t iel)
{
    
    //  std::cout<<"Computing average for iel "<<iel<<"\n";
    TPZGeoMesh *gmesh = pressuremesh->Reference();
    int dim = gmesh->Dimension();
    TPZCompEl *cel = pressuremesh->Element(iel);
    
    if (!cel || !cel->Reference() ) {
        DebugStop();
    }
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
    TPZGeoEl *gel = cel->Reference();
    int target_dim = gel->Dimension();
    if (target_dim == dim - 1 && gel->MaterialId() != fPressureSkeletonMatId) {
        DebugStop();
        
    }
    if (!intel) {
        DebugStop();
    }
    int nc = cel->NConnects();
    int order = cel->Connect(nc - 1).Order();
    TPZGeoElSide gelside(gel, gel->NSides() - 1);
    TPZStack<TPZCompElSide> celstack;
    gelside.EqualLevelCompElementList(celstack, 1, 0);
    
    
    int nequal = celstack.size();
    TPZManVector<TPZTransform<REAL>, 4> tr(nequal);
    for (int ieq = 0; ieq < nequal; ieq++) {
        // the transformation between the sides
        tr[ieq] = gelside.NeighbourSideTransform(celstack[ieq].Reference());
        // add the transformation between the side and volume of the element
        TPZGeoEl *right = celstack[ieq].Element()->Reference();
        TPZTransform<REAL> tmp = right->SideToSideTransform(celstack[ieq].Side(), right->NSides() - 1);
        tr[ieq] = tmp.Multiply(tr[ieq]);
    }
    if (celstack.size() == 1 && target_dim == dim - 1) {
        TPZCompElSide lowlevel = gelside.LowerLevelCompElementList2(1);
        if (!lowlevel) {
            DebugStop();
        }
        celstack.Push(lowlevel);
        tr.Resize(2);
        tr[1] = TPZTransform<REAL>(gelside.Dimension());
        gelside.SideTransform3(lowlevel.Reference(), tr[1]);
    } else if (celstack.size() != 2 && target_dim == dim - 1) {
        DebugStop();
    }
    
    TPZManVector<REAL> el_weights(nequal);
    REAL sum_weights = 0.;
    REAL max_weight = 0.;
    for (int ieq = 0; ieq < nequal; ieq++) {
        REAL index = celstack[ieq].Element()->Index();
        el_weights[ieq] = fPressureweights[index];
        sum_weights += el_weights[ieq];
        max_weight = max_weight < el_weights[ieq] ? el_weights[ieq] : max_weight;
        if(el_weights[ieq] <= 0.) DebugStop();
    }
    fPressureweights[iel] = max_weight;
    
    std::unique_ptr<TPZIntPoints> intp(gel->CreateSideIntegrationRule(gel->NSides() - 1, 2 * order));
    int nshape = intel->NShapeF();
    TPZFNMatrix<20, REAL> L2Mat(nshape, nshape, 0.), L2Rhs(nshape, 1, 0.);
    TPZFNMatrix<220, REAL> phi(nshape, 1, 0.), dshape(dim, nshape);
    int64_t npoints = intp->NPoints();
    for (int64_t ip = 0; ip < npoints; ip++) {
        TPZManVector<REAL, 3> pt(target_dim, 0.), pt1(target_dim + 1, 0.), sol1(1);
        REAL weight;
        intp->Point(ip, pt, weight);
        intel->Shape(pt, phi, dshape);
        TPZManVector<REAL,3> xref(3);
        gel->X(pt, xref);
        //           std::cout << "Values " << sol1 << " " << sol2 << std::endl;
        //projecao L2 da media das soluceos no espaco Lh, do esqueleto da malha
        for (int ieq = 0; ieq < nequal; ieq++) {
            tr[ieq].Apply(pt, pt1);
            TPZManVector<REAL,3> xeq(3,0.);
            // celstack[ieq].Element()->Print();
            celstack[ieq].Element()->Reference()->X(pt1,xeq);
            celstack[ieq].Element()->Solution(pt1, 0, sol1);//solucao a esquerda
                                                            //      std::cout << "xref " << xref << " ieq " << ieq << " xeq " << xeq << " sol " << sol1 << std::endl;
            for (int ishape = 0; ishape < nshape; ishape++)
            {
                L2Rhs(ishape, 0) += weight * phi(ishape, 0) * sol1[0] * el_weights[ieq] / (sum_weights);
            }
        }
        for (int ishape = 0; ishape < nshape; ishape++) {
            for (int jshape = 0; jshape < nshape; jshape++)
            {
                L2Mat(ishape, jshape) += weight * phi(ishape, 0) * phi(jshape, 0);
            }
        }
    }
    L2Mat.SolveDirect(L2Rhs, ECholesky);
    //apos este passo temos uma pressao que é continua ao longo das interfaces dos elementos, nos esqueletos. Falta suavizar nos vértices
    // L2Rhs.Print("Average pressure");
    
    //    std::cout << "average ";
    //    for (int i=0; i<nshape; i++) {
    //        std::cout << L2Rhs(i,0) << " ";
    //    }
    //    std::cout << std::endl;
    //store the average on solution of pressure mesh
    int count = 0;
    
    for (int ic = 0; ic < nc; ic++) {
        TPZConnect &c = cel->Connect(ic);
        int64_t seqnum = c.SequenceNumber();
        int64_t pos = pressuremesh->Block().Position(seqnum);
        int ndof = c.NShape() * c.NState();
        for (int idf = 0; idf < ndof; idf++) {
            pressuremesh->Solution()(pos + idf, 0) = L2Rhs(count++);
        }
    }
    
}


/// transfer the solution of the edge functions to the face functions
void TPZHybridHDivErrorEstimator::TransferEdgeSolution() {
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
                    pressureHybrid->Block()(neighblock, 0, i, 0) = pressureHybrid->Block()(edge_seqnum, 0, i, 0);
                }
            }
        }
    }
}


/// set the cornernode values equal to the averages
void TPZHybridHDivErrorEstimator::ComputeNodalAverages() {
    TPZCompMesh *pressureHybrid = PressureMesh();
    int fInterfaceMatid = fHybridizer.fLagrangeInterface;
    TPZGeoMesh *gmesh = pressureHybrid->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    pressureHybrid->LoadReferences();
    TPZMaterial *mat = pressureHybrid->FindMaterial(fInterfaceMatid);
    if (!mat) DebugStop();
    int nstate = mat->NStateVariables();
    int64_t nel = pressureHybrid->NElements();
    // load the pressure elements of dimension dim-1
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressureHybrid->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() != dim - 1) continue;
        gel->SetReference(cel);
    }
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressureHybrid->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        
        if (gel->Dimension() != dim - 1 || gel->MaterialId() != fInterfaceMatid) {
            // MatId of element not important for nodal average computation
            continue;
        }
        
        //percorre cada no do elemento de interface
        int ncorners = gel->NCornerNodes();
        for (int side = 0; side < ncorners; side++) {
            TPZGeoElSide gelside(gel, side);
            TPZCompElSide celside(intel,side);
            ComputeNodalAverage(celside);
        }
    }
}

/// compute the nodal average of all elements that share a point
void TPZHybridHDivErrorEstimator::ComputeNodalAverage(TPZCompElSide &celside)
{
    
    int lagrangematid = fHybridizer.fLagrangeInterface;
    TPZCompMesh *pressureHybrid = celside.Element()->Mesh();
    int dim = pressureHybrid->Dimension();
    TPZMaterial *mat = pressureHybrid->FindMaterial(lagrangematid);
    if (!mat) DebugStop();
    
    int nstate = mat->NStateVariables();
    TPZGeoElSide gelside(celside.Reference());
    TPZGeoEl *gel = gelside.Element();
    int side = gelside.Side();
    // celstack will contain all zero dimensional sides connected to the side
    TPZStack<TPZCompElSide> celstack;
    int onlyinterpolated = 1;
    int removeduplicates = 0;
    
    gelside.ConnectedCompElementList(celstack, onlyinterpolated, removeduplicates);
    celstack.Push(gelside.Reference());
    // This map stores the connects, the weight associated with the element
    // and the solution of that connect. The weight of Dirichlet condition is
    // higher and will be used later to impose the value of the BC in the
    // connects when needed
    std::map<int64_t, std::pair<REAL, TPZVec<STATE>>> connects;
    //para cada elemento que tem este noh procede como segue
    for (int elc = 0; elc < celstack.size(); elc++) {
        TPZCompElSide celside = celstack[elc];
        TPZGeoElSide gelside0 = celside.Reference();
        
        if (gelside0.Element()->Dimension() != dim - 1) {
            continue;
        }
        TPZInterpolatedElement *intel1 = dynamic_cast<TPZInterpolatedElement *>(celside.Element());
        if (!intel1) DebugStop();
        int64_t index = intel1->Index();
        REAL weight = fPressureweights[index];
        
        if(IsZero(weight)) {
            TPZMaterial *mat = intel1->Material();
            TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
            if (!bc) DebugStop();
            else{
                
                continue;//duvida
            }
            
        }
        
        int64_t conindex = intel1->ConnectIndex(celside.Side());
        
        if (connects.find(conindex) != connects.end()) {
            
            DebugStop();
        }
        
        TPZConnect &c = intel1->Connect(celside.Side());
        int64_t seqnum = c.SequenceNumber();
        if (c.NState() != nstate || c.NShape() != 1) DebugStop();
        TPZManVector<REAL,3> pt(0), x(3);
        gelside0.X(pt, x);
        TPZManVector<STATE, 3> sol(nstate, 0.);
        for (int istate = 0; istate < nstate; istate++) {
            
            sol[istate] = pressureHybrid->Block().Get(seqnum, 0, istate, 0);
            //      std::cout<<"sol "<<sol[istate]<<"\n";
        }
        //  std::cout<<"conindex "<<conindex << " weight "<<weight<<"\n";
        connects.insert({conindex, {weight, sol}});
    }
    
    TPZManVector<STATE, 3> averageSol(nstate, 0);
    int nconnects = 0;
    
    REAL sum_weight = 0.;
    for (auto it = connects.begin(); it != connects.end(); it++) {
        REAL weight = it->second.first;
        sum_weight += weight;
        // If the node is located on the Dirichlet BC, internal sides don't
        // contribute to the average of the connects
        for (int istate = 0; istate < nstate; istate++) {
            averageSol[istate] += it->second.second[istate]*weight;
            //   std::cout<<"averageSol += "<<it->second.second[istate] <<" * weight "<<weight<<"\n";
            nconnects++;
        }
    }
    
    for (int istate = 0; istate < nstate; istate++) {
        averageSol[istate] /= sum_weight;
        // std::cout<<"averasol /= "<< averageSol[istate] <<" sum_weight "<<sum_weight<<"\n";
    }
    
    for (auto it = connects.begin(); it != connects.end(); it++) {
        int64_t conindex = it->first;
        TPZConnect &c = pressureHybrid->ConnectVec()[conindex];
        int64_t seqnum = c.SequenceNumber();
        
        if (c.NState() != nstate || c.NShape() != 1) DebugStop();
        for (int istate = 0; istate < nstate; istate++) {
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "value before " << pressureHybrid->Block()(seqnum, 0, istate, 0) <<
                " value after " << averageSol[istate] << " diff "
                << pressureHybrid->Block()(seqnum, 0, istate, 0) - averageSol[istate] << " ncontributing connects "
                << nconnects;
                //                        res2.Print("Residual",sout);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            pressureHybrid->Block()(seqnum, 0, istate, 0) = averageSol[istate];
        }
    }
    
}


/// clone the meshes into the post processing mesh
void TPZHybridHDivErrorEstimator::CloneMeshVec() {
    
    for (int i = 0; i < fOriginal->MeshVector().size(); i++) {
        fPostProcMesh.MeshVector()[i] = fOriginal->MeshVector()[i]->Clone();
    }
    
}

/// compute the effectivity indices of the pressure error and flux error and store in the element solution
void TPZHybridHDivErrorEstimator::ComputeEffectivityIndices(TPZSubCompMesh *subcmesh) {
    int64_t nrows = subcmesh->ElementSolution().Rows();
    int64_t ncols = 5; // subcmesh->ElementSolution().Cols();

    if (subcmesh->ElementSolution().Cols() != 7) {
        // TODO I made some changes to be able to run the code again.
        //  Sometimes subcmesh->ElementSolution().Cols() equals 5  sometimes it's already 7.
        //  I'm not sure if this behaviour is expected.
        subcmesh->ElementSolution().Resize(nrows, 7);
    }

    int64_t nel = subcmesh->NElements();
    TPZFMatrix<STATE> &elsol = subcmesh->ElementSolution();
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

            REAL oscillatorytherm = 0;
            if (i == 2) {
                oscillatorytherm = subcmesh->ElementSolution()(el, i + 2);
                oscillatorytherm *= (hk / M_PI);
            }
            oscillatorytherm = 0;

            if (abs(ErrorEstimate) < tol) {
                subcmesh->ElementSolution()(el, ncols + i / 2) = 1.;

            } else {
                REAL EfIndex = (ErrorEstimate + oscillatorytherm) / ErrorExact;
                subcmesh->ElementSolution()(el, ncols + i / 2) = EfIndex;
            }
        }
    }
}

/// compute the effectivity indices of the pressure error and flux error and store in the element solution
void TPZHybridHDivErrorEstimator::ComputeEffectivityIndices() {
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
    int64_t nrows = cmesh->ElementSolution().Rows();
    int64_t ncols = cmesh->ElementSolution().Cols();
    
    //std::ostream &out;
    //    cmesh->ElementSolution().Print("ElSolution",std::cout);
    
    
    TPZFMatrix<REAL> dataIeff(nrows,1);
    
    //    REAL oscilatorytherm = 0;
    int dim = cmesh->Dimension();
    cmesh->ElementSolution().Resize(nrows, ncols+2);

    for (int64_t el = 0; el < nrows; el++) {
        
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if(subcmesh)
        {
            ComputeEffectivityIndices(subcmesh);
        }
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        if(gel->Dimension() != dim) continue;
        int nsides = gel->NSides();
        for(int is = 0; is<nsides; is++)
        {
            if(gel->SideDimension(is) != dim-1) continue;
            TPZStack<TPZCompElSide> equal;
            TPZGeoElSide gelside(gel,is);
            gelside.EqualLevelCompElementList(equal, 0, 0);
//            if(equal.size() != 1){
//                std::cout<<"Number of neighbour "<<equal.size()<<"\n";
//                DebugStop();
//            }
            TPZGeoElSide neighbour;
            TPZCompElSide selected;
            for(int i=0; i<equal.size(); i++)
            {
                TPZGeoEl *gequal = equal[i].Element()->Reference();
                int eldim = gequal->Dimension();
                if(eldim != dim-1) continue;
                int elmatid = gequal->MaterialId();
                if(fProblemConfig.bcmaterialids.find(elmatid) != fProblemConfig.bcmaterialids.end())
                {
                    neighbour = equal[i].Reference();
                    selected = equal[i];
                    break;
                }
            }

            if(!neighbour) continue;
            if(neighbour.Element()->Dimension() != dim-1) DebugStop();
            int64_t neighindex = selected.Element()->Index();
            for (int i = 0; i < 3; i += 2) {
                

                std::cout << "linha = " << el << " col = " << 4 + i / 2 << " neinEl " << neighindex << std::endl;
                
                if(neighindex > nrows){
                    std::cout<<" neighindex= "<< neighindex<<" nrows "<<nrows<<"\n";
                    DebugStop();
                }
                
                REAL NeighbourErrorEstimate = cmesh->ElementSolution()(neighindex, i + 1);
                REAL NeighbourErrorExact = cmesh->ElementSolution()(neighindex, i);
                REAL ErrorEstimate = cmesh->ElementSolution()(el, i + 1);
                REAL ErrorExact = cmesh->ElementSolution()(el, i);
                REAL sumErrorExact = sqrt(NeighbourErrorExact*NeighbourErrorExact+ErrorExact*ErrorExact);
                REAL sumErrorEstimate = sqrt(NeighbourErrorEstimate*NeighbourErrorEstimate+ErrorEstimate*ErrorEstimate);
                cmesh->ElementSolution()(neighindex,i+1) = 0.;
                cmesh->ElementSolution()(neighindex,i) = 0.;
                cmesh->ElementSolution()(el,i) = sumErrorExact;
                cmesh->ElementSolution()(el,i+1) = sumErrorEstimate;
            }
        }
    }
    for (int64_t el = 0; el < nrows; el++) {
        
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if(subcmesh)
        {
            ComputeEffectivityIndices(subcmesh);
        }
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        REAL hk = gel->CharacteristicSize();
        
        
        for (int i = 0; i < 3; i += 2) {
            
            //  std::cout<<"linha = "<<el<< "col = "<<4 + i / 2<<std::endl;
            
            REAL tol = 1.e-10;
            REAL ErrorEstimate = cmesh->ElementSolution()(el, i + 1);
            REAL ErrorExact = cmesh->ElementSolution()(el, i);
            
            TPZGeoEl *gel = cel->Reference();
            
            REAL hk = gel->CharacteristicSize();

            REAL oscilatorytherm = 0;
            if (i == 2) {
                oscilatorytherm = cmesh->ElementSolution()(el, i + 2);
                oscilatorytherm *= (hk / M_PI);
            }

            if (abs(ErrorEstimate) < tol) {
                cmesh->ElementSolution()(el, ncols + i / 2) = 1.;
                dataIeff(el,0)=1.;
            }
            else {
                REAL EfIndex = (ErrorEstimate + oscilatorytherm)/ ErrorExact;
                dataIeff(el,0)= EfIndex;
                cmesh->ElementSolution()(el, ncols + i / 2) = EfIndex;
            }
        }
        
    }
    
    //  cmesh->ElementSolution().Print("ElSolution",std::cout);
    //    ofstream out("IeffPerElement3DEx.nb");
    //    dataIeff.Print("Ieff = ",out,EMathematicaInput);


}

/// returns true if the material associated with the element is a boundary condition
/// and if the boundary condition is dirichlet type
bool TPZHybridHDivErrorEstimator::IsDirichletCondition(TPZGeoElSide gelside) {
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

/// return the value of the Dirichlet condition
void TPZHybridHDivErrorEstimator::GetDirichletValue(TPZGeoElSide gelside, TPZVec<STATE> &vals) {
    TPZGeoEl *gel = gelside.Element();
    int matid = gel->MaterialId();
    TPZMaterial *mat = fPostProcMesh.FindMaterial(matid);
    TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
    if (!bc) DebugStop();
    int typ = bc->Type();
    if (typ != 0) DebugStop();
    //TPZManVector<REAL,3> xco(3,0.);
    TPZVec<REAL> xco;
    xco.resize(3);
    if (bc->HasForcingFunction()) {
        gel->NodePtr(gelside.Side())->GetCoordinates(xco);
        bc->ForcingFunction()->Execute(xco, vals);
    } else {
        int nv = vals.size();
        for (int iv = 0; iv < nv; iv++) {
            vals[iv] = bc->Val2()(iv, 0);
        }
    }
}

void TPZHybridHDivErrorEstimator::PotentialReconstruction() {
    
    if (fPostProcMesh.MeshVector().size()) {
        DebugStop();
    }

    // Create the post processing mesh (hybridized H(div) mesh) with increased approximation order
    // for the border fluxes
    //    {
    //
    //        std::ofstream out("PressureOriginalBeforeProcessing.txt");
    //        fOriginal->MeshVector()[1]->Print(out);
    //
    //    }

    CreatePostProcessingMesh();
    
    //    {
    //        std::ofstream out("PressureOriginalPostProcessing.txt");
    //        fOriginal->MeshVector()[1]->Print(out);
    //        std::ofstream out2("PressureRecPostProcessing.txt");
    //        fPostProcMesh.MeshVector()[1]->Print(out2);
    //        std::ofstream outgvtk("gmesh.vtk");
    //        TPZVTKGeoMesh::PrintGMeshVTK(fPostProcMesh.Reference(), outgvtk);
    //    }
    
    // L2 projection for Dirichlet and Robin boundary condition for H1 reconstruction
    if (!fPostProcesswithHDiv) {
        TPZCompMesh *pressuremesh = PressureMesh();
        int target_dim = 1;
        // TODO ver se fica igual para dimensao maior
        ComputeBoundaryL2Projection(pressuremesh, target_dim);
        //BoundaryPressureProjection(pressuremesh, target_dim);
    }

    {
        std::ofstream out("PressureAverageMesh.txt");
        fPostProcMesh.MeshVector()[1]->Print(out);
        PlotLagrangeMultiplier("BeforeAverage");
    }
    
    // Calculates average pressure on interface edges and vertices
    if (fProblemConfig.makepressurecontinuous) {
        int dim = fPostProcMesh.Dimension();
        ComputeAveragePressures(dim - 1);
        // in three dimensions make the one-d polynoms compatible
        if (dim == 3) {
            ComputeAveragePressures(1);
        }
    }
    
    {
        std::ofstream out("PressureAverageMesh.txt");
        fPostProcMesh.MeshVector()[1]->Print(out);
        PlotLagrangeMultiplier("BeforeNodalAverage");
    }
    
    if(fProblemConfig.makepressurecontinuous)
    {
        ComputeNodalAverages();
    }
    {
        std::ofstream out("PressureNodalMesh.txt");
        fPostProcMesh.MeshVector()[1]->Print(out);
        PlotLagrangeMultiplier("AfterNodalAverage");
    }
    
    
    // in the case of hybrid hdiv, computing the error using h(div) spaces, nothing will be done
    if (!fPostProcesswithHDiv) {
        // std::ofstream file("PressureMeshToCopySkeleton.vtk");
        // TPZVTKGeoMesh::PrintGMeshVTK(fPostProcMesh.Reference(), file);
        CopySolutionFromSkeleton();
    }
    // transfer the continuous pressures to the multiphysics space
    {
        TPZManVector<TPZCompMesh *, 2> meshvec(2);
        meshvec[0] = fPostProcMesh.MeshVector()[0];//flux
        meshvec[1] = fPostProcMesh.MeshVector()[1];//pressure
        
        
        
        {
            std::ofstream out("MultiphysicsBeforeTransfer.txt");
            fPostProcMesh.Print(out);
            //  PlotState("MultiphysicsBeforeTransfer2D.vtk", 2, &fPostProcMesh);
        }
        
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, &fPostProcMesh);
        
        {
            std::ofstream out("MultiphysicsAfterTransfer.txt");
            fPostProcMesh.Print(out);
            PlotState("MultiphysicsAfterTransfer2D.vtk", fPostProcMesh.Dimension(), &fPostProcMesh);
        }
        
        
    }
    
    
    
    
    
    
#ifdef PZDEBUG
    {
        std::ofstream out("MeshWithSmoothPressure.txt");
        fPostProcMesh.Print(out);
        std::ofstream out2("PressureMeshSmooth.txt");
        fPostProcMesh.MeshVector()[1]->Print(out2);
        //        std::ofstream out3("FluxMeshSmooth.txt");
        //        fPostProcMesh.MeshVector()[0]->Print(out3);
    }
#endif
    
    
    //Resolver problema local com potencial continuo como condicao de Dirichlet
    //    {
    //        std::ofstream file("MeshToComputeStiff.vtk");
    //        TPZVTKGeoMesh::PrintGMeshVTK(fPostProcMesh.Reference(), file);
    //
    //
    //    }
    
    ComputeElementStiffnesses();
    
#ifdef PZDEBUG2
    {
        std::ofstream out("MeshBeforeLoadSol.txt");
        fPostProcMesh.Print(out);
        ofstream out2("SolBeforeLoadSolution.nb");
        fPostProcMesh.Solution().Print("SolBeforeLoadSolution=",out2,EMathematicaInput);
    }
#endif
    
    fPostProcMesh.LoadSolution(fPostProcMesh.Solution());
    
#ifdef PZDEBUG2
    {
        std::ofstream out("MeshAfterLoadSol.txt");
        fPostProcMesh.Print(out);
        //fPostProcMesh.Solution().Print("SolAfterLoadSolution");
        ofstream out2("SolAfterLoadSolution.nb");
        fPostProcMesh.Solution().Print("SolAfterLoadSolution=",out2,EMathematicaInput);
    }
#endif
    
    {
        TPZManVector<TPZCompMesh *, 2> meshvec(2);
        // fPostProcMesh[0] is the H1 or Hdiv mesh
        // fPostProcMesh[1] is the L2 mesh
        
        meshvec[0] = fPostProcMesh.MeshVector()[0];
        meshvec[1] = fPostProcMesh.MeshVector()[1];
        
        //        {
        //            //fPostProcMesh.Solution().Print("SolutionBeforeTranfer");
        //
        //            ofstream out2("SolutionBeforeTranfer.nb");
        //            fPostProcMesh.Solution().Print("SolutionBeforeTranfer=",out2,EMathematicaInput);
        //
        //        }
        
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, &fPostProcMesh);
        
        //        {
        //
        //       //  fPostProcMesh.ElementSolution().Print("SolutionAfterTranfer");
        //            ofstream out2("SolutionAfterTranfer.nb");
        //            fPostProcMesh.Solution().Print("SolutionAfterTranfer=",out2,EMathematicaInput);
        //
        //        }
        
#ifdef PZDEBUG2
        {
            std::ofstream out("PressureAfterTransfer.txt");
            fPostProcMesh.MeshVector()[1]->Print(out);
            std::ofstream out2("FluxAfterTransfer.txt");
            //fPostProcMesh.MeshVector()[0]->Print(out2);
            
        }
        VerifySolutionConsistency(PressureMesh());
#endif
    }
}

void TPZHybridHDivErrorEstimator::PlotLagrangeMultiplier(const std::string &filename, bool reconstructed) {
    
    TPZCompMesh *pressure = nullptr;

    if (!reconstructed) {
        pressure = fOriginal->MeshVector()[1];
    } else {
        pressure = PressureMesh();
    }

    std::ofstream out2("PressuretoStateGraph.txt");
    pressure->Print(out2);
    
    {
        TPZAnalysis an(pressure, false);
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
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(2, dim);
    }
    
}

static TPZMultiphysicsInterfaceElement *Extract(TPZElementGroup *cel)
{
    const TPZVec<TPZCompEl *> &elgr = cel->GetElGroup();
    for(int i=0; i<elgr.size(); i++)
    {
        TPZMultiphysicsInterfaceElement *interf = dynamic_cast<TPZMultiphysicsInterfaceElement *>(elgr[i]);
        if(interf) return interf;
    }
    return NULL;
}

static TPZMultiphysicsInterfaceElement *Extract(TPZCondensedCompEl *cond)
{
    TPZCompEl *cel = cond->ReferenceCompEl();
    TPZMultiphysicsInterfaceElement *interf = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cond);
    if(interf) return interf;
    TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
    if(elgr) return Extract(elgr);
    return NULL;
}

static TPZMultiphysicsInterfaceElement *Extract(TPZCompEl *cel)
{
    TPZMultiphysicsInterfaceElement *interf = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
    if(interf) return interf;
    TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>(cel);
    if(cond) return Extract(cond);
    return NULL;
}
/// identify the peripheral material objects and store the information in fHybridizer
void TPZHybridHDivErrorEstimator::IdentifyPeripheralMaterialIds() {
    int dim = fOriginal->Dimension();
    // identify the material id for interface elements
    int64_t nel = fOriginal->NElements();
    int numint_found = 0;
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fOriginal->Element(el);
        TPZMultiphysicsInterfaceElement *interf = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
        if(!interf)
        {
            TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>(cel);
            if(cond) interf = Extract(cond);
        }
        if (interf) {
            int matid = interf->Reference()->MaterialId();
            if(numint_found == 0)
            {
                fHybridizer.fInterfaceMatid.first = matid;
                numint_found++;
                break;
            }
            else if(numint_found == 1 && fHybridizer.fInterfaceMatid.first != matid) {
                fHybridizer.fInterfaceMatid.second = matid;
                numint_found++;
                break;
            }
        }
    }
    /// identify the material id of the pressure
    TPZCompMesh *pressure_mesh = fOriginal->MeshVector()[1];
    nel = pressure_mesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        if (cel) {
            TPZGeoEl *gel = cel->Reference();
            if (gel && gel->Dimension() == dim - 1) {
                fHybridizer.fLagrangeInterface = gel->MaterialId();
                break;
            }
        }
    }
    /// identify the material id of boundary elements
    TPZCompMesh *fluxmesh = fOriginal->MeshVector()[0];
    for (auto map_pair : fluxmesh->MaterialVec()) {
        TPZMaterial *mat = map_pair.second;
        TPZNullMaterial *nullmat = dynamic_cast<TPZNullMaterial *>(mat);
        if (nullmat) {
            fHybridizer.fHDivWrapMatid = nullmat->Id();
            break;
        }
    }
    
    fPressureSkeletonMatId = fHybridizer.fLagrangeInterface;
}

/// switch material object from mixed poisson to TPZMixedHdivErrorEstimate
void TPZHybridHDivErrorEstimator::SwitchMaterialObjects() {
    
    if(fPostProcesswithHDiv) {
        for (auto mat : fPostProcMesh.MaterialVec()) {
            TPZMixedPoisson *mixpoisson =
            dynamic_cast<TPZMixedPoisson *>(mat.second);
            if (mixpoisson) {
                TPZMixedHDivErrorEstimate<TPZMixedPoisson> *newmat = new TPZMixedHDivErrorEstimate<TPZMixedPoisson>(*mixpoisson);
                
                if (mixpoisson->HasForcingFunction()) {
                    newmat->SetForcingFunctionExact(mixpoisson->ForcingFunctionExact());
                    newmat->SetForcingFunction(mixpoisson->ForcingFunction());
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
    else {
        // switch the material of the HDiv approximation to a material for an H1 approximation
        for(auto mat : fPostProcMesh.MaterialVec())
        {
            TPZMixedPoisson *mixpoisson = dynamic_cast<TPZMixedPoisson *> (mat.second);
            if(mixpoisson)
            {
                TPZHDivErrorEstimateMaterial *newmat = new TPZHDivErrorEstimateMaterial(*mixpoisson);
                
                if (mixpoisson->HasForcingFunction()) {
                    newmat->SetForcingFunctionExact(mixpoisson->ForcingFunctionExact());
                    newmat->SetForcingFunction(mixpoisson->ForcingFunction());
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
}


void TPZHybridHDivErrorEstimator::VerifySolutionConsistency(TPZCompMesh *cmesh) {
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
                        std::stringstream sout;
                        sout << "\nSide solution =  " << sol0[0] << "\n";
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

void TPZHybridHDivErrorEstimator::PrepareElementsForH1Reconstruction() {
    
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
            if (!large) DebugStop();
            compNeighSides.Push(large);
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
                connectsToIncrement[seqNum] = 1;
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
        TPZMultiphysicsElement *mcel =
        dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mcel) DebugStop();
        
        std::set<int64_t> elementsToGroup;
        for (int iSide = 0; iSide < gel->NSides(); iSide++) {
            // AskPhil: is side creation expensive? Should we iterate from
            //  iSide = gel->NCornerNodes to gel->NSides-1 instead?
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
    
    for (auto matit : fPostProcMesh.MaterialVec()) {
        TPZMaterial *mat = matit.second;
        TPZHDivErrorEstimateMaterial *errormat =
        dynamic_cast<TPZHDivErrorEstimateMaterial *>(mat);
        if (errormat) {
            errormat->fNeumannLocalProblem = false;
        }
    }
    
    fPostProcMesh.CleanUpUnconnectedNodes();
    
#ifdef PZDEBUG2
    std::ofstream outTXT("PostProcessMeshAfterPreparingElements.txt");
    fPostProcMesh.Print(outTXT);
#endif
}

void TPZHybridHDivErrorEstimator::CopySolutionFromSkeleton() {
    
    TPZCompMesh *pressuremesh = PressureMesh();
    //        {
    //            std::ofstream out("MeshBeforeCopySkeleton.txt");
    //            pressuremesh->Print(out);
    //        }
    pressuremesh->Reference()->ResetReference();
    pressuremesh->LoadReferences();
    int dim = pressuremesh->Dimension();
    int64_t nel = pressuremesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressuremesh->Element(el);
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
            int nst = celstack.NElements();
            for (int ist = 0; ist < nst; ist++) {
                TPZCompElSide cneigh = celstack[ist];
                TPZGeoElSide gneigh = cneigh.Reference();
                if (gneigh.Element()->MaterialId() == this->fHybridizer.fLagrangeInterface || (gneigh.Element()->MaterialId() == this->fPressureSkeletonMatId )||
                    IsDirichletCondition(gneigh)) {
                    TPZInterpolatedElement *intelneigh = dynamic_cast<TPZInterpolatedElement *>(cneigh.Element());
                    if (!intelneigh) DebugStop();
                    TPZConnect &con_neigh = intelneigh->Connect(cneigh.Side());
                    int64_t con_seqnum = con_neigh.SequenceNumber();
                    int con_size = con_neigh.NState() * con_neigh.NShape();
                    if (con_size != c_blocksize) DebugStop();
                    for (int ibl = 0; ibl < con_size; ibl++) {
                        pressuremesh->Block()(c_seqnum, 0, ibl, 0) = pressuremesh->Block()(con_seqnum, 0, ibl, 0);
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
    //        {
    //            std::ofstream out("PrssureMeshAfterCopySkeleton.txt");
    //            pressuremesh->Print(out);
    //            std::ofstream out2("MeshAfterCopySkeleton.txt");
    //            fPostProcMesh.Print(out2);
    //        }
}


/// compute the pressure weights and material weights
// fills in the data structure pressureweights and matid_weights
void TPZHybridHDivErrorEstimator::ComputePressureWeights() {
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
        if (matid == fPressureSkeletonMatId ||
            matid == fHybridizer.fLagrangeInterface) {
            // fPressureweights[el] = 1.e12;
            continue;
        }
        if (!mat)  DebugStop();
        
        
        // if (gel->Dimension() != dim /*&& !bcmat*/) continue;
        TPZBndCond *bcmat = dynamic_cast<TPZBndCond *>(mat);
        
        if (bcmat) {
            if (bcmat->Type() == 0) {
                this->fPressureweights[el] = 1.e12;
                fMatid_weights[matid] = 1.e12;
                continue;
            }
            else if (bcmat->Type() == 4){
                this->fPressureweights[el] = bcmat->Val1()(0,0);//valor de Km
                fMatid_weights[matid] = bcmat->Val1()(0,0);
                
            }
            
            else {
                this->fPressureweights[el] = 0.;
                fMatid_weights[matid] = 0.;
                continue;
            }
        }
        
        else{
            
            
            TPZMixedPoisson *mixpoisson = dynamic_cast<TPZMixedPoisson *>(mat);
            if (!mixpoisson) DebugStop();
            
            REAL perm;
            mixpoisson->GetMaxPermeability(perm);
            if (IsZero(perm)) DebugStop();
            this->fPressureweights[el] = perm;
            fMatid_weights[matid] = perm;
        }
    }
}

void TPZHybridHDivErrorEstimator::PlotState(const std::string& filename, int targetDim, TPZCompMesh* cmesh) {
    
    std::ofstream out2("PressuretoStateGraph.txt");
    cmesh->Print(out2);
    
    {
        TPZAnalysis an(cmesh, false);
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("State");
        scalnames.Push("PressureReconstructed");
        
        std::string plotname;
        {
            std::stringstream out;
            out << filename << ".vtk";
            plotname = out.str();
        }
        an.DefineGraphMesh(targetDim, scalnames, vecnames, plotname);
        an.PostProcess(2, targetDim);
    }
}
