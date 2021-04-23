
//
//  TPZHDivErrorEstimator.cpp
//  ErrorEstimateHDiv
//
//  Created by Philippe Devloo on 10/06/18.
//

//#include "TPZHDivErrorEstimator.h" // DELETE ME
#include "TPZHybridH1ErrorEstimator.h"
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

//#include "TPZHDivErrorEstimateMaterial.h"

#include "TPZCompMeshTools.h"
#include "TPZPressureProjection.h"
#include "TPZMatLaplacianHybrid.h"
#include "TPZHybridH1ErrorEstimateMaterial.h"
#include "TPZGeoElSideAncestors.h"
#include "TPZGeoElSidePartition.h"

#include "TPZHDivErrorEstimateMaterial.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("HDivErrorEstimator"));
static LoggerPtr loggerF(Logger::getLogger("DebuggingF"));

#endif

TPZHybridH1ErrorEstimator::~TPZHybridH1ErrorEstimator() {
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
void TPZHybridH1ErrorEstimator::ComputeErrors(TPZVec<REAL>& errorVec, TPZVec<REAL> &elementerrors, bool store = true) {
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
    
    int64_t nErrorCols = 8;
    errorVec.resize(nErrorCols);
    errorVec.Fill(0);
    for (int64_t i = 0; i < nErrorCols; i++) {
        errorVec[i] = 0;
    }
    
    int64_t nelem = fPostProcMesh.NElements();
    fPostProcMesh.LoadSolution(fPostProcMesh.Solution());
    fPostProcMesh.ExpandSolution();
    fPostProcMesh.ElementSolution().Redim(nelem, nErrorCols-1);
    for(int64_t el = 0; el<nelem; el++)
    {
        TPZCompEl *cel = fPostProcMesh.Element(el);
        TPZSubCompMesh *subc = dynamic_cast<TPZSubCompMesh *>(cel);
        if(subc)
        {
            int64_t nelsub = subc->NElements();
            subc->ElementSolution().Redim(nelsub, 6);
        }
    }
    
#ifdef PZDEBUG1
    {
        std::ofstream out("MeshToComputeError2.txt");
        fPostProcMesh.Print(out);
        
    }
#endif

    an.PostProcessError(errorVec, store);//calculo do erro com sol exata e aprox e armazena no elementsolution
    
    std::cout << "\n\n############\n\n";
    
    TPZCompMeshTools::UnCondensedElements(&fPostProcMesh);
    TPZCompMeshTools::UnGroupElements(&fPostProcMesh);
    
    //Erro global
    std::ofstream myfile;
    myfile.open("ErrorsReconstruction.txt", std::ios::app);
    myfile << "\n\n Estimator errors for Problem " << fProblemConfig.problemname;
    myfile << "\n-------------------------------------------------- \n";
    myfile << "Ndiv = " << fProblemConfig.ndivisions <<" Order k= " << fProblemConfig.k << " Order n= "<< fProblemConfig.n<<"\n";
    //myfile << "DOF Total = " << fPostProcMesh.NEquations() << "\n";
    myfile << "||u-u_h|| = " << errorVec[0] << "\n";
    myfile << "||u_h-s_h|| = " << errorVec[1] << "\n";
    myfile << "e_{ex}: ||K^{0.5}.grad(u_h-u)|| = " << errorVec[2] << "\n";
    myfile << "n_{NC}: ||K^{0.5}.grad(u_h-s_h)|| = " << errorVec[3] << "\n";
    myfile << "n_{F} : ||K^{0.5}.[grad(u_h)-invK.T_h]|| = " << errorVec[5] << "\n";
    myfile << "||f-Proj(f)|| = " << errorVec[6] << "\n";
    myfile << "||f-Div(T_h)|| = " << errorVec[4] << "\n";
    //myfile <<"Residual ErrorL2= "<< errorVec[4] << "\n";
    //myfile <<"Global Index = "<< sqrt(errorVec[4] + errorVec[3]) / sqrt(errorVec[2]);
    
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
        TPZFMatrix<STATE> &elsol = fPostProcMesh.ElementSolution();
        elementerrors[gel->Index()] = elsol(i, 3);
    }
    
}

void TPZHybridH1ErrorEstimator::GlobalEffectivityIndex(){

    ////VO: Pesquisar/Perguntar pro Phil:{
    ////O material mixedpoisson têm 5 erros: [2
        ///error[0] = L2 for pressure
        ///error[1] = L2 for flux
        ///error[2] = L2 for div(flux)
        ///error[3] = Grad pressure (Semi H1)
        ///error[4] = Hdiv norm
    ////TPZMatLaplacianHybrid têm 4:
        ///error[0] = L2 for pressure
        ///error[1] = Grad pressure (Semi H1)
        ///error[2] = H1 norm
        ///error[3] = energy norm: Permeability*Semi_H1}


    
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
        cel->EvaluateError(elerror, false);
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
    
    std::ofstream myfile;
    myfile.open("ErrorEstimationResults.txt", std::ios::app);
    myfile << "\n\n Estimator errors for Problem "
    << fProblemConfig.problemname;
    myfile << "\n-------------------------------------------------- \n";
    myfile << "AdaptativStep " << fProblemConfig.adaptivityStep
    << " Order k= " << fProblemConfig.k
    << " Order n= " << fProblemConfig.n << "\n";
    myfile << "DOF Total = " << fPostProcMesh.NEquations() << "\n";
    myfile << "Global exact error = " << sqrt(globalerrors[2]) << "\n";
    myfile << "Global estimator = " << sqrt(globalerrors[3]) << "\n";
    myfile << "Global residual error = " << sqrt(globalerrors[4]) << "\n";
    myfile << "Ieff_global = " << Ieff_global << "\n";
    myfile.close();
}

void TPZHybridH1ErrorEstimator::PostProcessing(TPZAnalysis &an) {
    
    TPZMaterial *mat = fPostProcMesh.FindMaterial(*fProblemConfig.materialids.begin());
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
        vecnames.Push("FluxSigmaReconstructed");
        vecnames.Push("FluxReconstructed");
        scalnames.Push("POrder");
        
        int dim = fPostProcMesh.Reference()->Dimension();
        
        std::stringstream out;
        out << fProblemConfig.dir_name << "/" << fProblemConfig.problemname
        << "_k_" << fProblemConfig.k << "_n_"
        << fProblemConfig.n;
        if (fProblemConfig.ndivisions != -1) {
            out << "_Ndiv_" << fProblemConfig.ndivisions;
        }
        if (fProblemConfig.adaptivityStep != -1) {
            out << "_AdaptivityStep_" << fProblemConfig.adaptivityStep;
        }
        out << ".vtk";
        
        an.DefineGraphMesh(dim, scalnames, vecnames, out.str());
        an.PostProcess(2, dim);
    }
    else {
        std::cout << __PRETTY_FUNCTION__ << "\nPost Processing variable not found!\n";
    }
}

// a method for generating the HDiv mesh
TPZCompMesh *TPZHybridH1ErrorEstimator::CreateFluxMesh()
{
    TPZCompMesh *cmeshHdiv = fOriginal->MeshVector()[0]->Clone();

#ifdef PZDEBUG
    string command = "mkdir -p " + fDebugDirName;
    string dirPath = fDebugDirName + "/";
    system(command.c_str());
    {
        std::ofstream outCon(dirPath + "OriginalFluxConnects.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(cmeshHdiv, outCon, {}, false, true);
        std::ofstream outOriginalP(dirPath + "OriginalFlux.txt");
        cmeshHdiv->Print(outOriginalP);
        std::ofstream outGOriginalVTK(dirPath + "gFlux.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmeshHdiv, outGOriginalVTK);
    }
#endif

    VerifyBoundaryFluxConsistency(cmeshHdiv);

    int flux_order = fProblemConfig.k;
    cmeshHdiv->SetDefaultOrder(flux_order);
    cmeshHdiv->SetDimModel(fProblemConfig.gmesh->Dimension());

    cmeshHdiv->SetAllCreateFunctionsHDiv();

    for (auto mat : fOriginal->MaterialVec()) {
        TPZMatLaplacianHybrid *matlaplacian =
                dynamic_cast<TPZMatLaplacianHybrid *>(mat.second);
        if (matlaplacian) {
            TPZHybridH1ErrorEstimateMaterial *EEMat = new TPZHybridH1ErrorEstimateMaterial(*matlaplacian); //TPZHybridH1ErrorEstimateMaterial(mat.first,fProblemConfig.gmesh->Dimension());
            cmeshHdiv->MaterialVec()[mat.first] =  EEMat;

            for (auto bcmat : cmeshHdiv->MaterialVec()) {
                TPZBndCond *bc = dynamic_cast<TPZBndCond *>(bcmat.second);
                if (bc) {
                    bc->SetMaterial(EEMat);
                }
            }
        }
    }

    cmeshHdiv->AutoBuild();
    cmeshHdiv->InitializeBlock();

    TPZCompMeshTools::AdjustFluxPolynomialOrders(cmeshHdiv, fProblemConfig.n); //Increases internal flux order by "hdivmais"

#ifdef PZDEBUG
    {
        std::ofstream outCon(dirPath + "HdivFluxConnects.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(cmeshHdiv, outCon, {}, false, true);
        std::ofstream outOriginalP(dirPath + "HdivFlux.txt");
        cmeshHdiv->Print(outOriginalP);
        std::ofstream outGOriginalVTK(dirPath + "gHdivFlux.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmeshHdiv, outGOriginalVTK);
    }
#endif

#ifdef PZDEBUG
    {
        std::ofstream outOriginalP(dirPath + "HdivCondFlux.txt");
        cmeshHdiv->Print(outOriginalP);
    }
#endif

    TPZCompMesh *pressure = new TPZCompMesh(fProblemConfig.gmesh);
    //if(fisFluxFromMixedProblem){
    //active.Resize(3);
    int dimMesh = fProblemConfig.gmesh->Dimension();

    int potential_order = fProblemConfig.k+fProblemConfig.n;
    pressure->SetDefaultOrder(potential_order);
    pressure->SetDimModel(dimMesh);

    pressure->SetAllCreateFunctionsContinuous(); //H1 functions
    pressure->ApproxSpace().CreateDisconnectedElements(true);

    for(auto matid:fProblemConfig.materialids){
        TPZNullMaterial *material = new TPZNullMaterial(matid); material->SetDimension(dimMesh);
        pressure->InsertMaterialObject(material);
    }
    pressure->AutoBuild();
    pressure->ExpandSolution();

    TPZCompMesh *gspace = new TPZCompMesh(fProblemConfig.gmesh);
    {
        for (auto matid:fProblemConfig.materialids) {
            TPZNullMaterial *nullmat = new TPZNullMaterial(matid);
            nullmat->SetDimension(fProblemConfig.gmesh->Dimension());
            nullmat->SetNStateVariables(1);
            gspace->InsertMaterialObject(nullmat);
        }
        gspace->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
        gspace->SetDefaultOrder(0);//sao espacos de pressao media
        gspace->AutoBuild();
    }

    TPZMultiphysicsCompMesh HdivRecMesh;
    HdivRecMesh.SetReference(fOriginal->Reference());

    TPZManVector<TPZCompMesh *> mesh_vectors(4, 0);
    TPZManVector<int> active(4, 0);

    mesh_vectors[0] = cmeshHdiv;// flux
    mesh_vectors[1] = pressure->Clone();
    mesh_vectors[2] = gspace->Clone(); // g-space
    mesh_vectors[3] = fOriginal->MeshVector()[3]->Clone(); // avg-space

    active[0] = 1;
    active[1] = 1;
    active[2] = 1;

    for (auto mat : mesh_vectors[0]->MaterialVec()) {
        HdivRecMesh.MaterialVec()[mat.first] = mat.second;
    }

    for (auto mat : mesh_vectors[1]->MaterialVec()) {
        int matID = mat.first;
        if(mesh_vectors[0]->MaterialVec().find(matID) ==  mesh_vectors[0]->MaterialVec().end()){
            HdivRecMesh.MaterialVec()[matID] = mat.second;
        }
    }

    HdivRecMesh.BuildMultiphysicsSpace(active,mesh_vectors);

    bool keeponelagrangian = true, keepmatrix = false;
    TPZCompMeshTools::CreatedCondensedElements(&HdivRecMesh, keeponelagrangian, keepmatrix);

#ifdef PZDEBUG
    {
        std::ofstream outMultF(dirPath + "HdivMultMesh.txt");
        HdivRecMesh.Print(outMultF);
    }
#endif

    ComputeElementStiffnesses(HdivRecMesh);

    HdivRecMesh.LoadSolution(HdivRecMesh.Solution());

    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mesh_vectors, &HdivRecMesh);

#ifdef PZDEBUG
    {
        std::ofstream outFCon(dirPath + "HdivFluxAfterLoadSolConnects.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(cmeshHdiv, outFCon, {}, false, true);
        std::ofstream outF(dirPath + "HdivFluxAfterLoadSol.txt");
        cmeshHdiv->Print(outF);
        std::ofstream outMultF(dirPath + "HdivMultFluxAfterLoadSol.txt");
        HdivRecMesh.Print(outMultF);
        std::ofstream outvtk(dirPath + "fluxmesh.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmeshHdiv, outvtk);
    }
#endif
#ifdef PZDEBUG
    TPZAnalysis an(&HdivRecMesh, false);
    TPZStack<std::string> scalnames, vecnames;
    int dim = HdivRecMesh.Reference()->Dimension();

    vecnames.Push("FluxReconstructed");
    vecnames.Push("FluxExact");
    std::stringstream out;
    out << dirPath << "fluxmesh.vtk";

    std::set<int> matids ={1};
    an.DefineGraphMesh(dim,matids, scalnames, vecnames, out.str());
    an.PostProcess(2, dim);
#endif
    return cmeshHdiv->Clone();
}

// a method for creating the pressure mesh
TPZCompMesh *TPZHybridH1ErrorEstimator::CreatePressureMesh()
{
    // For H1 reconstruction, we need to build BC materials
    TPZCompMesh *pressureMesh = fOriginal->MeshVector()[1]->Clone();
    TPZGeoMesh *gmesh = pressureMesh->Reference();

#ifdef PZDEBUG
    {
        std::string dirPath = fDebugDirName + "/";
        std::ofstream outTXT(dirPath + "OriginalPressureMesh.txt");
        std::ofstream outVTK(dirPath + "OriginalPressureMesh.vtk");
        pressureMesh->Print(outTXT);
        TPZVTKGeoMesh::PrintCMeshVTK(pressureMesh, outVTK);
    }
#endif

    gmesh->ResetReference();
    int dim = gmesh->Dimension();

    // Delete compels of dimension dim - 1
    // Case HybridH1 : Delete all MatWrapId computational elements
    for (int64_t el = 0; el < pressureMesh->NElements(); el++) {
        TPZCompEl *cel = pressureMesh->Element(el);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (!gel) DebugStop();
        if (gel->Dimension() == dim - 1) {
            delete cel;
        }
    }
    //RemoveNullCompEl(pressureMesh);

    pressureMesh->ComputeNodElCon();
    pressureMesh->CleanUpUnconnectedNodes();
    pressureMesh->LoadReferences();
    pressureMesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    pressureMesh->ApproxSpace().CreateDisconnectedElements(true);

    //  Insert BC material into the pressure mesh material vector,
    //  Create computational element on BC elements
    AddBC2PressureMesh(pressureMesh);

    // Create internal skeleton elements for pressure reconstruction
    CreateSkeletonElements(pressureMesh);

    // TODO : Make it work after 2D works
    if (dim == 3) {
        CreateEdgeSkeletonMesh(pressureMesh);
    }

#ifdef PZDEBUG
    {
        string dirPath = fDebugDirName + "/";
        std::ofstream outTXT(dirPath + "PostProcPressureMesh.txt");
        std::ofstream outVTK(dirPath +"PostProcPressureMesh.vtk");
        pressureMesh->Print(outTXT);
        TPZVTKGeoMesh::PrintCMeshVTK(pressureMesh, outVTK);
    }
#endif

    return pressureMesh;
}

///  Insert BC material into the pressure mesh material vector,
///  Create computational element on BC elements
void TPZHybridH1ErrorEstimator::AddBC2PressureMesh(TPZCompMesh *pressureMesh){
    TPZCompMesh *mult = fOriginal;
    TPZGeoMesh *gmesh = pressureMesh->Reference();

    // Insert BC materials in pressure reconstruction mesh
    std::set<int> bcMatIDs = fProblemConfig.bcmaterialids;
    for (auto bcID : bcMatIDs) {
        TPZMaterial *mat = mult->FindMaterial(bcID);
        TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
        if (!bc) DebugStop();

        int volumetricMatId = bc->Material()->Id();
        TPZMaterial *pressuremat = pressureMesh->FindMaterial(volumetricMatId);
        if (!pressuremat) DebugStop();

        TPZMaterial *newbc = pressuremat->CreateBC(pressuremat, bc->Id(), bc->Type(), bc->Val1(), bc->Val2());
        if (bc->HasForcingFunction()) {
            newbc->SetForcingFunction(bc->ForcingFunction());
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

    // Create BC elements. We reset references at each step to allow for a discontinuous mesh globally
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
}

void TPZHybridH1ErrorEstimator::CreateSkeletonElements(TPZCompMesh *pressure_mesh) {

    TPZCompMesh* cmesh = fOriginal;
    TPZGeoMesh* gmesh = fOriginal->Reference();
    gmesh->ResetReference();
    cmesh->LoadReferences();
    int dim = gmesh->Dimension();

#ifdef PZDEBUG
    {
        std::string dirPath = fDebugDirName + "/";
        std::ofstream fileVTK(dirPath + "GeoMeshBeforePressureSkeleton.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, fileVTK);
        std::ofstream fileTXT(dirPath + "GeoMeshBeforePressureSkeleton.txt");
        gmesh->Print(fileTXT);
    }
#endif

    // Creation of geometric elements
    const int nel = gmesh->NElements();
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

            //Create Geometric element if there is no boundary neighbour and no skeleton neighbour were created.
            std::set<int> matIDs = fProblemConfig.bcmaterialids;
            matIDs.insert(fPressureSkeletonMatId);
            TPZGeoElSide neighSide = gelside.HasNeighbour(matIDs);
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
    fPostProcMesh.InsertMaterialObject(skeletonMat->NewMaterial());

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

/// create the post processed multiphysics mesh (which is necessarily
/// hybridized)
void TPZHybridH1ErrorEstimator::CreatePostProcessingMesh() {
    
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

    SwitchMaterialObjects();

    TPZManVector<TPZCompMesh *> mesh_vectors(5, 0);
    TPZManVector<int> active(5, 0);

    mesh_vectors[0] = CreateFluxMesh();
    mesh_vectors[1] = CreatePressureMesh();
    mesh_vectors[2] = fOriginal->MeshVector()[0];// flux
    mesh_vectors[3] = fOriginal->MeshVector()[1];// potential
    mesh_vectors[4] = ForceProjectionMesh();

    active[1] = 1;

    fPostProcMesh.BuildMultiphysicsSpace(active, mesh_vectors);

#ifdef PZDEBUG
    {
        std::string dirPath = fDebugDirName + "/";
        /*std::ofstream out(dirPath + "EnrichedFluxBorder.txt");
        mesh_vectors[0]->Print(out);*/
        std::ofstream out2(dirPath + "EnrichedPressure.txt");
        mesh_vectors[1]->Print(out2);
    }
#endif
}

void TPZHybridH1ErrorEstimator::GetPressureMatIDs(std::set<int> &matIDs){
    TPZCompMesh *pressureMesh = PressureMesh();
    int elMatID;
    for(int iel = 0 ; iel < pressureMesh->NElements() ; iel++){
        TPZCompEl *cel = pressureMesh->Element(iel);
        if(!cel) continue;
        elMatID = cel->Reference()->MaterialId();
        if (matIDs.find(elMatID) == matIDs.end()){
            matIDs.insert(elMatID);
        }
    }
}

TPZCompMesh *TPZHybridH1ErrorEstimator::ForceProjectionMesh(){
//L2 projection for forcing f
    TPZCompMesh *forceProj = new TPZCompMesh(fProblemConfig.gmesh);
    int dimMesh = fProblemConfig.gmesh->Dimension();

    int potential_order = fProblemConfig.k+fProblemConfig.n;
    forceProj->SetDefaultOrder(potential_order);
    forceProj->SetDimModel(dimMesh);

    forceProj->SetAllCreateFunctionsContinuous(); //H1 functions
    forceProj->ApproxSpace().CreateDisconnectedElements(true);

    for(auto matid:fProblemConfig.materialids){
        TPZNullMaterial *material = new TPZNullMaterial(matid);
        material->SetDimension(dimMesh);
        forceProj->InsertMaterialObject(material);
        material->SetForcingFunction(fProblemConfig.exact.operator*().ForcingFunction());
    }
    forceProj->AutoBuild();
    forceProj->ExpandSolution();

    for(int iel = 0 ; iel < forceProj->NElements(); iel++) {
        TPZGeoMesh *gmesh = forceProj->Reference();
        int dim = gmesh->Dimension();
        gmesh->ResetReference();
        forceProj->LoadReferences();
        TPZCompEl *cel = forceProj->Element(iel);
        if (!cel) DebugStop();
        TPZGeoEl *gel = cel->Reference();
        if (!gel) DebugStop();
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        int nc = cel->NConnects();
        TPZMaterial *material = cel->Material();
        int polyOrder = cel->Connect(nc - 1).Order();
        int order = material->IntegrationRuleOrder(polyOrder);

        int nshape = intel->NShapeF();
        TPZFNMatrix<20, REAL> L2Mat(nshape, nshape, 0.), L2Rhs(nshape, 1, 0.);
        TPZFNMatrix<220, REAL> phi(nshape, 1, 0.), dphi(dim, nshape);
        REAL weight = 0.;
        TPZAutoPointer<TPZIntPoints> intrule;
        intrule = gel->CreateSideIntegrationRule(gel->NSides() - 1, order);

        TPZManVector<int, 4> intorder(dim, order);
        intrule->SetOrder(intorder);
        int intrulepoints = intrule->NPoints();
        if (intrulepoints > 1000) {
            DebugStop();
        }

        TPZFMatrix<REAL> jac, axe, jacInv;
        REAL detJac;
        TPZManVector<REAL, 4> intpointtemp(3, 0.),x(3, 0.);
        for (int int_ind = 0; int_ind < intrulepoints; ++int_ind) {
            intrule->Point(int_ind, intpointtemp, weight);
            gel->Jacobian(intpointtemp, jac, axe, detJac, jacInv);
            weight *= fabs(detJac);

            TPZMaterialData data;
            intel->Shape(intpointtemp, phi, dphi);
            gel->X(intpointtemp,x);

            if(!material->ForcingFunction()) DebugStop();
            STATE force;
            TPZManVector<STATE> res(3);
            material->ForcingFunction()->Execute(x,res);
            force = res[0];

            //std::cout << "[" << intpointtemp[0] << ", " << intpointtemp[1] << ", " << intpointtemp[2] << "]\n";
            //std::cout<< "f: " << force <<"\n\n";
            //phi.Print(std::cout);



            for (int ishape = 0; ishape < nshape; ishape++) {
                L2Rhs(ishape, 0) += weight * phi(ishape, 0)*force;
                //std::cout<< "ishape/weight/phi/force/result: " << ishape << "/" << weight << "/" << phi(ishape, 0) << "/" << force << "/" << weight*phi(ishape, 0)*force << "\n\n";
            }
            for (int ishape = 0; ishape < nshape; ishape++) {
                for (int jshape = 0; jshape < nshape; jshape++) {
                    L2Mat(ishape, jshape) += weight * phi(ishape, 0) * phi(jshape, 0);
                }
            }
        }

        L2Mat.SolveDirect(L2Rhs, ECholesky);

        // Stores solution in the computational mesh
        int count = 0;
        TPZFMatrix<STATE> &sol = forceProj->Solution();
        for (int ic = 0; ic < nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            int64_t seqnum = c.SequenceNumber();
            int64_t pos = forceProj->Block().Position(seqnum);
            int ndof = c.NShape() * c.NState();
            for (int idf = 0; idf < ndof; idf++) {
                sol(pos + idf, 0) = L2Rhs(count++);
            }
        }
    }

    {
        std::string dirPath = fDebugDirName + "/";
        std::ofstream outCon(dirPath + "forceProj.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(forceProj, outCon, {1}, false, true);
    }

    return forceProj->Clone();
}

/// computing the element stifnesses will "automatically" compute the condensed form of the matrices
void TPZHybridH1ErrorEstimator::ComputeElementStiffnesses() {
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
}

/// computing the element stifnesses will "automatically" compute the condensed form of the matrices
void TPZHybridH1ErrorEstimator::ComputeElementStiffnesses(TPZCompMesh &cmesh) {
    std::cout << "Computing flux stiffness matrix " << std::endl;

    for (auto cel:cmesh.ElementVec()) {
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
void TPZHybridH1ErrorEstimator::IncreaseSideOrders(TPZCompMesh *mesh) {
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
        int order = cel->Connect(nc - 1).Order(); //VO: Setting the flux order equals the last element's order. Does it really increases the interpolation order? Why?
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

void TPZHybridH1ErrorEstimator::IncreasePressureSideOrders(TPZCompMesh *cmesh) {
    
    
    TPZGeoMesh *gmesh = cmesh->Reference();
    
    gmesh->ResetReference();
    cmesh->LoadReferences();

#ifdef PZDEBUG
    std::string dirPath = fDebugDirName + "/";
    std::set<int> matIDs = fProblemConfig.materialids;
    matIDs.insert(fProblemConfig.bcmaterialids.begin(),fProblemConfig.bcmaterialids.end());
    matIDs.insert(fPressureSkeletonMatId);
    {
        std::ofstream outCon(dirPath + "PressureConnectsB4IncreaseSideOrder.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(cmesh, outCon, matIDs, false, true);
    }
#endif
    
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
            if (!celside) {
                continue;//DebugStop();/// para nao incremenentar ordem na condicao de contorno
            }
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

#ifdef PZDEBUG
    {
        std::ofstream outCon(dirPath + "PressureConnectsAFTERIncreaseSideOrder.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(cmesh, outCon, matIDs, false, true);
    }
#endif

}

/// Find free matID number
void TPZHybridH1ErrorEstimator::FindFreeMatID(int &matID){
    TPZGeoMesh* gmesh = fOriginal->Reference();

    int maxMatId = std::numeric_limits<int>::min();
    const int nel = gmesh->NElements();

    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->Element(iel);
        if(gel) maxMatId = std::max(maxMatId, gel->MaterialId());
    }

    if (maxMatId == std::numeric_limits<int>::min()) maxMatId = 0;

    matID = maxMatId + 1;
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
void TPZHybridH1ErrorEstimator::CreateEdgeSkeletonMesh(TPZCompMesh *pressuremesh) {
    
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

void TPZHybridH1ErrorEstimator::RestrainSkeletonSides(TPZCompMesh *pressure_mesh) {

    TPZGeoMesh *gmesh = pressure_mesh->Reference();
    gmesh->ResetReference();
    pressure_mesh->LoadReferences();

#ifdef PZDEBUG
    {
        std::string dirPath = fDebugDirName + '/';
        std::ofstream out(dirPath + "MeshBeforeRestrainSkeleton.txt");
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
        //std::cout << "Restriction @ [" << xcenter << "]:"
        //          << "  Small El: " << small.Element()->Index() << ", Side: " << small.Side()
        //          << "  Large El: " << largerNeigh.Element()->Index() << ", Side: " << largerNeigh.Side() << "\n";
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
            //std::cout << "SubRestriction @ [" << xcenter << "]:"
            //          << "  Small El: " << small.Element()->Index() << ", Side: " << subsmall.Side()
             //         << "  Large El: " << largerNeigh.Element()->Index() << ", Side: " << subLargerNeigh.Side() << "\n";
            smallIntel->RestrainSide(subsmall.Side(), largeIntel, subLargerNeigh.Side());
        }
    }

        pressure_mesh->CleanUpUnconnectedNodes();
#ifdef PZDEBUG
    {
        std::string dirPath = fDebugDirName + '/';
        std::ofstream out(dirPath + "MeshAfterRestrainSkeleton.txt");
        pressure_mesh->Print(out);
    }
#endif
}

/// restrain the edge elements that have larger elements as neighbours
void TPZHybridH1ErrorEstimator::RestrainSmallEdges(TPZCompMesh *pressuremesh) {
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
void TPZHybridH1ErrorEstimator::AdjustNeighbourPolynomialOrders(TPZCompMesh *pressure_mesh) {
    //    TPZCompMesh *pressure_mesh = fPostProcMesh.MeshVector()[1];
    TPZGeoMesh *gmesh = pressure_mesh->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    // load the elements of lower dimension than dim
    int64_t nel = pressure_mesh->NElements();
    std::map<std::pair<int64_t, int>, int> polynomialorders;
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
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
            TPZCompEl *cel = pressure_mesh->Element(el);
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
        TPZCompEl *cel = pressure_mesh->Element(index);
        if (!cel) DebugStop();
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        intel->SetSideOrder(side, porder);
    }
}

/// return a pointer to the pressure mesh
TPZCompMesh *TPZHybridH1ErrorEstimator::PressureMesh() {
    return fPostProcMesh.MeshVector()[1];
}

TPZCompMesh *TPZHybridH1ErrorEstimator::FluxMesh() {
    return fPostProcMesh.MeshVector()[0];
}

/// compute the average pressures of the hybridized form of the H(div) mesh
void TPZHybridH1ErrorEstimator::ComputeAverageFacePressures() {
    DebugStop();
    /*TPZCompMesh *pressure = fOriginal->MeshVector()[1];
    TPZCompMesh *pressure_mesh = fPostProcMesh.MeshVector()[1];
    int fInterfaceMatid = fHybridizer.fLagrangeInterface;
    TPZGeoMesh *gmesh = pressure->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    pressure->LoadReferences();
    int64_t nel = pressure_mesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
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
            int64_t pos = pressure_mesh->Block().Position(seqnum);
            int ndof = c.NShape() * c.NState();
            for (int idf = 0; idf < ndof; idf++) {
                pressure_mesh->Solution()(pos + idf, 0) = L2Rhs(count++);
            }
        }
    }
    TPZManVector<TPZCompMesh *, 2> meshvec(2);
    meshvec[0] = fPostProcMesh.MeshVector()[0];
    meshvec[1] = fPostProcMesh.MeshVector()[1];
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, &fPostProcMesh);
*/}

/// compute the average pressures of across edges of the H(div) mesh
void TPZHybridH1ErrorEstimator::ComputeAveragePressures(int target_dim) {

    TPZCompMesh *pressure_mesh = PressureMesh();
    TPZGeoMesh *gmesh = pressure_mesh->Reference();

    //    std::ofstream out("PressureToAverage.txt");
    //    pressure_mesh->Print(out);

    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    int64_t nel = pressure_mesh->NElements();
    // load the pressure elements of dimension target_dim+1
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
        if(!cel) {
            continue;
        }

        TPZGeoEl *gel = cel->Reference();

        if (!gel || gel->Dimension() != target_dim) {
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

        {
            std::string dirName = "DebuggingConsistency";
            std::string command = "mkdir -p " + dirName;
            system(command.c_str());
        }

        ComputeAverage(pressure_mesh, el);
    }

    {
        std::string dirName = "DebuggingConsistency";
        std::string dirPath = dirName + "/";
        std::ofstream outCon(dirPath + "AverageB4LoadSolution.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(pressure_mesh, outCon, {1,2,3,fPressureSkeletonMatId}, false, true);
    }

    // Loads solution into the connects of the smaller skeletons
    pressure_mesh->LoadSolution(pressure_mesh->Solution());

    {
        std::string dirPath = "DebuggingConsistency/";
        std::ofstream outCon(dirPath + "AverageAfterLoadSolution.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(pressure_mesh, outCon, {1,2,3,fPressureSkeletonMatId}, false, true);
    }

    // apply the restraints to the edge connects
    if (target_dim == dim - 2) {
        pressure_mesh->LoadSolution(pressure_mesh->Solution());
        TransferEdgeSolution();
    }
}
//compute de L2 projection of Dirichlet boundary condition for Hdi-H1 reconstruction
void TPZHybridH1ErrorEstimator::ComputeBoundaryL2Projection(TPZCompMesh *pressuremesh, int target_dim){
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
        TPZFMatrix<STATE> &mesh_sol = pressuremesh->Solution();
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
}

void TPZHybridH1ErrorEstimator::BoundaryPressureProjection(TPZCompMesh *pressuremesh, int target_dim){
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
            
            
            
            
            
            if(ekbc.HasDependency()) DebugStop();
            ekbc.fMat.SolveDirect(efbc.fMat, ECholesky);
            //   efbc.fMat.Print("Solution",std::cout);
            
            
            TPZCompEl *celpressure = mphys->Element(1);
            int nc = celpressure->NConnects();
            int count = 0;
            TPZFMatrix<STATE> &mesh_sol = pressuremesh->Solution();
            for (int ic = 0; ic < nc; ic++) {
                TPZConnect &c = celpressure->Connect(ic);
                int64_t seqnum = c.SequenceNumber();
                int64_t pos = pressuremesh->Block().Position(seqnum);
                int ndof = c.NShape() * c.NState();
                for (int idf = 0; idf < ndof; idf++) {
                    mesh_sol(pos + idf, 0) = efbc.fMat(count++);
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



void TPZHybridH1ErrorEstimator::NewComputeBoundaryL2Projection(
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
        TPZFMatrix<STATE> &mesh_sol = pressuremesh->Solution();
        for (int ic = 0; ic < nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            int64_t seqnum = c.SequenceNumber();
            int64_t pos = pressuremesh->Block().Position(seqnum);
            int ndof = c.NShape() * c.NState();
            for (int idf = 0; idf < ndof; idf++) {
                mesh_sol(pos + idf, 0) = efbc(count++);
            }
        }
    }
    
    {
        std::ofstream out("PressureAfterL2Projection_2.txt");
        pressuremesh->Print(out);
    }
}

void TPZHybridH1ErrorEstimator::AdditionalAverageWeights(TPZGeoEl *large, TPZGeoEl *small, REAL &large_weight, REAL &small_weight, REAL sum){
    REAL diff_weight = large_weight+small_weight - 1.;
    if(!IsZero(diff_weight)) DebugStop();
    switch(fAverageMode){
        case 0:
            return;
        case 1:
            large_weight = 0.;
            small_weight = 1.;
            break;
        case 2: // w ~ 1/h^p
            REAL large_size, small_size;
            int large_order,small_order;

            //Compute h
            large_size = large->CharacteristicSize();
            small_size = small->CharacteristicSize();

            //Compute k
            large_order = large->Reference()->Connect(large->NSides()-1).Order();
            small_order = small->Reference()->Connect(small->NSides()-1).Order();

            //Compute 1/h^(k)
            REAL large_precision = 1.;
            REAL small_precision = 1.;
            for(int pindex = 0; pindex < large_order; pindex++){
                large_precision *= large_size;
            }
            for(int pindex = 0; pindex < small_order; pindex++){
                small_precision *= small_size;
            }
            large_precision = 1./large_precision;
            small_precision = 1./small_precision;

            REAL new_large_weight = large_weight*large_precision;
            REAL new_small_weight = small_weight*small_precision;
            REAL new_sum = new_large_weight+new_small_weight;

            //std::cout << "large0: " << large_weight << "; large_h: " << large_size << "; large_p: " <<  large_order << "; large_precision: " << large_precision << "; large1: " << new_large_weight/new_sum <<"\n";
            //std::cout << "small0: " << small_weight << "; small_h: " << small_size << "; small_p: " <<  small_order << "; small_precision: " << small_precision << "; small1: " << new_small_weight/new_sum <<"\n";

            large_weight = new_large_weight/new_sum;
            small_weight = new_small_weight/new_sum;

            if(!IsZero(large_weight + small_weight - 1.)) DebugStop();
    }
}

// compute the average of an element iel in the pressure mesh looking at its neighbours
void TPZHybridH1ErrorEstimator::ComputeAverage(TPZCompMesh *pressuremesh, int64_t iel)
{
    TPZGeoMesh *gmesh = pressuremesh->Reference();
    int dim = gmesh->Dimension();
    gmesh->ResetReference();
    pressuremesh->LoadReferences();
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
    {
        std::ofstream out("PressureB4Debug.txt");
        pressuremesh->Print(out);
    }
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
        if (!partition.HigherLevelNeighbours(smallerSkelSides, fPressureSkeletonMatId)) {
            DebugStop();
        }
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

        if(!noHangingSide){
            AdditionalAverageWeights(largeSkeletonSide.Element(), integrationGeoElSide.Element(), left_weight, right_weight,sum_weights);
        }

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

#ifdef LOG4CXX
            if(logger->isDebugEnabled()) {
                std::stringstream ss;
                phi.Print(ss);
                ss << detjac << '\n';
                LOGPZ_DEBUG(logger, ss.str())
            }
#endif
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

#ifdef LOG4CXX
    if(logger->isDebugEnabled()) {
        std::stringstream ss;
        L2Rhs.Print("Rhs =", ss, EMathematicaInput);
        L2Mat.Print("Stiffness =", ss, EMathematicaInput);
        LOGPZ_DEBUG(logger, ss.str())
    }
#endif

    L2Mat.SolveDirect(L2Rhs, ECholesky);
    // Stores solution in the computational mesh

    TPZFMatrix<STATE> &sol = pressuremesh->Solution();
    int count = 0;
    for (int ic = 0; ic < nc; ic++) {
        TPZConnect &c = cel->Connect(ic);
        int64_t seqnum = c.SequenceNumber();
        int64_t pos = pressuremesh->Block().Position(seqnum);
        int ndof = c.NShape() * c.NState();
        for (int idf = 0; idf < ndof; idf++) {
            sol(pos + idf, 0) = L2Rhs(count++);
        }
    }
}


/// transfer the solution of the edge functions to the face functions
void TPZHybridH1ErrorEstimator::TransferEdgeSolution() {
    // copy the solution associated with one-d edge connect to the corresponding side connect of the face mesh
    DebugStop();
    /*TPZCompMesh *pressure_mesh = PressureMesh();
    TPZGeoMesh *gmesh = pressure_mesh->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    if (dim != 3) {
        std::cout << __PRETTY_FUNCTION__ << " should not be called for mesh dimension " << dim << std::endl;
        return;
    }
    int lagrangematid = fHybridizer.fLagrangeInterface;
    TPZMaterial *mat = pressure_mesh->FindMaterial(lagrangematid);
    if (!mat) DebugStop();
    int nstate = mat->NStateVariables();
    int64_t nel = pressure_mesh->NElements();
    // load the pressure elements of dimension 1 and 2
    pressure_mesh->Reference()->ResetReference();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() > 2) continue;
        gel->SetReference(cel);
    }
    
    // loop over the edge elements
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
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
            int nshape_edge = pressure_mesh->Block().Size(edge_seqnum);
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
                int nshape_neigh = pressure_mesh->Block().Size(neighblock);
                if (nshape_edge != nshape_neigh) DebugStop();
                for (int i = 0; i < nshape_neigh; i++) {
                    pressure_mesh->Block()(neighblock, 0, i, 0) = pressure_mesh->Block()(edge_seqnum, 0, i, 0);
                }
            }
        }
    }*/
}


/// set the cornernode values equal to the averages
void TPZHybridH1ErrorEstimator::ComputeNodalAverages() {
    TPZCompMesh *pressure_mesh = PressureMesh();
    TPZGeoMesh *gmesh = pressure_mesh->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    pressure_mesh->LoadReferences();
    TPZMaterial *mat = pressure_mesh->FindMaterial(fPressureSkeletonMatId);
    if (!mat) DebugStop();
    int nstate = mat->NStateVariables();
    int64_t nel = pressure_mesh->NElements();
    // load the pressure elements of dimension dim-1
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
    }
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
        
        //percorre cada no do elemento de interface
        int ncorners = gel->NCornerNodes();
        for (int side = 0; side < ncorners; side++) {
            TPZGeoElSide gelside(gel, side);
            TPZCompElSide celside(intel,side);

            int64_t conindex = intel->ConnectIndex(side);
            TPZConnect &c = pressure_mesh->ConnectVec()[conindex];
            if (!c.HasDependency()) {
                if(LiesOnHangingSide(celside)) continue;
                ComputeNodalAverage(celside);
            }
        }
    }

    // Load the solution into connects with dependency
    pressure_mesh->LoadSolution(pressure_mesh->Solution());

    // Impose the solution on nodes which lies on hanging sides
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

        //percorre cada no do elemento de interface
        int ncorners = gel->NCornerNodes();
        for (int side = 0; side < ncorners; side++) {
            TPZGeoElSide gelside(gel, side);
            TPZCompElSide celside(intel,side);

            int64_t conindex = intel->ConnectIndex(side);
            TPZConnect &c = pressure_mesh->ConnectVec()[conindex];
            if (c.HasDependency()) continue;
            if(LiesOnHangingSide(celside)){
                ImposeHangingNodeSolution(celside);
            }
        }
    }
}

// If all nodal (0-dimensional) connects linked to this node have dependency, it lies* on a hanging side and
// instead of calculating the average we should impose the solution of the other connects in it.
// *We are not sure if this logic will work for every case/ref. pattern, specially in 3D. So I'm leaving this
// disclaimer. (Gustavo Batistela, 14/10/2020)
bool TPZHybridH1ErrorEstimator::LiesOnHangingSide(TPZCompElSide &node_celside){

    TPZCompMesh *pressure_mesh = PressureMesh();
    int dim = pressure_mesh->Dimension();
    TPZGeoElSide node_gelside(node_celside.Reference());

    // celstack will contain all zero dimensional sides connected to the side
    TPZStack<TPZCompElSide> celstack;
    int onlyinterpolated = 1;
    int removeduplicates = 0;

    node_gelside.ConnectedCompElementList(celstack, onlyinterpolated, removeduplicates);

    bool allNodeConnectsHaveDependency = true;
    for (int elc = 0; elc < celstack.size(); elc++) {
        TPZCompElSide celside = celstack[elc];
        TPZGeoElSide gelside = celside.Reference();
        if (gelside.Element()->Dimension() != dim - 1) {
            continue;
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(celside.Element());
        if (!intel) DebugStop();

        int64_t conindex = intel->ConnectIndex(celside.Side());
        TPZConnect &c = pressure_mesh->ConnectVec()[conindex];
        if (gelside.Dimension() == 0 && !c.HasDependency()) allNodeConnectsHaveDependency = false;
    }
    return allNodeConnectsHaveDependency;
}

void TPZHybridH1ErrorEstimator::ImposeHangingNodeSolution(TPZCompElSide &node_celside){

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

    // Impose solution of the first node connect in the celstack we can find

    celstack.Push(node_celside);

    TPZBlock &block =  pressure_mesh->Block();
    TPZFMatrix<STATE> &sol = pressure_mesh->Solution();
    for (int elc = 0; elc < celstack.size(); elc++) {
        TPZCompElSide celside = celstack[elc];
        if (celside.Reference().Dimension() != 0) continue;

        // Get solution of the neighbour
        TPZInterpolatedElement *neigh_intel = dynamic_cast<TPZInterpolatedElement *>(celside.Element());
        if (!neigh_intel) DebugStop();

        int64_t neigh_conindex = neigh_intel->ConnectIndex(celside.Side());
        TPZConnect &neigh_c = pressure_mesh->ConnectVec()[neigh_conindex];

        int64_t neigh_seqnum = neigh_c.SequenceNumber();
        if (neigh_c.NState() != nstate || neigh_c.NShape() != 1) DebugStop();
        TPZManVector<STATE, 3> neigh_sol(nstate, 0.);
        for (int istate = 0; istate < nstate; istate++) {
            neigh_sol[istate] = sol.at(block.at(neigh_seqnum, 0, istate, 0));
        }

        // Set solution to given connect
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(node_celside.Element());
        if (!intel) DebugStop();

        int64_t conindex = intel->ConnectIndex(side);
        TPZConnect &c = pressure_mesh->ConnectVec()[conindex];

        int64_t seqnum = c.SequenceNumber();
        if (c.NState() != nstate || c.NShape() != 1) DebugStop();
        for (int istate = 0; istate < nstate; istate++) {
            sol.at(block.at(seqnum, 0, istate, 0)) = neigh_sol[istate];
        }

        // TODO remove this
        {
            TPZManVector<REAL, 3> xicenter(gel->Dimension(), 0.);
            TPZManVector<REAL, 3> x_side(3, 0.);
            gel->CenterPoint(side, xicenter);
            gel->X(xicenter, x_side);
            //std::cout << __PRETTY_FUNCTION__ << "\nAll " << celstack.size() << " Connects Have Dependency!\nGel "
             //         << gel->Index() << ", Side " << side << " @ [" << x_side[0] << ", " << x_side[1] << ", "
             //         << x_side[2] << "]\n";
        }
        break;
    }
    return;
}

/// compute the nodal average of all elements that share a point
void TPZHybridH1ErrorEstimator::ComputeNodalAverage(TPZCompElSide &node_celside)
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

    // Impose solution of the first node connect in the celstack we can find

    celstack.Push(node_celside);

    // This map stores the connects, the weight associated with the element
    // and the solution of that connect. The weight of Dirichlet condition is
    // higher and will be used later to impose the value of the BC in the
    // connects when needed
    std::map<int64_t, std::pair<REAL, TPZVec<STATE>>> connects;
    TPZBlock &block =  pressure_mesh->Block();
    TPZFMatrix<STATE> &solMatrix = pressure_mesh->Solution();
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
            DebugStop();
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
            solMatrix.at(block.at(seqnum, 0, istate, 0))= averageSol[istate];
        }
    }
}


/// clone the meshes into the post processing mesh
void TPZHybridH1ErrorEstimator::CloneMeshVec() {
    
    for (int i = 0; i < fOriginal->MeshVector().size(); i++) {
        fPostProcMesh.MeshVector()[i] = fOriginal->MeshVector()[i]->Clone();
    }
    
}

/// compute the effectivity indices of the pressure error and flux error and store in the element solution
void TPZHybridH1ErrorEstimator::ComputeEffectivityIndices(TPZSubCompMesh *subcmesh) {
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

            REAL oscillatorytherm = 0;
            if (i == 2) {
                oscillatorytherm = elsol(el, i + 2);
                oscillatorytherm *= (hk / M_PI);
            }

            if (abs(ErrorEstimate) < tol) {
                elsol(el, ncols + i / 2) = 1.;

            } else {
                REAL EfIndex = (ErrorEstimate + oscillatorytherm) / ErrorExact;
                elsol(el, ncols + i / 2) = EfIndex;
            }
        }
    }
}

/// TODO: Suport flux' effectivity indices
/// compute the effectivity indices of the pressure error and flux error and store in the element solution
void TPZHybridH1ErrorEstimator::ComputeEffectivityIndices() {
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
    
    //std::ostream &out;
    //    cmesh->ElementSolution().Print("ElSolution",std::cout);
    
    
    TPZFMatrix<REAL> dataIeff(nrows,1);
    
    //    REAL oscilatorytherm = 0;
    int dim = cmesh->Dimension();
    elsol.Resize(nrows, ncols+2);
    REAL oscilatorytherm = 0;
    REAL fluxestimator = 0;
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

                //std::cout << "linha = " << el << " col = " << 4 + i / 2 << " neinEl " << neighindex << std::endl;
                
                if(neighindex > nrows){
                    std::cout<<" neighindex= "<< neighindex<<" nrows "<<nrows<<"\n";
                    DebugStop();
                }
                
                REAL NeighbourErrorEstimate = elsol(neighindex, i + 1);
                REAL NeighbourErrorExact = elsol(neighindex, i);
                REAL ErrorEstimate = elsol(el, i + 1);
                REAL ErrorExact = elsol(el, i);
                REAL sumErrorExact = sqrt(NeighbourErrorExact*NeighbourErrorExact+ErrorExact*ErrorExact);
                REAL sumErrorEstimate = sqrt(NeighbourErrorEstimate*NeighbourErrorEstimate+ErrorEstimate*ErrorEstimate);
                elsol(neighindex,i+1) = 0.;
                elsol(neighindex,i) = 0.;
                elsol(el,i) = sumErrorExact;
                elsol(el,i+1) = sumErrorEstimate;
            }
        }
    }
    REAL globalResidual =0., globalProjResidual =0.;;
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

        REAL elementResidual,elementProjResidual;
        elementResidual = (hk/M_PI)*elsol(el, 4);
        elementProjResidual = (hk/M_PI)*elsol(el, 6);
        globalResidual += elementResidual*elementResidual;
        globalProjResidual += elementProjResidual*elementProjResidual;

        for (int i = 0; i < 3; i += 2) {
            
            //  std::cout<<"linha = "<<el<< "col = "<<4 + i / 2<<std::endl;
            
            REAL tol = 1.e-10;
            oscilatorytherm = 0.;
            fluxestimator = 0.;
            REAL ErrorEstimate = elsol(el, i + 1);
            REAL ErrorExact = elsol(el, i);
           // cmesh->ElementSolution().Print(std::cout);

            TPZGeoEl *gel = cel->Reference();
            
            REAL hk = gel->CharacteristicSize();
            
            if(i==2){
                oscilatorytherm = elsol(el, i + 2);
                oscilatorytherm *= (hk/M_PI);
                fluxestimator = elsol(el, i + 3);
            }
            
            
            if (abs(ErrorEstimate) < tol) {
                elsol(el, ncols + i / 2) = 1.;
                dataIeff(el,0)=1.;
                
            }
            else {
                REAL EfIndex = sqrt(ErrorEstimate*ErrorEstimate +(oscilatorytherm+fluxestimator)*(oscilatorytherm+fluxestimator))/ErrorExact;
                dataIeff(el,0)= EfIndex;
                
                elsol(el, ncols + i / 2) = EfIndex;
            }
        }
    }

    globalResidual = sqrt(globalResidual);
    globalProjResidual = sqrt(globalProjResidual);
    std::cout << "\n\n";
    std::cout << "Residual = " << globalResidual << "\n";
    std::cout << "ProjResidual = " << globalProjResidual << "\n";
    std::cout << "\n\n";

    
    //  cmesh->ElementSolution().Print("ElSolution",std::cout);
    //    ofstream out("IeffPerElement3DEx.nb");
    //    dataIeff.Print("Ieff = ",out,EMathematicaInput);


}

/// returns true if the material associated with the element is a boundary condition
/// and if the boundary condition is dirichlet type
bool TPZHybridH1ErrorEstimator::IsDirichletCondition(TPZGeoElSide gelside) {
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
void TPZHybridH1ErrorEstimator::GetDirichletValue(TPZGeoElSide gelside, TPZVec<STATE> &vals) {
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

void TPZHybridH1ErrorEstimator::FluxReconstruction() {
    DebugStop();
}

void TPZHybridH1ErrorEstimator::CreateReconstructionSpaces() {
    if (fPostProcMesh.MeshVector().size()) {
        DebugStop();
    }

#ifdef PZDEBUG
    std::string command = "mkdir -p " + fDebugDirName;
    std::string dirPath = fDebugDirName + "/";
    system(command.c_str());
    {
        std::ofstream outCon(dirPath + "OriginalPressureConnects.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(fOriginal->MeshVector()[1], outCon, {1,2,3}, false, true);
        std::ofstream outMultCon(dirPath + "OriginalConnects.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(fOriginal, outMultCon, {-9,-8,-6,-5,1,2,3}, false, true);
        std::ofstream outOriginalP(dirPath + "OriginalPressure.txt");
        fOriginal->MeshVector()[1]->Print(outOriginalP);
        std::ofstream outGOriginalVTK(dirPath + "gOriginal.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fOriginal->Reference(), outGOriginalVTK);
        std::ofstream outGOriginalTXT(dirPath + "gOriginal.txt");
        fOriginal->Reference()->Print(outGOriginalTXT);
    }
#endif

    CreatePostProcessingMesh();

#ifdef PZDEBUG
    {
        std::ofstream out(dirPath + "MultiphysicsMeshInPotentialReconstruction.txt");
        fPostProcMesh.Print(out);
        std::ofstream outOrig(dirPath + "PressureConnectsBeforeReconstruction.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(fPostProcMesh.MeshVector()[1], outOrig, {1,2,3,fPressureSkeletonMatId}, false, true);
        std::ofstream outGReconstVTK(dirPath + "gReconstruct.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fPostProcMesh.Reference(), outGReconstVTK);
        std::ofstream outGReconstTXT(dirPath + "gReconstruct.txt");
        fOriginal->Reference()->Print(outGReconstTXT);
    }
#endif

    PotentialReconstruction();
};

void TPZHybridH1ErrorEstimator::PotentialReconstruction() {

    if (!fPostProcMesh.MeshVector().size()) {
        DebugStop();
    }

    PrepareElementsForH1Reconstruction();

    // Compute continuos pressure on the skeleton;
    MakeSkeletonContinuous();

    VerifySkeletonContinuity(PressureMesh());

#ifdef PZDEBUG
    {
        std::ofstream out("MeshWithSmoothPressure.txt");
        fPostProcMesh.Print(out);
        std::ofstream out2("PressureMeshSmooth.txt");
        fPostProcMesh.MeshVector()[1]->Print(out2);
    }
#endif

    ComputeElementStiffnesses();
    
#ifdef PZDEBUG
    {
        std::string command = "mkdir -p " + fDebugDirName +  "/DebuggingLoadSol";
        system(command.c_str());
        std::string dirPath = fDebugDirName + "/DebuggingLoadSol/";
        std::ofstream out(dirPath + "MeshBeforeLoadSol.txt");
        fPostProcMesh.Print(out);
        ofstream out2(dirPath + "SolBeforeLoadSolution.nb");
        fPostProcMesh.Solution().Print("SolBeforeLoadSolution=",out2,EMathematicaInput);
        ofstream out3(dirPath + "PressureWithAverage.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(fPostProcMesh.MeshVector()[1], out3, {1,2,3}, false, true);
        ofstream out4(dirPath + "PressureWithAverage.txt");
    }
#endif
    
    fPostProcMesh.LoadSolution(fPostProcMesh.Solution());
    
#ifdef PZDEBUG
    {
        std::string dirPath = fDebugDirName + "/DebuggingLoadSol/";
        std::ofstream out(dirPath + "MeshAfterLoadSol.txt");
        fPostProcMesh.Print(out);
        std::ofstream outP(dirPath + "PotentialAfterLoadSol.txt");
        fPostProcMesh.MeshVector()[1]->Print(outP);
        //fPostProcMesh.Solution().Print("SolAfterLoadSolution");
        ofstream out2(dirPath + "SolAfterLoadSolution.nb");
        fPostProcMesh.Solution().Print("SolAfterLoadSolution=",out2,EMathematicaInput);
        ofstream out3(dirPath + "PressureAfterLoadSolution.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(fPostProcMesh.MeshVector()[1], out3, {1,2,3}, false, true);
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
        
#ifdef PZDEBUG
        {
            string dirPath = fDebugDirName + "/";
            std::ofstream outCon(dirPath + "PressureConnectsAfterLoadSolution.txt");
            TPZCompMeshTools::PrintConnectInfoByGeoElement(fPostProcMesh.MeshVector()[1], outCon, {1,2,3}, false, true);
            //fPostProcMesh.MeshVector()[0]->Print(out2);
            
        }
        VerifySolutionConsistency(PressureMesh());
#endif
    }
}

void TPZHybridH1ErrorEstimator::MakeSkeletonContinuous(){

    ComputePressureWeights();

    // L2 Projection
    TPZCompMesh *pressuremesh = PressureMesh();
    int target_dim = 1;
    // TODO ver se fica igual para dimensao maior
    ComputeBoundaryL2Projection(pressuremesh, target_dim);
    //BoundaryPressureProjection(pressuremesh, target_dim);

    {
        std::string dirPath = fDebugDirName + "/";
        std::ofstream out(dirPath + "PressureAverageMesh.txt");
        fPostProcMesh.MeshVector()[1]->Print(out);
        bool reconstructed = true;
        PlotLagrangeMultiplier(dirPath + "BeforeAverage",reconstructed);
    }

    {
        std::string dirPath = "DebuggingConsistency/";
        std::ofstream outCon(dirPath + "AverageB4Average.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(fPostProcMesh.MeshVector()[1], outCon, {fPressureSkeletonMatId}, false, true);
    }

    // Calculates average pressure on interface edges and vertices
    int dim = fPostProcMesh.Dimension();
    ComputeAveragePressures(dim - 1);
    // in three dimensions make the one-d polynoms compatible
    if (dim == 3) {
        ComputeAveragePressures(1);
    }

    //VerifyAverage(dim-1);


    {
        std::string dirPath = fDebugDirName + "/";
        std::ofstream out(dirPath+"PressureAverageMesh.txt");
        fPostProcMesh.MeshVector()[1]->Print(out);
        PlotLagrangeMultiplier(dirPath+"BeforeNodalAverage");
    }

    {
        std::string dirPath = "DebuggingConsistency/";
        std::ofstream outCon(dirPath + "AverageB4NodalAverage.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(fPostProcMesh.MeshVector()[1], outCon, {fPressureSkeletonMatId}, false, true);
    }

    ComputeNodalAverages();

    {
        std::string dirPath = "DebuggingConsistency/";
        std::ofstream outCon(dirPath + "AverageAfterNodalAverage.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(fPostProcMesh.MeshVector()[1], outCon, {fPressureSkeletonMatId}, false, true);
    }

    {
        std::string dirPath = fDebugDirName + "/";
        std::ofstream out(dirPath + "PressureNodalMesh.txt");
        fPostProcMesh.MeshVector()[1]->Print(out);
        PlotLagrangeMultiplier("AfterNodalAverage");
    }


    {
        std::string dirPath = "DebuggingConsistency/";
        std::ofstream outCon(dirPath + "AverageAfterLoadSolAverage.txt");
        TPZCompMeshTools::PrintConnectInfoByGeoElement(fPostProcMesh.MeshVector()[1], outCon, {1,2,3,fPressureSkeletonMatId}, false, true);
    }

    CopySolutionFromSkeleton();

    // transfer the continuous pressures to the multiphysics space
    {
        TPZManVector<TPZCompMesh *, 2> meshvec(2);
        meshvec[0] = fPostProcMesh.MeshVector()[0];//flux
        meshvec[1] = fPostProcMesh.MeshVector()[1];//pressure
        fPostProcMesh.CleanUpUnconnectedNodes();


        std::string command = "mkdir -p " + fDebugDirName +  "/DebuggingTransfer";
        system(command.c_str());

        {
            std::string dirPath = fDebugDirName +  "/DebuggingTransfer/";
            std::ofstream out(dirPath + "PressureBeforeTransferFromMeshes.txt");
            TPZCompMeshTools::PrintConnectInfoByGeoElement(fPostProcMesh.MeshVector()[1], out);

            std::ofstream outMultiphysics(dirPath + "MultiphysicsBeforeTransferFromMeshes.txt");
            std::set<int> matIDs;
            GetPressureMatIDs(matIDs);
            TPZCompMeshTools::PrintConnectInfoByGeoElement(&fPostProcMesh, outMultiphysics,matIDs);
        }
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, &fPostProcMesh);
        {
            std::string dirPath = fDebugDirName +  "/DebuggingTransfer/";
            std::ofstream out(dirPath + "PressureAfterTransferFromMeshes.txt");
            TPZCompMeshTools::PrintConnectInfoByGeoElement(fPostProcMesh.MeshVector()[1], out);
            std::ofstream outMultiphysics(dirPath + "MultiphysicsAfterTransferFromMeshes.txt");
            std::set<int> matIDs;
            GetPressureMatIDs(matIDs);
            TPZCompMeshTools::PrintConnectInfoByGeoElement(&fPostProcMesh, outMultiphysics,matIDs);
        }
    }
}

void TPZHybridH1ErrorEstimator::PlotLagrangeMultiplier(const std::string &filename, bool reconstructed) {
    
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
        std::set<int> matids ={fPressureSkeletonMatId};
        an.DefineGraphMesh(dim,matids, scalnames, vecnames, plotname);
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
void TPZHybridH1ErrorEstimator::IdentifyPeripheralMaterialIds() {
    DebugStop();
    /*int dim = fOriginal->Dimension();
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
    
    fPressureSkeletonMatId = fHybridizer.fLagrangeInterface;*/
}

/// switch material object from mixed poisson to TPZMixedHdivErrorEstimate
void TPZHybridH1ErrorEstimator:: SwitchMaterialObjects() {
    for (auto mat : fPostProcMesh.MaterialVec()) {
        TPZMatLaplacianHybrid *matlaplacian=
                dynamic_cast<TPZMatLaplacianHybrid *>(mat.second);
        if (matlaplacian) {
            TPZHybridH1ErrorEstimateMaterial *newmat = new TPZHybridH1ErrorEstimateMaterial(*matlaplacian);

            for (auto bcmat : fPostProcMesh.MaterialVec()) {
                TPZBndCond *bc = dynamic_cast<TPZBndCond *>(bcmat.second);
                if (bc) {
                    bc->SetMaterial(newmat);
                }
            }
            fPostProcMesh.MaterialVec()[newmat->Id()] = newmat;
            delete matlaplacian;
        }
    }
}

void TPZHybridH1ErrorEstimator::VerifyBoundaryFluxConsistency(TPZCompMesh* fluxmesh){
    TPZGeoMesh *gmesh = fluxmesh->Reference();
    gmesh->ResetReference();
    fluxmesh->LoadReferences();
    int meshdim = gmesh->Dimension();
    int nel = fluxmesh->NElements(), nEffectiveFluxes = 0;
    for(int iel = 0; iel < nel ;iel++) {
        TPZCompEl *cel = fluxmesh->Element(iel);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        int matID = gel->MaterialId();
        if(meshdim != fluxmesh->Dimension()) DebugStop();
        bool IsBC = (fProblemConfig.bcmaterialids.find(matID) != fProblemConfig.bcmaterialids.end());
        if(!IsBC && matID != fHybridizer.fFluxMatId) {
            DebugStop();
        }
        nEffectiveFluxes++;
        int sideminus = gel->NSides(meshdim-2);
        // Check if every corner
        for(int iside = 0 ; iside < sideminus ; iside++){
            TPZGeoElSide minusSide(gel,iside);
            if (minusSide.Dimension() != 0) DebugStop();
            TPZStack<TPZCompElSide> minusNeigh;
            minusSide.ConnectedCompElementList(minusNeigh, 1, 0);
            if(minusNeigh.size() == 0) DebugStop();
        }
    }
    if(nEffectiveFluxes == 0) DebugStop();
}

void TPZHybridH1ErrorEstimator::VerifySkeletonContinuity(TPZCompMesh* pressuremesh){
    TPZGeoMesh *gmesh = pressuremesh->Reference();
    gmesh->ResetReference();
    int nel = pressuremesh->NElements();
    for(int iel = 0; iel < nel ;iel++) {
        TPZCompEl *cel = pressuremesh->Element(iel);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if (gel->MaterialId() == fPressureSkeletonMatId) {
            gel->SetReference(cel);
        }
    }

    TPZBlock &block =  pressuremesh->Block();
    TPZFMatrix<STATE> &sol = pressuremesh->Solution();
    for(int iel = 0; iel < nel ; iel++) {
        TPZCompEl *cel = pressuremesh->Element(iel);
        if (!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if (gel->MaterialId() != fPressureSkeletonMatId) continue;
        for (int iside = 0; iside < gel->NSides(); iside++) {
            TPZGeoElSide gelside(gel, iside);
            if (gelside.Dimension() > 1) DebugStop();  // 3D not supported yet
            if (gelside.Dimension() == 0) {
                TPZCompElSide celside(cel, gelside.Side());
                TPZConnect &elcon = pressuremesh->ConnectVec()[celside.ConnectIndex()];
                if (elcon.NShape() != 1 || elcon.NState() != 1) DebugStop();
                STATE elsol = sol.at(block.at(elcon.SequenceNumber(), 0, 0, 0));
                TPZStack<TPZCompElSide> neighSides;
                gelside.ConnectedCompElementList(neighSides, 1, 0);

                for (int ind = 0; ind < neighSides.size(); ind++) {
                    TPZCompElSide neighCelSide = neighSides[ind];
                    if (!(neighCelSide.Element())) continue;
                    TPZGeoElSide neighGelSide = neighCelSide.Reference();
                    if (!(neighGelSide.Element()) || neighGelSide.Element()->MaterialId() != fPressureSkeletonMatId)
                        DebugStop();
                    if (neighGelSide.Dimension() != 0) continue;

                    TPZConnect &c1 = pressuremesh->ConnectVec()[neighCelSide.ConnectIndex()];
                    if (c1.NShape() != 1 || c1.NState() != 1) DebugStop();
                    STATE c1sol = sol.at(block.at(c1.SequenceNumber(), 0, 0, 0));

#ifdef LOG4CXX
                    if (logger->isDebugEnabled()) {
                        if (!IsZero(elsol - c1sol)) {
                            std::stringstream sout;

                            TPZManVector<REAL> x(3);
                            TPZManVector<REAL> pt(gelside.Dimension(), 0);
                            gelside.X(pt, x);

                            TPZManVector<REAL> xcenter(3, 0);
                            TPZGeoElSide gelHigherSide(gelside.Element());
                            TPZManVector<REAL> ptcenter(gelHigherSide.Dimension(), 0);
                            gelHigherSide.X(ptcenter, xcenter);

                            TPZManVector<REAL> xneigh(3);
                            TPZManVector<REAL> ptneigh(neighGelSide.Dimension(), 0);
                            neighGelSide.X(ptneigh, xneigh);

                            TPZManVector<REAL> xneighcenter(3, 0);
                            TPZGeoElSide neighHigherSide(neighGelSide.Element());
                            TPZManVector<REAL> ptneighcenter(neighHigherSide.Dimension(), 0);
                            neighHigherSide.X(ptneighcenter, xneighcenter);

                            sout << "\ngel/side =  " << gelside.Id() << "/" << gelside.Side();
                            sout << "; gel center coordinates: [ " << xcenter[0] << ", " << xcenter[1] << ", "
                                 << xcenter[2] << "]\n";
                            sout << "neigh/side =  " << neighGelSide.Id() << "/" << neighGelSide.Side();
                            sout << "; neigh center coordinates: [ " << xneighcenter[0] << ", " << xneighcenter[1]
                                 << ", "
                                 << xneighcenter[2] << "]\n";
                            sout << "Side solution =  " << elsol << "\n";
                            sout << "Neigh solution = " << c1sol << "\n";
                            sout << "Diff = " << elsol - c1sol << "\n";
                            sout << "Side coord:  [" << x[0] << ", " << x[1] << ", " << x[2] << "]\n";
                            sout << "Neigh coord: [" << xneigh[0] << ", " << xneigh[1] << ", " << xneigh[2] << "]\n";

                            c1.Print(*pressuremesh, sout);

                            //std::cout << sout.str(); // TODO remove
                            LOGPZ_DEBUG(logger, sout.str())
                        }
                    }
#endif
                }
            }
            if (gelside.Dimension() == 1) {
                TPZStack<TPZCompElSide> celstack;
                gelside.EqualLevelCompElementList(celstack, 1, 0);

                TPZCompElSide large = gelside.LowerLevelCompElementList2(1);
                if (large) celstack.Push(large);

                if (celstack.size() == 0) continue;

                int intOrder = 3;

                TPZIntPoints *intRule = gelside.CreateIntegrationRule(intOrder);

                // Iterates through the comp sides connected to the reference gelside
                int nstack = celstack.size();
                for (int ist = 0; ist < nstack; ist++) {
                    TPZCompElSide cneighbour = celstack[ist];
                    if (!cneighbour) continue;
                    TPZGeoElSide neighbour = cneighbour.Reference();
                    if (!(neighbour.Element()) || neighbour.Element()->MaterialId() != fPressureSkeletonMatId)
                        DebugStop();

                    // Filters comp sides in elements of highest dimension (2 or 3)
                    if (neighbour.Element()->Dimension() > 1) DebugStop(); // 3D not supported
                    if (neighbour.Element()->Dimension() == 0) DebugStop();

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

                        TPZManVector<STATE> sol0(1), sol1(1);
                        cel->Solution(pt0, 0, sol0);
                        cneighbour.Element()->Solution(pt1, 0, sol1);

#ifdef LOG4CXX
                        if (logger->isDebugEnabled()) {
                            if (!IsZero(sol1[0] - sol0[0])) {
                                std::stringstream sout;
                                sout << "\ngel/side =  " << gelside.Id() << "/" << gelside.Side() << "\n";
                                sout << "neigh/side =  " << neighbour.Id() << "/" << neighbour.Side() << "\n";
                                sout << "Side solution =  " << sol0[0] << "\n";
                                sout << "Neigh solution = " << sol1[0] << "\n";
                                sout << "Diff = " << sol1[0] - sol0[0] << "\n";
                                sout << "Side coord:  [" << x0[0] << ", " << x0[1] << ", " << x0[2] << "]\n";
                                sout << "Neigh coord: [" << x1[0] << ", " << x1[1] << ", " << x1[2] << "]\n";

                                //std::cout << sout.str(); // TODO remove
                                LOGPZ_DEBUG(logger, sout.str())
                            }
                        }
#endif
                    }
                }
                delete intRule;
            }
        }
    }
}

void TPZHybridH1ErrorEstimator::VerifySolutionConsistency(TPZCompMesh *cmesh) {
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
            
            int intOrder = 3;
            
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
                            sout << "\ngel/side =  " << gelside.Id() <<"/" << gelside.Side() << "\n";
                            sout << "neigh/side =  " << neighbour.Id() <<"/" << neighbour.Side() << "\n";
                            sout << "Side solution =  " << sol0[0] << "\n";
                            sout << "Neigh solution = " << sol1[0] << "\n";
                            sout << "Diff = " << sol1[0] - sol0[0] << "\n";
                            sout << "Side coord:  [" << x0[0] << ", " << x0[1] << ", " << x0[2] << "]\n";
                            sout << "Neigh coord: [" << x1[0] << ", " << x1[1] << ", " << x1[2] << "]\n";
                            //std::cout << sout.str(); // TODO remove
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

void TPZHybridH1ErrorEstimator::PrepareElementsForH1Reconstruction() {

    // This vector stores the connects from elements which have a neighbour of
    // an internal boundary material. We don't want to condense these connects,
    // so we are later incrementing the number of elements connected to them.
    // Then we compute the stiffness matrix and load the solution of the
    // internal degrees of freedom.
    TPZManVector<int64_t> connectsToIncrement(fPostProcMesh.NConnects(), -1);
    fPostProcMesh.ComputeNodElCon();

    TPZCompMesh *pressureMesh = fPostProcMesh.MeshVector()[1];
    pressureMesh->LoadReferences();

#ifdef PZDEBUG
    {
        string dirPath = fDebugDirName + "/";
        std::ofstream txtPostProcMesh(dirPath + "PostProcMeshB4IncrementingConnects.txt");
        fPostProcMesh.ShortPrint(txtPostProcMesh);
        std::ofstream txtPressureMesh(dirPath + "PressureMeshB4IncrementingConnects.txt");
        pressureMesh->Print(txtPressureMesh);
    }
#endif

    for (int64_t el = 0; el < pressureMesh->NElements(); el++) {
        TPZCompEl *cel = pressureMesh->Element(el);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (!gel) continue;

        if (gel->MaterialId() != fPressureSkeletonMatId) continue;

        TPZGeoElSide skelSide(gel, gel->NSides() - 1);
        //if (skelSide.HasNeighbour(fProblemConfig.bcmaterialids)) continue;

        TPZStack<TPZCompElSide> compNeighSides;
        skelSide.EqualLevelCompElementList(compNeighSides, 1, 0);
        if (compNeighSides.size() ==0) DebugStop();

        for (int i = 0; i < compNeighSides.size(); i++) {
            TPZCompEl *neighCel = compNeighSides[i].Element();
            TPZInterpolatedElement *neighIntEl = dynamic_cast<TPZInterpolatedElement *>(neighCel);
            if (!neighIntEl) DebugStop();

            int sideNum = compNeighSides[i].Side();
            int nCon = neighIntEl->NSideConnects(sideNum);

            for (int iCon = 0; iCon < nCon; iCon++) {
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
    TPZCompMesh *fluxMesh = fPostProcMesh.MeshVector()[0];
    int nfc = 0;
    if(fluxMesh && fPostProcMesh.GetActiveApproximationSpaces()[0]){
        nfc = fluxMesh->NConnects();
    }
    for (int64_t i = 0; i < fPostProcMesh.NConnects(); i++) {
        if (connectsToIncrement[i] == 1) {
            fPostProcMesh.ConnectVec()[nfc + i].IncrementElConnected();
        }
    }
#ifdef PZDEBUG
    {
        std::string dirPath = fDebugDirName + "/";
        std::ofstream txtPostProcMesh(dirPath + "PostProcMeshAfterIncrementingConnects.txt");
        fPostProcMesh.ShortPrint(txtPostProcMesh);
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
        TPZHybridH1ErrorEstimateMaterial *errormat = dynamic_cast<TPZHybridH1ErrorEstimateMaterial *>(mat);
        if (errormat) {
            errormat->fNeumannLocalProblem = false;
        }
    }

    fPostProcMesh.CleanUpUnconnectedNodes();
    pressureMesh->CleanUpUnconnectedNodes();

#ifdef PZDEBUG
    string dirPath = fDebugDirName + "/";
    std::ofstream outTXT(dirPath + "PostProcessMeshAfterPreparingElements.txt");
    fPostProcMesh.Print(outTXT);
#endif
}

void TPZHybridH1ErrorEstimator::CopySolutionFromSkeleton() {
    
    TPZCompMesh *pressuremesh = PressureMesh();
    //        {
    //            std::ofstream out("MeshBeforeCopySkeleton.txt");
    //            pressuremesh->Print(out);
    //        }
    pressuremesh->Reference()->ResetReference();
    pressuremesh->LoadReferences();
    int dim = pressuremesh->Dimension();
    int64_t nel = pressuremesh->NElements();

    TPZBlock &block = pressuremesh->Block();
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
            //
            TPZGeoElSide gelside(gel, is);
            TPZConnect &c = intel->Connect(is);
            int64_t c_seqnum = c.SequenceNumber();
            int c_blocksize = c.NShape() * c.NState();
            //TPZGeoElSide gelside(gel,is);
            TPZStack<TPZCompElSide> celstack;
            gelside.EqualLevelCompElementList(celstack, 1, 0);
            //std::cout << "celstack size: "<< celstack.size() <<"\n";
            int nst = celstack.NElements();
            for (int ist = 0; ist < nst; ist++) {
                TPZCompElSide cneigh = celstack[ist];
                TPZGeoElSide gneigh = cneigh.Reference();
                //std::cout << "Gel/side: " << intel->Index() << "/" << gelside.Side() << "\nNeighbour Gel/side: " <<cneigh.Element()->Index() << "/" << cneigh.Side() << "\n\n";
                if ((gneigh.Element()->MaterialId() == this->fPressureSkeletonMatId)||
                    IsDirichletCondition(gneigh)) {
                    TPZInterpolatedElement *intelneigh = dynamic_cast<TPZInterpolatedElement *>(cneigh.Element());
                    if (!intelneigh) DebugStop();
                    TPZConnect &con_neigh = intelneigh->Connect(cneigh.Side());
                    int64_t con_seqnum = con_neigh.SequenceNumber();
                    int con_size = con_neigh.NState() * con_neigh.NShape();
                    if (con_size != c_blocksize) DebugStop();
                    for (int ibl = 0; ibl < con_size; ibl++) {
                        sol.at(block.at(c_seqnum,0,ibl,0)) = sol.at(block.at(con_seqnum,0,ibl,0));
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
}


/// compute the pressure weights and material weights
// fills in the data structure pressureweights and matid_weights
void TPZHybridH1ErrorEstimator::ComputePressureWeights() {
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
        // TODO: Case HybridSquared: insert if( ... || matid == fLagrangeInterface)
        if (matid == fPressureSkeletonMatId ) {
            continue;
        }
        if (!mat)  DebugStop();

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
            // Neumann
            else {
                this->fPressureweights[el] = 0.;
                fMatid_weights[matid] = 0.;
                continue;
            }
        }
        
        else{


            TPZMatLaplacianHybrid *matlaplacian = dynamic_cast<TPZMatLaplacianHybrid *>(mat);
            if (!matlaplacian) DebugStop();
            
            REAL perm;
            matlaplacian->GetMaxPermeability(perm);
            if (IsZero(perm)) DebugStop();
            this->fPressureweights[el] = perm;
            fMatid_weights[matid] = perm;
        }
    }
}

void TPZHybridH1ErrorEstimator::PlotState(const std::string& filename, int targetDim, TPZCompMesh* cmesh) {
    
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

void TPZHybridH1ErrorEstimator::VerifyAverage(int target_dim) {


    TPZCompMesh *pressure_mesh = PressureMesh();
    TPZGeoMesh *gmesh = pressure_mesh->Reference();

    //    std::ofstream out("PressureToAverage.txt");
    //    pressure_mesh->Print(out);

    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    int64_t nel = pressure_mesh->NElements();
    // load the pressure elements of dimension target_dim+1
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
        if (!cel) {
            continue;
        }

        TPZGeoEl *gel = cel->Reference();

        if (!gel || gel->Dimension() != target_dim) {
            continue;
        }

        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);

        // Skip calculation if the element is a boundary condition
        int matid = gel->MaterialId();
        TPZMaterial *mat = pressure_mesh->FindMaterial(matid);
        // TODO change this. Look for matIDs in bcMatIds instead. Only cast in debug mode for further checks
#ifdef PZDEBUG
        TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
        if (bc) continue;

#endif
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

        // Start veryfication

        TPZGeoMesh *gmesh = pressure_mesh->Reference();
        int dim = gmesh->Dimension();
        gmesh->ResetReference();
        pressure_mesh->LoadReferences();
        int nc = cel->NConnects();
        int order = cel->Connect(nc - 1).Order();

#ifdef PZDEBUG
        for (int ic = 0; ic < nc; ic++) {
            if (cel->Connect(ic).HasDependency()) DebugStop();
        }
#endif

        std::cout << "Computing average for compel " << el << ", gel: " << cel->Reference()->Index() << "\n";

        int target_dim = gel->Dimension();
        if (target_dim == dim - 1 && gel->MaterialId() != fPressureSkeletonMatId) {
            DebugStop();
        }

        // We denote large the skeleton element passed to this method, which is whether the skeleton between two elements
        // of the same level or the larger skeleton in the region of a hanging node.
        TPZGeoElSide largeSkeletonSide(gel, gel->NSides() - 1);

        // Stack of TPZCompElSides containing the volumetric neighbours from which we obtain the pressure solution to
        // be used in the average computation.
        TPZCompElSide SkelSide = largeSkeletonSide.Reference();
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
                TPZCompEl *smallSkelCel = smallerSkelSides[i].Element()->Reference();
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
        TPZTransform<REAL> tmp = volumeNeigh->SideToSideTransform(volumeNeighSides[0].Side(),
                                                                  volumeNeigh->NSides() - 1);
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
            std::unique_ptr<TPZIntPoints> intpoints(
                    gel->CreateSideIntegrationRule(integrationGeoElSide.Side(), 2 * order));

            REAL left_weight = fMatid_weights[leftVolumeGel->MaterialId()] / sum_weights;
            rightVolumeGel = volumeNeighSides[iskel + 1].Element()->Reference();
            REAL right_weight = fMatid_weights[rightVolumeGel->MaterialId()] / sum_weights;

            if (!noHangingSide) {
                AdditionalAverageWeights(largeSkeletonSide.Element(), integrationGeoElSide.Element(), left_weight,
                                         right_weight, sum_weights);
            }

            int64_t nintpoints = intpoints->NPoints();
            //std::cout <<"########\n\n skelId = " << gel->Index() << "; leftVolId: " <<  leftVolumeGel->Index() <<"; rightVolId: " <<  rightVolumeGel->Index() << "\n";
            for (int64_t ip = 0; ip < nintpoints; ip++) {
                TPZManVector<REAL, 3> pt_left_skel(target_dim, 0.), pt_right_skel(target_dim, 0.);
                TPZManVector<REAL, 3> pt_left_vol(target_dim + 1, 0.), pt_right_vol(target_dim + 1, 0.);
                TPZManVector<REAL, 3> x_left(3, 0.);

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
                TPZVec<STATE> left_sol, right_sol,skel_sol;
                volumeNeighSides[0].Element()->Solution(pt_left_vol, 0, left_sol);
                volumeNeighSides[iskel + 1].Element()->Solution(pt_right_vol, 0, right_sol);

                STATE average_sol = left_weight * left_sol[0] + right_weight * right_sol[0];

                TPZFNMatrix<9, REAL> jac(dim, dim), jacinv(dim, dim), axes(dim, 3);
                REAL detjac;
                integrationGeoElSide.Jacobian(pt_right_skel, jac, axes, detjac, jacinv);

                SkelSide.Element()->Solution(pt_left_skel,0,skel_sol);

                leftVolumeGel->X(pt_left_vol,x_left);

                //std::cout <<"X: [" << x_left[0] << "," <<x_left[1] << "," <<x_left[2] << "]: skelSol: " << skel_sol[0] << "; leftSol/weight: " << left_sol[0] << "/" << left_weight << "; rightSol/weight: " << right_sol[0] << "/" << right_weight <<"\n";

                if(!IsZero(skel_sol[0]- (left_sol[0]+right_sol[0])/2)) {
                    DebugStop();
                }
            }
            //std::cout <<"\n";
        }
    }
}
