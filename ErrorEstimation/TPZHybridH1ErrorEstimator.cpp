
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
#include "TPZElementMatrixT.h"
#include "TPZMixedHdivErrorEstimate.h"
#include "TPZAnalysis.h"
#include "TPZNullMaterial.h"
#include "TPZInterfaceEl.h"

#include "TPZHybridH1CreateHDivReconstruction.h"
#include "TPZHybridH1CreateH1Reconstruction.h"

#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"

#include "TPZMultiphysicsCompMesh.h"
#include "pzmultiphysicscompel.h"

#include "TPZVTKGeoMesh.h"

#include "TPZCompMeshTools.h"
#include "TPZMatLaplacianHybrid.h"
#include "TPZHybridH1ErrorEstimateMaterial.h"
#include "TPZGeoElSideAncestors.h"

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
    TPZLinearAnalysis an(&fPostProcMesh, false);
    
    if (fExact) {
        an.SetExact(fExact->ExactSolution());
    }
    
#ifdef ERRORESTIMATION_DEBUG2
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
    
#ifdef ERRORESTIMATION_DEBUG1
    {
        std::ofstream out("MeshToComputeError2.txt");
        fPostProcMesh.Print(out);
        
    }
#endif

    an.PostProcessError(errorVec, store);//calculo do erro com sol exata e aprox e armazena no elementsolution
    
    std::cout << "\n############\n\n";
    
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
    //myfile <<"Residual ErrorL2= "<< errorVec[4] << "\n";
    //myfile <<"Global Index = "<< sqrt(errorVec[4] + errorVec[3]) / sqrt(errorVec[2]);

    myfile.close();

    double globalIndex;
    ComputeEffectivityIndices(globalIndex);
    
    //GlobalEffectivityIndex();
    
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

void TPZHybridH1ErrorEstimator::FillVTKoutputVariables(TPZStack<std::string> &scalnames,TPZStack<std::string> &vecnames){
    DebugStop();
}

void TPZHybridH1ErrorEstimator::PostProcessing(TPZAnalysis &an) {
    
    TPZMaterial *mat = fPostProcMesh.FindMaterial(*fProblemConfig.materialids.begin());
    int varindex = -1;
    if (mat) varindex = mat->VariableIndex("PressureFEM");
    if (varindex != -1) {
        TPZStack<std::string> scalnames, vecnames;
        if (fExact) {
            scalnames.Push("PressureExact");
            scalnames.Push("PressureErrorExact");
            scalnames.Push("EnergyErrorExact");
            scalnames.Push("PressureEffectivityIndex");
            scalnames.Push("EnergyEffectivityIndex");
            scalnames.Push("EnergyErrorEstimated");
            vecnames.Push("FluxExact");
        }
        scalnames.Push("PressureFEM");
        scalnames.Push("PressureReconstructed");
        scalnames.Push("PressureErrorEstimate");
        scalnames.Push("NCIndex");
        scalnames.Push("NRIndex");
        scalnames.Push("NFIndex");
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

/// create the post processed multiphysics mesh (which is necessarily
/// hybridized)
void TPZHybridH1ErrorEstimator::CreatePostProcessingMesh() {

if (fPostProcMesh.MeshVector().size()) {
        DebugStop();
}

#ifdef ERRORESTIMATION_DEBUG666
    {
        std::ofstream out("OriginalFlux.txt");
        fOriginal->MeshVector()[0]->Print(out);
        std::ofstream out2("OriginalPotential.txt");
        fOriginal->MeshVector()[1]->Print(out2);
        std::ofstream out3("fPostProcMeshMeshHybrid.txt");
        fPostProcMesh.Print(out3);
        std::ofstream out4("OriginalMeshHybrid.txt");
        fOriginal->Print(out4);
    }
#endif

    // initialize the post processing mesh
    fPostProcMesh.SetReference(fOriginal->Reference());
    int dim = fOriginal->Dimension(); 
    fOriginal->CopyMaterials(fPostProcMesh);

    SwitchMaterialObjects();

    TPZManVector<TPZCompMesh *> mesh_vectors(5, 0);
    TPZManVector<int> active(5, 0);

    auto myHdivMeshCreator = new TPZHybridH1CreateHDivReconstruction();
    auto myHdivMesh = myHdivMeshCreator->CreateFluxReconstructionMesh();
    mesh_vectors[0] = myHdivMesh->MeshVector()[0]; //CreateFluxReconstructionMesh()->MeshVector()[0];
    //myHdivMeshCreator->PostProcess(myHdivMesh);

    auto myH1MeshCreator = new TPZHybridH1CreateH1Reconstruction();
    TPZMultiphysicsCompMesh* myH1Mesh = myH1MeshCreator->CreateH1ReconstructionMesh();
    myH1MeshCreator->PostProcess(myH1Mesh);

    mesh_vectors[1] = myH1Mesh->MeshVector()[1]; //sh
    mesh_vectors[2] = fOriginal->MeshVector()[0];// flux
    mesh_vectors[3] = fOriginal->MeshVector()[1];// potential
    mesh_vectors[4] = ForceProjectionMesh();

    active[1] = 1;

    fPostProcMesh.SetAllCreateFunctionsMultiphysicElem();

    fPostProcMesh.BuildMultiphysicsSpace(active, mesh_vectors);

#ifdef ERRORESTIMATION_DEBUG
    {
        std::string dirPath = fDebugDirName + "/";
        /*std::ofstream out(dirPath + "EnrichedFluxBorder.txt");
        mesh_vectors[0]->Print(out);*/
        std::ofstream out2(dirPath + "EnrichedPressure.txt");
        mesh_vectors[1]->Print(out2);
    }
#endif
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
        TPZNullMaterial<> *material = new TPZNullMaterial<>(matid);
        material->SetDimension(dimMesh);
        forceProj->InsertMaterialObject(material);
        int porder = 5;
        material->SetForcingFunction(fProblemConfig.exact->ForceFunc(),porder);
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
        TPZMaterial *matabstr = cel->Material();
        TPZMatSingleSpace *material = dynamic_cast<TPZMatSingleSpace*>(matabstr);
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
        TPZManVector<REAL, 4> intpointtemp(dim, 0.),x(3, 0.);
        for (int int_ind = 0; int_ind < intrulepoints; ++int_ind) {
            intrule->Point(int_ind, intpointtemp, weight);
            gel->Jacobian(intpointtemp, jac, axe, detJac, jacInv);
            weight *= fabs(detJac);

            TPZMaterialDataT<STATE> data;
            intel->Shape(intpointtemp, phi, dphi);
            gel->X(intpointtemp,x);

            TPZMaterialT<STATE> *material = dynamic_cast<TPZMaterialT<STATE>*>(matabstr);
            if(!material->HasForcingFunction()) DebugStop();
            STATE force;
            TPZManVector<STATE> res(3);
            material->ForcingFunction()(x,res);
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
        TPZCompMeshTools::PrintConnectInfoByGeoElement(forceProj, outCon, {}, false, true);
    }

    return forceProj->Clone();
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

/// return a pointer to the pressure mesh
TPZCompMesh *TPZHybridH1ErrorEstimator::PressureMesh() {
    return fPostProcMesh.MeshVector()[1];
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
void TPZHybridH1ErrorEstimator::ComputeEffectivityIndices(double &globalIndex) {
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
    elsol.Resize(nrows, ncols+3);
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
    REAL n1 = 0.,n2n3 = 0.,ex = 0.;
    for (int64_t el = 0; el < nrows; el++) {
        
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if(subcmesh)
        {
            std::cout << "Computing submesh effectivity indices\n";
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

                ex += elsol(el, 2)*elsol(el, 2);
                n1 += elsol(el, 3)*elsol(el, 3);
                n2n3 += (oscilatorytherm + elsol(el, 5))*(oscilatorytherm + elsol(el, 5));
            }
            
            
            if (abs(ErrorEstimate) < tol) {
                elsol(el, ncols + i / 2) = 1.;
                dataIeff(el,0)=1.;
                
            }
            else {
                REAL EfIndex = sqrt(ErrorEstimate*ErrorEstimate +(oscilatorytherm+fluxestimator)*(oscilatorytherm+fluxestimator))/ErrorExact;
                dataIeff(el,0)= EfIndex;
                
                elsol(el, ncols + i / 2) = EfIndex;
                if(i == 2){
                    elsol(el, ncols + 2) = sqrt(ErrorEstimate*ErrorEstimate +(oscilatorytherm+fluxestimator)*(oscilatorytherm+fluxestimator));
                }

            }
        }
    }

    globalResidual = sqrt(globalResidual);
    globalProjResidual = sqrt(globalProjResidual);
    globalIndex = sqrt((n1+n2n3)/ex);
    std::cout << "\n\n";
    std::cout << "Residual = " << globalResidual << "\n";
    std::cout << "ProjResidual = " << globalProjResidual << "\n";
    std::cout << "I_{EF} = " << globalIndex << "\n";
    std::cout << "Est. Error = " << sqrt(n1+n2n3) << "\n";
    std::cout << "\n\n";

    {
        std::ofstream myfile;
        myfile.open("ErrorsReconstruction.txt", std::ios::app);
        myfile << "||f-Proj(f)|| = " << globalProjResidual << "\n";
        myfile << "||f-Div(T_h)|| = " << globalResidual << "\n";
        myfile << "I_{EF} = " << globalIndex << "\n";
        myfile << "Est. Error = " << sqrt(n1+n2n3) << "\n";
        myfile.close();
    }


    //  cmesh->ElementSolution().Print("ElSolution",std::cout);
    //    ofstream out("IeffPerElement3DEx.nb");
    //    dataIeff.Print("Ieff = ",out,EMathematicaInput);


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

void TPZHybridH1ErrorEstimator::PlotState(const std::string& filename, int targetDim, TPZCompMesh* cmesh) {
    
    std::ofstream out2("PressuretoStateGraph.txt");
    cmesh->Print(out2);
    
    {
        TPZLinearAnalysis an(cmesh, false);
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

