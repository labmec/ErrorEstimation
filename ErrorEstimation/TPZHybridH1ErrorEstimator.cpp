
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
#include "TPZNullMaterialCS.h"
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
    TPZVec<TPZCompMesh *> &meshvec = fMultiphysicsReconstructionMesh->MeshVector();
    fMultiphysicsReconstructionMesh->Reference()->ResetReference();
    for (int i = 0; i < 4; i++) {
        if (!meshvec[i]) continue;
//        TPZCompMesh *cmesh = meshvec[i];
//        cmesh->SetReference(nullptr);
//        for (int64_t ic = 0; ic < meshvec[i]->NConnects(); ic++) {
//            TPZConnect &c = meshvec[i]->ConnectVec()[ic];
//            c.RemoveDepend();
//        }
        meshvec[i] = 0;
    }
    for (int64_t ic = 0; ic < fMultiphysicsReconstructionMesh->NConnects(); ic++) {
        TPZConnect &c = fMultiphysicsReconstructionMesh->ConnectVec()[ic];
        c.RemoveDepend();
    }
}

TPZVec<REAL> TPZHybridH1ErrorEstimator::PostProcess() {
    TPZLinearAnalysis an(fMultiphysicsReconstructionMesh, false);
    
    // The solution is expanded to store errors,
    // Therefore it is required to account for the original solution and the errors.
    int numErrors = 7;
    numErrors++; 
    TPZVec<REAL> errorVec = ComputeErrors(&an,numErrors);
    
    std::cout << "\n############\n";
    std::cout << "Computing Error H1 reconstruction\n";
    std::cout <<       "||u_h-u||:                \t" << 
    (errorVec)[0]<< "\n||u_h-s_h||:              \t" <<
    (errorVec)[1]<< "\n||Grad(u_h)-Grad(u)||:    \t" <<
    (errorVec)[2]<< "\n||Grad(u_h)-Grad(s_h)||:  \t"   <<
    (errorVec)[3]<< "\n||Grad(u_h)+t_h||:        \t"   <<
    (errorVec)[4]<< "\n||f-div(t_h)||:           \t"   <<
    (errorVec)[5]<< "\n||f-projL2(f)||:          \t" << (errorVec)[6] <<"\n";

    
    TPZCompMeshTools::UnCondensedElements(fMultiphysicsReconstructionMesh);
    TPZCompMeshTools::UnGroupElements(fMultiphysicsReconstructionMesh);
    
    //Erro global
    std::ofstream myfile;
    myfile.open("ErrorsReconstruction.txt", std::ios::app);
    myfile << "\n\n Estimator errors for Problem " << *fproblemname;
    myfile << "\n-------------------------------------------------- \n";
    myfile << "Ndiv = " << fnDivisions <<" Order k= " << forderFEM_k << " Order n= "<< forderFEM_n<<"\n";
    //myfile << "DOF Total = " << fMultiphysicsReconstructionMesh->NEquations() << "\n";
    myfile << "||u-u_h|| = " << (errorVec)[0] << "\n";
    myfile << "||u_h-s_h|| = " << (errorVec)[1] << "\n";
    myfile << "e_{ex}: ||K^{0.5}.grad(u_h-u)|| = " << (errorVec)[2] << "\n";
    myfile << "n_{NC}: ||K^{0.5}.grad(u_h-s_h)|| = " << (errorVec)[3] << "\n";
    myfile << "n_{F} : ||K^{0.5}.[grad(u_h)-invK.T_h]|| = " << (errorVec)[4] << "\n";
    myfile <<"Residual ErrorL2= "<< (errorVec)[5] << "\n";
    //myfile <<"Global Index = "<< sqrt(errorVec[4] + errorVec[3]) / sqrt(errorVec[2]);

    myfile.close();

    double globalIndex;
    ComputeEffectivityIndices(globalIndex);
    
    PrintSolutionVTK(an);

    return errorVec;
}

/// compute the element errors comparing the reconstructed solution based on average pressures
/// with the original solution
void TPZHybridH1ErrorEstimator::PostProcess(REAL threshold, std::set<int64_t> &geltodivide) {
    
    TPZVec<REAL> errorVec = PostProcess();
    
    errorVec.resize(fMultiphysicsReconstructionMesh->Reference()->NElements());
    for (REAL & elementerror : errorVec) {
        elementerror = 0;
    }
    STATE maxerror = 0.;
    for (int64_t i = 0; i < fMultiphysicsReconstructionMesh->NElements(); i++) {
        TPZCompEl *cel = fMultiphysicsReconstructionMesh->Element(i);
        if (!cel) continue;
        TPZGeoEl* gel = fMultiphysicsReconstructionMesh->Element(i)->Reference();
        if (!gel) continue;
        TPZFMatrix<STATE> &elsol = fMultiphysicsReconstructionMesh->ElementSolution();
        STATE error = elsol(i,3);
        (errorVec)[gel->Index()] = error;
        if(maxerror < error) maxerror = error;
    }
    geltodivide.clear();
    for (int64_t i = 0; i<errorVec.size(); i++) {
        REAL elementerror = errorVec[i];
        if(elementerror > threshold*maxerror) geltodivide.insert(i);
    }

}

void TPZHybridH1ErrorEstimator::FillVTKoutputVariables(TPZStack<std::string> &scalnames,TPZStack<std::string> &vecnames){
     if (fExact) {
            scalnames.Push("u");
            scalnames.Push("uh_minus_u");
            scalnames.Push("KGradUh_minus_KGradU");
            scalnames.Push("L2EffectivityIndex");
            scalnames.Push("EffectivityIndex");
            vecnames.Push("minus_KGradU");
        }
        scalnames.Push("uh");
        scalnames.Push("sh");
        scalnames.Push("uh_minus_sh");
        scalnames.Push("KGradSh_minus_KGradUh");
        scalnames.Push("residual");
        scalnames.Push("th_plus_KGradUh");
        vecnames.Push("minus_KGradUh");
        vecnames.Push("th");
        vecnames.Push("minus_KGradSh");
        scalnames.Push("EstimatedError");
        scalnames.Push("POrder");
}

/// create the post processed multiphysics mesh (which is necessarily
/// hybridized)
void TPZHybridH1ErrorEstimator::CreatePostProcessingMesh() {

if (fMultiphysicsReconstructionMesh->MeshVector().size()) {
        DebugStop();
}

if(fH1conformMesh == NULL || fHDivconformMesh == NULL){
    DebugStop();
}

#ifdef ERRORESTIMATION_DEBUG666
    {
        std::ofstream out("OriginalFlux.txt");
        fOriginal->MeshVector()[0]->Print(out);
        std::ofstream out2("OriginalPotential.txt");
        fOriginal->MeshVector()[1]->Print(out2);
        std::ofstream out3("fPostProcMeshMeshHybrid.txt");
        fMultiphysicsReconstructionMesh->Print(out3);
        std::ofstream out4("OriginalMeshHybrid.txt");
        fOriginal->Print(out4);
    }
#endif

    // initialize the post processing mesh
    fMultiphysicsReconstructionMesh->SetReference(fOriginal->Reference());
    int dim = fOriginal->Dimension(); 
    fOriginal->CopyMaterials(*fMultiphysicsReconstructionMesh);
    
    auto skelmat = new TPZNullMaterialCS<>(fSkeletonMatId,dim,1);
    fMultiphysicsReconstructionMesh->InsertMaterialObject(skelmat);

    SwitchMaterialObjects();

    TPZManVector<TPZCompMesh *> mesh_vectors(5, 0);
    TPZManVector<int> active(5, 0);

    mesh_vectors[0] = fHDivconformMesh;
    mesh_vectors[1] = fH1conformMesh; //sh
    mesh_vectors[2] = fOriginal->MeshVector()[0];// flux
    mesh_vectors[3] = fOriginal->MeshVector()[1];// potential
    mesh_vectors[4] = ForceProjectionMesh();

    active[1] = 1;

    fMultiphysicsReconstructionMesh->SetAllCreateFunctionsMultiphysicElem();

    fMultiphysicsReconstructionMesh->BuildMultiphysicsSpace(active, mesh_vectors);

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
    TPZGeoMesh *gmesh = fOriginal->Reference();
    TPZCompMesh *forceProj = new TPZCompMesh(gmesh);
    int dimMesh = gmesh->Dimension();

    int potential_order = forderFEM_k+forderFEM_n;
    forceProj->SetDefaultOrder(potential_order);
    forceProj->SetDimModel(dimMesh);

    forceProj->SetAllCreateFunctionsContinuous(); //H1 functions
    forceProj->ApproxSpace().CreateDisconnectedElements(true);

    for(auto matid:fmaterialids){
        TPZNullMaterial<> *material = new TPZNullMaterial<>(matid);
        material->SetDimension(dimMesh);
        forceProj->InsertMaterialObject(material);
        int porder = 5;
        material->SetForcingFunction(fExact->ForceFunc(),porder);
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
        std::ofstream outCon(fFolderOutput + "forceProj.txt");
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
    return fMultiphysicsReconstructionMesh->MeshVector()[1];
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
    
    TPZCompMesh *cmesh = fMultiphysicsReconstructionMesh;
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
                if(fbcmaterialids.find(elmatid) != fbcmaterialids.end())
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
    REAL globalResidual =0., globalProjResidual =0.;
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
        elementResidual = (hk/M_PI)*elsol(el, 5);
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
                oscilatorytherm = elsol(el, i + 3);
                oscilatorytherm *= (hk/M_PI);
                fluxestimator = elsol(el, i + 2);

                ex += elsol(el, 2)*elsol(el, 2);
                n1 += elsol(el, 3)*elsol(el, 3);
                n2n3 += (oscilatorytherm + elsol(el, 4))*(oscilatorytherm + elsol(el, 4));
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

    std::stringstream ss;
    ss << "\n\n";
    ss << "frac{hk}{pi}||f-Div(T_h)|| =\t" << globalResidual << "\n";
    ss << "frac{hk}{pi}||f-Proj(f)||  =\t" << globalProjResidual << "\n";
    ss << "I_{EF}                       =\t" << globalIndex << "\n";
    ss << "Estimated Error              =\t" << sqrt(n1+n2n3) << "\n";
    ss << "\n\n";
    std::cout << ss.str();

    {
        std::ofstream myfile;
        myfile.open(fFolderOutput + "ErrorsReconstruction.txt", std::ios::app);
        myfile << ss.str();
        myfile.close();
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

/// switch material object from mixed poisson to TPZMixedHdivErrorEstimate
void TPZHybridH1ErrorEstimator:: SwitchMaterialObjects() {
    for (auto mat : fMultiphysicsReconstructionMesh->MaterialVec()) {
        TPZMatLaplacianHybrid *matlaplacian=
                dynamic_cast<TPZMatLaplacianHybrid *>(mat.second);
        if (matlaplacian) {
            TPZHybridH1ErrorEstimateMaterial *newmat = new TPZHybridH1ErrorEstimateMaterial(*matlaplacian);

            for (auto bcmat : fMultiphysicsReconstructionMesh->MaterialVec()) {
                TPZBndCond *bc = dynamic_cast<TPZBndCond *>(bcmat.second);
                if (bc) {
                    bc->SetMaterial(newmat);
                }
            }
            fMultiphysicsReconstructionMesh->MaterialVec()[newmat->Id()] = newmat;
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

