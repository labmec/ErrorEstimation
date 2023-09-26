//
// Created by victor on 16/03/2021.
//

#include "Solver.h"
#include <TPZMultiphysicsCompMesh.h>
#include "DataStructure.h"
#include "MeshInit.h"
#include "TPZCompMeshTools.h"
#include "TPZCreateMultiphysicsSpace.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "Tools.h"
#include "Output.h"
#include "pzvisualmatrix.h"
#include "MeshInit.h"
#include "TPZHybridH1ErrorEstimator.h"
#include "TPZHybridH1CreateHDivReconstruction.h"
#include "TPZHybridH1CreateH1Reconstruction.h"
#include "InputTreatment.h"
#include "ForcingFunction.h"
#include "TPZFrontSym.h"
#include "DataStructure.h"
#include "TPZPostProcessError.h"

using namespace std;
void Solve(ProblemConfig &config, PreConfig &preConfig){

    TPZCompMesh *cmesh = 0;
//    TPZCompMesh *cmesh = InsertCMeshH1(config,preConfig);
    TPZMultiphysicsCompMesh *multiCmesh = new TPZMultiphysicsCompMesh(config.gmesh);
    int interfaceMatID = -10;
    int fluxMatID = -10;
    int hybridLevel = 1;
    const clock_t start = clock();

    switch(preConfig.mode){
        case 0: //H1
            config.gmesh->ResetReference();
            cmesh = InsertCMeshH1(config,preConfig);
            config.PorderIncrement();
            //TPZCompMeshTools::CreatedCondensedElements(cmesh, false, false);
            SolveH1Problem(cmesh, config, preConfig);
            if (preConfig.estimateError){
                EstimateError(config, preConfig, fluxMatID, multiCmesh);
            }
            
            break;
        case 1: //Hybrid
            CreateHybridH1ComputationalMesh(multiCmesh, interfaceMatID, fluxMatID,preConfig, config,hybridLevel);
            SolveHybridH1Problem(multiCmesh, interfaceMatID, config, preConfig,hybridLevel);
            if (preConfig.estimateError){
                EstimateError(config, preConfig, fluxMatID, multiCmesh);
            }
            break;
        case 2: //Mixed
            CreateMixedComputationalMesh(multiCmesh, preConfig, config);
            SolveMixedProblem(multiCmesh, config, preConfig);
            if (preConfig.estimateError) EstimateError(config, preConfig, fluxMatID, multiCmesh);
            break;
        default:
            DebugStop();
            break;
    }
    FlushTime(preConfig,start);

    if(preConfig.debugger) Tools::DrawCompMesh(config,preConfig,cmesh,multiCmesh);
    delete multiCmesh;
}

void SolveDiff(PreConfig &hybConfig, PreConfig &mixConfig,char *argv[]){
    if (hybConfig.mode != 1 || mixConfig.mode != 2) DebugStop();

    ProblemConfig conf;
    Configure(conf,hybConfig.refLevel,hybConfig,argv);

    TPZMultiphysicsCompMesh *multiHyb = new TPZMultiphysicsCompMesh(conf.gmesh);
    TPZMultiphysicsCompMesh *multiMix = new TPZMultiphysicsCompMesh(conf.gmesh);

    int interfaceMatID = -10;
    int fluxMatID = -10;
    int hybridLevel = 1;

    CreateHybridH1ComputationalMesh(multiHyb, interfaceMatID, fluxMatID,hybConfig, conf,hybridLevel);
    SolveHybridH1Problem(multiHyb, interfaceMatID, conf, hybConfig,hybridLevel);

    CreateMixedComputationalMesh(multiMix, mixConfig, conf);
    SolveMixedProblem(multiMix, conf, mixConfig);

    if(hybConfig.debugger){
        TPZCompMesh *cmesh;
        Tools::DrawCompMesh(conf,hybConfig,cmesh,multiHyb);
        conf.problemname = mixConfig.problem;
        Tools::DrawCompMesh(conf,hybConfig,cmesh,multiHyb);
    }

    TPZMultiphysicsCompMesh *multiHybMix = new TPZMultiphysicsCompMesh(conf.gmesh);
    CreateHybMixCompMesh(multiHyb, multiMix, multiHybMix,hybConfig,conf);

    PostProcessHybMix(multiHybMix,hybConfig,conf);
}

void PostProcessHybMix(TPZMultiphysicsCompMesh *multHybMix,PreConfig &pConfig, ProblemConfig &config){

    TPZLinearAnalysis an(multHybMix);

    std::cout << "Post Processing ""Hyb - Mix"" difference " << std::endl;
    an.SetExact(config.exact.operator*().ExactSolution());

    TPZVec<REAL> errorVec; 
    //The following line suddenly breaks. NEvalErrors no longer belongs to the Material class.
    int64_t nErrorCols = -1;DebugStop();
    errorVec.resize(nErrorCols);
    for (int64_t i = 0; i < nErrorCols; i++) {
        errorVec[i] = 0;
    }

    int64_t nelem = multHybMix->NElements();
    multHybMix->LoadSolution(multHybMix->Solution());
    multHybMix->ExpandSolution();
    multHybMix->ElementSolution().Redim(nelem, nErrorCols - 1);

    an.PostProcessError(errorVec, true);
    PrintErrorsDiff(errorVec,config);

    ////PostProcess
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("hybP");
    scalnames.Push("mixP");
    scalnames.Push("exactP");
    scalnames.Push("DiffP");

    vecnames.Push("hybF");
    vecnames.Push("mixF");
    vecnames.Push("exactF");
    vecnames.Push("DiffF");

    int dim = pConfig.dim;
    std::string plotname;
    {
        std::stringstream out;
        out << pConfig.plotfile << "/" << config.problemname <<"_" << pConfig.topologyFileName << "_k-" << config.k << "_n-" << config.n;

        if(dim == 2) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";
        if(dim == 3) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";

        plotname = out.str();
    }
    int resolution = 3;
    an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
    an.PostProcess(resolution, dim);

}
void PrintErrorsDiff(TPZVec<REAL> errorVec, ProblemConfig &config){
    std::cout << "\nPotential error for HybH1 approx.: " << errorVec[0] << std::endl;
    std::cout << "Potential error for Mixed approx.: " << errorVec[1] << std::endl;
    std::cout << "Potential difference bet. HybH1 & Mix approx.: " << errorVec[2] << std::endl;

    std::cout << "Flux error for HybH1 approx.: " << errorVec[3] << std::endl;
    std::cout << "Flux error for Mixed approx.: " << errorVec[4] << std::endl;
    std::cout << "Flux difference bet. HybH1 & Mix approx.: " << errorVec[5] << std::endl;

    //Erro global
    std::ofstream myfile;
    myfile.open("HybMixDiff.txt", std::ios::app);
    myfile << "\n\n Estimator errors for Problem " << config.problemname;
    myfile << "\n-------------------------------------------------- \n";
    myfile << "Ndiv = " << config.ndivisions <<" Order k= " << config.k << " Order n= "<< config.n<<"\n";
    myfile << "\nPotential error for HybH1 approx.: " << errorVec[0] << std::endl;
    myfile << "Potential error for Mixed approx.: " << errorVec[1] << std::endl;
    myfile << "Potential difference bet. HybH1 & Mix approx.: " << errorVec[2] << std::endl;
    myfile << "Flux error for HybH1 approx.: " << errorVec[3] << std::endl;
    myfile << "Flux error for Mixed approx.: " << errorVec[4] << std::endl;
    myfile << "Flux difference bet. HybH1 & Mix approx.: " << errorVec[5] << std::endl;
}
void EstimateError(ProblemConfig &config, PreConfig &preConfig, int fluxMatID, TPZMultiphysicsCompMesh *multiCmesh){

    if(preConfig.mode == 0){
        TPZCompMesh* cmeshH1 = config.gmesh->Reference();
        
        TPZPostProcessError error(cmeshH1);
        
        TPZVec<STATE> estimatedelementerror;
        error.ComputeElementErrors(estimatedelementerror);
        
        TPZFMatrix<STATE> true_elerror(cmeshH1->ElementSolution());
        TPZFMatrix<STATE> estimate_elerror(error.MultiPhysicsMesh()->ElementSolution());
        
        STATE maxerror = 0.;
        int64_t nel = estimate_elerror.Rows();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = cmeshH1->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            TPZCompEl *mphys = gel->Reference();
            int64_t elindex2 = mphys->Index();
            STATE erro = estimate_elerror(elindex2,2);
            if(maxerror < erro){
                maxerror = erro;
            }
        }
        
        std::cout << "max estimated error within elements = " << maxerror <<"\n";
        
        // Select elements to hp refinement
        if(0){
            std::set<int64_t> gelstohref;
            std::map<int64_t,int> gelstoPplus;
            gelstohref.clear();
            gelstoPplus.clear();
            REAL threshold = config.division_threshold;
            for (int64_t i = 0; i < estimate_elerror.Rows(); i++) {
                REAL elementerror = estimate_elerror(i,2);
                if(elementerror > threshold*maxerror){
                    TPZCompEl* cel = cmeshH1->Element(i);
                    TPZGeoEl* gel = cel->Reference();
                    gelstohref.insert(gel->Index());
                } else if (threshold*maxerror >= elementerror && elementerror >(1-threshold)*maxerror){
                    TPZCompEl* cel = cmeshH1->Element(i);
                    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(cel);
                    int porder = intel->GetPreferredOrder();
                    TPZGeoEl* gel = cel->Reference();
                    if (gel->Dimension() != cmeshH1->Dimension()) continue;
                    gelstoPplus[gel->Index()] = porder+1;
                }
                
                //Selects elements adjacent to the selected node
//                TPZGeoNode node = config.gmesh->NodeVec()[0];
//                TPZCompEl* cel = cmeshH1->Element(i);
//                TPZGeoEl* gel = cel->Reference();
//                const int nnodes = gel->NNodes();
//                std::set<int64_t> nodeindices;
//                gel->GetNodeIndices(nodeindices);
//                if (nodeindices.count(node.Id())){
//                    gelstohref.insert(gel->Index());
//                }
            }
            config.fElIndexDivide.push_back(gelstohref);
            config.fElIndexPplus.push_back(gelstoPplus);
        }
        else{// h-ref to singular node adjacent elements and posible p-ref in ohers
            std::set<int64_t> gelstohref;
            std::map<int64_t,int> gelstoPplus;
            gelstohref.clear();
            gelstoPplus.clear();
            REAL threshold = config.division_threshold;

            for (int64_t i = 0; i < estimate_elerror.Rows(); i++) {
                REAL elementerror = estimate_elerror(i,2);
                
                //Selects elements adjacent to the selected node
                TPZGeoNode node = config.gmesh->NodeVec()[0];
                TPZCompEl* cel = cmeshH1->Element(i);
                TPZGeoEl* gel = cel->Reference();
                const int nnodes = gel->NNodes();
                std::set<int64_t> nodeindices;
                gel->GetNodeIndices(nodeindices);
                if (nodeindices.count(node.Id())){
                    gelstohref.insert(gel->Index());
                }
                else if (elementerror > (1-threshold)*maxerror){
                    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(cel);
                    int porder = intel->GetPreferredOrder();
                    gelstoPplus[gel->Index()] = porder+1;
                }
                
                config.fElIndexDivide.push_back(gelstohref);
                config.fElIndexPplus.push_back(gelstoPplus);
            }
        }
        
        PostProcessing(cmeshH1,true_elerror, estimate_elerror, config);
        
    }
        
    //if(preConfig.topologyMode != 2) DebugStop();
    if(preConfig.mode == 1){
        EstimatorConfig *estimatorConfig = new EstimatorConfig(multiCmesh,config,fluxMatID);

        // the empty multiphysics mesh is created here
        auto myHdivMeshCreator = new TPZHybridH1CreateHDivReconstruction(estimatorConfig);
        // here we do not only create the mesh, we also project the solution onto de hdiv space
        auto myHdivMesh = myHdivMeshCreator->CreateFluxReconstructionMesh();
        //myHdivMeshCreator->PostProcess();

        auto myH1MeshCreator = new TPZHybridH1CreateH1Reconstruction(estimatorConfig);
        TPZMultiphysicsCompMesh* myH1Mesh = myH1MeshCreator->CreateH1ReconstructionMesh();
        //myH1MeshCreator->PostProcess();

        TPZHybridH1ErrorEstimator *HybridH1Estimate = new TPZHybridH1ErrorEstimator(estimatorConfig);
        HybridH1Estimate->SetH1conformMesh(myH1MeshCreator->GetReconstructionMesh());
        HybridH1Estimate->SetHDivConformMesh(myHdivMeshCreator->GetReconstructionMesh());
        HybridH1Estimate->CreatePostProcessingMesh();
        std::set<int64_t> geltodivide;
        HybridH1Estimate->PostProcess(config.division_threshold,geltodivide);
        config.fElIndexDivide.push_back(geltodivide);
        
        delete HybridH1Estimate;
        delete myH1MeshCreator;
        delete myHdivMesh;
    }

    if (preConfig.mode == 2){
        TPZHDivErrorEstimator test(*multiCmesh);
        test.SetAnalyticSolution(config.exact);
        test.PotentialReconstruction();

        TPZManVector<REAL> elementerrors;
        TPZVec<REAL> errorVec;
        std::string outVTK = config.dir_name + "/out.vtk";
        test.ComputeErrors(errorVec, elementerrors,outVTK);
    }
}


void DrawMesh(ProblemConfig &config, PreConfig &preConfig, TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *multiCmesh) {

    std::stringstream ref;
    ref << "_ref-" << 1/preConfig.h <<" x " << 1/preConfig.h;
    std::string refinement =  ref.str();

    std::ofstream out(preConfig.plotfile + "/gmesh"+ refinement + ".vtk");
    std::ofstream out2(preConfig.plotfile + "/gmesh"+ refinement + "txt");
    std::ofstream out3(preConfig.plotfile + "/cmesh.txt");

    TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, out);
    config.gmesh->Print(out2);

    if (preConfig.mode == 0) cmesh->Print(out3);
    else multiCmesh->Print(out3);
}

void CreateHybMixCompMesh(TPZMultiphysicsCompMesh *multiHyb, TPZMultiphysicsCompMesh *multiMix, TPZMultiphysicsCompMesh *multHybMix,PreConfig &hybConfig, ProblemConfig &ConfHyb){

    InsertMaterialMixHyb(multHybMix,hybConfig,ConfHyb);

    TPZManVector<TPZCompMesh *> mesh_vectors(3, 0);
    TPZManVector<int> active(3, 1);

    mesh_vectors[0] = multiHyb->MeshVector()[1]->Clone();      // HybH1 potential
    mesh_vectors[1] = multiMix->MeshVector()[1]->Clone();      // Mixed potential
    mesh_vectors[2] = multiMix->MeshVector()[0]->Clone();      // HDiv mixed flux

    multHybMix->BuildMultiphysicsSpace(active,mesh_vectors);
    multHybMix->LoadReferences();
    multHybMix->InitializeBlock();

    //int64_t nelem = multHybMix->NElements();
    //multHybMix->LoadSolution(multHybMix->Solution());
    //multHybMix->ExpandSolution();
    //multHybMix->ElementSolution().Redim(nelem, 8 - 1);

    std::ofstream outHybP("hybP.txt");
    std::ofstream outMixP("mixP.txt");
    std::ofstream outMixF("mixF.txt");
    multHybMix->MeshVector()[0]->Print(outHybP);
    multHybMix->MeshVector()[1]->Print(outMixP);
    multHybMix->MeshVector()[2]->Print(outMixF);

    std::ofstream outOrigHybP("OrigHybP.txt");
    std::ofstream outOrigMixP("OrigMixP.txt");
    std::ofstream outOrigMixF("OrigMixF.txt");
    multiHyb->MeshVector()[1]->Print(outOrigHybP);
    multiMix->MeshVector()[1]->Print(outOrigMixP);
    multiMix->MeshVector()[0]->Print(outOrigMixF);

    std::ofstream outMult("multHybMix.txt");
    multHybMix->Print(outMult);
    std::ofstream outGeo("geomeshTest.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(multHybMix->Reference(), outGeo);
}

void CreateMixedComputationalMesh(TPZMultiphysicsCompMesh *cmesh_Mixed, PreConfig &pConfig, ProblemConfig &config){
    InsertMaterialMixed(cmesh_Mixed, config,pConfig);
    TPZManVector<TPZCompMesh *, 4> meshvector(4);
    CreateMixedAtomicMeshes(meshvector,pConfig, config);

    TPZManVector<int> active(4, 1);

    cmesh_Mixed->BuildMultiphysicsSpace(active,meshvector);
    CreateCondensedMixedElements(cmesh_Mixed);
    cmesh_Mixed->LoadReferences();
    cmesh_Mixed->InitializeBlock();
}
void CreateCondensedMixedElements(TPZMultiphysicsCompMesh *cmesh_Mixed){

    int numActiveSpaces = cmesh_Mixed->MeshVector().size();
    if(numActiveSpaces == 4){
        int64_t nconnects = cmesh_Mixed->NConnects();
        for (int64_t ic = 0; ic<nconnects; ic++) {
            TPZConnect &c = cmesh_Mixed->ConnectVec()[ic];
            if(c.LagrangeMultiplier() == 4) c.IncrementElConnected();
        }
    }
    bool keeponelagrangian = true, keepmatrix = false;
    TPZCompMeshTools::CreatedCondensedElements(cmesh_Mixed, keeponelagrangian, keepmatrix);
}

void CreateHybridH1ComputationalMesh(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,int &interFaceMatID,int &fluxMatID , PreConfig &pConfig, ProblemConfig &config,int hybridLevel){
    auto spaceType = TPZCreateMultiphysicsSpace::EH1Hybrid;
    cmesh_H1Hybrid->SetAllCreateFunctionsMultiphysicElem();
    if (hybridLevel == 2) {
        spaceType = TPZCreateMultiphysicsSpace::EH1HybridSquared;
    }
    else if(hybridLevel != 1) {
        DebugStop();
    }

    TPZCreateMultiphysicsSpace createspace(config.gmesh, spaceType);
    //TPZCreateMultiphysicsSpace createspace(config.gmesh);
    std::cout << cmesh_H1Hybrid->NEquations();
    (pConfig.type == 2) ? 
        createspace.SetMaterialIds({2,3}, {-6,-5}) :
        createspace.SetMaterialIds({1,}, {-2,-1});

    createspace.fH1Hybrid.fHybridizeBCLevel = 1;//opcao de hibridizar o contorno
    createspace.ComputePeriferalMaterialIds();

    TPZManVector<TPZCompMesh *> meshvec;

    int pOrder = config.n+config.k;
    createspace.CreateAtomicMeshes(meshvec,pOrder,config.k);

    InsertMaterialHybrid(cmesh_H1Hybrid, config,pConfig);
    createspace.InsertPeriferalMaterialObjects(cmesh_H1Hybrid);
    cmesh_H1Hybrid->BuildMultiphysicsSpace(meshvec);
    createspace.InsertLagranceMaterialObjects(cmesh_H1Hybrid);

    {
        int64_t nel = cmesh_H1Hybrid->NElements();
        for(int64_t el = 0; el<nel; el++) {
            auto cel = dynamic_cast<TPZMultiphysicsElement *>(cmesh_H1Hybrid->Element(el));
            cel->InitializeIntegrationRule();
            auto gel = cel->Reference();
            int dim = gel->Dimension();
            TPZManVector<int,3> order(dim);
            cel->GetIntegrationRule().GetOrder(order);
            for(int i=0; i<dim; i++) order[i]+=1;
            cel->GetIntegrationRule().SetOrder(order);
        }
    }
    
    createspace.AddInterfaceElements(cmesh_H1Hybrid);
    
#ifdef PZDEBUG
    {
        std::ofstream out("mphysics.txt");
        cmesh_H1Hybrid->Print(out);
    }
#endif
    createspace.GroupandCondenseElements(cmesh_H1Hybrid);

    cmesh_H1Hybrid->InitializeBlock();
    cmesh_H1Hybrid->ComputeNodElCon();

    interFaceMatID = createspace.fH1Hybrid.fLagrangeMatid.first;
    fluxMatID = createspace.fH1Hybrid.fFluxMatId;
}

void SolveH1Problem(TPZCompMesh *cmeshH1,struct ProblemConfig &config, struct PreConfig &pConfig){

    config.exact.operator*().fSignConvention = 1;

    std::cout << "Solving H1 " << std::endl;

    TPZLinearAnalysis an(cmeshH1);
    //an.SetExact(config.exact->ExactSolution());
    
#ifdef PZ_USING_MKL
     TPZSSpStructMatrix<STATE> strmat(cmeshH1);
    strmat.SetNumThreads(0);
    //        strmat.SetDecomposeType(ELDLt);
#else
//    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmeshH1);
//    strmat.SetNumThreads(0);
            TPZSkylineStructMatrix<STATE> strmat(cmeshH1);
            strmat.SetNumThreads(0);
#endif

    std::set<int> matids;
    for (auto matid : config.materialids) matids.insert(matid);

    for(auto mat:config.bcmaterialids){
        matids.insert(mat);
    }

    strmat.SetMaterialIds(matids);
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();//resolve o problema ate aqui

    int64_t nelem = cmeshH1->NElements();
    cmeshH1->LoadSolution(cmeshH1->Solution());
    cmeshH1->ExpandSolution();
    cmeshH1->ElementSolution().Redim(nelem, 10);

    ////Calculo do erro
    std::cout << "Computing Error H1 " << std::endl;

    an.SetExact(config.exact.operator*().ExactSolution());

    StockErrorsH1(an,cmeshH1,pConfig.Erro,pConfig.Log,pConfig);

    ////PostProcess
    if(pConfig.debugger) {
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Solution");
        vecnames.Push("Derivative");
        scalnames.Push("ExactSolution");
        vecnames.Push("ExactFlux");
        int dim = cmeshH1->Reference()->Dimension();

        std::string plotname;
        {
            std::stringstream out;
            out << pConfig.plotfile << "/" << config.problemname <<"_" << pConfig.topologyFileName << "_k-" << config.k << "_n-" << config.n;

            if(dim == 2) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";
            if(dim == 3) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";

            plotname = out.str();
        }
        int resolution=3;
        
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(resolution,dim);
    }
    std::cout << "FINISHED!" << std::endl;
}

void SolveHybridH1Problem(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,int InterfaceMatId, struct ProblemConfig &config,struct PreConfig &pConfig,int hybridLevel){

    config.exact.operator*().fSignConvention = 1;
    std::cout << "\nStarting Primal Hybrid FEM solving step" << std::endl;

    TPZLinearAnalysis an(cmesh_H1Hybrid);

#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> strmat(cmesh_H1Hybrid);
    strmat.SetNumThreads(0);
    //        strmat.SetDecomposeType(ELDLt);
#else
    TPZFStructMatrix<STATE> strmat(cmesh_H1Hybrid);
    
//    TPZSpStructMatrix<STATE> strmat(cmesh_H1Hybrid);
//    strmat.SetNumThreads(0);
//        TPZFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmesh_H1Hybrid);
//        strmat.SetNumThreads(2);
//        strmat.SetDecomposeType(ELDLt);
//    TPZSkylineStructMatrix<STATE> strmat(cmesh_H1Hybrid);
//    strmat.SetNumThreads(0);
#endif
    std::set<int> matIds;
    for (auto matid : config.materialids) matIds.insert(matid);
    for (auto matidbc : config.bcmaterialids) matIds.insert(matidbc);

    matIds.insert(InterfaceMatId);
    strmat.SetMaterialIds(matIds);
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE>* direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();

    int64_t nelem = cmesh_H1Hybrid->NElements();
    cmesh_H1Hybrid->LoadSolution(cmesh_H1Hybrid->Solution());
    cmesh_H1Hybrid->ExpandSolution();
    cmesh_H1Hybrid->ElementSolution().Redim(nelem, 5);

    std::cout << "DOF = " << cmesh_H1Hybrid->NEquations() << std::endl;
    if(pConfig.debugger) {
        std::cout << "\n############\n";
        std::cout << "Computing Error HYBRID_H1" << std::endl;
        an.SetExact(config.exact.operator*().ExactSolution());

        ////Calculo do erro
        StockErrors(an,cmesh_H1Hybrid,pConfig.Erro,pConfig.Log,pConfig);
        ////PostProcess
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Pressure");
        scalnames.Push("PressureExact");
        scalnames.Push("Kperm");
        vecnames.Push("Flux");
        vecnames.Push("ExactFlux");

        int dim = pConfig.dim;
        std::string plotname;
        {
            std::stringstream out;
            out << pConfig.plotfile << "/" << config.problemname <<"_" << pConfig.topologyFileName << "_k-" << config.k << "_n-" << config.n;

            if(dim == 2) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";
            if(dim == 3) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";

            plotname = out.str();
        }
        int resolution = 3;
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(resolution, dim);
    }
}

void SolveMixedProblem(TPZMultiphysicsCompMesh *cmesh_Mixed,struct ProblemConfig config,struct PreConfig &pConfig) {

    config.exact.operator*().fSignConvention = 1;
    RenumType optBW = RenumType::ESloan;

    std::cout << "Solving Mixed " << std::endl;
    TPZLinearAnalysis an(cmesh_Mixed, optBW); //Cria objeto de análise que gerenciará a analise do problema
    if(false){
        cout<<"Total ecuaciones:"<<an.Solution().Rows()<<endl;
    }
    //MKL solver
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> strmat(cmesh_Mixed);
    //strmat.SetNumThreads(8);
    strmat.SetNumThreads(0);
#else
    TPZSkylineStructMatrix<STATE> strmat(cmesh_Mixed);
    strmat.SetNumThreads(0);
#endif
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE>* direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
#ifdef ERRORESTIMATION_DEBUG2
    const string matrixNamevtk("matrixRigidezMixedProblem.vtk");
    TPZMatrix<REAL> * matrizRigidez = an.Solver().Matrix().operator->();
    //VisualMatrixVTK((TPZFMatrix<REAL>&)(*matrizRigidez),matrixNamevtk);
#endif
    an.Solve();

    ////PostProcess
    if(pConfig.debugger) {

        ////Calculo do erro
        std::cout << "Computing Error MIXED " << std::endl;

        an.SetExact(config.exact.operator*().ExactSolution());

        std::cout << "DOF = " << cmesh_Mixed->NEquations() << std::endl;

        StockErrors(an,cmesh_Mixed,pConfig.Erro,pConfig.Log,pConfig);

        int dim = config.gmesh->Dimension();
        std::string plotname;
        {
            std::stringstream out;
            out << pConfig.plotfile << "/" << config.problemname <<"_" << pConfig.topologyFileName << "_k-" << config.k << "_n-" << config.n;

            if(dim == 2) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";
            if(dim == 3) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";

            plotname = out.str();
        }

        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Pressure");
        scalnames.Push("ExactPressure");
        vecnames.Push("Flux");
        vecnames.Push("ExactFlux");

        int resolution = 0;
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(resolution, dim);
    }
}

void StockErrorsH1(TPZAnalysis &an,TPZCompMesh *cmesh, ofstream &Erro, TPZVec<REAL> *Log,PreConfig &pConfig){

    TPZManVector<REAL,6> Errors;
    Errors.resize(pConfig.numErrors);
    bool store_errors = true;

    an.PostProcessError(Errors, store_errors, Erro);

    if ((*Log)[0] != -1) {
        for (int j = 0; j < 3; j++) {
            (*pConfig.rate)[j] =
                    (log10(Errors[j]) - log10((*Log)[j])) /
                    (log10(pConfig.h) - log10(pConfig.hLog));
            Erro << "rate " << j << ": " << (*pConfig.rate)[j] << std::endl;
        }
    }

    Erro << "h = " << pConfig.h << std::endl;
    Erro << "DOF = " << cmesh->NEquations() << std::endl;
    for (int i = 0; i < pConfig.numErrors; i++)
        (*Log)[i] = Errors[i];
    Errors.clear();
}

void StockErrors(TPZAnalysis &an,TPZMultiphysicsCompMesh *cmesh, ofstream &Erro, TPZVec<REAL> *Log,PreConfig &pConfig){

    TPZManVector<REAL,6> Errors;
    Errors.Fill(0);
    Errors.resize(pConfig.numErrors);
    bool store_errors = false;

    an.PostProcessError(Errors, store_errors, Erro);

    std::cout << "||u_h-u||:                   \t" <<Errors[0] << 
               "\n||Grad(u_h)-Grad(u)||:       \t" << Errors[1]<< 
               "\nError[0]+Error[1]:           \t"<< Errors[2]<<
               "\n||K^0.5(Grad(u_h)-Grad(u))||:\t"<< Errors[3]<<"\n\n";

    if ((*Log)[0] != -1) {
        for (int j = 0; j < 3; j++) {
            (*pConfig.rate)[j] =
                    (log10(Errors[j]) - log10((*Log)[j])) /
                    (log10(pConfig.h) - log10(pConfig.hLog));
            Erro << "rate " << j << ": " << (*pConfig.rate)[j] << std::endl;
        }
    }

    Erro << "h = " << pConfig.h << std::endl;
    Erro << "DOF = " << cmesh->NEquations() << std::endl;
    for (int i = 0; i < pConfig.numErrors; i++)
        (*Log)[i] = Errors[i];
    Errors.clear();
}

void FluxErrorCreateCompMesh(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,int &interFaceMatID,int &fluxMatID , PreConfig &pConfig, ProblemConfig &config){
    auto spaceType = TPZCreateMultiphysicsSpace::EH1Hybrid;

    TPZCreateMultiphysicsSpace createspace(config.gmesh, spaceType);

    std::cout << cmesh_H1Hybrid->NEquations();

    createspace.SetMaterialIds({1}, {-2,-1});
    createspace.fH1Hybrid.fHybridizeBCLevel = 1;//opcao de hibridizar o contorno
    createspace.ComputePeriferalMaterialIds();

    TPZManVector<TPZCompMesh *> meshvec;

    int pOrder = config.n+config.k;
    createspace.CreateAtomicMeshes(meshvec,pOrder,config.k);

    FluxErrorInsertMaterial(cmesh_H1Hybrid, config,pConfig);
    createspace.InsertPeriferalMaterialObjects(cmesh_H1Hybrid);
    cmesh_H1Hybrid->BuildMultiphysicsSpace(meshvec);
    createspace.InsertLagranceMaterialObjects(cmesh_H1Hybrid);

    createspace.AddInterfaceElements(cmesh_H1Hybrid);
    createspace.GroupandCondenseElements(cmesh_H1Hybrid);

    cmesh_H1Hybrid->InitializeBlock();
    cmesh_H1Hybrid->ComputeNodElCon();

    interFaceMatID = createspace.fH1Hybrid.fLagrangeMatid.first;
    fluxMatID = createspace.fH1Hybrid.fFluxMatId;
}

bool PostProcessing(TPZCompMesh * cmeshH1, TPZFMatrix<STATE> true_elerror, TPZFMatrix<STATE> estimate_elerror, ProblemConfig &config) {
    
    TPZLinearAnalysis an(cmeshH1);
    
    int64_t nels = cmeshH1->ElementVec().NElements();
    cmeshH1->ElementSolution().Redim(nels, 6);
    
    //Compute the effectivity index
    std::cout<<"***** Computing effectivity index *****"<<std::endl;
    
    {
        REAL sum = 0.;
        int64_t nel = true_elerror.Rows();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = cmeshH1->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            TPZCompEl *mphys = gel->Reference();
            int64_t elindex2 = mphys->Index();
            REAL aux =  estimate_elerror(elindex2,2) + estimate_elerror(elindex2,3);
            sum += aux*aux; // To compute global effectivity index
            true_elerror(el,0) = estimate_elerror(elindex2,2);
            true_elerror(el,1) = true_elerror(el,2);
            if (true_elerror(el,1) > 1.e-10) {
                true_elerror(el,2) = true_elerror(el,0)/true_elerror(el,1);
            }
        }
        
        cmeshH1->ElementSolution() = true_elerror;
        
        TPZManVector<REAL> errorsum(5, 0.);
        cmeshH1->EvaluateError(false, errorsum);
        REAL globeffind = sqrt(sum)/errorsum[2];
        std::cout << "Global full error estimate: "<<sqrt(sum)<<std::endl;
        std::cout << "Global effectivity index: " << globeffind << std::endl;
        
        std::ofstream fileouput;
        fileouput.open("adaptivityresults.txt",ios::app);
        fileouput << std::setw(15) << cmeshH1->NEquations();
        fileouput << std::setw(15) << errorsum[2];
        fileouput << std::setw(15) << sqrt(sum);
        fileouput << std::setw(15) << globeffind << std::endl;
        fileouput.close();
    }
    
    {
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Solution");
        scalnames.Push("ExactSolution");
        scalnames.Push("EstimatedError");
        scalnames.Push("TrueError");
        scalnames.Push("POrder");
        scalnames.Push("EffectivityIndex");
        vecnames.Push("Flux");
        vecnames.Push("ExactFlux");

        std::string plotname;
        {   //Creating the directories for the output
            std::string foldername = "ErrorEstimate/";
            std::string command = "mkdir -p " + foldername;
            system(command.c_str());
            
            std::stringstream ss;
            ss << "__p-" << config.k;
            REAL threshold = config.division_threshold;
            int th = (int)(100.*threshold);
            if(config.division_threshold!=-1)
                 ss << "__thr-"<< th;
            std::string problemName = config.problemname;
            problemName += ss.str();
            std::string command2 = "mkdir -p " + foldername + problemName;
            system(command2.c_str());
            
            // Enabling file name
            std::stringstream out;
            out << "ErrorEstimate/" + problemName + "/";
            out <<"ErrorEstimationH1_" << cmeshH1->GetDefaultOrder() << "_" << cmeshH1->Reference()->Dimension()
            <<  "Ndofs " << cmeshH1->NEquations() << ".vtk";
            plotname = out.str();
        }
        
        //an.SetExact(config.exact->ExactSolution());
        //an.SetExact(config.exact.operator*().ExactSolution());
        an.DefineGraphMesh(cmeshH1->Dimension(), scalnames, vecnames, plotname);
        int resolution = 0;
        an.PostProcess(resolution);
    }

    return true;
}
