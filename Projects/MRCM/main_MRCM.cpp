//
//  New_main.cpp
//  ErrorEstimateHDiv
//
//  Created by Denise De Siqueira on 01/04/19.
//

#include "pzlog.h"
#include "TPZRefPatternDataBase.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"

#include "tpzquadraticquad.h"
#include "tpzquadratictrig.h"
#include "tpzgeoelrefpattern.h"
#include "tpzarc3d.h"


#include "ProblemConfig.h"

#include "mixedpoisson.h"
#include "TPZMatLaplacianHybrid.h"
#include "TPZVecL2.h"
#include "pzbndcond.h"
#include "TPZNullMaterial.h"
#include "TPZLagrangeMultiplier.h"

#include "pzintel.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZHybridizeHDiv.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZMultiphysicsCompMesh.h"

#include "TPZHDivErrorEstimatorH1.h"

#include "Tools.h"

#include "TPZBFileStream.h"

#include "TPZCreateMultiphysicsSpace.h"
#include <tuple>
#include <memory>




bool neumann = true;
bool h1solution = false;
bool hybridh1 = false;
bool mcrn = true;

void InsertFluxMultipliers(TPZHybridizeHDiv &hybridize, TPZVec<TPZCompMesh *> &meshvec,
                           int lower_order);

void InsertMaterialObjectsH1Hybrid(TPZMultiphysicsCompMesh *cmesh, ProblemConfig &config);
void SolveH1Problem(TPZCompMesh *cmeshH1,struct ProblemConfig &config);
//TPZCompMesh *CMeshH1(const ProblemConfig &problem);
void SolveHybridH1Problem(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,int InterfaceMatId,struct ProblemConfig config);
void SolveMixedProblem(TPZMultiphysicsCompMesh *cmesh_Mixed, ProblemConfig &config);

void CreateMixedComputationalMesh(TPZMultiphysicsCompMesh *cmesh_Mixed, ProblemConfig &config);
void SwitchLagrangeMultipliers(TPZCompMesh *cmesh, ProblemConfig &config);


int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    for (int ndiv = 0; ndiv < 3; ndiv++) {
        
        ProblemConfig config;
        // Internal hDiv order
        config.porder = 4;
        config.ndivisions = ndiv;
        config.dimension = 2;
        config.prefine = false;
        // Penalty term to make the normal component of the fluxes compatible
        config.fMRCMBeta = 10;
        config.hdivmais = 0;
        
        int orderlagrange = 1;
        config.H1Hybridminus = config.porder-orderlagrange;
        
        config.exact = new TLaplaceExample1;
        config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
        
        config.problemname = "ESinSin k=1 e lagrange order 1";
        
        config.dir_name = "MRCM_SINSIN";
        std::string command = "mkdir " + config.dir_name;
        system(command.c_str());
        
        // geometric mesh
        TPZManVector<int, 4> bcids(4, -1);
        TPZGeoMesh *gmesh = CreateGeoMesh(2, bcids);
        //    if(1)
        //    {
        //        TPZManVector<TPZGeoEl *> sub;
        //        gmesh->Element(0)->Divide(sub);
        //        DivideLowerDimensionalElements(gmesh);
        //    }
        
        config.gmesh = gmesh;
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        
        
        
        UniformRefinement(config.ndivisions, gmesh);
        //        int refinement_depth = 2;
        //        RandomRefine(config, 1, refinement_depth);
        
#ifdef PZDEBUG
        {
            std::ofstream out("gmesh.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
            std::ofstream out2("gmeshInitial.txt");
            gmesh->Print(out2);
            
        }
#endif
        
        //problema H1
        if(h1solution) {
            
            config.exact.operator*().fSignConvention = -1;
            TPZCompMesh *cmeshH1 = CMeshH1(config);
            {
                ofstream arg1("CompMeshH1.txt");
                cmeshH1->Print(arg1);
            }
            SolveH1Problem(cmeshH1, config);
        }
        
        if(hybridh1){
            config.exact.operator*().fSignConvention = 1;
            TPZCreateMultiphysicsSpace createspace(gmesh,TPZCreateMultiphysicsSpace::EH1HybridSquared);
            
            createspace.SetMaterialIds({1}, {-2,-1});
            createspace.fH1Hybrid.fHybridizeBCLevel = 2;//opcao de hibridizar o contorno
            createspace.ComputePeriferalMaterialIds();
            
            
            std::cout<<"---Original PerifericalMaterialId --- "<<std::endl;
            std::cout <<" fMatWrapId + = "<<createspace.fH1Hybrid.fMatWrapId<<std::endl;
            std::cout <<" fLagrangeMatid + = "<<createspace.fH1Hybrid.fLagrangeMatid.first<<std::endl;
            std::cout <<" fLagrangeMatid - = "<<createspace.fH1Hybrid.fLagrangeMatid.second<<std::endl;
            std::cout <<" fFluxMatId = "<<createspace.fH1Hybrid.fFluxMatId<<std::endl;
            std::cout << "fSecond Lagrange MatID = " <<createspace.fH1Hybrid.fSecondLagrangeMatid<<std::endl;
            std::cout << "fInterfacePressure = " <<createspace.fH1Hybrid.fInterfacePressure<<std::endl;
            std::cout << "fIBC hybridization level = " <<createspace.fH1Hybrid.fHybridizeBCLevel<<std::endl;
            
            TPZManVector<TPZCompMesh *> meshvec;
            
            createspace.CreateAtomicMeshes(meshvec,config.porder,orderlagrange);
            
            TPZMultiphysicsCompMesh *cmesh_H1Hybrid = new TPZMultiphysicsCompMesh(gmesh);
            InsertMaterialObjectsH1Hybrid(cmesh_H1Hybrid, config);
            createspace.InsertPeriferalMaterialObjects(cmesh_H1Hybrid);
            
#ifdef PZDEBUG
            {
                std::ofstream out3("pressure.txt");
                meshvec[0]->Print(out3);
                std::ofstream out4("flux.txt");
                meshvec[1]->Print(out4);
                std::ofstream out("gmeshIncremented.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
                std::ofstream out2("gmeshIncremented.txt");
                gmesh->Print(out2);
                
            }
#endif
            cmesh_H1Hybrid->BuildMultiphysicsSpace(meshvec);
            createspace.InsertLagranceMaterialObjects(cmesh_H1Hybrid);
            createspace.AddInterfaceElements(cmesh_H1Hybrid);
#ifdef PZDEBUG
            {
                std::map<int,int> matelem;
                int64_t nel = cmesh_H1Hybrid->NElements();
                for (int64_t el = 0; el<nel; el++) {
                    TPZCompEl *cel = cmesh_H1Hybrid->Element(el);
                    TPZGeoEl *gel = cel->Reference();
                    matelem[gel->MaterialId()]++;
                }
                std::cout << __PRETTY_FUNCTION__ << " number of computational elements by material \n";
                for (auto it : matelem) {
                    std::cout << "Material id " << it.first << " number of elements " << it.second << std::endl;
                }
            }
#endif
            cmesh_H1Hybrid->ComputeNodElCon();
#ifdef PZDEBUG
            {
                std::ofstream out("mphysicsmeshBeforeCondense.txt");
                cmesh_H1Hybrid->Print(out);
            }
#endif
            createspace.GroupandCondenseElements(cmesh_H1Hybrid);
            
            cmesh_H1Hybrid->InitializeBlock();
            cmesh_H1Hybrid->ComputeNodElCon();
            
#ifdef PZDEBUG
            {
                std::ofstream out("mphysicsmesh.txt");
                cmesh_H1Hybrid->Print(out);
            }
#endif
            //Solve Hybrid problem
            
            SolveHybridH1Problem(cmesh_H1Hybrid,createspace.fH1Hybrid.fLagrangeMatid.first,config);
            
            //Post Processing for Lagrange Multiplier
            //    {
            //        TPZAnalysis an(meshvec[1],false);
            //
            //        TPZStack<std::string> scalnames, vecnames;
            //        scalnames.Push("State");
            //
            //        int dim = 1;
            //        std::string plotname("LagrangeMultiplier.vtk");
            //        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
            //        an.PostProcess(0, dim);
            //
            //
            //    }
            
#ifdef PZDEBUG
            {
                
                std::ofstream out2("H1HybridMesh.txt");
                cmesh_H1Hybrid->Print(out2);
                std::ofstream out3("gmeshHybridH1.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out3);
                
            }
#endif
        }
        
        if(mcrn)
        {
            TPZMultiphysicsCompMesh *cmesh_MCRN = new TPZMultiphysicsCompMesh(gmesh);
            config.exact.operator*().fSignConvention = 1;
            CreateMixedComputationalMesh(cmesh_MCRN, config);
            TPZVec<TPZCompMesh *> meshvec = cmesh_MCRN->MeshVector();
//            TPZCreateMultiphysicsSpace createspace(gmesh,TPZCreateMultiphysicsSpace::EH1HybridSquared);
//
//            createspace.SetMaterialIds({1}, {-2,-1});
//            createspace.fH1Hybrid.fHybridizeBCLevel = 2;//opcao de hibridizar o contorno
//            createspace.ComputePeriferalMaterialIds();
//
//
//            std::cout<<"---Original PerifericalMaterialId --- "<<std::endl;
//            std::cout <<" fMatWrapId + = "<<createspace.fH1Hybrid.fMatWrapId<<std::endl;
//            std::cout <<" fLagrangeMatid + = "<<createspace.fH1Hybrid.fLagrangeMatid.first<<std::endl;
//            std::cout <<" fLagrangeMatid - = "<<createspace.fH1Hybrid.fLagrangeMatid.second<<std::endl;
//            std::cout <<" fFluxMatId = "<<createspace.fH1Hybrid.fFluxMatId<<std::endl;
//            std::cout << "fSecond Lagrange MatID = " <<createspace.fH1Hybrid.fSecondLagrangeMatid<<std::endl;
//            std::cout << "fInterfacePressure = " <<createspace.fH1Hybrid.fInterfacePressure<<std::endl;
//            std::cout << "fIBC hybridization level = " <<createspace.fH1Hybrid.fHybridizeBCLevel<<std::endl;
//
//            TPZManVector<TPZCompMesh *> meshvec;
//
//            createspace.CreateAtomicMeshes(meshvec,config.porder,orderlagrange);
//
//            InsertMaterialObjectsH1Hybrid(cmesh_H1Hybrid, config);
//            createspace.InsertPeriferalMaterialObjects(cmesh_H1Hybrid);
            
#ifdef PZDEBUG
            {
                meshvec[0]->ComputeNodElCon();
                std::ofstream out3("flux.txt");
                meshvec[0]->Print(out3);
                std::ofstream out4("pressure.txt");
                meshvec[1]->ComputeNodElCon();
                meshvec[1]->Print(out4);
                std::ofstream out("gmeshIncremented.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
                std::ofstream out2("gmeshIncremented.txt");
                gmesh->Print(out2);
                
            }
#endif
//            cmesh_H1Hybrid->BuildMultiphysicsSpace(meshvec);
//            createspace.InsertLagranceMaterialObjects(cmesh_H1Hybrid);
//            createspace.AddInterfaceElements(cmesh_H1Hybrid);
#ifdef PZDEBUG
            {
                std::map<int,int> matelem;
                int64_t nel = cmesh_MCRN->NElements();
                for (int64_t el = 0; el<nel; el++) {
                    TPZCompEl *cel = cmesh_MCRN->Element(el);
                    if(!cel) continue;
                    TPZGeoEl *gel = cel->Reference();
                    if(!gel) continue;
                    matelem[gel->MaterialId()]++;
                }
                std::cout << __PRETTY_FUNCTION__ << " number of computational elements by material \n";
                for (auto it : matelem) {
                    std::cout << "Material id " << it.first << " number of elements " << it.second << std::endl;
                }
            }
#endif
            cmesh_MCRN->ComputeNodElCon();
//            createspace.GroupandCondenseElements(cmesh_H1Hybrid);
            
            cmesh_MCRN->InitializeBlock();
            cmesh_MCRN->ComputeNodElCon();
            std::cout << "Beta " << config.fMRCMBeta << std::endl;
            std::cout << "order " << config.porder << std::endl;
            std::cout << "lagrange order " << orderlagrange << std::endl;
            std::cout << "HDivMais " << config.hdivmais << std::endl;
            SolveMixedProblem(cmesh_MCRN, config);
        }
        
    }
    
    return 0.;
    
}

void InsertMaterialObjectsH1Hybrid(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, ProblemConfig &config)
{
    TPZGeoMesh *gmesh = cmesh_H1Hybrid->Reference();
    int dim = gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;
    
    // Creates Poisson material
    TPZMatLaplacianHybrid *material = new TPZMatLaplacianHybrid(matID, dim);
    
    cmesh_H1Hybrid->InsertMaterialObject(material);
    if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
        material->SetForcingFunction(
                                     config.exact.operator*().ForcingFunction());
        material->SetForcingFunctionExact(config.exact.operator*().Exact());
    }
    //    TPZMaterial * mat(material);
    //    cmesh->InsertMaterialObject(mat);
    
    // Inserts boundary conditions
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 1.);
    TPZMaterial *BCond0 =
    material->CreateBC(material, -1, dirichlet, val1, val2);
    if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
        BCond0->SetForcingFunction(config.exact.operator*().Exact());
    }
    val2.Zero();
    TPZMaterial *BCond1 = material->CreateBC(material, -2, neumann, val1, val2);
    
    cmesh_H1Hybrid->InsertMaterialObject(BCond0);
    cmesh_H1Hybrid->InsertMaterialObject(BCond1);
}

//TPZCompMesh *CMeshH1(const ProblemConfig &problem) {
//    
//    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
//    TPZMaterial *mat = 0;
//    
//    
//    for (auto matid : problem.materialids) {
//        TPZMatPoisson3d *mix = new TPZMatPoisson3d(matid, cmesh->Dimension());
//        mix->SetForcingFunctionExact(problem.exact.Exact());
//        mix->SetForcingFunction(problem.exact.ForcingFunction());
//        
//        if (!mat) mat = mix;
//        cmesh->InsertMaterialObject(mix);
//        
//    }
//    
//    for (auto matid : problem.bcmaterialids) {
//        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
//        int bctype = 0;
//        val2.Zero();
//        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
//        bc->TPZMaterial::SetForcingFunction(problem.exact.Exact());
//        
//        cmesh->InsertMaterialObject(bc);
//    }
//    
//    cmesh->SetDefaultOrder(problem.porder);//ordem
//    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
//    
//    cmesh->AutoBuild();
//    
//
//    return cmesh;
//}

void SolveH1Problem(TPZCompMesh *cmeshH1,struct ProblemConfig &config){
    
    TPZAnalysis an(cmeshH1);
    
    
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmeshH1);
    strmat.SetNumThreads(0);
    //        strmat.SetDecomposeType(ELDLt);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmeshH1);
    strmat.SetNumThreads(0);
    //        TPZSkylineStructMatrix strmat3(cmesh_HDiv);
    //        strmat3.SetNumThreads(8);
#endif
    
    std::set<int> matids;
    matids.insert(1);
    
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
    an.Solve();//resolve o problema misto ate aqui
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Solution");
    vecnames.Push("Derivative");
    scalnames.Push("ExactSolution");
    
    
    
    
    
    int dim = cmeshH1->Reference()->Dimension();
    
    std::string plotname;
    {
        std::stringstream out;
        out << config.dir_name << "/" << "H1_Problem" << config.porder << "_" << dim
        << "D_" << config.problemname << "Ndiv_ " << config.ndivisions << ".vtk";
        plotname = out.str();
    }
    int resolution=0;
    an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
    an.PostProcess(resolution,dim);
    
    an.SetExact(config.exact.operator*().ExactSolution());
    
    TPZManVector<REAL> errorvec(10, 0.);
    int64_t nelem = cmeshH1->NElements();
    cmeshH1->LoadSolution(cmeshH1->Solution());
    cmeshH1->ExpandSolution();
    cmeshH1->ElementSolution().Redim(nelem, 10);
    
    an.PostProcessError(errorvec);//calculo do erro com sol exata e aprox
    
    std::cout << "Computed errors " << errorvec << std::endl;
    
    
    //Erro
    
    ofstream myfile;
    myfile.open("ArquivosErrosH1.txt", ios::app);
    myfile << "\n\n Error for H1 formulation " << config.problemname;
    myfile << "\n-------------------------------------------------- \n";
    myfile << "Ndiv = " << config.ndivisions << " Order = " << config.porder << "\n";
    myfile << "DOF Total = " << cmeshH1->NEquations() << "\n";
    myfile << "Energy norm = " << errorvec[0] << "\n";//norma energia
    myfile << "error norm L2 = " << errorvec[1] << "\n";//norma L2
    myfile << "Semi norm H1 = " << errorvec[2] << "\n";//norma L2
    myfile.close();
    
    
}

void SolveHybridH1Problem(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,int InterfaceMatId,struct ProblemConfig config){
    
    TPZAnalysis an(cmesh_H1Hybrid);
    
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh_H1Hybrid);
    strmat.SetNumThreads(0);
    //        strmat.SetDecomposeType(ELDLt);
#else
    //    TPZFrontStructMatrix<TPZFrontSym<STATE> > strmat(Hybridmesh);
    //    strmat.SetNumThreads(2);
    //    strmat.SetDecomposeType(ELDLt);
    TPZSkylineStructMatrix strmat(cmesh_H1Hybrid);
    strmat.SetNumThreads(0);
#endif
    
    
    std::set<int> matIds;
    
    
    for (auto matid : config.materialids) {
        
        matIds.insert(matid);
    }
    
    
    for (auto matidbc : config.bcmaterialids) {
        
        matIds.insert(matidbc);
    }
    
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
    
    //Pos processamento
    
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    scalnames.Push("PressureExact");
    
    int dim = 2;
    std::string plotname;
    {
        std::stringstream out;
        out << config.dir_name << "/" << "HybridH1" << config.porder << "_" << dim
        << "D_" << config.problemname << "Ndiv_ " << config.ndivisions << "HdivMais"
        << config.hdivmais << ".vtk";
        plotname = out.str();
    }
    int resolution=0;
    an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
    an.PostProcess(resolution, dim);
    
    
    //Calculo do Erro
    
    an.SetExact(config.exact.operator*().ExactSolution());
    
    TPZManVector<REAL> errorvec(5, 0.);
    int64_t nelem = cmesh_H1Hybrid->NElements();
    cmesh_H1Hybrid->LoadSolution(cmesh_H1Hybrid->Solution());
    cmesh_H1Hybrid->ExpandSolution();
    cmesh_H1Hybrid->ElementSolution().Redim(nelem, 5);
    
    an.PostProcessError(errorvec);//calculo do erro com sol exata e aprox
    
    std::cout << "Computed errors " << errorvec << std::endl;
    
    
    //Erro
    
    ofstream myfile;
    myfile.open("ArquivosErrosH1Hibrido.txt", ios::app);
    myfile << "\n\n Estimator errors for Problem " << config.problemname;
    myfile << "\n-------------------------------------------------- \n";
    myfile << "Ndiv = " << config.ndivisions << " Order = " << config.porder << "\n";
    myfile << "DOF Total = " << cmesh_H1Hybrid->NEquations() << "\n";
    myfile << "error norm L2 = " << errorvec[0] << "\n";
    myfile << "semi norm H1 = " << errorvec[1] << "\n";
    myfile << "H1 norm = " << errorvec[2] << "\n";
    myfile << "energy norm = " << errorvec[3] << "\n";
    myfile.close();
    
}

void BuildFluxMesh(TPZCompMesh *cmesh_flux, ProblemConfig &config){
    
    int dim = config.gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;
    
    int flux_order = config.porder;
    cmesh_flux->SetDefaultOrder(flux_order);
    cmesh_flux->SetDimModel(config.gmesh->Dimension());
    
    cmesh_flux->SetAllCreateFunctionsHDiv();
    
    TPZNullMaterial *material = new TPZNullMaterial(matID);
    material->SetDimension(dim);
    cmesh_flux->InsertMaterialObject(material);
    
    //(EE)Create Boundary conditions
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    //(PROPRIETARY)
    TPZMaterial *BCond0 = material->CreateBC(material,-1,dirichlet,val1,val2);
    cmesh_flux->InsertMaterialObject(BCond0);
    
    TPZMaterial *BCond1 = material->CreateBC(material,-2,neumann,val1,val2);
    cmesh_flux->InsertMaterialObject(BCond1);
    
    cmesh_flux->AutoBuild();
    cmesh_flux->InitializeBlock();
}

void BuildPotentialMesh(TPZCompMesh *cmesh_p, ProblemConfig &config){
    int matID = 1;
    int dim = config.gmesh->Dimension();
    
    int potential_order = config.porder;
    cmesh_p->SetDefaultOrder(potential_order);
    //cmesh_p->SetDefaultOrder(config.porder);
    //cmesh_p->SetDefaultOrder(config.porder + config.hdivmais);
    cmesh_p->SetDimModel(dim);
    
    cmesh_p->SetAllCreateFunctionsContinuous(); //H1 functions
    cmesh_p->ApproxSpace().CreateDisconnectedElements(true);
    
    TPZNullMaterial *material = new TPZNullMaterial(matID);
    material->SetDimension(dim);
    cmesh_p->InsertMaterialObject(material);
    
    cmesh_p->AutoBuild();
    cmesh_p->ExpandSolution();
    
    TPZAdmChunkVector<TPZConnect> &nodeIt = cmesh_p->ConnectVec();
    for(auto &nodo : nodeIt){
        nodo.SetLagrangeMultiplier(1);
    }
}

void InsertMaterialMixed(TPZMultiphysicsCompMesh *cmesh_mixed, ProblemConfig config){
    
    int dim = config.gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;
    
    cmesh_mixed -> SetDefaultOrder(config.porder);
    cmesh_mixed -> SetDimModel(dim);
    cmesh_mixed->SetAllCreateFunctionsMultiphysicElem();
    
    TPZMixedPoisson *material = new TPZMixedPoisson(matID,dim); //Using standard PermealityTensor = Identity.
    material->SetForcingFunction(config.exact.operator*().ForcingFunction());
    material->SetForcingFunctionExact(config.exact.operator*().Exact());
    cmesh_mixed->InsertMaterialObject(material);
    
    //Boundary Conditions
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    TPZMaterial *BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2);
    BCond0->SetForcingFunction(config.exact.operator*().Exact());
    
    TPZMaterial *BCond1 = material->CreateBC(material, -2, neumann, val1, val2);
    
    cmesh_mixed->InsertMaterialObject(BCond0);
    cmesh_mixed->InsertMaterialObject(BCond1);
}



void CreateMixedComputationalMesh(TPZMultiphysicsCompMesh *cmesh_Mixed, ProblemConfig &config){
    
    int matID = 1;
    int dim = config.gmesh->Dimension();
    
    //Flux mesh creation
    TPZCompMesh *cmesh_flux = new TPZCompMesh(config.gmesh);
    BuildFluxMesh(cmesh_flux, config);
    
    //Potential mesh creation
    TPZCompMesh *cmesh_p = new TPZCompMesh(config.gmesh);
    BuildPotentialMesh(cmesh_p, config);
    
    //Multiphysics mesh build
    InsertMaterialMixed(cmesh_Mixed, config);
    TPZManVector<int> active(2, 1);
    TPZManVector<TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmesh_flux;
    meshvector[1] = cmesh_p;
    TPZCompMeshTools::AdjustFluxPolynomialOrders(meshvector[0], config.hdivmais); //Increases internal flux order by "hdivmais"
    TPZCompMeshTools::SetPressureOrders(meshvector[0], meshvector[1]);//Set the pressure order the same as the internal flux
    
    
    {
        TPZHybridizeHDiv hybridize;
        hybridize.fInterfaceMatid.second = -11;
        REAL Lagrange_term_multiplier = -1.;
        TPZManVector<TPZCompMesh *, 3> meshvec_Hybrid = meshvector;
        hybridize.InsertPeriferalMaterialObjects(cmesh_Mixed, Lagrange_term_multiplier);
        int lagrangematid = hybridize.fInterfaceMatid.second;
        TPZMaterial *mat = cmesh_Mixed->FindMaterial(lagrangematid);
        TPZLagrangeMultiplier *lagrange = dynamic_cast<TPZLagrangeMultiplier *>(mat);
        lagrange->SetMultiplier(1.);
//        hybridize.InsertPeriferalMaterialObjects(meshvec_Hybrid);
        //        TPZManVector<int> active = cmesh_Mixed->GetActiveApproximationSpaces();
        hybridize.HybridizeInternalSides(meshvec_Hybrid);
        // this is where we insert the dim-1 flux elements on top of the pressure interface
        // elements
        InsertFluxMultipliers(hybridize, meshvec_Hybrid,config.H1Hybridminus);
        cmesh_Mixed->BuildMultiphysicsSpace(active,meshvec_Hybrid);
        
        hybridize.CreateInterfaceElements(cmesh_Mixed);
        SwitchLagrangeMultipliers(cmesh_Mixed,config);
    
#ifdef PZDEBUG
    {
        std::ofstream out("mphysicsmeshBeforeCondense.txt");
        cmesh_Mixed->Print(out);
    }
#endif
        int keeponelagrangian = 1;
        cmesh_Mixed->ComputeNodElCon();
        hybridize.GroupandCondenseElements(cmesh_Mixed,keeponelagrangian);
    }
    cmesh_Mixed->LoadReferences();
    cmesh_Mixed->InitializeBlock();
#ifdef PZDEBUG
    {
        std::ofstream out("mphysicsmeshAfterCondense.txt");
        cmesh_Mixed->Print(out);
    }
#endif
}

void InsertFluxMultipliers(TPZHybridizeHDiv &hybrid, TPZVec<TPZCompMesh *> &meshvec, int lower_order)
{
    TPZCompMesh *flux = meshvec[0];
    TPZGeoMesh *gmesh = flux->Reference();
    TPZNullMaterial *nullmat = new TPZNullMaterial(hybrid.fLagrangeInterface);
    flux->InsertMaterialObject(nullmat);
    std::set<int> matids;
    matids.insert(hybrid.fLagrangeInterface);
    gmesh->ResetReference();
    TPZCompMesh *pressure = meshvec[1];
    int64_t nel = pressure->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = pressure->Element(el);
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        if(!intel) continue;
        TPZGeoEl *gel = intel->Reference();
        if(gel->MaterialId() != hybrid.fLagrangeInterface) continue;
        int porder = intel->GetPreferredOrder();
        intel->PRefine(porder-lower_order);
        flux->SetDefaultOrder(porder-lower_order);
        int64_t index;
        TPZCompEl *celflux = flux->ApproxSpace().CreateCompEl(gel, *flux, index);
        celflux->Reference()->ResetReference();
    }
    flux->ExpandSolution();
}

#include "TPZMRCMLagrangeMultiplier.h"

void SwitchLagrangeMultipliers(TPZCompMesh *cmesh, ProblemConfig &config)
{
    for (auto it : cmesh->MaterialVec()) {
        TPZMaterial *mat = it.second;
        TPZLagrangeMultiplier *lagrange = dynamic_cast<TPZLagrangeMultiplier *>(mat);
        if(!lagrange) continue;
        int dim = lagrange->Dimension();
        int nstate = lagrange->NStateVariables();
        int matid = it.first;
        STATE multiplier = lagrange->Multiplier();
        TPZMRCMLagrangeMultiplier *mrcmlagrange = new TPZMRCMLagrangeMultiplier(matid,dim,nstate);
        mrcmlagrange->SetMultiplier(multiplier);
        mrcmlagrange->SetBeta(config.fMRCMBeta);
        delete it.second;
        cmesh->MaterialVec()[it.first] = mrcmlagrange;
    }
}

void SolveMixedProblem(TPZMultiphysicsCompMesh *cmesh_Mixed,struct ProblemConfig &config) {

    config.exact.operator*().fSignConvention = 1;
    bool optBW = true;

    std::cout << "Solving Mixed " << std::endl;
    TPZAnalysis an(cmesh_Mixed, optBW); //Cria objeto de análise que gerenciará a analise do problema

    //MKL solver
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh_Mixed);
    strmat.SetNumThreads(0);
    //        strmat.SetDecomposeType(ELDLt);
#else
    //    TPZFrontStructMatrix<TPZFrontSym<STATE> > strmat(Hybridmesh);
        //    strmat.SetNumThreads(2);
        //    strmat.SetDecomposeType(ELDLt);
        TPZSkylineStructMatrix strmat(cmesh_H1Hybrid);
        strmat.SetNumThreads(0);
#endif
    //std::set<int> matIds;
    //for (auto matid : config.materialids) matIds.insert(matid);
    //for (auto matidbc : config.bcmaterialids) matIds.insert(matidbc);

    //matIds.insert(InterfaceMatId);
    //strmat.SetMaterialIds(matIds);
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE>* direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();

    //Previous solver
    /*TPZSkylineStructMatrix matskl(cmesh_Mixed); //caso simetrico ***
    int numthreads = 0;
    matskl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Assemble();//Assembla a matriz de rigidez (e o vetor de carga) global
    an.Solve();*/

    ////Calculo do erro
    std::cout << "Computing Error MIXED " << std::endl;

    an.SetExact(config.exact.operator*().ExactSolution());

    std::cout << "DOF = " << cmesh_Mixed->NEquations() << std::endl;

    ////PostProcess
    {

        int dim = config.gmesh->Dimension();
        std::string plotname;
        {
            std::stringstream out;
            out << config.dir_name << "/"
                << config.problemname << "_Mixed_k-" << config.porder
                << "_Mais-" << config.hdivmais << "_ref-" << config.gmesh->NElements()
            << "_beta-" << config.fMRCMBeta << ".vtk";
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
    TPZManVector<REAL,10> errorvec;
    int64_t nel = cmesh_Mixed->NElements();
    cmesh_Mixed->ElementSolution().Redim(nel, 5);
    an.PostProcessError(errorvec);//calculo do erro com sol exata e aprox
    
    std::cout << "Computed errors " << errorvec << std::endl;

}

