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

#include "TPZMatLaplacianHybrid.h"
#include "TPZHybridElasticity2D.h"

#include "pzintel.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZHybridizeHDiv.h"

#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontSym.h"
#include "pzskylstrmatrix.h"
#include "TPZMultiphysicsCompMesh.h"

#include "Tools.h"

#include "TPZBFileStream.h"

#include "TPZCreateHybridH1Space.h"
#include <tuple>
#include <memory>


bool neumann = true;
bool h1solution = true;
bool hybridh1 = true;

std::complex<STATE> integrateF = 0;

void InsertMaterialObjectsH1Hybrid(TPZMultiphysicsCompMesh *cmesh, ProblemConfig &config);
void SolveH1Problem(TPZCompMesh *cmeshH1,struct ProblemConfig &config);
//TPZCompMesh *CMeshH1(const ProblemConfig &problem);
void SolveHybridH1Problem(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,int InterfaceMatId,struct ProblemConfig config);
using namespace std;
int main(int argc, char *argv[]) {
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    int nstate = 1;
    for (int ndiv = 0; ndiv < 1; ndiv++) {

        ProblemConfig config;

        config.porder = 3;
        config.ndivisions = ndiv;
        config.dimension = 2;
        config.prefine = false;

        int orderlagrange = 1;

//        config.exact = new TLaplaceExample1;
//        config.exact.operator*().fExact = TLaplaceExample1::EConst;
        config.exactelast = new TElasticity2DAnalytic;
        config.exactelast.operator*().fProblemType = TElasticity2DAnalytic::EBend;
        nstate = 2;

        config.problemname = "Bend k=1 e lagrange order 2";

        config.dir_name = "HybridH1_Bend";
        std::string command = "mkdir -p " + config.dir_name;
        system(command.c_str());

        // geometric mesh
        TPZManVector<int, 4> bcids(4, -1);
        TPZGeoMesh *gmesh = Tools::CreateGeoMesh(2, bcids);
        //    if(1)
        //    {
        //        TPZManVector<TPZGeoEl *> sub;
        //        gmesh->Element(0)->Divide(sub);
        //        DivideLowerDimensionalElements(gmesh);
        //    }

        config.gmesh = gmesh;
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        
        
    Tools::UniformRefinement(config.ndivisions, gmesh);
    int refinement_depth = 3;
    Tools::RandomRefinement(config.gmesh, 5, refinement_depth);
    
#ifdef ERRORESTIMATION_DEBUG
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmeshInitial.txt");
        gmesh->Print(out2);
        
    }
#endif
        
        //problema H1
    if(h1solution) {

        if(config.exact)
        {
            config.exact.operator*().fSignConvention = -1;
        }
        TPZCompMesh *cmeshH1 = Tools::CMeshH1(config);
        {
            ofstream arg1("CompMeshH1.txt");
            cmeshH1->Print(arg1);
        }
        SolveH1Problem(cmeshH1, config);
    }

    if(hybridh1){
        if(config.exact) {
            config.exact.operator*().fSignConvention = 1;
        }
        TPZCreateHybridH1Space createspace(gmesh,TPZCreateHybridH1Space::EH1Hybrid);
        createspace.fH1Hybrid.fNState = nstate;
        createspace.SetMaterialIds({1}, {-2,-1});
        createspace.fH1Hybrid.fHybridizeBCLevel = 1;//opcao de hibridizar o contorno
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
        cmesh_H1Hybrid->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;
        InsertMaterialObjectsH1Hybrid(cmesh_H1Hybrid, config);
        createspace.InsertPeriferalMaterialObjects(cmesh_H1Hybrid);

#ifdef ERRORESTIMATION_DEBUG
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

#ifdef ERRORESTIMATION_DEBUG
        {
            std::map<int,int> matelem;
            int64_t nel = cmesh_H1Hybrid->NElements();
            for (int64_t el = 0; el<nel; el++) {
                TPZCompEl *cel = cmesh_H1Hybrid->Element(el);
                TPZGeoEl *gel = cel->Reference();
//                TPZManVector<REAL,3> center(3);
//                TPZGeoElSide gelside(gel);
//                gelside.CenterX(center);
//                std::cout << "Matid " << gel->MaterialId() << " center " << center << std::endl;
                matelem[gel->MaterialId()]++;
            }
            std::cout << __PRETTY_FUNCTION__ << " number of computational elements by material \n";
            for (auto it : matelem) {
                std::cout << "Material id " << it.first << " number of elements " << it.second << std::endl;
            }
        }
#endif


        createspace.InsertLagranceMaterialObjects(cmesh_H1Hybrid);
        createspace.AddInterfaceElements(cmesh_H1Hybrid);
#ifdef ERRORESTIMATION_DEBUG
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
#ifdef ERRORESTIMATION_DEBUG
        {
            std::ofstream out("mphysicsmeshBeforeCondense.txt");
            cmesh_H1Hybrid->Print(out);
        }
#endif
        createspace.GroupandCondenseElements(cmesh_H1Hybrid);

        cmesh_H1Hybrid->InitializeBlock();
        cmesh_H1Hybrid->ComputeNodElCon();
        
#ifdef ERRORESTIMATION_DEBUG
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
        
#ifdef ERRORESTIMATION_DEBUG
        {
            
            std::ofstream out2("H1HybridMesh.txt");
            cmesh_H1Hybrid->Print(out2);
            std::ofstream out3("gmeshHybridH1.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out3);
            
        }
#endif
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
    TPZMaterialT<STATE> *matref = 0;
    int nstate = -1;

    if(config.exact) {
        // Creates Poisson material
        TPZMatLaplacianHybrid *material = new TPZMatLaplacianHybrid(matID, dim);
        if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
            material->SetForcingFunction(config.exact->ForceFunc(), 5);
            material->SetExactSol(config.exact->ExactSolution(), 5);
        }
        nstate = 1;
        cmesh_H1Hybrid->InsertMaterialObject(material);
        matref = material;
    } else if(config.exactelast) {
        REAL elast = 1.;
        REAL poisson = 0.3;
        REAL fx(0.), fy(0.);
        auto *matelast = new TPZHybridElasticity2D(matID,elast,poisson,fx,fy);
        matelast->SetForcingFunction(config.exactelast->ForceFunc(), 5);
        matelast->SetExactSol(config.exactelast->ExactSolution(), 5);
        cmesh_H1Hybrid->InsertMaterialObject(matelast);
        matref = matelast;
        nstate = 2;


    }
    //    TPZMaterial * mat(material);
    //    cmesh->InsertMaterialObject(mat);

    // Inserts boundary conditions
    TPZFNMatrix<4,STATE> val1(nstate, nstate, 0.);
    TPZManVector<STATE, 2> val2(nstate, 1.);
    auto *BCond0 = matref->CreateBC(matref, -1, dirichlet, val1, val2);
    if (config.exact && config.exact.operator*().fExact != TLaplaceExample1::ENone) {
        BCond0->SetForcingFunctionBC(config.exact->ExactSolution(),5);
    }
    if (config.exactelast && config.exactelast.operator*().fProblemType != TElasticity2DAnalytic::ENone) {
        BCond0->SetForcingFunctionBC(config.exactelast->ExactSolution(),5);
    }
    auto *BCond1 = matref->CreateBC(matref, -2, neumann, val1, val2);
    if (config.exactelast && config.exactelast.operator*().fProblemType != TElasticity2DAnalytic::ENone) {
        BCond1->SetForcingFunctionBC(config.exactelast->ExactSolution(),5);
    }
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

    TPZLinearAnalysis an(cmeshH1);


#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> strmat(cmeshH1);
    strmat.SetNumThreads(0);
    //        strmat.SetDecomposeType(ELDLt);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmeshH1);
    strmat.SetNumThreads(0);
    strmat.SetDecomposeType(ELDLt);
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
    int nstate = cmeshH1->FindMaterial(1)->NStateVariables();
    if(nstate == 1) {
        scalnames.Push("Solution");
        vecnames.Push("Derivative");
        scalnames.Push("ExactSolution");
    } else {
        vecnames.Push("displacement");
        scalnames.Push("sig_x");
        scalnames.Push("sig_y");
        scalnames.Push("tau_xy");
    }




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

    if(nstate == 1) {
        an.SetExact(config.exact.operator*().ExactSolution());
    } else {
        an.SetExact(config.exactelast.operator*().ExactSolution());
    }

    TPZManVector<REAL> errorvec(10, 0.);
    int64_t nelem = cmeshH1->NElements();
    cmeshH1->LoadSolution(cmeshH1->Solution());
    cmeshH1->ExpandSolution();
    if(nstate == 1) {
        cmeshH1->ElementSolution().Redim(nelem, 10);
    } else {
        cmeshH1->ElementSolution().Redim(nelem, 6);
        errorvec.resize(6);
    }
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
    
    TPZLinearAnalysis an(cmesh_H1Hybrid);

#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> strmat(cmesh_H1Hybrid);
    strmat.SetNumThreads(0);
    //        strmat.SetDecomposeType(ELDLt);
#else
    //    TPZFrontStructMatrix<TPZFrontSym<STATE> > strmat(Hybridmesh);
    //    strmat.SetNumThreads(2);
    //    strmat.SetDecomposeType(ELDLt);
    TPZSkylineStructMatrix<STATE> strmat(cmesh_H1Hybrid);
    strmat.SetNumThreads(0);
#endif
        
    
    std::set<int> matIds;
    
    
    for (auto matid : config.materialids) {
        
        matIds.insert(matid);
    }
    TPZMaterial *matref = cmesh_H1Hybrid->FindMaterial(*matIds.begin());
    int nstate = matref->NStateVariables();
        
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
    int nerrors = 5;
    if(nstate == 1) {
        scalnames.Push("Pressure");
        scalnames.Push("PressureExact");
        nerrors = 5;
    } else if (nstate == 2) {
        vecnames.Push("displacement");
        scalnames.Push("SigmaX");
        scalnames.Push("SigmaY");
        scalnames.Push("TauXY");
        nerrors = 6;
    }
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
    if(config.exact) {
        an.SetExact(config.exact.operator*().ExactSolution());
    } else if(config.exactelast) {
        an.SetExact(config.exactelast.operator*().ExactSolution());
    }

        TPZManVector<REAL> errorvec(nerrors, 0.);
        int64_t nelem = cmesh_H1Hybrid->NElements();
        cmesh_H1Hybrid->LoadSolution(cmesh_H1Hybrid->Solution());
        cmesh_H1Hybrid->ExpandSolution();
        cmesh_H1Hybrid->ElementSolution().Redim(nelem, nerrors);

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
