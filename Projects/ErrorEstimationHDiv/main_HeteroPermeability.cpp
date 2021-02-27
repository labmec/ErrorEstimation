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
#include "TPZVecL2.h"
#include "pzbndcond.h"

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

#include "Tools.h"
#include "TPZHDivErrorEstimatorH1.h"
#include "TPZHybridHDivErrorEstimator.h"

#include "TPZBFileStream.h"
#include <tuple>
#include <memory>


TPZMultiphysicsCompMesh *CreateNeumannHDivMesh(const ProblemConfig &problem);
TPZCompMesh *CreateNeumannFluxHDivMesh(const ProblemConfig &problem);
TPZMultiphysicsCompMesh *CreateBoundaryLayerHDivMesh(const ProblemConfig &problem);
void Neumann1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff);
void Neumann2(const TPZVec<REAL> &pt, TPZVec<STATE> &ff);
void Neumann3(const TPZVec<REAL> &pt, TPZVec<STATE> &ff);
void Neumann4(const TPZVec<REAL> &pt, TPZVec<STATE> &ff);
void Dirichlet1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff);
void Dirichlet2(const TPZVec<REAL> &pt, TPZVec<STATE> &ff);
void Dirichlet3(const TPZVec<REAL> &pt, TPZVec<STATE> &ff);
void Dirichlet4(const TPZVec<REAL> &pt, TPZVec<STATE> &ff);


void ExataOmega1(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
void ExataOmega2(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
void ExataOmega3(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
void ExataOmega4(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &ff);


bool mixedsolution = true;


REAL  alpha=1.;//sqrt(0.1);

bool IsreadMesh = true;

int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
   for (int i = 0; i < 2; i++) {
    
    ProblemConfig config;
   
    
    config.porder = 1;
    config.hdivmais = 1;
    config.ndivisions = 1;
    config.alpha=alpha;
    
    config.prefine=false;
    config.makepressurecontinuous = true;
    
    config.steklovexample = true;
    
    config.GalvisExample = false;
    
       
        config.dir_name = "HPermeability";
        //config.dir_name= "ESinSin";
        std::string command = "mkdir " + config.dir_name;
        system(command.c_str());
    
    
    int dim=3;
    
    //malha geometrica
    TPZGeoMesh *gmesh = nullptr;
    
    if(IsreadMesh){
        
    
    if(config.steklovexample) {

            config.exact.operator*().fExact =
                TLaplaceExample1::ESteklovNonConst;
            config.problemname = "ESteklovNonConst k=1 n=3";

            TPZGmshReader gmsh;
            std::string meshfilename = "MalhaSteklov.msh";
#ifdef MACOSX
            meshfilename = "../" + meshfilename;
#endif
            gmsh.GetDimNamePhysical()[2]["Omega1"] = 1;
            gmsh.GetDimNamePhysical()[2]["Omega2"] = 2;
            gmsh.GetDimNamePhysical()[2]["Omega3"] = 3;
            gmsh.GetDimNamePhysical()[2]["Omega4"] = 4;

            gmsh.GetDimNamePhysical()[1]["boundary1"] =5;
        gmsh.GetDimNamePhysical()[1]["boundary2"] =6;
        gmsh.GetDimNamePhysical()[1]["boundary3"] =7;
        gmsh.GetDimNamePhysical()[1]["boundary4"] =8;
        
        config.materialids.insert(1);
        config.materialids.insert(2);
        config.materialids.insert(3);
        config.materialids.insert(4);
        
        config.bcmaterialids.insert(5);
        config.bcmaterialids.insert(6);
        config.bcmaterialids.insert(7);
        config.bcmaterialids.insert(8);
        
        
        gmsh.SetFormatVersion("4.1");
        gmesh = gmsh.GeometricGmshMesh(meshfilename);
        gmsh.PrintPartitionSummary(std::cout);
        gmesh->SetDimension(dim);
        config.gmesh = gmesh;
    
    }
    
    else if(config.GalvisExample) {

        config.exact.operator*().fExact = TLaplaceExample1::EGalvisNonConst;
        config.problemname = "EGalvisNonConst k=1 n=3 Up=1";
        std::string meshfilename = "Galvismesh.msh";
#ifdef MACOSX
        meshfilename = "../" + meshfilename;
#endif
        TPZGmshReader gmsh;

        gmsh.GetDimNamePhysical()[2]["Omega1"] = 1;
        gmsh.GetDimNamePhysical()[2]["Omega2"] = 2;
        gmsh.GetDimNamePhysical()[2]["Omega3"] = 3;
        gmsh.GetDimNamePhysical()[2]["Omega4"] = 4;
        
        gmsh.GetDimNamePhysical()[1]["boundary1"] =5;
        gmsh.GetDimNamePhysical()[1]["boundary2"] =6;
        gmsh.GetDimNamePhysical()[1]["boundary3"] =7;
        gmsh.GetDimNamePhysical()[1]["boundary4"] =8;
        gmsh.GetDimNamePhysical()[1]["boundary5"] =9;
        gmsh.GetDimNamePhysical()[1]["boundary6"] =10;
        gmsh.GetDimNamePhysical()[1]["boundary7"] =11;
        gmsh.GetDimNamePhysical()[1]["boundary8"] =12;
        
        config.materialids.insert(1);
        config.materialids.insert(2);
        config.materialids.insert(3);
        config.materialids.insert(4);
        
        config.bcmaterialids.insert(5);
        config.bcmaterialids.insert(6);
        config.bcmaterialids.insert(7);
        config.bcmaterialids.insert(8);
        config.bcmaterialids.insert(9);
        config.bcmaterialids.insert(10);
        config.bcmaterialids.insert(11);
        config.bcmaterialids.insert(12);
        
        
        gmsh.SetFormatVersion("4.1");
        gmesh = gmsh.GeometricGmshMesh(meshfilename);
        gmsh.PrintPartitionSummary(std::cout);
        gmesh->SetDimension(dim);
        config.gmesh = gmesh;
        
    }
    
    
    else{

            std::string meshfilename = "../MeshHetero.msh";
            TPZGmshReader gmsh;  
            
            gmsh.GetDimNamePhysical()[2]["Omega1"] = 1;
            gmsh.GetDimNamePhysical()[2]["Omega2"] = 2;
            gmsh.GetDimNamePhysical()[2]["Omega3"] = 3;
            gmsh.GetDimNamePhysical()[2]["Omega4"] = 4;
            
            gmsh.GetDimNamePhysical()[1]["neumann1"] =5;
            gmsh.GetDimNamePhysical()[1]["neumann2"] =6;
            gmsh.GetDimNamePhysical()[1]["neumann3"] =7;
            gmsh.GetDimNamePhysical()[1]["neumann4"] =8;
            
            config.materialids.insert(1);
            config.materialids.insert(2);
            config.materialids.insert(3);
            config.materialids.insert(4);
            
            config.bcmaterialids.insert(5);
            config.bcmaterialids.insert(6);
            config.bcmaterialids.insert(7);
            config.bcmaterialids.insert(8);
            
            
            gmsh.SetFormatVersion("4.1");
            gmesh = gmsh.GeometricGmshMesh(meshfilename);
            gmsh.PrintPartitionSummary(std::cout);
            gmesh->SetDimension(dim);
            config.gmesh = gmesh;
    
    }
    }
    else {

        config.dir_name = "BoundaryLayer";
        config.exact = new TLaplaceExample1;
        config.exact.operator*().fExact = TLaplaceExample1::EBoundaryLayer;

        TPZManVector<int, 4> bcids(4, -1);
        gmesh = Tools::CreateGeoMesh(1, bcids);
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        config.gmesh = gmesh;
    }

       Tools::UniformRefinement(config.ndivisions, gmesh);
    
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmeshInitial.txt");
        gmesh->Print(out2);

    }
        
        TPZGeoMesh *markEstimatorMesh = new TPZGeoMesh();
        *markEstimatorMesh = *gmesh;
    
    TPZManVector<TPZCompMesh*, 2> meshvec_HDiv(2, 0);
    
    
    TPZMultiphysicsCompMesh *cmesh_HDiv=CreateNeumannHDivMesh(config);//CreateBoundaryLayerHDivMesh(config);//Hdiv x L2
    
    cmesh_HDiv->InitializeBlock();
    
       
       if(mixedsolution) Tools::SolveMixedProblem(cmesh_HDiv,config);
    
//    if(mixedsolution)
//    {
//        
//        // Creating the directory
//        std::string command = "mkdir " + config.dir_name;
//        system(command.c_str());
//        
//        TPZAnalysis an(cmesh_HDiv);
//        
//        {
//        TPZStack<std::string> scalnames, vecnames;
//        scalnames.Push("ExactPressure");
//        vecnames.Push("ExactFlux");
//
//        std::stringstream sout;
//        sout << config.dir_name << "/" <<  "ExactSolution_"<<config.problemname<<"Nref_"<<config.ndivisions<<".vtk";
//
//        an.DefineGraphMesh(2, scalnames, vecnames, sout.str());
//        an.PostProcess(2,2);
//        }
//        
//        
//#ifdef USING_MKL
//        TPZSymetricSpStructMatrix strmat(cmesh_HDiv);
//        strmat.SetNumThreads(0);
//        //        strmat.SetDecomposeType(ELDLt);
//        an.SetStructuralMatrix(strmat);
//#else
//        TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmesh_HDiv);
//        strmat.SetNumThreads(0);
//        //        TPZSkylineStructMatrix strmat3(cmesh_HDiv);
//        //        strmat3.SetNumThreads(8);
//#endif
//        
//        TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
//        direct->SetDirect(ELDLt);
//        an.SetSolver(*direct);
//        delete direct;
//        direct = 0;
//        an.Assemble();
//        an.Solve();//resolve o problema misto ate aqui
//        TPZStack<std::string> scalnames, vecnames;
//        scalnames.Push("Pressure");
//        scalnames.Push("ExactPressure");
//        vecnames.Push("Flux");
//        vecnames.Push("ExactFlux");
//        
//        std::stringstream sout;
//        sout << config.dir_name << "/" <<  "OriginalMixed_Order_"<<config.porder<<"Nref_"<<config.ndivisions<<".vtk";
// 
//        //an.DefineGraphMesh(2, scalnames, vecnames, "Original_Misto.vtk");
//        an.DefineGraphMesh(2, scalnames, vecnames, sout.str());
//        an.PostProcess(2,2);
//
//        
//    }
     
        
    
    
    meshvec_HDiv = cmesh_HDiv->MeshVector();
    
    //cria malha hibrida
    
    TPZHybridizeHDiv hybrid;
    auto HybridMesh = hybrid.Hybridize(cmesh_HDiv);
    HybridMesh->CleanUpUnconnectedNodes();//enumerar adequadamente os connects
    HybridMesh->AdjustBoundaryElements();
    delete cmesh_HDiv;
    delete meshvec_HDiv[0];
    delete meshvec_HDiv[1];
    
    int n=hybrid.fLagrangeInterface;
    int n1=hybrid.fHDivWrapMatid;
    std::pair<int,int> n2=hybrid.fInterfaceMatid;
    
    std::cout<<"---Original PerifericalMaterialId --- "<<std::endl;
    std::cout <<" LagrangeInterface = "<<n<<std::endl;
    std::cout <<" HDivWrapMatid = "<<n1<<std::endl;
    std::cout <<" InterfaceMatid = "<<n2<<std::endl;
    
    
    cmesh_HDiv=(HybridMesh);//malha hribrida
    meshvec_HDiv[0] = (HybridMesh)->MeshVector()[0];//malha Hdiv
    meshvec_HDiv[1] = (HybridMesh)->MeshVector()[1];//malha L2
    
//    {
//
//        //                std::ofstream outgeo("HrybridGeometria.txt");
//        //                std::get<0>(HybridMesh)->Reference()->Print(outgeo);
//                        std::ofstream out("OriginalHybridMesh.txt");
//                        (HybridMesh)->Print(out);
//        //
//        std::ofstream out2("OriginalFluxMesh.txt");
//        meshvec_HDiv[0]->Print(out2);
//
//        std::ofstream out3("OriginalPotentialMesh.txt");
//        meshvec_HDiv[1]->Print(out3);
//
//
//    }


       Tools::SolveHybridProblem(cmesh_HDiv,n2,config,false);
    

    
//    {
//        std::ofstream out("OriginalHybridMesh.txt");
//        (HybridMesh)->Print(out);
//    }

  //  PlotLagrangreMultiplier(meshvec_HDiv[1],config);
    
    
    //reconstroi potencial e calcula o erro
    {
        
       TPZHDivErrorEstimatorH1 HDivEstimate(*cmesh_HDiv);
        HDivEstimate.SetProblemConfig(config);
        HDivEstimate.SetPostProcUpliftOrder(config.hdivmais);
        HDivEstimate.SetAnalyticSolution(config.exact);
        
        HDivEstimate.fperformUplift = true;
        HDivEstimate.fUpliftOrder = 1;
        
 
        /*
      
        TPZHybridHDivErrorEstimator HDivEstimate(*cmesh_HDiv);
        
        HDivEstimate.fProblemConfig = config;
        HDivEstimate.fUpliftPostProcessOrder = config.hdivmais;
        
        if(config.GalvisExample || config.steklovexample){
            
            HDivEstimate.SetAnalyticSolution(config.exact);
        }
        */
        
        HDivEstimate.PotentialReconstruction();
        
        TPZManVector<REAL> elementerrors,errvec;
        std::string vtkPath = "hetero_permeability_error_results.vtk";
        HDivEstimate.ComputeErrors(errvec,elementerrors, vtkPath);

        Tools::hAdaptivity(HDivEstimate.PostProcMesh(), markEstimatorMesh, config);
    }
    
//    delete cmesh_HDiv;
//    delete meshvec_HDiv[0];
//    delete meshvec_HDiv[1];
   // return 0;
       
   }
    
    
    
}


void Exata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    solp[0] = pt[0]*pt[1];
    flux(0,0) = pt[1];
    flux(1,0) = pt[0];
    flux(2,0) = 0.;
}

void ExataOmega1(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
    flux.Resize(3, 1);
    
    REAL x,y,z;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    solp[0] = x*x-y*y;
    flux(0,0) = -2.*x;
    flux(1,0) = 2*y;
    flux(2,0) = 0.;
    //Exata(pt,solp,flux);
    //   std::cout<<"funcao em omega1 pt = " << pt << " sol " << solp[0] << std::endl;
}

void ExataOmega2(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
    flux.Resize(3, 1);
    
    REAL x,y,z;
    
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    REAL alpha2 = pow(alpha, 2);
    
    solp[0] = (x*x-y*y)/alpha2 ;
    
    flux(0,0) = -(2.)*x;
    flux(1,0) = (2.)*y;
    flux(2,0) = 0.;
    
//    flux(0,0) = -(2.)*x/alpha2;
//    flux(1,0) = (2.)*y/alpha2;
//    flux(2,0) = 0.;
    //    Exata(pt,solp,flux);
    
    //  std::cout<<"funcao em omega2 pt " << pt << " sol " << solp[0] << std::endl;
    
}

void ExataOmega3(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    solp.resize(1);
    flux.Resize(3, 1);
    
    ExataOmega1(pt, solp, flux);
    //    Exata(pt,solp,flux);
    
    //  std::cout<<"funcao em omega3 pt = " << pt << " sol " << solp[0] << std::endl;
    
    
}

void ExataOmega4(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
    flux.Resize(3, 1);
    
    REAL x,y,z;
    
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    REAL alpha4 = pow(alpha, 4);
    
    solp[0] = (x*x-y*y)/(alpha4);
    
    flux(0,0) = -2.*x;
    flux(1,0) = (2.)*y;
    flux(2,0) = 0.;
    
    
//    flux(0,0) = -2.*x/(alpha4);
//    flux(1,0) = (2.)*y/(alpha4);
//    flux(2,0) = 0.;
    //    Exata(pt,solp,flux);
    
    // std::cout<<"funcao em omega4 pt " << pt << " sol " << solp[0] << " flux " << flux << std::endl;
    
}

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &ff)
{
    ff[0]=0.;
}

void Neumann1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff)
{
    
//    TPZFMatrix<STATE> flux;
//    ExataOmega1(pt, ff, flux);
//    ff[0]=0.;
    
   
    
    TPZFMatrix<STATE> normal(3,1);
    normal.Zero();
    
    normal(0,0)=1.;
    REAL x,y,z;
    x = pt[0];
    y = pt[1];
    z = pt[2];
    //-K grad u. eta
    ff[0] = -2*x*normal(0,0) + 2*y*normal(1,0);//flux(0,0)*normal(0,0)+flux(1,0)*normal(1,0)+flux(2,0)*normal(2,0);
   //   std::cout<<"derivada normal em omega1 pt " << pt << " normal " << ff[0] << std::endl;
    
}

void Neumann2(const TPZVec<REAL> &pt, TPZVec<STATE> &ff)
{
    
//    TPZFMatrix<STATE> flux;
//    ExataOmega2(pt, ff, flux);
//    ff[0]=0.;
 //    REAL alpha2 = pow(alpha, 2);
    
    TPZFMatrix<STATE> normal(3,1);
    normal.Zero();
    
    normal(1,0) = 1.;
   
    
    REAL x,y,z;
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    ff[0] = -2*x*normal(0,0) + 2*y*normal(1,0);
    
    //ff[0] = (flux(0,0)*normal(0,0)+flux(1,0)*normal(1,0)+flux(2,0)*normal(2,0))* alpha2;
   //    std::cout<<"derivada normal em omega2 pt " << pt << " normal " << ff[0] << std::endl;//
    
}

void Neumann3(const TPZVec<REAL> &pt, TPZVec<STATE> &ff)
{
    
//    TPZFMatrix<STATE> flux;
//
//    ExataOmega3(pt, ff, flux);
//    ff[0]=0.;
    
    TPZFMatrix<STATE> normal(3,1);
    normal.Zero();
    
    normal(0,0) = (-1.);
    
    REAL x,y,z;
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    ff[0] = -2*x*normal(0,0) + 2*y*normal(1,0);
    
    //ff[0] = flux(0,0)*normal(0,0)+flux(1,0)*normal(1,0)+flux(2,0)*normal(2,0);
  //    std::cout<<"derivada normal em omega3 pt " << pt << " normal " << ff[0] << std::endl;
    
}

void Neumann4(const TPZVec<REAL> &pt, TPZVec<STATE> &ff)
{
    
//    TPZFMatrix<STATE> flux;
//
//    ExataOmega4(pt, ff, flux);
//    ff[0]=0.;
//    REAL alpha4 = pow(alpha, 4);
    
    TPZFMatrix<STATE> normal(3,1);
    normal.Zero();
    
    normal(1,0) = (-1.);
    
    REAL x,y,z;
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    ff[0] = -2*x*normal(0,0) + 2*y*normal(1,0);
    
    
  //  ff[0] = ((flux(0,0)*normal(0,0)+flux(1,0)*normal(1,0)+flux(2,0)*normal(2,0)))*alpha4;
  //    std::cout<<"derivada normal em omega4 pt " << pt << " normal " << ff[0] << std::endl;
    
}

void Dirichlet1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff)
{
    
    TPZFMatrix<STATE> flux;
    TPZVec<STATE> sol;
    
    ExataOmega1(pt, sol, flux);
    ff[0]=sol[0];
    
   //  std::cout<<"Dirichlet em omega1 pt " << pt << " normal " << ff[0] << std::endl;//
    
}

void Dirichlet2(const TPZVec<REAL> &pt, TPZVec<STATE> &ff)
{
    
    TPZFMatrix<STATE> flux;
    TPZVec<STATE> sol;
    
    ExataOmega2(pt, sol, flux);
    ff[0]=sol[0];
    
    
}
void Dirichlet3(const TPZVec<REAL> &pt, TPZVec<STATE> &ff)
{
    
    TPZFMatrix<STATE> flux;
    TPZVec<STATE> sol;
    
    ExataOmega3(pt, sol, flux);
    ff[0]=sol[0];
    
    
}

void Dirichlet4(const TPZVec<REAL> &pt, TPZVec<STATE> &ff)
{
    
    TPZFMatrix<STATE> flux;
    TPZVec<STATE> sol;
    
    ExataOmega4(pt, sol, flux);
    ff[0]=sol[0];
    
    
}




TPZMultiphysicsCompMesh *CreateNeumannHDivMesh(const ProblemConfig &problem) {
    
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(problem.gmesh);
    TPZMaterial *mat = NULL;
    for (auto matid : problem.materialids) {
        TPZMixedPoisson *mix =NULL;

        TPZAutoPointer<TPZFunction<STATE> > solexata;
        TPZFMatrix<REAL> K(3,3,0),invK(3,3,0);
        K.Identity();
        invK.Identity();
        REAL k1=2;
        REAL k2=5;
        
     
        switch (matid) {
            case 1:
            {
                mix = new TPZMixedPoisson(matid, cmesh->Dimension());
                mix->SetInternalFlux(0);
                
                if(problem.steklovexample){

                    mix->SetForcingFunctionExact(
                        problem.exact.operator*().Exact());
                    K(0, 0) = 5.;
                    K(1,1)=5.;
                    K(2,2)=5.;
                    invK(0,0)=1./5.;
                    invK(1,1)=1./5.;
                    invK(2,2)=1./5.;
                    mix->SetPermeabilityTensor(K, invK);
                    if (!mat) mat = mix;
                    cmesh->InsertMaterialObject(mix);
                    
                }
                
                else if(problem.GalvisExample) {
                    mix->SetForcingFunctionExact(
                        problem.exact.operator*().Exact());
                    mix->SetForcingFunction(
                        problem.exact.operator*().ForcingFunction());
                    K(0, 0) = k1 * k2;
                    K(1, 1) = k1 * k2;
                    K(2, 2) = k1 * k2;
                    invK(0, 0) = 1. / k1 * k2;
                    invK(1, 1) = 1. / k1 * k2;
                    invK(2, 2) = 1. / k1 * k2;
                    mix->SetPermeabilityTensor(K, invK);
                    if (!mat) mat = mix;
                    cmesh->InsertMaterialObject(mix);

                }
                else{
                    
                    solexata = new TPZDummyFunction<STATE>(ExataOmega1,10);
                    TPZAutoPointer<TPZFunction<STATE>> sol(solexata);
                    mix->SetForcingFunctionExact(sol);
                    mix->SetPermeabilityTensor(K, invK);
                    if (!mat) mat = mix;
                    cmesh->InsertMaterialObject(mix);
                    
                }
                
                
//                mix->SetPermeabilityTensor(K, invK);
//                if (!mat) mat = mix;
//                cmesh->InsertMaterialObject(mix);
            }
                break;
                case 2:
            {
                mix = new TPZMixedPoisson(matid, cmesh->Dimension());
                mix->SetInternalFlux(0);
                
                if(problem.steklovexample){

                    mix->SetForcingFunctionExact(
                        problem.exact.operator*().Exact());
                    mix->SetPermeabilityTensor(K, invK);
                    if (!mat) mat = mix;
                    cmesh->InsertMaterialObject(mix);
                    
                }
                
                else if(problem.GalvisExample) {
                    mix->SetForcingFunctionExact(
                        problem.exact.operator*().Exact());
                    mix->SetForcingFunction(
                        problem.exact.operator*().ForcingFunction());
                    K(0, 0) = k2;
                    K(1, 1) = k2;
                    K(2, 2) = k2;
                    invK(0, 0) = 1. / k2;
                    invK(1, 1) = 1. / k2;
                    invK(2, 2) = 1. / k2;
                    mix->SetPermeabilityTensor(K, invK);
                    if (!mat) mat = mix;
                    cmesh->InsertMaterialObject(mix);

                }
                
                
                else{
                    
                    solexata = new TPZDummyFunction<STATE>(ExataOmega2,10);
                    TPZAutoPointer<TPZFunction<STATE>> sol(solexata);
                    mix->SetForcingFunctionExact(sol);
                    
                    STATE alpha2=(problem.alpha)*(problem.alpha);
                    
                    K(0,0) = alpha2;
                    K(1,1) = alpha2;
                    K(2,2) = alpha2;
                    
                    invK(0,0) = 1./alpha2;
                    invK(1,1) = 1./alpha2;
                    invK(2,2) = 1./alpha2;
                    mix->SetPermeabilityTensor(K, invK);
                    if (!mat) mat = mix;
                    cmesh->InsertMaterialObject(mix);
                }
                
//                mix->SetPermeabilityTensor(K, invK);
//                if (!mat) mat = mix;
//                cmesh->InsertMaterialObject(mix);
                
                
                
            }
                 break;
            case 3:{
                
                mix = new TPZMixedPoisson(matid, cmesh->Dimension());
                
                if(problem.steklovexample){

                    mix->SetForcingFunctionExact(
                        problem.exact.operator*().Exact());
                    K(0, 0) = 5.;
                    K(1,1)=5.;
                    K(2,2)=5.;
                    invK(0,0)=1./5.;
                    invK(1,1)=1./5.;
                    invK(2,2)=1./5.;
                    mix->SetPermeabilityTensor(K, invK);
                    if (!mat) mat = mix;
                    cmesh->InsertMaterialObject(mix);
                    
                }
                
               else if(problem.GalvisExample) {
                    mix->SetForcingFunctionExact(
                        problem.exact.operator*().Exact());
                    mix->SetForcingFunction(
                        problem.exact.operator*().ForcingFunction());
                    K(0, 0) = k1;
                    K(1, 1) = k1;
                    K(2, 2) = k1;
                    invK(0, 0) = 1. / k1;
                    invK(1, 1) = 1. / k1;
                    invK(2, 2) = 1. / k1;
                    mix->SetPermeabilityTensor(K, invK);
                    if (!mat) mat = mix;
                    cmesh->InsertMaterialObject(mix);

                }
                
                
                else{
                    
                    solexata = new TPZDummyFunction<STATE>(ExataOmega3,10);
                    TPZAutoPointer<TPZFunction<STATE>> sol(solexata);
                    mix->SetForcingFunctionExact(sol);
                    mix->SetPermeabilityTensor(K, invK);
                    if (!mat) mat = mix;
                    cmesh->InsertMaterialObject(mix);
                }
                
                
//                mix->SetPermeabilityTensor(K, invK);
//                if (!mat) mat = mix;
//                cmesh->InsertMaterialObject(mix);
                
            }
                 break;
            case 4:{
                
                mix = new TPZMixedPoisson(matid, cmesh->Dimension());
                mix->SetInternalFlux(0);
                
                if(problem.steklovexample){

                    mix->SetForcingFunctionExact(
                        problem.exact.operator*().Exact());
                    K(0, 0) = 1.;
                    K(1,1)=1.;
                    K(2,2)=1.;
                    invK(0,0)=1./1.;
                    invK(1,1)=1./1.;
                    invK(2,2)=1./1.;
                    mix->SetPermeabilityTensor(K, invK);
                    if (!mat) mat = mix;
                    cmesh->InsertMaterialObject(mix);
                    
                }
               else if(problem.GalvisExample) {
                    mix->SetForcingFunctionExact(
                        problem.exact.operator*().Exact());
                    mix->SetForcingFunction(
                        problem.exact.operator*().ForcingFunction());

                    K.Identity();
                    invK.Identity();

                    mix->SetPermeabilityTensor(K, invK);
                    if (!mat) mat = mix;
                    cmesh->InsertMaterialObject(mix);

                }

                else{
                    
                    
                    solexata = new TPZDummyFunction<STATE>(ExataOmega4,10);
                    TPZAutoPointer<TPZFunction<STATE>> sol(solexata);
                    mix->SetForcingFunctionExact(sol);
                    
                    STATE alpha4=pow(problem.alpha, 4);
                
                    
                    K(0,0) = alpha4;
                    K(1,1) = alpha4;
                    K(2,2) = alpha4;
                    
                    invK(0,0) = 1./alpha4;
                    invK(1,1) = 1./alpha4;
                    invK(2,2) = 1./alpha4;
                    mix->SetPermeabilityTensor(K, invK);
                    if (!mat) mat = mix;
                    cmesh->InsertMaterialObject(mix);
                    
                }
//                mix->SetPermeabilityTensor(K, invK);
//                if (!mat) mat = mix;
//                cmesh->InsertMaterialObject(mix);
                
                
                
            }
                
            default:
                break;
        }
   
    }
    
    
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 1;
        TPZAutoPointer<TPZFunction<STATE> > bcfunction;
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        
        if(problem.steklovexample)
        {

                bctype=0;
                val2.Zero();
                TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
                bc->TPZMaterial::SetForcingFunction(
                    problem.exact.operator*().Exact());
                cmesh->InsertMaterialObject(bc);

        }
        
        else if(problem.GalvisExample){
            bctype=0;
            val2.Zero();
            TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
            cmesh-> InsertMaterialObject(bc);
            
            
        }
        
        else{
            
        switch(matid){

        case 5:
        {
            bcfunction=new TPZDummyFunction<STATE>(Dirichlet1,5);
         //   bcfunction=new TPZDummyFunction<STATE>(Neumann1,5);
            int bctype=0;
            bc->SetType(bctype);
            TPZAutoPointer<TPZFunction<STATE>> func(bcfunction);
            bc->TPZMaterial::SetForcingFunction(func);
            cmesh->InsertMaterialObject(bc);
        }
            
            break;

        case 6:
        {
            bcfunction=new TPZDummyFunction<STATE>(Neumann2,5);
            int bctype = 1;
            bc->SetType(bctype);

            TPZAutoPointer<TPZFunction<STATE>> func(bcfunction);
            bc->TPZMaterial::SetForcingFunction(func);
            cmesh->InsertMaterialObject(bc);
        }
            break;
        case 7:
        {
            //bcfunction=new TPZDummyFunction<STATE>(Neumann3,5);
            bcfunction=new TPZDummyFunction<STATE>(Dirichlet3,5);
            int  bctype=0;
            bc->SetType(bctype);
            TPZAutoPointer<TPZFunction<STATE>> func(bcfunction);
            bc->TPZMaterial::SetForcingFunction(func);
            cmesh->InsertMaterialObject(bc);
        }
            break;

        case 8:
        {
            bcfunction=new TPZDummyFunction<STATE>(Neumann4,5);
            int bctype = 1;
            bc->SetType(bctype);
            TPZAutoPointer<TPZFunction<STATE>> func(bcfunction);
            bc->TPZMaterial::SetForcingFunction(func);
            cmesh->InsertMaterialObject(bc);
        }
    }
}
        
}
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    
    TPZManVector<int> active(2,1);
    TPZManVector<TPZCompMesh *> meshvector(2,0);
    
    meshvector[0] = CreateNeumannFluxHDivMesh(problem);
    meshvector[1] = Tools::CreatePressureMesh(problem);
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    bool keepmatrix = false;
    TPZCompMeshTools::CreatedCondensedElements(cmesh, true, keepmatrix);
    
    {
        std::ofstream out("MixedMesh.txt");
        cmesh->Print(out);
    }
    
    return cmesh;
}

TPZCompMesh *CreateNeumannFluxHDivMesh(const ProblemConfig &problem) {
    int dim = problem.gmesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial *mat = NULL;
    problem.gmesh->ResetReference();
    for (auto matid : problem.materialids) {
        TPZVecL2 *mix = new TPZVecL2(matid);
        mix->SetDimension(dim);
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 1;
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        cmesh->InsertMaterialObject(bc);
        
    }
    
    cmesh->SetDefaultOrder(problem.porder);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmesh->AutoBuild();
    if (problem.hdivmais) {
        int64_t nel = cmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if (gel->Dimension() == dim) {
                int side = gel->NSides() - 1;
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
                intel->SetSideOrder(side, problem.porder + problem.hdivmais);//seta ordem +hdivmais
                intel->SetPreferredOrder(problem.porder+problem.hdivmais);
            }
        }
        
        if(problem.prefine){
            Tools::Prefinamento(cmesh, problem.ndivisions, problem.porder);
        }
        
    }
    cmesh->InitializeBlock();
    return cmesh;
    
}

TPZMultiphysicsCompMesh *CreateBoundaryLayerHDivMesh(const ProblemConfig &problem) {
    
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(problem.gmesh);
    TPZMaterial *mat = NULL;
    for (auto matid : problem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        
        TPZFMatrix<REAL> K(3,3,0),invK(3,3,0);
        K.Identity();
        invK.Identity();

        mix->SetForcingFunctionExact(problem.exact.operator*().Exact());
        mix->SetPermeabilityTensor(K, invK);
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);


    }
    
    
    
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        TPZAutoPointer<TPZFunction<STATE> > bcfunction;
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);

        bc->TPZMaterial::SetForcingFunction(problem.exact.operator*().Exact());
        cmesh->InsertMaterialObject(bc);
    }
    
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    
    TPZManVector<int> active(2,1);
    TPZManVector<TPZCompMesh *> meshvector(2,0);
    
    meshvector[0] = CreateNeumannFluxHDivMesh(problem);
    meshvector[1] = Tools::CreatePressureMesh(problem);
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    bool keepmatrix = false;
    TPZCompMeshTools::CreatedCondensedElements(cmesh, true, keepmatrix);
    
    {
        std::ofstream out("MixedMesh.txt");
        cmesh->Print(out);
    }
    
    return cmesh;
}
