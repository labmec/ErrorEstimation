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

#include "TPZHybridHDivErrorEstimator.h"
#include "Tools.h"

#include "TPZBFileStream.h"
#include <tuple>
#include <memory>


bool IsgmeshReader=true;
bool neumann=true;

bool mixedsolution=true;


int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    

    ProblemConfig config;
        
    config.ndivisions= 0;
        
    config.porder = 1;
    config.hdivmais = 1;
    config.ndivisions=0;
    config.prefine=false;
    config.makepressurecontinuous = true;
    
    config.exact.fExact = TLaplaceExample1::EBubble;//EConst;//ESinSin;//ESinMark;//EArcTanSingular;//EArcTan;//
    config.problemname ="Bubble";//"EConst";//"ESinSin";//" ESinMark";////"EArcTanSingular_PRef";//""ArcTang";//
    
    config.dir_name= "Bubble";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    
 //   FunctionTest();
    
    int dim=3;
    
    //malha geometrica
    TPZGeoMesh *gmesh = nullptr;
    
    if(IsgmeshReader){
        
        
        std::string meshfilename = "../LCircle.msh";
    
        if(dim==3)
        {
            meshfilename = "../Cube.msh";
        }
        TPZGmshReader gmsh;
      //  gmsh.GetDimNamePhysical().resize(4);
      //  gmsh.GetDimPhysicalTagName().resize(4);
        if(dim==2)
        {
            gmsh.GetDimNamePhysical()[1]["dirichlet"] =2;
            gmsh.GetDimNamePhysical()[2]["domain"] = 1;
        }
        else
        {
            gmsh.GetDimNamePhysical()[2]["dirichlet"] =2;
            gmsh.GetDimNamePhysical()[3]["domain"] = 1;
        }
        config.materialids.insert(1);
        config.bcmaterialids.insert(2);
    
        
        gmsh.SetFormatVersion("4.1");
        gmesh = gmsh.GeometricGmshMesh(meshfilename);
        gmsh.PrintPartitionSummary(std::cout);
        gmesh->SetDimension(dim);
        config.gmesh = gmesh;

    }
    
    else{
    
        gmesh = CreateGeoMesh(2);
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        config.bcmaterialids.insert(-2);
        config.gmesh = gmesh;
        gmesh->SetDimension(dim);
        
        
        
    }
    
    UniformRefinement(config.ndivisions, gmesh);

   #ifdef PZDEBUG
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmeshInitial.txt");
        gmesh->Print(out2);
        
    }
#endif
        
    TPZManVector<TPZCompMesh*, 2> meshvec_HDiv(2, 0);
    
    TPZMultiphysicsCompMesh *cmesh_HDiv=nullptr;
    

    cmesh_HDiv = CreateHDivMesh(config);//Hdiv x L2
    cmesh_HDiv->InitializeBlock();
    
#ifdef PZDEBUG
    {
        
        std::ofstream out2("MalhaMista.txt");
        cmesh_HDiv->Print(out2);
   
    }
#endif
    
    
    
    if(mixedsolution)
    {

        
        TPZAnalysis an(cmesh_HDiv);
        
        
#ifdef USING_MKL
        TPZSymetricSpStructMatrix strmat(cmesh_HDiv);
        strmat.SetNumThreads(0);
        //        strmat.SetDecomposeType(ELDLt);
        an.SetStructuralMatrix(strmat);
#else
        TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmesh_HDiv);
        strmat.SetNumThreads(0);
        //        TPZSkylineStructMatrix strmat3(cmesh_HDiv);
        //        strmat3.SetNumThreads(8);
#endif
        
        TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
        direct->SetDirect(ELDLt);
        an.SetSolver(*direct);
        delete direct;
        direct = 0;
        an.Assemble();
        an.Solve();//resolve o problema misto ate aqui
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Pressure");
        scalnames.Push("ExactPressure");
        vecnames.Push("Flux");
        vecnames.Push("ExactFlux");
        
        std::stringstream sout;
       
        sout << config.dir_name << "/"  "OriginalMixed_Order_"<<config.problemname<<"Order"<< config.porder<<"Nref_"<<config.ndivisions<<".vtk";
        
        //an.DefineGraphMesh(2, scalnames, vecnames, "Original_Misto.vtk");
        an.DefineGraphMesh(dim, scalnames, vecnames, sout.str());
        int resolution=2;
        an.PostProcess(resolution,dim);
        
        
    }
    
    
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
    int n2=hybrid.fInterfaceMatid;
    
    std::cout<<"---Original PerifericalMaterialId --- "<<std::endl;
    std::cout <<" LagrangeInterface = "<<n<<std::endl;
    std::cout <<" HDivWrapMatid = "<<n1<<std::endl;
    std::cout <<" InterfaceMatid = "<<n2<<std::endl;
    
    
    cmesh_HDiv=(HybridMesh);//malha hribrida
    meshvec_HDiv[0] = (HybridMesh)->MeshVector()[0];//malha Hdiv
    meshvec_HDiv[1] = (HybridMesh)->MeshVector()[1];//malha L2

   #ifdef PZDEBUG
    {
        
        std::ofstream out2("OriginalFluxMesh.txt");
        meshvec_HDiv[0]->Print(out2);

        std::ofstream out3("OriginalPotentialMesh.txt");
        meshvec_HDiv[1]->Print(out3);
       
    }
#endif

    
    SolveHybridProblem(cmesh_HDiv,n2,config);

#ifdef PZDEBUG
    {
        std::ofstream out("OriginalHybridMesh.txt");
        (HybridMesh)->Print(out);
    }
#endif
 //   PlotLagrangreMultiplier(meshvec_HDiv[1],config);

    
    
    //reconstroi potencial e calcula o erro
    {
        TPZHybridHDivErrorEstimator HDivEstimate(*cmesh_HDiv);
        
        HDivEstimate.fProblemConfig = config;
        HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
        HDivEstimate.SetAnalyticSolution(config.exact);
        
        HDivEstimate.PotentialReconstruction();
        
        TPZManVector<REAL> elementerrors;
        HDivEstimate.ComputeErrors(elementerrors);
        
    }
    
    delete cmesh_HDiv;
    delete meshvec_HDiv[0];
    delete meshvec_HDiv[1];
    return 0;
}

