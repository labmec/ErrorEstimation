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
#include "TPZHDivErrorEstimatorH1.h"

#include "Tools.h"

#include "TPZBFileStream.h"
#include <tuple>
#include <memory>




bool IsgmeshReader = true;//para ler a malha
bool neumann = true; //para o problema local de neumann da forlmulacao Mark

bool mixedsolution = true;//se quiser rodar o prolbema misto

bool PostProcessingFEM = true;//para graficos da solucao FEM

int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(ECube);
    
    
    ///
    
    
    ///
    
    
    for(int ndiv=0; ndiv<5; ndiv++){
        ProblemConfig config;
        
        config.porder = 1;
        config.hdivmais = 1;
        config.ndivisions = ndiv;
        config.dimension = 3;
        config.prefine=false;
        config.makepressurecontinuous = true;
        config.adaptivityStep = ndiv;
        config.TensorNonConst = true;//para problem 3d com tensor nao constante
        
        config.exact.fExact = TLaplaceExample1::EBubble;//ESinSinDirNonHom;//EArcTan;//ESinSin;
        config.problemname = "EBubblePermNonConst";
        
        bool RunMark = false;
        
        config.dir_name= "TestePaper";//"ReconstructionH1";//
        std::string command = "mkdir " + config.dir_name;
        system(command.c_str());
        
        int dim = config.dimension;
        
        //malha geometrica
        TPZGeoMesh *gmesh = nullptr;
        
        if(IsgmeshReader){
            
            
            std::string meshfilename = "../LCircle.msh";
            //std::string meshfilename = "../esfera2.msh";
            
            if(dim==3)
            {
                meshfilename = "Cube.msh";
                
            }
            TPZGmshReader gmsh;
            //  gmsh.GetDimNamePhysical().resize(4);
            //  gmsh.GetDimPhysicalTagName().resize(4);
            if(dim==2)
            {
                gmsh.GetDimNamePhysical()[1]["dirichlet"] = 2;
                // gmsh.GetDimNamePhysical()[1]["boundary"] =2;
                gmsh.GetDimNamePhysical()[2]["domain"] = 1;
                //  gmsh.GetDimNamePhysical()[1]["boundary2"] =3;
                
            }
            else
            {
                gmsh.GetDimNamePhysical()[2]["dirichlet"] =2;
                gmsh.GetDimNamePhysical()[3]["domain"] = 1;
            }
            config.materialids.insert(1);
            config.bcmaterialids.insert(2);
            //config.bcmaterialids.insert(3);
            
            
            gmsh.SetFormatVersion("4.1");
            gmesh = gmsh.GeometricGmshMesh(meshfilename);
            gmsh.PrintPartitionSummary(std::cout);
            gmesh->SetDimension(dim);
            config.gmesh = gmesh;
            
        }
        
        else{
            TPZManVector<int,4> bcids(4,-1);
            gmesh = CreateGeoMesh(2, bcids);
            config.materialids.insert(1);
            config.bcmaterialids.insert(-1);
            config.gmesh = gmesh;
            gmesh->SetDimension(dim);
            
            
            
        }
        
#ifdef PZDEBUG
        {
            std::ofstream out("gmesh.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
            std::ofstream out2("gmeshInitial.txt");
            gmesh->Print(out2);
            
        }
#endif
        
        TPZGeoMesh *hybridEstimatorMesh = new TPZGeoMesh();
        *hybridEstimatorMesh = *gmesh;
        
        
        
        UniformRefinement(config.ndivisions, gmesh);
        RandomRefine(config, config.ndivisions);
        
#ifdef PZDEBUG2
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
        
#ifdef PZDEBUG2
        {
            
            std::ofstream out2("MalhaMista.txt");
            cmesh_HDiv->Print(out2);
            
        }
#endif
        
        
        
        if(mixedsolution) SolveMixedProblem(cmesh_HDiv,config);
        
        
        meshvec_HDiv = cmesh_HDiv->MeshVector();
        
        //cria malha hibrida
        
        TPZHybridizeHDiv hybrid;
        auto HybridMesh = hybrid.Hybridize(cmesh_HDiv);
        HybridMesh->CleanUpUnconnectedNodes();//enumerar adequadamente os connects
        HybridMesh->AdjustBoundaryElements();
        delete cmesh_HDiv;
        delete meshvec_HDiv[0];
        delete meshvec_HDiv[1];
        
        
        std::cout<<"---Original PerifericalMaterialId --- "<<std::endl;
        std::cout <<" LagrangeInterface = "<<hybrid.fLagrangeInterface<<std::endl;
        std::cout <<" HDivWrapMatid = "<<hybrid.fHDivWrapMatid<<std::endl;
        std::cout <<" InterfaceMatid = "<<hybrid.fInterfaceMatid<<std::endl;
        
        
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
        
        
        SolveHybridProblem(cmesh_HDiv,hybrid.fInterfaceMatid,config);
        
#ifdef PZDEBUG2
        {
            std::ofstream out("OriginalHybridMesh.txt");
            (HybridMesh)->Print(out);
        }
#endif
        //   PlotLagrangreMultiplier(meshvec_HDiv[1],config);
        
        
        
        //reconstroi potencial e calcula o erro
        {
            
            if(RunMark){
                TPZHDivErrorEstimatorH1 HDivEstimate(*cmesh_HDiv);
                HDivEstimate.fProblemConfig = config;
                HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
                HDivEstimate.SetAnalyticSolution(config.exact);
                HDivEstimate.fperformUplift = true;
                HDivEstimate.fUpliftOrder = 1;
                
                HDivEstimate.PotentialReconstruction();
                
                TPZManVector<REAL> elementerrors;
                HDivEstimate.ComputeErrors(elementerrors);
                
            }
            
            else{
                
                TPZHybridHDivErrorEstimator HDivEstimate(*cmesh_HDiv);
                HDivEstimate.fProblemConfig = config;
                HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
                HDivEstimate.SetAnalyticSolution(config.exact);
                HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
                
                HDivEstimate.fPostProcesswithHDiv = true;
                
                HDivEstimate.PotentialReconstruction();
                
                TPZManVector<REAL> elementerrors;
                HDivEstimate.ComputeErrors(elementerrors);
                // hAdaptivity(&HDivEstimate.fPostProcMesh, hybridEstimatorMesh);
                
            }
            
            
        }
        
        //    delete cmesh_HDiv;
        //    delete meshvec_HDiv[0];
        //    delete meshvec_HDiv[1];
        //return 0;
    }
}

