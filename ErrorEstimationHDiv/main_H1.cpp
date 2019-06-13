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




bool IsgmeshReader = false;
bool neumann = true;

bool mixedsolution = false;


int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    

    ProblemConfig config;
    
    config.porder = 1;
    config.hdivmais = 1;
    config.ndivisions = 0;
    config.dimension = 2;
    config.prefine=false;
    config.makepressurecontinuous = true;
    
    config.exact.fExact = TLaplaceExample1::ESinSinDirNonHom;//;//EConst;//ESinSin;//ESinMark;//EArcTanSingular;//EArcTan;//
    config.problemname ="ESinSinDirNonHom";//"EConst";//"ESinSin";//" ESinMark";////"EArcTanSingular_PRef";//""ArcTang";//
    
    config.dir_name= "ESinSinDirNonHom";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    
    //geometric mesh

     TPZGeoMesh *gmesh = ReadGeometricMesh(config, IsgmeshReader);
    
#ifdef PZDEBUG
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmeshInitial.txt");
        gmesh->Print(out2);
        
    }
#endif
    
 
    UniformRefinement(config.ndivisions, gmesh);
   // RandomRefine(config, config.ndivisions);


        
    TPZManVector<TPZCompMesh*, 2> meshvec_HDiv(2, 0);
    
    TPZMultiphysicsCompMesh *cmesh_HDiv = CreateHDivMesh(config);//Hdiv x L2
    cmesh_HDiv->InitializeBlock();
    
#ifdef PZDEBUG
    {
        
        std::ofstream out2("MalhaMista.txt");
        cmesh_HDiv->Print(out2);
        
    }
#endif
    
    TPZMultiphysicsCompMesh *hybridmesh= HybridSolveProblem(cmesh_HDiv,meshvec_HDiv ,config);
  

//
//    
//    
//    if(mixedsolution) SolveMixedProblem(cmesh_HDiv,config);
//
//    
//    meshvec_HDiv = cmesh_HDiv->MeshVector();
//    
//    //cria malha hibrida
//    
//    TPZHybridizeHDiv hybrid;
//    auto HybridMesh = hybrid.Hybridize(cmesh_HDiv);
//    HybridMesh->CleanUpUnconnectedNodes();//enumerar adequadamente os connects
//    HybridMesh->AdjustBoundaryElements();
//    delete cmesh_HDiv;
//    delete meshvec_HDiv[0];
//    delete meshvec_HDiv[1];
//    
//    
//    std::cout<<"---Original PerifericalMaterialId --- "<<std::endl;
//    std::cout <<" LagrangeInterface = "<<hybrid.fLagrangeInterface<<std::endl;
//    std::cout <<" HDivWrapMatid = "<<hybrid.fHDivWrapMatid<<std::endl;
//    std::cout <<" InterfaceMatid = "<<hybrid.fInterfaceMatid<<std::endl;
//    
//    
//    cmesh_HDiv=(HybridMesh);//malha hribrida
//    meshvec_HDiv[0] = (HybridMesh)->MeshVector()[0];//malha Hdiv
//    meshvec_HDiv[1] = (HybridMesh)->MeshVector()[1];//malha L2

//   #ifdef PZDEBUG
//    {
//        
//        std::ofstream out2("OriginalFluxMesh.txt");
//        meshvec_HDiv[0]->Print(out2);
//
//        std::ofstream out3("OriginalPotentialMesh.txt");
//        meshvec_HDiv[1]->Print(out3);
//
//    }
//#endif
//
//
//    SolveHybridProblem(cmesh_HDiv,hybrid.fInterfaceMatid,config);
//
//#ifdef PZDEBUG
//    {
//        std::ofstream out("OriginalHybridMesh.txt");
//        (HybridMesh)->Print(out);
//    }
//#endif
 //   PlotLagrangreMultiplier(meshvec_HDiv[1],config);

    
    
    //reconstroi potencial e calcula o erro
    {
        //TPZHybridHDivErrorEstimator HDivEstimate(*cmesh_HDiv);
        TPZHybridHDivErrorEstimator HDivEstimate(*hybridmesh);
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

