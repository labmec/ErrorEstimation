//
// Created by Gustavo Batistela on 12/07/19.
//

#include "TPZRefPatternDataBase.h"
#include "TPZGmshReader.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "TPZGeoLinear.h"

#include <cmath>
#include <tuple>

#include "Tools.h"
#include "ProblemConfig.h"
#include "TPZHDivErrorEstimatorH1.h"
#include "TPZHybridHDivErrorEstimator.h"

//#include "pzelchdiv.h"

bool readGeoMeshFromFile = false;
bool postProcessWithHDiv = false;
int refinementSteps = 13;

void TracingTriangleBug(TPZMultiphysicsCompMesh* multiphysics);

int main(int argc, char* argv[]) {

#ifdef LOG4CXX
    InitializePZLOG();
#endif

    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);

    // Creates geometric mesh

    TPZGeoMesh *gmeshOriginal =
        nullptr; // CreateLShapeMesh(bcIDs);//CreateLShapeMesh(bcIDs);

    ProblemConfig config;

    config.porder = 1;
    config.hdivmais = 1;

    config.dimension = 2;
    config.makepressurecontinuous = true;

    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;

    config.dir_name = "TriangularLShapeMesh";
    config.problemname = "ESinSinMark_UniRef";

    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    if (readGeoMeshFromFile) {
        std::string meshfilename = "../LMesh.msh"; //"../LMesh3.msh";

        TPZGmshReader gmsh;
        gmsh.GetDimNamePhysical()[1]["dirichlet"] = 2;
        gmsh.GetDimNamePhysical()[2]["domain"] = 1;
        gmsh.SetFormatVersion("4.1");
        gmeshOriginal = gmsh.GeometricGmshMesh(meshfilename);
        gmsh.PrintPartitionSummary(std::cout);
        config.materialids.insert(1);
        config.bcmaterialids.insert(2);

    }
else

    {
        TPZManVector<int, 4> bcIDs(8, -1);
        gmeshOriginal = CreateLShapeMesh(bcIDs);//CreateGeoMesh(1, bcIDs);//
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        //UniformRefinement(3, gmeshOriginal);
        
    }
    
    gmeshOriginal->SetDimension(config.dimension);
    config.gmesh = gmeshOriginal;
    
    
        
        
    TPZManVector<TPZCompMesh*, 2> meshvec_HDiv(2, 0);
        
    TPZMultiphysicsCompMesh* cmesh_HDiv = CreateHDivMesh(config); //Hdiv x L2
    cmesh_HDiv->InitializeBlock();
    SolveMixedProblem(cmesh_HDiv, config);
    
    meshvec_HDiv = cmesh_HDiv->MeshVector();
    
    
    for (int iSteps = 1; iSteps < refinementSteps; iSteps++) {
    
    
        config.adaptivityStep = iSteps;
        
        UniformRefinement(iSteps, gmeshOriginal);
        
    
        

              TPZManVector<TPZCompMesh*, 2> meshvec_HDiv(2, 0);
              
              TPZMultiphysicsCompMesh* cmesh_HDiv = nullptr;
              
              
              cmesh_HDiv = CreateHDivMesh(config);//Hdiv x L2
              cmesh_HDiv->InitializeBlock();
               #ifdef PZDEBUG2
              {
                  std::ofstream out("MultiPhysicsMesh.txt");
                  cmesh_HDiv->Print(out);
                  std::ofstream outvtk("MultiPhysicsMesh.vtk");
                  
                  TPZVTKGeoMesh::PrintCMeshVTK(cmesh_HDiv,outvtk);
        
                  
              }
              #endif
              
              meshvec_HDiv = cmesh_HDiv->MeshVector();
              
              //cria malha hibrida
              TPZHybridizeHDiv hybrid;
              auto HybridMesh = hybrid.Hybridize(cmesh_HDiv);
              HybridMesh->CleanUpUnconnectedNodes();
              HybridMesh->AdjustBoundaryElements();
              
              delete cmesh_HDiv;
              delete meshvec_HDiv[0];
              delete meshvec_HDiv[1];
              
              cmesh_HDiv = (HybridMesh);//malha hribrida
              meshvec_HDiv[0] = (HybridMesh)->MeshVector()[0];//malha Hdiv
              meshvec_HDiv[1] = (HybridMesh)->MeshVector()[1];//malha L2
              
              SolveHybridProblem(cmesh_HDiv, hybrid.fInterfaceMatid, config,false);
    
    
   
        //reconstroi potencial e calcula o erro
        {
            TPZHybridHDivErrorEstimator HDivEstimate(*cmesh_HDiv);
            HDivEstimate.fProblemConfig = config;
            HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
            HDivEstimate.SetAnalyticSolution(config.exact);
            HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
            HDivEstimate.fPostProcesswithHDiv = postProcessWithHDiv;
        
            HDivEstimate.PotentialReconstruction();
            TPZManVector<REAL> elementerrors;
            HDivEstimate.ComputeErrors(elementerrors);
            hAdaptivity(&HDivEstimate.fPostProcMesh, gmeshOriginal, config);
    }
    
   // return 0;
        
    }
}
