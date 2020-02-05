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

#include "ProblemConfig.h"
#include "Tools.h"
#include "TPZHybridHDivErrorEstimator.h"
#include "TPZHDivErrorEstimatorH1.h"

bool readGeoMeshFromFile = false;
bool postProcessWithHDiv = false;
int refinementSteps = 5;


int main(int argc, char* argv[]) {

#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    // Creates geometric mesh
    TPZManVector<int, 4> bcids(4, -1);
    TPZGeoMesh* gmeshOriginal = CreateGeoMesh(2, bcids);//CreateLCircleGeoMesh();//
    
#ifdef PZDEBUG2
    {
        std::ofstream out("AdaptivityInitialGmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshOriginal, out);
    }
#endif
    
    for (int iSteps = 1; iSteps < refinementSteps; iSteps++) {
        
        ProblemConfig config;
        
        config.porder = 1;
        config.hdivmais = 1;
        config.dimension = 2;
        config.makepressurecontinuous = true;
        config.adaptivityStep = iSteps;
        
        TLaplaceExample1 example;
        config.exact.fExact = example.ESinSinDirNonHom;
        config.dir_name = "TestAdaptivityH1";
        config.problemname = "ESinSinNonHom";
        
        std::string command = "mkdir " + config.dir_name;
        system(command.c_str());
        
        UniformRefinement(iSteps, gmeshOriginal);
        
        //malha geometrica
        TPZGeoMesh* gmesh = new TPZGeoMesh();
        *gmesh = *gmeshOriginal;
        
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        config.gmesh = gmesh;
        
        TPZManVector<TPZCompMesh*, 2> meshvec_HDiv(2, 0);
        
        TPZMultiphysicsCompMesh* cmesh_HDiv = nullptr;
        
        
        cmesh_HDiv = CreateHDivMesh(config);//Hdiv x L2
        cmesh_HDiv->InitializeBlock();
        
       // SolveMixedProblem(cmesh_HDiv, config);
        
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
           // hAdaptivity(&HDivEstimate.fPostProcMesh, gmeshOriginal);
        }
    }
    
    return 0;
}

