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
int refinementSteps = 3;


int main(int argc, char* argv[]) {

#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    // Creates geometric mesh
    TPZManVector<int, 4> bcIDs(8, -1);
    TPZGeoMesh* gmeshOriginal = CreateQuadLShapeMesh(bcIDs);//CreateLShapeMesh(bcIDs);
    
    ProblemConfig config;
    config.porder = 1;
    config.hdivmais = 1;
    config.dimension = 2;
    config.makepressurecontinuous = true;
    
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    //UniformRefinement(1, gmeshOriginal);
    
    gmeshOriginal->SetDimension(config.dimension);
    config.gmesh = gmeshOriginal;

#ifdef PZDEBUG
    {
        std::ofstream out("InitialGmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshOriginal, out);
    }
#endif
    
    
    TLaplaceExample1 example;
    config.exact.fExact = example.EConst;
    config.dir_name = "TestDebugMonday";
    config.problemname = "Debug";
    
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());
    
    
    TPZManVector<TPZCompMesh*, 2> meshvec_HDiv(2, 0);
    
    TPZMultiphysicsCompMesh* cmesh_HDiv = CreateHDivMesh(config); //Hdiv x L2
    cmesh_HDiv->InitializeBlock();
    
    {
        std::ofstream out("MultiphysicsMesh.txt");
        cmesh_HDiv->Print(out);
        std::ofstream outvtk("MultiphysicsMesh.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh_HDiv, outvtk);
    }
    
    SolveMixedProblem(cmesh_HDiv, config);
    
    DebugStop();
    
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
    
    SolveHybridProblem(cmesh_HDiv, hybrid.fInterfaceMatid, config, true);
    
    
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
    
    return 0;
}

