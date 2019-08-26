//
// Created by Gustavo Batistela on 12/07/19.
//

#include "TPZRefPatternDataBase.h"
#include "TPZGmshReader.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "TPZGeoLinear.h"

#include <math.h>
#include <tuple>

#include "ProblemConfig.h"
#include "Tools.h"
#include "TPZHybridHDivErrorEstimator.h"
#include "TPZHDivErrorEstimatorH1.h"

bool readGeoMeshFromFile = false;

TPZGeoMesh *CreateGeoMesh();
TPZGeoMesh *CreateLCircleGeoMesh();


int main(int argc, char *argv[]) {

#ifdef LOG4CXX
    InitializePZLOG();
#endif

    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);

    // Creates geometric mesh
    TPZGeoMesh *gmeshOriginal = CreateGeoMesh();
    
#ifdef PZDEBUG
    {
        std::ofstream out("Gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshOriginal, out);
        
    }
#endif
    
    

    int refinementSteps = 5;

    // Copies meshes to be used with each proposal
    TPZGeoMesh *hybridEstimatorMesh = new TPZGeoMesh();
    *hybridEstimatorMesh = *gmeshOriginal;
    TPZGeoMesh *markEstimatorMesh = new TPZGeoMesh();
    *markEstimatorMesh = *gmeshOriginal;
    
    
    
#ifdef PZDEBUG
    {
        std::ofstream out("GmeshPosCopy.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(hybridEstimatorMesh, out);
        
    }
#endif

    // Run tests with hdiv++ proposal
    
    for (int i = 0; i < refinementSteps; i++) {
        ProblemConfig config;
        config.dir_name = "AdaptivityHybridSinSin";
        config.adaptivityStep = i;

//        config.gmesh = new TPZGeoMesh();
        config.gmesh = hybridEstimatorMesh;

        config.materialids.insert(1);
        config.bcmaterialids.insert(2);
        config.bcmaterialids.insert(3);

        config.porder = 1;
        config.hdivmais = 1;
        config.dimension = 2;
        config.prefine = false;
        config.makepressurecontinuous = true;

        config.exact.fExact = TLaplaceExample1::ESinMark;
        config.problemname = "AdaptivityTestProposal";

        std::string command = "mkdir " + config.dir_name;
        system(command.c_str());

        TPZMultiphysicsCompMesh *mixedMesh = nullptr;

        mixedMesh = CreateHDivMesh(config); // H(div) x L2
        mixedMesh->InitializeBlock();

        TPZManVector<TPZCompMesh *, 2> mixedMeshVector(2, 0);
        mixedMeshVector = mixedMesh->MeshVector();

        // Hybridizes mixed mesh
        TPZHybridizeHDiv hybrid;
        auto HybridMesh = hybrid.Hybridize(mixedMesh);
        HybridMesh->CleanUpUnconnectedNodes();
        HybridMesh->AdjustBoundaryElements();
        delete mixedMesh;
        delete mixedMeshVector[0];
        delete mixedMeshVector[1];

        mixedMesh = (HybridMesh); // Substitute mixed by hybrid mesh
        mixedMeshVector[0] = (HybridMesh)->MeshVector()[0]; // Hdiv mesh
        mixedMeshVector[1] = (HybridMesh)->MeshVector()[1]; // L2 mesh

        // Solves problem
        SolveHybridProblem(mixedMesh, hybrid.fInterfaceMatid, config);

        // Reconstructs pressure and calculates error
        TPZHybridHDivErrorEstimator HDivEstimate(*mixedMesh);
        HDivEstimate.fProblemConfig = config;
        HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
        HDivEstimate.SetAnalyticSolution(config.exact);

        HDivEstimate.PotentialReconstruction();

        TPZManVector<REAL> elementErrors;
        HDivEstimate.ComputeErrors(elementErrors);

        delete mixedMesh;
        delete mixedMeshVector[0];
        delete mixedMeshVector[1];

        hAdaptivity(&HDivEstimate.fPostProcMesh, hybridEstimatorMesh);
    }

    
    // Run tests with Ainsworth's proposal
    for (int i = 0; i < refinementSteps; i++) {

        ProblemConfig config;
        config.dir_name = "AdaptivityMarkSin";
        config.adaptivityStep = i;

//        config.gmesh = new TPZGeoMesh();
        config.gmesh = markEstimatorMesh;

        config.materialids.insert(1);
        config.bcmaterialids.insert(2);
        config.bcmaterialids.insert(3);

        config.porder = 1;
        config.hdivmais = 1;
        config.dimension = 2;
        config.prefine = false;
        config.makepressurecontinuous = true;

        config.exact.fExact = TLaplaceExample1::ESinMark;
        config.problemname = "AdaptivityTestMark";

        std::string command = "mkdir " + config.dir_name;
        system(command.c_str());
        
//        std::string plotname;
//        {
//            std::stringstream out;
//            out << config.dir_name << "/" << "GMeshT_" << i<<".vtk";
//            plotname = out.str();
//            TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, out);
//
//        }
        TPZMultiphysicsCompMesh *mixedMesh = nullptr;

        mixedMesh = CreateHDivMesh(config); // H(div) x L2
        mixedMesh->InitializeBlock();

        TPZMultiphysicsCompMesh *hybridmesh = HybridSolveProblem(mixedMesh, config);

        // Reconstructs pressure and calculates error
        TPZHDivErrorEstimatorH1 HDivEstimate(*hybridmesh);
        HDivEstimate.fProblemConfig = config;
        HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;

        HDivEstimate.SetAnalyticSolution(config.exact);

        HDivEstimate.fperformUplift = true;
        HDivEstimate.fUpliftOrder = 1;

        HDivEstimate.PotentialReconstruction();

        TPZManVector<REAL> elementErrors;
        HDivEstimate.ComputeErrors(elementErrors);

        hAdaptivity(&HDivEstimate.fPostProcMesh, markEstimatorMesh);
    }

    return 0;
}



TPZGeoMesh *CreateGeoMesh() {
    TPZGeoMesh * gmesh = nullptr;
    if (readGeoMeshFromFile) {
        TPZGmshReader gmsh;
        std::string meshfilename = "LCircle.msh";

        gmsh.GetDimNamePhysical()[1]["dirichlet"] = 2;
        gmsh.GetDimNamePhysical()[2]["domain"] = 1;

        gmsh.SetFormatVersion("4.1");
        gmsh.PrintPartitionSummary(std::cout);

        gmesh = gmsh.GeometricGmshMesh(meshfilename);
        gmesh->SetDimension(2);
    } else {
        gmesh = CreateLCircleGeoMesh();
    }
    int initialRefinement = 1;
    UniformRefinement(initialRefinement, gmesh);

#ifdef PZDEBUG
    {
        std::ofstream out("OriginalGeometricMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
#endif
    return gmesh;
}
