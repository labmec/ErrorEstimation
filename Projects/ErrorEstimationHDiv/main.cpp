//
//  main.cpp
//  ErrorEstimateHDiv
//
//  Created by Denise De Siqueira on 01/04/19.
//

#include "ProblemConfig.h"
#include "TPZBFileStream.h"
#include "TPZGmshReader.h"
#include "TPZHDivErrorEstimator.h"
#include "TPZHybridizeHDiv.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZRefPatternDataBase.h"
#include "Tools.h"
#include "pzbndcond.h"
#include "pzlog.h"
#include "tpzgeoelrefpattern.h"
#include <tuple>

TPZGeoMesh *CreateLShapeGeoMesh(int nCoarseRef);

void EstimateError(ProblemConfig &config, TPZMultiphysicsCompMesh *cmesh_HDiv, TPZHybridizeHDiv &hybrid);

TPZMultiphysicsCompMesh *CreateHybridCompMesh(const ProblemConfig &config, TPZHybridizeHDiv &hybridizer);

void RunSingularProblemHDiv();
void RunHPQuadProblemHDiv();
void RunHPCubeProblemHDiv();

int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeAllUniformRefPatterns();

    RunSingularProblemHDiv();
    //RunHPQuadProblemHDiv();
    //RunHPCubeProblemHDiv();

    return 0;
}

void RunSingularProblemHDiv() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::EArcTan;
    config.problemname = "SinMarkLShapeCircle";
    config.dir_name = "ErrorResults";
    config.porder = 2;
    config.hdivmais = 3;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nRef = 2;

    config.ndivisions = 0;
    config.gmesh = CreateLShapeGeoMesh(nRef);

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    TPZHybridizeHDiv hybridizer;
    TPZMultiphysicsCompMesh *hybridMesh = CreateHybridCompMesh(config, hybridizer);

    Tools::SolveHybridProblem(hybridMesh, hybridizer.fInterfaceMatid, config, false);

    // Reconstruct potential and estimate error
    EstimateError(config, hybridMesh, hybridizer);

    delete hybridMesh->MeshVector()[0];
    delete hybridMesh->MeshVector()[1];
    delete hybridMesh;
}

TPZMultiphysicsCompMesh *CreateHybridCompMesh(const ProblemConfig &config, TPZHybridizeHDiv &hybridizer) {

    TPZMultiphysicsCompMesh *cmesh_HDiv = Tools::CreateHDivMesh(config); // Hdiv x L2

#ifdef PZDEBUG
    {
        ofstream out("MixedMesh.txt");
        cmesh_HDiv->Print(out);
    }
#endif

    // Hybridize mesh
    TPZMultiphysicsCompMesh* hybridMesh = hybridizer.Hybridize(cmesh_HDiv);
    hybridMesh->CleanUpUnconnectedNodes(); // Enumereate connects correctly
    hybridMesh->AdjustBoundaryElements();

    delete cmesh_HDiv;

    std::cout << "---Original PerifericalMaterialId --- " << std::endl;
    std::cout << " LagrangeInterface = " << hybridizer.fLagrangeInterface << std::endl;
    std::cout << " HDivWrapMatid = " << hybridizer.fHDivWrapMatid << std::endl;
    std::cout << " InterfaceMatid = " << hybridizer.fInterfaceMatid << std::endl;
    return hybridMesh;
}

void EstimateError(ProblemConfig &config, TPZMultiphysicsCompMesh *cmesh_HDiv, TPZHybridizeHDiv &hybrid) {
    bool postProcWithHDiv = false;
    TPZHDivErrorEstimator HDivEstimate(*cmesh_HDiv, postProcWithHDiv);
    //HDivEstimate.SetHybridizer(hybrid);
    HDivEstimate.SetProblemConfig(config);

    HDivEstimate.SetPostProcUpliftOrder(config.hdivmais);
    HDivEstimate.SetAnalyticSolution(config.exact);

    HDivEstimate.PotentialReconstruction();

    TPZManVector<REAL> elementerrors;
    TPZVec<REAL> errorVec;
    std::string outVTK = "postProcErrors.vtk";
    HDivEstimate.ComputeErrors(errorVec, elementerrors, outVTK);
}

TPZGeoMesh *CreateLShapeGeoMesh(int nCoarseRef) {

    TPZVec<int> bcIDs(8, -1);
    TPZGeoMesh *gmesh = Tools::CreateQuadLShapeMesh(bcIDs);
    gmesh->SetDimension(2);
    gmesh->BuildConnectivity();

    Tools::UniformRefinement(nCoarseRef, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    return gmesh;
}

void RunHPQuadProblemHDiv() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
    config.problemname = "ESinSin";
    config.dir_name = "HP-ESinSin";
    config.porder = 1;
    config.hdivmais = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nElem = 3;
    config.ndivisions = nElem;

    TPZManVector<int, 4> bcIDs(4, -1);
    config.gmesh = Tools::CreateGeoMesh(nElem, bcIDs);

    //Tools::RefineElements(config.gmesh, {1, 3});
    //Tools::RefineElements(config.gmesh, {12});

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    TPZHybridizeHDiv hybridizer;
    TPZMultiphysicsCompMesh *hybridMesh = CreateHybridCompMesh(config, hybridizer);

    Tools::SolveHybridProblem(hybridMesh, hybridizer.fInterfaceMatid, config, false);

    // Reconstruct potential and estimate error
    EstimateError(config, hybridMesh, hybridizer);

    delete hybridMesh->MeshVector()[0];
    delete hybridMesh->MeshVector()[1];
    delete hybridMesh;
}

void RunHPCubeProblemHDiv() {
    ProblemConfig config;
    config.dimension = 3;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
    config.problemname = "SinSinHP";
    config.dir_name = "HP_3D";
    config.porder = 1;
    config.hdivmais = 3;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nElem = 2;
    config.ndivisions = nElem;

    TPZManVector<int, 4> bcIDs(6, -1);
    TPZManVector<int, 3> nelDiv(3, nElem);
    nelDiv[2] = 1;
    config.gmesh = Tools::CreateCubeGeoMesh(nelDiv, bcIDs);

    Tools::RefineElements(config.gmesh, {1, 3});

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    TPZHybridizeHDiv hybridizer;
    TPZMultiphysicsCompMesh *hybridMesh = CreateHybridCompMesh(config, hybridizer);

    Tools::SolveHybridProblem(hybridMesh, hybridizer.fInterfaceMatid, config, false);

    // Reconstruct potential and estimate error
    EstimateError(config, hybridMesh, hybridizer);

    delete hybridMesh->MeshVector()[0];
    delete hybridMesh->MeshVector()[1];
    delete hybridMesh;
}
