//
//  main.cpp
//  ErrorEstimateHDiv
//
//  Created by Denise De Siqueira on 01/04/19.
//

#include "TPZGmshReader.h"
#include "TPZRefPatternDataBase.h"
#include "pzlog.h"
#include "tpzgeoelrefpattern.h"
#include "ProblemConfig.h"
#include "pzbndcond.h"
#include "TPZHybridizeHDiv.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZHybridHDivErrorEstimator.h"
#include "Tools.h"
#include "TPZBFileStream.h"
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

    //RunSingularProblemHDiv();
    RunHPQuadProblemHDiv();
    //RunHPCubeProblemHDiv();

    return 0;
}

void RunSingularProblemHDiv() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinMark;
    config.problemname = "SinMarkLShapeCircle";
    config.dir_name = "ErrorResults";
    config.porder = 1;
    config.hdivmais = 3;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nRef = 0;

    config.ndivisions = nRef;
    config.gmesh = CreateLShapeGeoMesh(nRef);

    std::string command = "mkdir " + config.dir_name;
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

    cout << "---Original PerifericalMaterialId --- " << endl;
    cout << " LagrangeInterface = " << hybridizer.fLagrangeInterface << endl;
    cout << " HDivWrapMatid = " << hybridizer.fHDivWrapMatid << endl;
    cout << " InterfaceMatid = " << hybridizer.fInterfaceMatid << endl;
    return hybridMesh;
}

void EstimateError(ProblemConfig &config, TPZMultiphysicsCompMesh *cmesh_HDiv, TPZHybridizeHDiv &hybrid) {

    TPZHybridHDivErrorEstimator HDivEstimate(*cmesh_HDiv);
    HDivEstimate.SetHybridizer(hybrid);
    HDivEstimate.fProblemConfig = config;
    HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
    HDivEstimate.SetAnalyticSolution(config.exact);

    HDivEstimate.fPostProcesswithHDiv = false;

    HDivEstimate.PotentialReconstruction();

    TPZManVector<REAL> elementerrors;
    TPZVec<REAL> errorVec;
    HDivEstimate.ComputeErrors(errorVec, elementerrors, true);
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
    config.hdivmais = 3;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nElem = 2;
    config.ndivisions = nElem;

    TPZManVector<int, 4> bcIDs(4, -1);
    config.gmesh = Tools::CreateGeoMesh(nElem, bcIDs);

   // Tools::RefineElements(config.gmesh, {1, 3});

    string command = "mkdir " + config.dir_name;
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

    string command = "mkdir " + config.dir_name;
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
