//
// mainQuadMesh
// ErrorEstimationHDiv
//
// This file performs different simulations and error estimations on a basic 2D
// quadrilateral mesh. The purpose of this target is to test different boundary
// conditions, polynomial orders and exact functions with different error
// estimators. It shouldn't perform adaptivity at this point.
//
// Created by Gustavo Batistela on 28/04/20.
//

#include "ProblemConfig.h"
#include "TPZBFileStream.h"
#include "TPZHDivErrorEstimatorH1.h"
#include "TPZHybridHDivErrorEstimator.h"
#include "TPZHybridizeHDiv.h"
#include "TPZVTKGeoMesh.h"
#include "Tools.h"
#include "pzlog.h"
#include "tpzgeoelrefpattern.h"
#include <tuple>

// TODO: expand to include homogeneous BCs, mixed cases, etc.
enum class BCType { Dirichlet, Neumann, Robin };
enum class ErrorEstimator { HybridHDiv, MarkHDiv };

void RunProblem(int nDivisions, int pOrder, int hDivPlusPlus,
                TLaplaceExample1::EExactSol exactSol, BCType bcType,
                ErrorEstimator estimator);

std::string ToString(const BCType &bcType);
std::string ToString(const ErrorEstimator &errorEstimator);
std::string ToString(const TLaplaceExample1::EExactSol &exactSol);

void ConfigureGeoMesh(int nDivisions, const BCType &bcType,
                      ProblemConfig &config);

void RunHybridHDivEstimation(ProblemConfig &config,
                             TPZMultiphysicsCompMesh *mixedMesh,
                             TPZHybridizeHDiv &hybrid);

void RunMarkHDivEstimation(ProblemConfig &config,
                           TPZMultiphysicsCompMesh *mixedMesh);
int main(int argc, char *argv[]) {

    const int nDivisions = 4; // Creates a 'nDivisions' x 'nDivisions' mesh
    const int pOrder = 1;
    const int hDivPlusPlus = 1;
    const TLaplaceExample1::EExactSol exactSol =
        // TLaplaceExample1::ENone
        // TLaplaceExample1::EConst
        // TLaplaceExample1::EX
        TLaplaceExample1::ESinSin
        // TLaplaceExample1::ECosCos
        // TLaplaceExample1::EArcTan
        // TLaplaceExample1::EArcTanSingular
        // TLaplaceExample1::ESinDist
        // TLaplaceExample1::E10SinSin
        // TLaplaceExample1::E2SinSin
        // TLaplaceExample1::ESinSinDirNonHom
        // TLaplaceExample1::ESinMark
        // TLaplaceExample1::ESteklovNonConst
        // TLaplaceExample1::EGalvisNonConst
        // TLaplaceExample1::EBoundaryLayer
        // TLaplaceExample1::EBubble
        // TLaplaceExample1::ESinCosCircle
        ;

    const BCType bcType =
        BCType::Dirichlet
        // BCType::Neumann
        // BCType::Robin
        ;

    const ErrorEstimator estimator =
        ErrorEstimator::HybridHDiv
        // ErrorEstimator::MarkHDiv
        ;

    RunProblem(nDivisions, pOrder, hDivPlusPlus, exactSol, bcType, estimator);

    return 0;
}

void RunProblem(const int nDivisions, const int pOrder, const int hDivPlusPlus,
                const TLaplaceExample1::EExactSol exactSol, const BCType bcType,
                const ErrorEstimator estimator) {

#ifdef LOG4CXX
    InitializePZLOG();
#endif

    ProblemConfig config;
    config.porder = pOrder;
    config.hdivmais = hDivPlusPlus;

    config.dimension = 2;
    config.makepressurecontinuous = true;

    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = exactSol;

    config.problemname = ToString(exactSol) + "_" + ToString(bcType);

    config.dir_name = "Results_" + ToString(estimator);
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    ConfigureGeoMesh(nDivisions, bcType, config);

    TPZMultiphysicsCompMesh *mixedMesh = CreateHDivMesh(config);
    mixedMesh->InitializeBlock();

#ifdef PZDEBUG
    {
        std::string fileName = config.dir_name + "/MixedMesh.txt";
        std::ofstream outTXT(fileName);
        mixedMesh->Print(outTXT);
    }
#endif

    TPZManVector<TPZCompMesh *, 2> meshVector(2, 0);
    meshVector = mixedMesh->MeshVector();

    // Hybridizes mixed mesh
    TPZHybridizeHDiv hybrid;
    auto HybridMesh = hybrid.Hybridize(mixedMesh);
    HybridMesh->CleanUpUnconnectedNodes();
    HybridMesh->AdjustBoundaryElements();

    delete mixedMesh;
    mixedMesh = HybridMesh;

    meshVector[0] = HybridMesh->MeshVector()[0]; // Flux mesh
    meshVector[1] = HybridMesh->MeshVector()[1]; // Pressure mesh

#ifdef PZDEBUG
    std::cout << ":: Original Periferical Material Id" << '\n';
    std::cout << "  LagrangeInterface = " << hybrid.fLagrangeInterface << '\n';
    std::cout << "  HDivWrap = " << hybrid.fHDivWrapMatid << '\n';
    std::cout << "  Interface = " << hybrid.fInterfaceMatid << std::endl;
#endif

    // Solves finite element problem
    SolveHybridProblem(mixedMesh, hybrid.fInterfaceMatid, config, true);

#ifdef PZDEBUG
    {
        std::string fileName = config.dir_name + "/HybridizedMixedMesh.txt";
        std::ofstream outTXT(fileName);
        mixedMesh->Print(outTXT);
    }
#endif

    switch (estimator) {
    case ErrorEstimator::HybridHDiv:
        RunHybridHDivEstimation(config, mixedMesh, hybrid);
        break;
    case ErrorEstimator::MarkHDiv:
        RunMarkHDivEstimation(config, mixedMesh);
        break;
    default:
        DebugStop();
    }
}

void RunMarkHDivEstimation(ProblemConfig &config,
                           TPZMultiphysicsCompMesh *mixedMesh) {
    TPZHDivErrorEstimatorH1 HDivEstimate(*mixedMesh);
    HDivEstimate.fProblemConfig = config;
    HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
    HDivEstimate.SetAnalyticSolution(config.exact);
    HDivEstimate.fperformUplift = true;
    HDivEstimate.fUpliftOrder = config.hdivmais;

    HDivEstimate.PotentialReconstruction();

    TPZManVector<REAL> elementerrors;
    HDivEstimate.ComputeErrors(elementerrors);
}

void RunHybridHDivEstimation(ProblemConfig &config,
                             TPZMultiphysicsCompMesh *mixedMesh,
                             TPZHybridizeHDiv &hybrid) {

    TPZHybridHDivErrorEstimator HDivEstimate(*mixedMesh);
    HDivEstimate.SetHybridizer(hybrid);
    HDivEstimate.fProblemConfig = config;
    HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
    HDivEstimate.SetAnalyticSolution(config.exact);
    HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;

    HDivEstimate.fPostProcesswithHDiv = false;

    HDivEstimate.PotentialReconstruction();

    TPZManVector<REAL> elementerrors;
    HDivEstimate.ComputeErrors(elementerrors);
}

void ConfigureGeoMesh(const int nDivisions, const BCType &bcType,
                      ProblemConfig &config) { // Creates geometric mesh
    int bcId = 0;
    switch (bcType) {
    case BCType::Dirichlet:
        bcId = -1;
        break;
    case BCType::Neumann:
        bcId = -2;
        break;
    case BCType::Robin:
        bcId = -3;
        break;
    default:
        DebugStop();
    }

    TPZManVector<int, 4> bcids(4, bcId);

    config.gmesh = CreateGeoMesh(nDivisions, bcids);
    config.materialids.insert(1);
    config.bcmaterialids.insert(bcId);

#ifdef PZDEBUG
    {
        std::string fileName = config.dir_name + "/gmesh.vtk";
        ofstream outVTK(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, outVTK);
        fileName = config.dir_name + "/gmesh.txt";
        ofstream outTXT(fileName);
        config.gmesh->Print(outTXT);
    }
#endif
}

std::string ToString(const BCType &bcType) {
    switch (bcType) {
    case BCType::Dirichlet:
        return "Dirichlet";
    case BCType::Neumann:
        return "Neumann";
    case BCType::Robin:
        return "Robin";
    default:
        return "[Unknown BCType]";
    }
}

std::string ToString(const ErrorEstimator &errorEstimator) {
    switch (errorEstimator) {
    case ErrorEstimator::HybridHDiv:
        return "HybridHDiv";
    case ErrorEstimator::MarkHDiv:
        return "MarkHDiv";
    default:
        return "[Unknown Error Estimator]";
    }
}

std::string ToString(const TLaplaceExample1::EExactSol &exactSol) {
    switch (exactSol) {
    case TLaplaceExample1::ENone:
        return "ENone";
    case TLaplaceExample1::EConst:
        return "EConst";
    case TLaplaceExample1::EX:
        return "EX";
    case TLaplaceExample1::ESinSin:
        return "ESinSin";
    case TLaplaceExample1::ECosCos:
        return "ECosCos";
    case TLaplaceExample1::EArcTan:
        return "EArcTan";
    case TLaplaceExample1::EArcTanSingular:
        return "EArcTanSingular";
    case TLaplaceExample1::ESinDist:
        return "ESinDist";
    case TLaplaceExample1::E10SinSin:
        return "E10SinSin";
    case TLaplaceExample1::E2SinSin:
        return "E2SinSin";
    case TLaplaceExample1::ESinSinDirNonHom:
        return "ESinSinDirNonHom";
    case TLaplaceExample1::ESinMark:
        return "ESinMark";
    case TLaplaceExample1::ESteklovNonConst:
        return "ESteklovNonConst";
    case TLaplaceExample1::EGalvisNonConst:
        return "EGalvisNonConst";
    case TLaplaceExample1::EBoundaryLayer:
        return "EBoundaryLayer";
    case TLaplaceExample1::EBubble:
        return "EBubble";
    case TLaplaceExample1::ESinCosCircle:
        return "ESinCosCircle";
    default:
        return "[Unknown Exact Solution]";
    }
}
