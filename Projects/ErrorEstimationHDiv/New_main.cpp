//
//  New_main.cpp
//  ErrorEstimateHDiv
//
//  Created by Denise De Siqueira on 01/04/19.
//

#include "TPZGmshReader.h"
#include "TPZRefPatternDataBase.h"
#include "TPZVTKGeoMesh.h"
#include "pzlog.h"

#include "tpzarc3d.h"
#include "tpzgeoelrefpattern.h"
#include "tpzquadraticquad.h"
#include "tpzquadratictrig.h"

#include "ProblemConfig.h"

#include "TPZVecL2.h"
#include "mixedpoisson.h"
#include "pzbndcond.h"

#include "pzintel.h"

#include "TPZCompMeshTools.h"
#include "TPZHybridizeHDiv.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZMultiphysicsCompMesh.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "TPZHDivErrorEstimatorH1.h"
#include "TPZHybridHDivErrorEstimator.h"

#include "Tools.h"

#include "TPZBFileStream.h"
#include <memory>
#include <tuple>

bool IsgmeshReader = false; // para ler a malha
bool neumann = true;        // para o problema local de neumann da forlmulacao Mark

bool mixedsolution = true; // se quiser rodar o problema misto

bool PostProcessingFEM = true; // para graficos da solucao FEM

TPZGeoMesh *CreateLShapeGeoMesh(int nCoarseRef, int nInternalRef, TPZStack<int64_t> &mhmIndexes);
int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(ECube);

    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinMark;
    config.problemname = "SinMarkLShapeCircle";
    config.dir_name = "HDIV_CILAMCE";
    config.porder = 1;
    config.hdivmais = 3;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nCoarseRef = 0;
    int nInternalRef = 0;

    config.ndivisions = nCoarseRef;
    TPZStack<int64_t> mhmIndexes;
    config.gmesh = CreateLShapeGeoMesh(nCoarseRef, nInternalRef, mhmIndexes);

    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    TPZManVector<TPZCompMesh *, 2> meshvec_HDiv(2, 0);

    TPZMultiphysicsCompMesh *cmesh_HDiv = nullptr;

    cmesh_HDiv = CreateHDivMesh(config); // Hdiv x L2
    cmesh_HDiv->InitializeBlock();

#ifdef PZDEBUG
    {

        std::ofstream out2("MalhaMista.txt");
        cmesh_HDiv->Print(out2);
    }
#endif

    //  SolveMixedProblem(cmesh_HDiv,config);

    meshvec_HDiv = cmesh_HDiv->MeshVector();

    // cria malha hibrida

    TPZHybridizeHDiv hybrid;
    auto HybridMesh = hybrid.Hybridize(cmesh_HDiv);
    HybridMesh->CleanUpUnconnectedNodes(); // enumerar adequadamente os connects
    HybridMesh->AdjustBoundaryElements();
    delete cmesh_HDiv;
    delete meshvec_HDiv[0];
    delete meshvec_HDiv[1];

    std::cout << "---Original PerifericalMaterialId --- " << std::endl;
    std::cout << " LagrangeInterface = " << hybrid.fLagrangeInterface << std::endl;
    std::cout << " HDivWrapMatid = " << hybrid.fHDivWrapMatid << std::endl;
    std::cout << " InterfaceMatid = " << hybrid.fInterfaceMatid << std::endl;

    cmesh_HDiv = (HybridMesh);                       // malha hribrida
    meshvec_HDiv[0] = (HybridMesh)->MeshVector()[0]; // malha Hdiv
    meshvec_HDiv[1] = (HybridMesh)->MeshVector()[1]; // malha L2

#ifdef PZDEBUG2
    {

        std::ofstream out2("OriginalFluxMesh.txt");
        meshvec_HDiv[0]->Print(out2);

        std::ofstream out3("OriginalPotentialMesh.txt");
        meshvec_HDiv[1]->Print(out3);
    }
#endif

    SolveHybridProblem(cmesh_HDiv, hybrid.fInterfaceMatid, config, false);


#ifdef PZDEBUG2
    {
        std::ofstream out("OriginalHybridMesh.txt");
        (HybridMesh)->Print(out);
    }
#endif
    //  PlotLagrangeMultiplier(meshvec_HDiv[1],config);

    // reconstroi potencial e calcula o erro
    {

        TPZHybridHDivErrorEstimator HDivEstimate(*cmesh_HDiv);
        HDivEstimate.SetHybridizer(hybrid);
        HDivEstimate.fProblemConfig = config;
        HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
        HDivEstimate.SetAnalyticSolution(config.exact);
        HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;

        HDivEstimate.fPostProcesswithHDiv = false;

        HDivEstimate.PotentialReconstruction();

        TPZManVector<REAL> elementerrors;
        TPZVec<REAL> errorVec;
        HDivEstimate.ComputeErrors(errorVec, elementerrors, true);
        //         hAdaptivity(&HDivEstimate.fPostProcMesh, hybridEstimatorMesh,config);
    }

    delete cmesh_HDiv;
    delete meshvec_HDiv[0];
    delete meshvec_HDiv[1];
    // return 0;
}

TPZGeoMesh *CreateLShapeGeoMesh(int nCoarseRef, int nInternalRef, TPZStack<int64_t> &mhmIndexes) {

    TPZVec<int> bcIDs(8, -1);
    TPZGeoMesh *gmesh = CreateQuadLShapeMesh(bcIDs);
    gmesh->SetDimension(2);
    gmesh->BuildConnectivity();

    UniformRefinement(nCoarseRef, gmesh);
    DivideLowerDimensionalElements(gmesh);

    int64_t nElem = gmesh->NElements();
    for (int64_t i = 0; i < nElem; i++) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (gel->Dimension() != gmesh->Dimension() || gel->NSubElements() > 0) continue;
        mhmIndexes.Push(i);
    }

    UniformRefinement(nInternalRef, gmesh);
    DivideLowerDimensionalElements(gmesh);

    for (int64_t i = 0; i < mhmIndexes.size(); i++) {
        std::cout << mhmIndexes[i] << '\n';
    }
    std::cout << '\n';

    return gmesh;
}
