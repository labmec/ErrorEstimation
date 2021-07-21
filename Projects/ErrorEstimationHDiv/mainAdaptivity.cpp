//
// Created by Gustavo Batistela on 12/07/19.
//

#include "TPZRefPatternDataBase.h"
#include "tpzarc3d.h"
#include <cmath>
#include "ProblemConfig.h"
#include "TPZHDivErrorEstimator.h"
#include "Tools.h"

constexpr bool postProcessWithHDiv = false;
constexpr int refinementSteps = 4;

TPZMultiphysicsCompMesh *CreateHybridCompMesh(const ProblemConfig &config, TPZHybridizeHDiv &hybridizer);
TPZGeoMesh *CreateLShapeGeoMesh(int nCoarseRef, int nInternalRef, TPZStack<int64_t> &mhmIndexes);

int main() {

    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);

    // Creates geometric mesh
    TPZGeoMesh *gmeshOriginal;

    ProblemConfig config;

    config.porder = 1;
    config.hdivmais = 1;

    config.dimension = 2;
    config.makepressurecontinuous = true;

    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinMark;

    config.dir_name = "AdaptivityLShape";
    config.problemname = "ESinSinMark";

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    TPZManVector<int, 4> bcids(8, -1);
    bcids[1] = -1;
    gmeshOriginal = Tools::CreateLShapeMesh(bcids);
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);

    Tools::UniformRefinement(2, 2 , gmeshOriginal) ;
    Tools::DivideLowerDimensionalElements(gmeshOriginal);

    for (int iSteps = 0; iSteps < refinementSteps; iSteps++) {

        config.gmesh = new TPZGeoMesh(*gmeshOriginal);
        config.adaptivityStep = iSteps;

        TPZMultiphysicsCompMesh *mixedCompMesh = Tools::CreateMixedMesh(config); // Hdiv x L2
        mixedCompMesh->InitializeBlock();
        Tools::SolveMixedProblem(mixedCompMesh, config);

        // Estimate error and run adaptive process
        {
            TPZHDivErrorEstimator HDivEstimate(*mixedCompMesh, postProcessWithHDiv);
            HDivEstimate.SetAnalyticSolution(config.exact);

            HDivEstimate.PotentialReconstruction();
            TPZManVector<REAL> elementerrors;
            TPZManVector<REAL> errorvec;
            std::string vtkPath = "adaptivity_error_results.vtk";
            HDivEstimate.ComputeErrors(errorvec, elementerrors, vtkPath);
            Tools::hAdaptivity(HDivEstimate.PostProcMesh(), gmeshOriginal, config);
        }
    }
    return 0;
}

TPZGeoMesh *CreateLShapeGeoMesh(int nCoarseRef, int nInternalRef, TPZStack<int64_t> &mhmIndexes) {

    TPZVec<int> bcIDs(8, -1);
    TPZGeoMesh *gmesh = Tools::CreateQuadLShapeMesh(bcIDs);
    gmesh->SetDimension(2);
    gmesh->BuildConnectivity();

    Tools::UniformRefinement(nCoarseRef, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    int64_t nElem = gmesh->NElements();
    for (int64_t i = 0; i < nElem; i++) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (gel->Dimension() != gmesh->Dimension() || gel->NSubElements() > 0) continue;
        mhmIndexes.Push(i);
    }

    Tools::UniformRefinement(nInternalRef, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    for (int64_t i = 0; i < mhmIndexes.size(); i++) {
        std::cout << mhmIndexes[i] << '\n';
    }
    std::cout << '\n';

    return gmesh;
}

TPZMultiphysicsCompMesh *CreateHybridCompMesh(const ProblemConfig &config, TPZHybridizeHDiv &hybridizer) {

    TPZMultiphysicsCompMesh *cmesh_HDiv = Tools::CreateMixedMesh(config); // Hdiv x L2

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
