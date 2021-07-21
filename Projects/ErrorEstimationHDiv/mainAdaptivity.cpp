//
// Created by Gustavo Batistela on 12/07/19.
//

#include "ProblemConfig.h"
#include "TPZHDivErrorEstimator.h"
#include "TPZRefPatternDataBase.h"
#include "Tools.h"

constexpr bool postProcessWithHDiv = false;
constexpr int refinementSteps = 7;

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