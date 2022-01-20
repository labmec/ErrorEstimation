//
// Created by Gustavo Batistela on 12/07/19.
//

#include "ProblemConfig.h"
#include "TPZHDivErrorEstimator.h"
#include "TPZRefPatternDataBase.h"
#include "Tools.h"

void RunAdaptivitySuite(int refinementSteps);
void RunUniformRefinementSuite(int refinementSteps);

int main() {

    constexpr int adaptivitySteps = 8;
    constexpr int uniformSteps = 6;
    RunAdaptivitySuite(adaptivitySteps);
    RunUniformRefinementSuite(uniformSteps);
}

void RunAdaptivitySuite(int refinementSteps) {

    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);

    // Creates geometric mesh
    TPZGeoMesh *gmeshOriginal;

    ProblemConfig config;

    config.porder = 1;
    config.hdivmais = 2;

    config.dimension = 2;
    config.makepressurecontinuous = true;

    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::EBoundaryLayer;

    config.dir_name = "HDivAdaptivityEquations";
    config.problemname = "Adaptivity";
    {
        std::string command = "mkdir -p " + config.dir_name;
        system(command.c_str());
    }

    TPZManVector<int, 4> bcids(4, -1);
    gmeshOriginal = Tools::CreateGeoMesh(3, bcids);
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);

    Tools::DivideLowerDimensionalElements(gmeshOriginal);

    for (int iStep = 0; iStep < refinementSteps; iStep++) {

        config.gmesh = new TPZGeoMesh(*gmeshOriginal);
        config.adaptivityStep = iStep;

        TPZMultiphysicsCompMesh *mixedCompMesh = Tools::CreateMixedMesh(config); // Hdiv x L2
        mixedCompMesh->InitializeBlock();
        Tools::SolveMixedProblem(mixedCompMesh, config);
        config.cmesh = mixedCompMesh;

        // Estimate error and run adaptive process
        {
            bool postProcWithHDiv = false;
            TPZHDivErrorEstimator HDivEstimate(*mixedCompMesh, postProcWithHDiv);
            HDivEstimate.SetAnalyticSolution(config.exact);
            HDivEstimate.SetAdaptivityStep(iStep);
            HDivEstimate.PotentialReconstruction();
            TPZManVector<REAL> elementerrors;
            TPZManVector<REAL> errorvec;
            std::stringstream outVTK;
            outVTK << config.dir_name << "/" << config.problemname << "-Errors.vtk";
            auto stringVTK = outVTK.str();
            HDivEstimate.ComputeErrors(errorvec, elementerrors, stringVTK);
            Tools::hAdaptivity(HDivEstimate.PostProcMesh(), gmeshOriginal, config);

            std::stringstream outTXT;
            outTXT << config.dir_name << "/" << config.problemname << "-Errors-Step" << config.adaptivityStep << ".txt";
            std::ofstream fileTXT(outTXT.str());
            Tools::PrintErrors(fileTXT, config, errorvec);
        }
    }
}

void RunUniformRefinementSuite(int refinementSteps) {

    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);

    ProblemConfig config;

    config.porder = 1;
    config.hdivmais = 2;

    config.dimension = 2;
    config.makepressurecontinuous = true;

    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::EBoundaryLayer;

    config.dir_name = "HDivAdaptivityEquations";
    config.problemname = "UnifRef";
    {
        std::string command = "mkdir -p " + config.dir_name;
        system(command.c_str());
    }

    for (int iStep = 0; iStep < refinementSteps; iStep++) {

        // Creates geometric mesh
        TPZGeoMesh *gmeshOriginal;

        TPZManVector<int, 4> bcids(4, -1);
        int nelem = 3 * static_cast<int>(pow(2.0, iStep));
        gmeshOriginal = Tools::CreateGeoMesh(nelem, bcids);
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);

        Tools::DivideLowerDimensionalElements(gmeshOriginal);

        config.gmesh = new TPZGeoMesh(*gmeshOriginal);
        config.adaptivityStep = iStep;

        TPZMultiphysicsCompMesh *mixedCompMesh = Tools::CreateMixedMesh(config); // Hdiv x L2
        mixedCompMesh->InitializeBlock();
        Tools::SolveMixedProblem(mixedCompMesh, config);
        config.cmesh = mixedCompMesh;

        {
            bool postProcessWithHDiv = false;
            TPZHDivErrorEstimator HDivEstimate(*mixedCompMesh, postProcessWithHDiv);
            HDivEstimate.SetAnalyticSolution(config.exact);
            HDivEstimate.SetAdaptivityStep(iStep);
            HDivEstimate.PotentialReconstruction();
            TPZManVector<REAL> elementerrors;
            TPZManVector<REAL> errorvec;
            std::stringstream outVTK;
            outVTK << config.dir_name << "/" << config.problemname << "-Errors.vtk";
            auto stringVTK = outVTK.str();
            HDivEstimate.ComputeErrors(errorvec, elementerrors, stringVTK);

            std::stringstream outTXT;
            outTXT << config.dir_name << "/" << config.problemname << "-Errors-Step" << config.adaptivityStep << ".txt";
            std::ofstream fileTXT(outTXT.str());
            Tools::PrintErrors(fileTXT, config, errorvec);
        }
    }
}
