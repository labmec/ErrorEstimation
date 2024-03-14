// Comparison between HybridH1 and Mixed Poisson
// Modified from main_HybridH1.cpp
//
// Created by victor on 28/04/2020.
//

#include "InputTreatment.h"
#include "Solver.h"
#include "Output.h"
#include "tpzgeoelrefpattern.h"
#include "DataStructure.h"
#include "Tools.h"
#include <tuple>

int main(int argc, char *argv[]) {

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    PreConfig pConfig;
    pConfig.k = 2; // Porder for H1
    pConfig.n = 2;
    pConfig.problem = "ESteepWave";         //// {"ESinSin","EArcTan",ESteklovNonConst", "EBubble2D", "ELaplace","ESing2D, "EProb","ESinMarkHom", "EBubble2DTemp"}
    pConfig.integrationorder = 11;
    pConfig.maxIter = 15;                     //// Maximum iterations for computing the exact solution (only for ELaplace)
    pConfig.approx = "H1";                 //// {"H1","Hybrid", "Mixed"}
    pConfig.topology = "Quadrilateral";        //// Triangular, Quadrilateral, LQuad, Tetrahedral, Hexahedral, Prism
    pConfig.refLevel = 5;                      //// How many uniform refinements
    pConfig.numberAdapativitySteps = 0;        //// Maximum number of adapativity refinement steps.
    pConfig.estimateError = true;              //// Wheater Error Estimation procedure is invoked
    pConfig.debugger = false;                   //// Print geometric and computational mesh for the simulation (Error estimate not involved).
    pConfig.vtkResolution = 3;                 //// Vtk resolution. Set 0 to see a paraview mesh equals the  simulation mesh.

    // this is where the type in pConfig is set
    EvaluateEntry(argc,argv,pConfig);
    InitializeOutstream(pConfig,argv);

    ProblemConfig config;
    config.division_threshold = 0.5;
    for (config.adaptivityStep = 0; config.adaptivityStep < pConfig.numberAdapativitySteps+1; config.adaptivityStep++) { //ndiv = 1 corresponds to a 2x2 mesh.
        pConfig.h = 1./pConfig.exp;
        
        Configure(config,pConfig.refLevel,pConfig,argv);
        
        Solve(config,pConfig);
    }
    
    std::string command = "cp Erro.txt " + pConfig.plotfile + "/Erro.txt";
    system(command.c_str());
    FlushTable(pConfig,argv);

    return 0.;
}
