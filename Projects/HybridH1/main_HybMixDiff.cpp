//
// Created by victor on 29/03/2021.
//
// Runs Mixed and Hybrid simulation and compute its difference

// Implementation plan:
// Run both simulations, Create material with both meshes and compute difference

#include "InputTreatment.h"
#include "Solver.h"
#include "Output.h"
#include "tpzgeoelrefpattern.h"
#include "DataStructure.h"
#include "Tools.h"
#include <tuple>

int main(int argc, char *argv[]) {

#ifdef LOG4CXX
    InitializePZLOG();
#endif
    PreConfig hybConfig, mixConfig;
    hybConfig.k = 1;
    hybConfig.n = 2;
    hybConfig.problem = "EBubble2D";                //// {"ESinSin","EArcTan",ESteklovNonConst"}
    hybConfig.approx = "Hybrid";                 //// {"H1","Hybrid", "Mixed"}
    hybConfig.topology = "Quadrilateral";        //// Triangular, Quadrilateral, Tetrahedral, Hexahedral, Prism

    hybConfig.refLevel = 3;                       //// How many refinements
    hybConfig.debugger = true;                    //// Print geometric and computational mesh

    // Copy the Hybrid set up to the Mixed problem
    CopyHybSetup(hybConfig,mixConfig);

    // Correspondence between strings to integers(ex. "Hybrid" -> 2), and definition of output names
    DataInitialization(argc,argv,hybConfig,mixConfig);

    // Solve hybtid and mixed problems, then compute solution difference
    SolveDiff(hybConfig,mixConfig,argv);

    return 0.;
}