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

#ifdef LOG4CXX
    InitializePZLOG();
#endif
    PreConfig pConfig;
    pConfig.k = 1;
    pConfig.n = 2;
    pConfig.problem = "ESinSin";               //// {"ESinSin","EArcTan",ESteklovNonConst"}
    pConfig.approx = "Hybrid";                 //// {"H1","Hybrid", "Mixed"}
    pConfig.topology = "Quadrilateral";        //// Triangular, Quadrilateral, Tetrahedral, Hexahedral, Prism
    pConfig.refLevel = 1;                      //// How many refinements
    pConfig.estimateError = true;              //// Wheater Error Estimation procedure is invoked
    pConfig.debugger = true;                  //// Print geometric and computational mesh

    EvaluateEntry(argc,argv,pConfig);
    InitializeOutstream(pConfig,argv);
//#define PZDEBUG

    for (int ndiv = 1; ndiv < /*pConfig.refLevel+1*/3; ndiv++) {     //ndiv = 1 corresponds to a 2x2 mesh.
        pConfig.h = 1./pConfig.exp;
        ProblemConfig config;
        Configure(config,ndiv,pConfig,argv);

        Solve(config,pConfig);

        pConfig.hLog = pConfig.h;
        pConfig.exp *=2;
    }
    std::string command = "cp Erro.txt " + pConfig.plotfile + "/Erro.txt";
    system(command.c_str());
    FlushTable(pConfig,argv);

    return 0.;
}