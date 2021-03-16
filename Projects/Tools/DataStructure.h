//
// Created by victor on 16/03/2021.
//

#ifndef ERRORESTIMATION_DATASTRUCTURE_H
#define ERRORESTIMATION_DATASTRUCTURE_H

#include <set>
#include "TPZAnalyticSolution.h"

struct PreConfig{
    std::ofstream Erro, timer;
    TPZVec<REAL> *rate, *Log;
    int refLevel = -1;

    int k = 1;
    int n = 1;
    int dim = 1;
    int topologyMode = -1;

    std::string problem;
    std::string approx;
    std::string topology;           //Topology' name typed as input
    std::string topologyFileName;   //Simplified name used for naming files/directories

    REAL perm_Q1 = 5;      /// Permeability coefficient of even quadrants (Steklov only)
    REAL perm_Q2 = 1;

    REAL hLog = -1, h = -1000;
    int numErrors = 4;

    std::string plotfile;
    int mode = -1;           // 0 = "H1"; 1 = "Hybrid"; 2 = "Mixed";
    int argc = 1;
    int type= -1;

    bool estimateError;
    bool debugger = true;
    int exp = 2; // Initial exponent of mesh refinement (numElem = 2*2^exp)
};

#endif //ERRORESTIMATION_DATASTRUCTURE_H
