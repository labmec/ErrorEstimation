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
    /// element type to be used
    int topologyMode = -1;

    /// determines the problem that will be run. This string will translate in a value of the "type" data structure
    std::string problem;
    /// determines the the approximation space, will translate in an integer value of mode data (should be an enum)
    std::string approx;
    /// indicates the element topology that should be used
    std::string topology;           //Topology' name typed as input will translate in a value of topologyMode (should be an enum)
    
    std::string topologyFileName;   //Simplified name used for naming files/directories

    REAL perm_Q1 = 5;      /// Permeability coefficient of even quadrants (Steklov only)
    REAL perm_Q2 = 1;

    REAL hLog = -1, h = -1000;
    int numErrors = 4;

    std::string plotfile;
    int mode = -1;           // 0 = "H1"; 1 = "Hybrid"; 2 = "Mixed";
    // argc is the number of arguments passed by the command line. Will be set to 5 if the problem is run from command line
    int argc = 1;
    // the type is an integer indicating which problem will be approximated (it should be an enum)
    int type= -1;
    int maxIter = 15;

    bool estimateError;
    bool debugger = true;
    int exp = 2; // Initial exponent of mesh refinement (numElem = 2*2^exp)
};

#endif //ERRORESTIMATION_DATASTRUCTURE_H
