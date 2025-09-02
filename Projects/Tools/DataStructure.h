//
// Created by victor on 16/03/2021.
//

#ifndef ERRORESTIMATION_DATASTRUCTURE_H
#define ERRORESTIMATION_DATASTRUCTURE_H

#include <set>
#include "TPZAnalyticSolution.h"
/// @brief Configuration of the error estimation process
struct PreConfig{
    std::ofstream Erro, timer;
    /// @brief vectors where the error values and estimated convergence rates are stored
    // the memory is allocated in InitializeOutstream
    TPZVec<REAL> *rate = 0, *Log = 0;

    /// @brief This variable seems to be dead weight? It is not used anywhere.
    int refLevel = -1;

    /// @brief Indicates a polynomial order. Not used anywhere.
    int k = 1;
    /// @brief Indicates the increment in polynomial order between the interface order and internal order. Not used anywhere
    int n = 1;
   /// @brief Indicates the dimension of the problem. Is copied to ProblemConfig datastructure at some point
   int dim = 1;

    /// determines the problem that will be run. This string will translate in a value of the "type" data structure. 
    std::string problem;
    // the type is an integer indicating which problem will be approximated (it should be an enum)
    // type is number from 0 to 18
    int type= -1;
    /// integral order of the forcing function. It is set once in the main program
    int integrationorder = -1;
    /// determines the the approximation space, will translate in an integer value of mode data (should be an enum) the values are "H1", "Hybrid", "Mixed"
    std::string approx;
    // this parameter indicates whether we are estimating the error of H1, hybrid o mixed approximation
    int mode = -1;           // 0 = "H1"; 1 = "Hybrid"; 2 = "Mixed";
    /// indicates the element topology that should be used
    std::string topology;           //Topology' name typed as input will translate in a value of topologyMode (should be an enum)
    /// element type to be used. It is translated from the value topology in the EvaluateEntry function
    int topologyMode = -1;
    
    /// @brief Filename initialized in the function InitializeOutputStream
    std::string topologyFileName;   //Simplified name used for naming files/directories
    /// @brief string added to the directory name determined by the variable plotfile
    std::string FEMoutput = "FEMsimulation/";
    /// @brief name determined by the problem, element type, approximation space, polynomial order
    std::string plotfile;

    /// @brief Number of adaptivity cycles that will be executed. Initialized in the main program
    int numberAdapativitySteps =-1;

    /// @brief Intended to determine the resolution of the vtk files. Never used
    int vtkResolution = -1;

    REAL perm_Q1 = 5;      /// Permeability coefficient of even quadrants (Steklov only)
    REAL perm_Q2 = 1;

    /// @brief never initialized. Mentioned in the function StockErrorsH1 in Solver.cpp
    REAL hLog = -1;
    /// @brief @brief Is set to 1/pConfig.exp in main function
    REAL h = -1000;
    /// @brief Number of error values that will be post processed
    int numErrors = 4;

    // argc is the number of arguments passed by the command line. Will be set to 5 if the problem is run from command line
    int argc = 1;

    /// @brief Number of terms for the exact solution Laplace2D data fMaxIter
    int maxIter = 15;

    /// @brief flag indicating whether to estimate the error
    bool estimateError;
    /// @brief variable controlling the debugging output
    bool debugger = true;
    /// @brief variable used to set the value of h at each adaptivity step. Its value is never updated
    int exp = 2; // Initial exponent of mesh refinement (numElem = 2*2^exp)

    PreConfig() {
        std::string FEMsimulationOutput = "FEMsimulation/";
        std::string temp = FEMsimulationOutput;
        /// remove trailing slash
        temp.pop_back();
        std::string command = "mkdir -p " + temp;
    }
};

#endif //ERRORESTIMATION_DATASTRUCTURE_H
