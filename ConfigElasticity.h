//
//  ConfigElasticity.hpp
//  ALL_BUILD
//
//  Created by Denise De Siqueira on 04/09/23.
//

#ifndef ConfigElasticity_hpp
#define ConfigElasticity_hpp

#include <stdio.h>

struct ConfigElasticity{
    
    
    const int xdiv = 1; //Number of elements in each direction
    const int pOrder = 1; // Polynomial degree
    // Family of HDiv approximation spaces.
    // The possible choices are HDivFamily::EHDivStandard, HDivFamily::EHDivConstant and HDivFamily::EHDivKernel
    HDivFamily hdivfam = HDivFamily::EHDivStandard;
    
    
    
    
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
    /// integral order of the forcing function
    int integrationorder = -1;
    /// determines the the approximation space, will translate in an integer value of mode data (should be an enum)
    std::string approx;
    /// indicates the element topology that should be used
    std::string topology;           //Topology' name typed as input will translate in a value of topologyMode (should be an enum)
    
    std::string topologyFileName;   //Simplified name used for naming files/directories
    std::string FEMoutput = "ElasticityFEMSilumation/";

    int numberAdapativitySteps =-1;

    int vtkResolution = -1;


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

    ConfigElasticity() {
        std::string FEMsimulationOutput = "ElasticityFEMSilumation/";
        std::string temp = FEMsimulationOutput;
        temp.pop_back();
        std::string command = "mkdir -p " + temp;
    }
};




#endif /* ConfigElasticity_hpp */
