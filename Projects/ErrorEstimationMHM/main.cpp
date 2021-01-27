#include "pzlog.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"

#include <string>
#include <cmath>
#include <set>

#include "Tools.h"
#include "meshgen.h"
#include "ToolsMHM.h"
#include "TPZMHMHDivErrorEstimator.h"

using namespace std;


void RunNonConvexProblem();
void RunHighGradientProblem();

int main() {
    InitializePZLOG();

    RunHighGradientProblem();
    
}

void RunNonConvexProblem(){
    ProblemConfig config;
    config.dimension = 2;
    config.prefine = false;
    config.TensorNonConst = false;
    config.MeshNonConvex = true;

    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinMark;

    config.problemname = "LShape";

    config.dir_name = "LShapeMHMSingular";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    for (int orderp = 1; orderp < 2; orderp++) {
        config.porder = orderp;
        for (int hdivmais = 3; hdivmais < 4; hdivmais++) {
            config.hdivmais = hdivmais;

            for (int ndiv = 0; ndiv < 1; ndiv++) {
                config.ndivisions = ndiv;
                config.materialids.insert(1);
                config.bcmaterialids.insert(-1);

                MHMTest(config);
            }
        }
    }
    
}

void RunHighGradientProblem(){
    ProblemConfig config;
    config.dimension = 2;
    config.prefine = false;
    config.TensorNonConst = false;
    config.MeshNonConvex = false;

    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::EBoundaryLayer;

    config.problemname = "BoundaryLeyer";

    config.dir_name = "HighGradient";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    for (int orderp = 1; orderp < 2; orderp++) {
        config.porder = orderp;
        for (int hdivmais = 3; hdivmais < 4; hdivmais++) {
            config.hdivmais = hdivmais;

            for (int ndiv = 3; ndiv < 4; ndiv++) {
                config.ndivisions = ndiv;
                config.materialids.insert(1);
                config.bcmaterialids.insert(-1);

                MHMTestDenise(config);
            }
        }
    }
}

