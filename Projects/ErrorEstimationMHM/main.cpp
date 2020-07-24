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

int main() {
    InitializePZLOG();

    ProblemConfig config;
    config.dimension = 2;
    config.prefine = false;
    config.TensorNonConst = false;
    config.MeshNonConvex = true;

    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;

    config.problemname = "NonConvexMesh";

    config.dir_name = "NonConvexMHM";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    for (int orderp = 1; orderp < 2; orderp++) {
        config.porder = orderp;
        for (int hdivmais = 1; hdivmais < 2; hdivmais++) {
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
