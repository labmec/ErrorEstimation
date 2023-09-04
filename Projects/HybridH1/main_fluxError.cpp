//
// Created by victor on 16/03/2021.
//
//// This target focuses on the p-convergence of the flux contribution (NF) to the estimated error;
//// As such, the laplacian equation is evaluated over a single quadrilateral,
//// with homogeneous Neumann BC over 3 sides and a linear Neumann condition over the last side.
//// For a given value of k, the flux index will be evaluated for different enrichment values.

#include "InputTreatment.h"
#include "Solver.h"
#include "Output.h"
#include "tpzgeoelrefpattern.h"
#include "DataStructure.h"
#include "Tools.h"
#include "TPZCreateMultiphysicsSpace.h"
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
    {
        pConfig.h = 1.;
        ProblemConfig config;
        FluxErrorConfigure(config,  pConfig);

        TPZMultiphysicsCompMesh *multiCmesh = new TPZMultiphysicsCompMesh(config.gmesh);
        int interfaceMatID = -10;
        int fluxMatID = -10;

        FluxErrorCreateCompMesh(multiCmesh, interfaceMatID, fluxMatID,pConfig, config);

    }
    return 0.;
}