//
//  ProblemConfig.h
//  ErrorEstimation
//
//  Created by Philippe Devloo on 09/05/18.
//

#ifndef ProblemConfig_h
#define ProblemConfig_h

#include <set>
#include "TPZAnalyticSolution.h"

/// class to guide the error estimator
struct ProblemConfig
{
    
    /// geometric mesh on which the computational meshes are based
    TPZGeoMesh *gmesh = 0;
    /// polynomial order of the original mesh
    int porder = 2;
    /// increment in internal order of flux and pressure
    int hdivmais = 1;

    /// Instead of specifying "porder" and "hdivmais", one may specify "k" and "n"
    /// Enrichment order
    int n = 1;
    /// Flux order for HDiv configuration or Lagrange coefficient order for Primal Hybrid.
    int k = 1;

    /// option to compute the error based on continuous pressures or not
    bool makepressurecontinuous = true;
    
    /// number of uniform refinements applied to the mesh
    int ndivisions = -1;
    int adaptivityStep = -1;
    int dimension = 0;
    bool prefine = false;
    bool steklovexample = false;
    bool GalvisExample = false;
    bool TensorNonConst = false;
    bool MeshNonConvex = false;
    
    STATE alpha=1;
    STATE Km = 0.;
    STATE coefG = 0.;
    
    /// directory where the files will be stored
    std::string dir_name = ".";
    /// name identifying the problem
    std::string problemname;
    /// set of materialids in the mesh
    std::set<int> materialids;
    /// set of boundary condition material ids
    std::set<int> bcmaterialids;
    /// exact solution
    TPZAutoPointer<TLaplaceExample1> exact;

    ProblemConfig() = default;

    ProblemConfig(const ProblemConfig &cp) = default;

    ProblemConfig &operator=(const ProblemConfig &cp) = default;
};

#endif /* ProblemConfig_h */
