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
#include "TPZMultiphysicsCompMesh.h"

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
    int ninternalref = -1;
    int adaptivityStep = -1; // Keep this variable for compatibility to maintain support for Gustavo's code. Useless for HybridH1 simulations.
    int maxPrefine = 4;
    
    int dimension = 0;
    bool prefine = false;
    bool steklovexample = false;
    bool GalvisExample = false;
    bool TensorNonConst = false;
    bool MeshNonConvex = false;
    
    STATE alpha=1;
    STATE Km = 0.;
    STATE coefG = 0.;
    
    /// the elements with error larger than 0.7*max_error will be divided
    REAL division_threshold = 0.7;
    
    /// directory where the files will be stored
    std::string dir_name = ".";
    /// name identifying the problem
    std::string problemname;
    /// set of materialids in the mesh
    std::set<int> materialids;
    /// set of boundary condition material ids
    std::set<int> bcmaterialids;

    int vtkResolution = -1;
    /// exact solution
    TPZAutoPointer<TLaplaceExample1> exact;
    
    /// set of elements to be divided after each adaptivity step
    std::list<std::set<int64_t> > fElIndexDivide;

    /// set of elements to be increment his p-order after each adaptivity step
    std::list<std::map<int64_t,int> > fElIndexPplus;
    
    ProblemConfig() = default;

    ProblemConfig(const ProblemConfig &cp) = default;

    ProblemConfig &operator=(const ProblemConfig &cp) = default;
    
    void ApplyDivision();
    
    void ApplyTwoonOneRestraint();
    
    void DivideEncircledElements();
    
    void DivideBoundaryElements();
    
    void PorderIncrement();
    
};

struct EstimatorConfig{
   TPZMultiphysicsCompMesh *fOriginal =NULL;

   /// name identifying the problem
   std::string *fproblemname;
   /// set of materialids in the mesh
   std::set<int> fmaterialids;
   /// set of boundary condition material ids
   std::set<int> fbcmaterialids;
   /// Material id of lagrange coefficents
   int fLagrangeMatId =-666;
    /// Mateiral for computing the continuous skeleton solution
    int fSkeletonMatId = -666;
   /// Polynomial orders of original problem
   int fk = -666;
   int fn = -666;
   /// Number of divisions of uniformly refined meshes
   int fnDivisions = -666;

   int fAdaptivityStep = -1;

    int fvtkResolution = -1;
    
    REAL fdivision_threshold;
   /// exact solution
   TPZAutoPointer<TLaplaceExample1> fExact;

   EstimatorConfig(TPZMultiphysicsCompMesh *multimesh,ProblemConfig pConfig,int lagrangeMatId){
       fOriginal = multimesh;
        fproblemname = &(pConfig.problemname);
        fmaterialids = pConfig.materialids;
        fbcmaterialids=pConfig.bcmaterialids;
        fnDivisions = pConfig.ndivisions;
        fExact = pConfig.exact;
        fAdaptivityStep = pConfig.adaptivityStep;
        fvtkResolution = pConfig.vtkResolution;
        fk = pConfig.k;
        fn = pConfig.n;
        fdivision_threshold = pConfig.division_threshold;

        fLagrangeMatId = lagrangeMatId;
   }
};

//EstimatorConfig::EstimatorConfig(TPZMultiphysicsCompMesh *multimesh,ProblemConfig pConfig,int lagrangeMatId){
//    
//}

#endif /* ProblemConfig_h */
