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
    int porder = 1;
    int hdivmais = 1;
    bool makepressurecontinuous = 0;
    
    int ndivisions=1;
    bool prefine=false;
    STATE alpha=1;
    std::string dir_name;
    
    
    std::string problemname;
    std::set<int> materialids;
    std::set<int> bcmaterialids;
    TLaplaceExample1 exact;
    
    ProblemConfig(){};
    ProblemConfig(const ProblemConfig &cp) : gmesh(cp.gmesh), porder(cp.porder), hdivmais(cp.hdivmais),
    makepressurecontinuous(cp.makepressurecontinuous),
    problemname(cp.problemname),
    materialids(cp.materialids), bcmaterialids(cp.bcmaterialids),exact(cp.exact),ndivisions(cp.ndivisions),prefine(cp.prefine),alpha(cp.alpha),dir_name(cp.dir_name)
    {
    }
    
    ProblemConfig &operator=(const ProblemConfig &cp)
    {
        gmesh = cp.gmesh;
        porder = cp.porder;
        hdivmais = cp.hdivmais;
        makepressurecontinuous = cp.makepressurecontinuous;
        problemname = cp.problemname;
        materialids = cp.materialids;
        bcmaterialids = cp.bcmaterialids;
        exact = cp.exact;
        
        ndivisions = cp.ndivisions;
        prefine = cp.prefine;
        alpha = cp.alpha;
        dir_name = cp.dir_name;
        
        return *this;
    }
};

#endif /* ProblemConfig_h */
