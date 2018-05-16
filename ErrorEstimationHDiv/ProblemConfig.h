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

struct ProblemConfig
{
    TPZGeoMesh *gmesh = 0;
    int porder = 1;
    bool hdivmais = false;
    std::set<int> materialids;
    std::set<int> bcmaterialids;
    TLaplaceExample1 exact;
    
    ProblemConfig(){};
    ProblemConfig(const ProblemConfig &cp) : gmesh(cp.gmesh), porder(cp.porder), hdivmais(cp.hdivmais),
    materialids(cp.materialids), bcmaterialids(cp.bcmaterialids),exact(cp.exact)
    {
    }
    
    ProblemConfig &operator=(const ProblemConfig &cp)
    {
        gmesh = cp.gmesh;
        porder = cp.porder;
        hdivmais = cp.hdivmais;
        materialids = cp.materialids;
        bcmaterialids = cp.bcmaterialids;
        exact = cp.exact;
        return *this;
    }
};

#endif /* ProblemConfig_h */
