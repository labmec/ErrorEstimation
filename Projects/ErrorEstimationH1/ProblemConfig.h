//
// Created by Gustavo on 15/07/18.
//

#include <iostream>
#include "pzgmesh.h"

#ifndef PROBLEMCONFIG_H
#define PROBLEMCONFIG_H

// Simulation Case
struct SimulationCase {

    int nthreads = 0;
    int numinitialrefine = 0;
    int porder = 1;
    std::string dir_name;
    TPZGeoMesh * gmesh = nullptr;
    
    /// name identifying the problem
    std::string problemname;
    
    std::set<int> materialids;
    std::set<int> bcmaterialids;
    TLaplaceExample1 exact;

    SimulationCase() {
        dir_name = "dump";
    }

    SimulationCase(const SimulationCase &cp)
            : nthreads(cp.nthreads),
              numinitialrefine(cp.numinitialrefine),
              porder(cp.porder),
              dir_name(cp.dir_name),
              gmesh(cp.gmesh),
    materialids(cp.materialids),
    bcmaterialids (cp.bcmaterialids),problemname(cp.problemname),
    exact(cp.exact)
    {
    }
};
#endif // PROBLEMCONFIG_H
