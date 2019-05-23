//
// Created by Gustavo on 21/05/19.
//
// This class configures problems to be used by the error estimator.
//

#ifndef PROBLEMCONFIG_H
#define PROBLEMCONFIG_H

#include <set>
#include "TPZAnalyticSolution.h"

struct ProblemConfig {

    // Geometric mesh on which the computational meshes are based
    TPZGeoMesh *gmesh = nullptr;

    // Polynomial order of the solution mesh
    int pOrder = 1;
    // Enrichment order
    int hdivplus = 1;
    bool makePressureContinuous = false;

    // Refinement-related variables
    int nDivisions = 1;
    bool pRefine = false;

    // Permeability tensor alpha parameter
    STATE alpha = 1;

    // Problem and directory names
    std::string dirName;
    std::string problemName;

    // Material handling
    std::set<int> matIDs;
    std::set<int> bcMatIDs;

    // Exact solution to validate the estimator
    TLaplaceExample1 exactSolution;

    // Default constructors and destructor
    ProblemConfig() = default;
    ProblemConfig(const ProblemConfig &cp) = default;
    ProblemConfig &operator=(const ProblemConfig &cp) = default;
};

#endif // PROBLEMCONFIG_H
