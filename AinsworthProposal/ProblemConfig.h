//
// Created by Gustavo on 21/05/19.
//
// This class is used to manage problems to be used with error estimators.
//

#ifndef PROBLEMCONFIG_H
#define PROBLEMCONFIG_H

#include <set>
#include <tuple>
#include "TPZAnalyticSolution.h"

using namespace std;
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
    string dirName;
    string problemName;
    string meshFileName;

    // Material handling
    set<tuple<string, int, int>> DomainMats;
    set<tuple<string, int, int>> BCMats;

    // Exact solution to validate the estimator
    TLaplaceExample1 exactSolution;

    // Default constructors and destructor
    ProblemConfig() = default;
    ProblemConfig(const ProblemConfig &cp) = default;
    ProblemConfig &operator=(const ProblemConfig &cp) = default;

    // Adds material to problem config
    void InsertDomainMat(string matName, int matID, int matDim) {
        DomainMats.insert(std::make_tuple(matName, matID, matDim));
    }

    void InsertBCMat(string bcName, int bcID, int bcDim) {
        BCMats.insert(std::make_tuple(bcName, bcID, bcDim));
    }

    // Reads Gmsh file
    TPZGeoMesh *CreateGeoMesh() {

        TPZGmshReader reader;
        for (const auto &DomainMat : DomainMats) {
            reader.GetDimNamePhysical()[get<2>(DomainMat)][get<0>(DomainMat)] = get<1>(DomainMat);
        }

        for (const auto &BCMat : BCMats) {
            reader.GetDimNamePhysical()[get<2>(BCMat)][get<0>(BCMat)] = get<1>(BCMat);
        }

        reader.SetFormatVersion("4.1");

#ifdef MACOSX
        meshFileName = "../" + meshFileName;
#endif

        gmesh = reader.GeometricGmshMesh(meshFileName);
        gmesh->SetDimension(reader.Dimension());

#ifdef PZDEBUG
        reader.PrintPartitionSummary(std::cout);
#endif
        return gmesh;
    }

};

#endif // PROBLEMCONFIG_H
