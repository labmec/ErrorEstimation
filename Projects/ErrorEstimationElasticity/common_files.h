#pragma once

#include <fstream>

#include "pzlog.h"
#include "pzgmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZMHMeshControl.h"
#include <iostream>
#include <string>

#include <cmath>
#include <set>

enum EConfigElastic {
    EThiago = 0, EAxiSymmetric = 1, EThiagoPlus = 2, EAxiSymmetricPlus = 3, EThiagoPlusPlus = 4
};

enum EConfigDarcy {
    EPhil = 0
};

static std::string ConfigRootnameElastic[5] = {
    "Mixed",
    "Mixed_AxiSymmetric",
    "MixedPlus",
    "Mixed_AxiSymmetricPlus",
    "MixedPlusPlus"
};

static std::string ConfigRootnameDarcy[1] = {
    "Primal"
};

// local enum for mesh types @ToDo these names might lead to confusion. We should consider changing.
enum EElementType {
    ETriangular = 0, ESquare = 1, ETrapezoidal = 2
};

const int dim = 2; // Dimension of the problem
const int matID = 1; // Material of the volumetric element
const int matLagrange = -10; // Material of the Lagrange multipliers
const int matBCbott = -1, matBCtop = -2, matBCleft = -3, matBCright = -4; // Materials of the boundary conditions
const int dirichlet = 0, neumann = 1, mixed = 2, dirichletvar = 4, pointtype = 5; // Boundary conditions of the problem ->default: Dirichlet on left and right

/**
 * @brief Funcao para criar a malha geometrica do problema a ser simulado
 * @note A malha sera unidim5ensional formada por nel elementos de tamanho elsize
 * @param uNDiv number of divisions ortogonal to the plates performed on the domain
 * @param vNDiv number of divisions parallel to the plates performed on the domain
 * @param nel numero de elementos
 * @param elsize tamanho dos elementos
 */
TPZGeoMesh *CreateGMesh(int nelx, int nely, double hx, double hy, double x0, double y0, EElementType meshType);

/**
* @brief Funcao para criar a malha geometrica do problema a ser simulado
* @note A malha sera tridimensional formada por nel elementos de tamanho elsize
* @param nelx, nely number of elements in the x and y direction (the number of elements in the z direction = nelx
* @param hx, hy size of the domain in x and y (size in z = hx)
* @param x0, y0 bottom left point coordinate (bottom left in z = 0)
* @param meshType = triangle, quadrilateral or trapeze
*/
TPZGeoMesh *CreateGMesh3D(int nelx, int nely, double hx, double hy, double x0, double y0, EElementType meshType);

/// Build the system of equations and post process the solution
void SolveProblem(TPZMultiphysicsCompMesh &cmesh, std::stringstream &rootname);

/// Compute the error norms
void ComputeError(TPZCompMesh &cmesh, std::stringstream &rootname, int href, int pref, TPZMHMeshControl &mhm);
