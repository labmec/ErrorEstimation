
#ifndef MAZE_COMMON_H
#define MAZE_COMMON_H


#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzlog.h"
#include "pzreal.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoelside.h"
#include "tpzgeoelrefpattern.h"
#include "tpzautopointer.h"
#include "TPZLinearAnalysis.h"

#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZStructMatrixT.h"

#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZVTKGeoMesh.h"
#include "TPZGenGrid2D.h"
#include "TPZMHMixedMeshControl.h"
#include "TPZMHMixedMeshChannelControl.h"
#include "TPZHybridizeHDiv.h"
#include "ConfigCasesMaze.h"
#include <ToolsMHM.h>
#include <cmath>
#include <opencv2/opencv.hpp>
#include <set>
#include <string>

#include "TPZPersistenceManager.h"

using namespace std;
using namespace cv;

// Creating the computational flux mesh
TPZCompMesh *CMeshFlux(TPZGeoMesh * gmesh,int pOrder);

// Creating the computational pressure mesh
TPZCompMesh *CMeshPressure(TPZGeoMesh * gmesh, int pOrder,ConfigCasesMaze Conf);

// Creating the computational multphysics mesh
TPZCompMesh *CMeshMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec,ConfigCasesMaze &Conf);

// Read a mesh from a png file. The size of the domain will be npix_x by npix_y (read from the image) . Return l=nx; h=ny.
TPZGeoMesh *GeoMeshFromPng(string name, double &l, double &h);

// Create a geometric mesh with the given parameters, nx and ny are the coarse elements number. The total number of elements are defined by the image read.
TPZGeoMesh *GenerateGeoMesh(string name, int nx, int ny);

// Compute the geometric mesh coarse indices
void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices);

// Insert the necessary objects material in the computational mesh
void InsertMaterialObjects(TPZMHMixedMeshControl &control);

// Solve the mixed problem with "Conf" configuration
// Conf contains the maze information and the problem boundary conditions
TPZCompMesh* MixedTest(ConfigCasesMaze &Conf, int nx, int ny);

// Solve the maze using MHM. By default (2x2 coarse elements)
// Conf contains the maze information and the problem boundary conditions
int MHMTest(ConfigCasesMaze &Conf);

void EstimateError(TPZMHMHDivErrorEstimator &errorEstimator, ProblemConfig &config);

void LocateElementsToAdapt(TPZMHMHDivErrorEstimator &errorEstimator, ProblemConfig &config);

#endif
