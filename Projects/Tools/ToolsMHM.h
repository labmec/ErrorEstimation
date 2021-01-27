//
//  ToolsMHM.hpp
//  AdaptivityTest
//
//  Created by Denise De Siqueira on 01/11/19.
//

#ifndef ToolsMHM_hpp
#define ToolsMHM_hpp

#include <stdio.h>

#endif /* ToolsMHM_hpp */


#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <fstream>
#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"
#include "pzreal.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZParSkylineStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "TPZGeoLinear.h"

#include "pzpoisson3d.h"
//#include "pzhybridpoisson.h"
#include "pzpoisson3dreferred.h"
#include "mixedpoisson.h"
#include "pzelasmat.h"
#include "pzelasthybrid.h"
#include "pzmat1dlin.h"
#include "TPZVecL2.h"
#include "TPZMatLaplacianHybrid.h"
#include "TPZLagrangeMultiplier.h"

#include "pzbuildmultiphysicsmesh.h"
#include "pzelementgroup.h"
#include "TPZCompMeshTools.h"
#include "pzcondensedcompel.h"
#include "pzfunction.h"
#include "pzgraphmesh.h"
#include "pzfmatrix.h"

#include "pzlog.h"

#include "TPZVTKGeoMesh.h"
#include "pzvisualmatrix.h"

//#include "pzgengrid.h"

#include "TPZGenGrid2D.h"
#include "TPZExtendGridDimension.h"
#include "pzcheckgeom.h"

#include "TPZMHMeshControl.h"
#include "TPZMHMixedMeshControl.h"
#include "TPZMHMixedMeshChannelControl.h"

#include "TPZMHMixedHybridMeshControl.h"
#include "TPZHybridizeHDiv.h"
#include "meshgen.h"
//#include "ConfigCasesMaze.h"
#include <iostream>
#include <string>
//#include <opencv2/opencv.hpp>
#include <math.h>
#include <set>
#include "pzsolve.h"

#include "TPZPersistenceManager.h"

#include "TPZMHMHDivErrorEstimator.h"
#include "Tools.h"
#include "TPZGenGrid2D.h"

using namespace std;

#define LMESH



// Creating the computational H1 mesh
TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int p_order, ProblemConfig &Conf);

TPZCompMesh *CMesh_H1(struct ProblemConfig &problem);

// Creating the computational flux mesh
TPZCompMesh *CMeshFlux(TPZGeoMesh * gmesh,int pOrder);

// Creating the computational pressure mesh
TPZCompMesh *CMeshPressure(TPZGeoMesh * gmesh, int pOrder,ProblemConfig &Conf);

// Creating the computational multphysics mesh
TPZCompMesh *CMeshMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec,ProblemConfig &Conf);

// Read a mesh from a png file. The size of the domain will be npix_x by npix_y (read from the image) . Return l=nx; h=ny.
// on entry image_size_x and image_size_y are equal to the number of MHM elements in x and y
TPZGeoMesh *GeoMeshFromPng(string name, int &image_size_x, int &image_size_y);

// Create a geometric mesh with the given parameters, nx and ny are the coarse elements number. The total number of elements are defined by the image read.
TPZGeoMesh *GenerateGeoMesh(string name, int nx, int ny);

// Create a geoElSide map to be open because it has flux.
// cmesh is the flux mesh obtained from a Hdiv previus simulation.
std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>> IdentifyChanel (TPZCompMesh *cmesh);


// Compute the geometric mesh coarse indices
void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices);

// Insert the necessary objects material in the computational mesh
void InsertMaterialObjects(TPZMHMixedMeshControl &control,const ProblemConfig &config);

// Solve the H1 problem with "Conf" configuration
// Conf contains the maze information and the problem boundary conditions
int H1Test(ProblemConfig &Conf);

// Solve the mixed problem with "Conf" configuration
// Conf contains the maze information and the problem boundary conditions
TPZCompMesh* MixedTest(ProblemConfig &Conf);

// Solve the maze using MHM. By default (2x2 coarse elements)
// Conf contains the maze information and the problem boundary conditions
int MHMTest(ProblemConfig &Conf);
int MHMTestDenise(ProblemConfig &Conf);

TPZGeoMesh *CreateCircleGeoMesh();
TPZGeoMesh *CreateLMHMMesh(int nDiv, TPZVec<int64_t>& coarseIndexes);

/// Solve the problem composed of a multiphysics mesh composed of compmeshes - applies to MHM and MHM-H(div)
void SolveProblem(TPZAutoPointer<TPZCompMesh> cmesh, const TPZVec<TPZAutoPointer<TPZCompMesh> >& compmeshes, TPZAnalyticSolution *analytic, const std::string& prefix, TRunConfig config);
//adpative refinement
void MHMAdaptivity(TPZMHMixedMeshControl *mhm, TPZGeoMesh* gmeshToRefine, ProblemConfig& config);
