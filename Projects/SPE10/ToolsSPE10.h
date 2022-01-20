//
// Created by Gustavo Batistela on 10/26/21.
//

#ifndef TOOLSSPE10_H
#define TOOLSSPE10_H

#include "pzgmesh.h"
#include <TPZMHMHDivErrorEstimator.h>
#include <TPZMHMixedMeshControl.h>

namespace SPE10 {

constexpr int nx = 220;
constexpr int ny = 60;
constexpr int layer = 36;
constexpr int n_cells = nx * ny;
static TPZManVector<REAL, n_cells> *perm_vec;

[[maybe_unused]] TPZGeoMesh *CreateFineGridGeoMesh();
[[maybe_unused]] TPZGeoMesh *CreateMHMGeoMesh();
[[maybe_unused]] TPZGeoMesh *CreateRefinementGeoMesh(int nx, int ny);
[[maybe_unused]] TPZGeoMesh *CreateLineRefinementGeoMesh(int nx);

STATE PermeabilityFunction(const TPZVec<REAL> &x);

void ReadSPE10CellPermeabilities();

void InsertMaterials(TPZCompMesh *cmesh);

void EstimateError(TPZMHMixedMeshControl &mhm, TPZMHMHDivErrorEstimator &estimator);
void EstimateError(TPZMHMixedMeshControl &mhm, TPZMHMHDivErrorEstimator &estimator, TPZMultiphysicsCompMesh * ref_sol);

void SolveMHMProblem(TPZMHMixedMeshControl &mhm, int adaptivity_step);

TPZGeoMesh *CreateSPE10CoarseGeoMesh();

void CreateSPE10MHMCompMesh(TPZMHMixedMeshControl &mhm, const std::vector<int64_t> &skelsToDivide, int nInternalRef);

} // namespace SPE10

#endif // TOOLSSPE10_H
