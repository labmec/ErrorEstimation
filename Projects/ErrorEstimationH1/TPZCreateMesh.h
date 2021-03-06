//
//  TPZCreateHDivMesh.hpp
//  PZ
//
//  Created by Philippe Devloo on 6/21/16.
//
//

#ifndef TPZCreateHDivMesh_hpp
#define TPZCreateHDivMesh_hpp

#include <stdio.h>
#include "pzcmesh.h"

TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim);
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim, bool disconnected, bool referred);
TPZCompMesh *CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> &meshvec);

void CreateAllMeshes(TPZGeoMesh *gmesh, TPZVec<TPZCompMesh *> &meshvec, int porder, int dim, int hdivplusplus);

/// Reconstruct the multi physics mesh from the fluxmesh
void ReconstructHDivMesh(TPZVec<TPZCompMesh *> &meshvec, int hdivplusplus);

/// adjust the polynomial orders of the hdiv elements such that the internal order is higher than the sideorders
void AdjustFluxPolynomialOrders(TPZCompMesh *fluxmesh, int hdivplusplus);

/// set the pressure order acording to the order of internal connect of the elements of the fluxmesh
void SetPressureOrders(TPZCompMesh *fluxmesh, TPZCompMesh *pressuremesh);

/// uncondense the elements unwrap the elements
void UnwrapMesh(TPZCompMesh *cmesh);

void UnitPressure(TPZCompMesh *PressMesh);


#endif /* TPZCreateHDivMesh_hpp */
