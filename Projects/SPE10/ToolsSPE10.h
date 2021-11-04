//
// Created by Gustavo Batistela on 10/26/21.
//

#ifndef TOOLSSPE10_H
#define TOOLSSPE10_H

#include "pzgmesh.h"

namespace SPE10 {
[[maybe_unused]] TPZGeoMesh *CreateFineGridGeoMesh();
[[maybe_unused]] TPZGeoMesh *CreateMHMGeoMesh();
[[maybe_unused]] TPZGeoMesh *CreateRefinementGeoMesh(int nx, int ny);
[[maybe_unused]] TPZGeoMesh *CreateLineRefinementGeoMesh(int nx);
}

#endif // TOOLSSPE10_H
