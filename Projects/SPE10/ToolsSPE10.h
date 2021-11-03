//
// Created by Gustavo Batistela on 10/26/21.
//

#ifndef TOOLSSPE10_H
#define TOOLSSPE10_H

#include "pzgmesh.h"

namespace SPE10 {
[[maybe_unused]] [[maybe_unused]] TPZGeoMesh *CreateFineGridGeoMesh();
[[maybe_unused]] TPZGeoMesh *CreateMHMGeoMesh();
[[maybe_unused]] TPZGeoMesh *CreateRefinementGeoMesh(int nx, int ny);

}

#endif // TOOLSSPE10_H
