//
//  ToolsMHM.hpp
//
//  Created by Denise De Siqueira on 01/11/19.
//

#ifndef TOOLSMHM_H
#define TOOLSMHM_H

#include "pzvec.h"
#include "pzgmesh.h"

// Compute the geometric mesh coarse indices
void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices);

TPZGeoMesh *CreateLMHMMesh(int nDiv, TPZVec<int64_t> &coarseIndexes);

#endif // TOOLSMHM_H
