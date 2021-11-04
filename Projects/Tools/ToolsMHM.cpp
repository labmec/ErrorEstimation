//
//  ToolsMHM.cpp
//
//  Created by Denise De Siqueira on 01/11/19.
//

#include "ToolsMHM.h"
#include "pzvec.h"
#include "pzstack.h"
#include "pzreal.h"
#include "pzgmesh.h"
#include "tpzgeoblend.h"
#include "TPZHybridizeHDiv.h"
#include "TPZGenGrid2D.h"

#include <cmath>
#include <set>

// compute the coarse indices of the geometric mesh
void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices)
{
    coarseindices.Resize(gmesh->NElements());
    int count = 0;
    for (int64_t el=0; el<gmesh->NElements(); el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel || gel->Dimension() != gmesh->Dimension()) continue;
        if(gel->Father()) continue;
        coarseindices[count] = el;
        count++;
    }
    coarseindices.Resize(count);
}

TPZGeoMesh *CreateLMHMMesh(int nDiv, TPZVec<int64_t>& coarseIndexes) {

    int factor = static_cast<int>(lround(pow(2, nDiv) + 0.5));
    
    auto xElements = factor * 6;
    auto yElements = factor * 4;
    
    TPZManVector<int, 2> nx(2);
    nx[0] = xElements;
    nx[1] = yElements;
    
    TPZManVector<REAL> x0(3,0.), x1(3,1.);
    x1[2] = 0.;
    
    TPZGenGrid2D gen(nx, x0, x1, 1, 0);
    
    auto *gmesh = new TPZGeoMesh();
    gen.Read(gmesh);
    gen.SetBC(gmesh, 4, -1);
    gen.SetBC(gmesh, 5, -1);
    gen.SetBC(gmesh, 6, -1);
    gen.SetBC(gmesh, 7, -1);
    
    gmesh->SetDimension(2);
    gmesh->BuildConnectivity();
    
    // Assigns matIDs to create L elements
    {
        int coarseIndex = 0;
        
        // Get number of 2D elements
        int64_t nelem = xElements * yElements;
        
        coarseIndexes.Resize(nelem, -1);
        for (int64_t elem = 0; elem < nelem; elem++) {
            TPZGeoEl *gel = gmesh->ElementVec()[elem];
            if (gel->Dimension() != 2) DebugStop();
            
            auto lineInPattern = elem / nx[0] % 4;
            auto colInPattern = elem % nx[0] % 6;
            
            // IDs of elements in the neighbourhood to which a coarse index has been already assigned
            auto leftEl = elem - 1;
            auto bottomEl = elem - nx[0];
            
            if (lineInPattern == 0) {
                if (colInPattern % 2 == 0) {
                    coarseIndexes[elem] = coarseIndex;
                    coarseIndex++;
                } else {
                    coarseIndexes[elem] = coarseIndexes[leftEl];
                }
            } else if (lineInPattern == 1) {
                if (colInPattern == 0 || colInPattern == 2 || colInPattern == 5) {
                    coarseIndexes[elem] = coarseIndexes[bottomEl];
                } else if (colInPattern == 4) {
                    coarseIndexes[elem] = coarseIndexes[leftEl];
                } else {
                    coarseIndexes[elem] = coarseIndex;
                    coarseIndex++;
                }
            } else if (lineInPattern == 2) {
                if (colInPattern == 1 || colInPattern == 4) {
                    coarseIndexes[elem] = coarseIndexes[bottomEl];
                } else if (colInPattern == 2) {
                    coarseIndexes[elem] = coarseIndexes[leftEl];
                } else {
                    coarseIndexes[elem] = coarseIndex;
                    coarseIndex++;
                }
            } else if (lineInPattern == 3) {
                if (colInPattern == 0) {
                    coarseIndexes[elem] = coarseIndexes[bottomEl];
                } else if (colInPattern % 2 == 1) {
                    coarseIndexes[elem] = coarseIndexes[leftEl];
                } else {
                    coarseIndexes[elem] = coarseIndexes[leftEl] + 1;
                }
            }
        }
    }

    return gmesh;
}