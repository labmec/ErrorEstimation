//
// Created by Gustavo Batistela on 11/8/21.
//

#ifndef TPZMULTISCALEGRIDGEN2D_H
#define TPZMULTISCALEGRIDGEN2D_H

#include <pzgmesh.h>

class TPZMultiscaleGridGen2D {

private:
    // Member variables
    const TPZManVector<int, 2> fNDivFineGrid{0};
    const int fNElemCoarseGrid{0};
    const TPZManVector<REAL, 3> fMinX{0};
    const TPZManVector<REAL, 3> fMaxX{0};

    std::map<int, std::map<int, TPZRefPattern>> fRefPatterns;

    TPZGeoMesh *fGeoMesh = nullptr;

};

#endif // TPZMULTISCALEGRIDGEN2D_H
