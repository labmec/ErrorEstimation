//
// Created by Gustavo Batistela on 11/8/21.
//

#ifndef TPZMULTISCALEGRIDGEN2D_H
#define TPZMULTISCALEGRIDGEN2D_H

#include <pzgmesh.h>

struct RefTree {
    int fActualNumber{0};
    RefTree *fChildLeft = nullptr;
    RefTree *fChildRight = nullptr;

    RefTree() = default;

    explicit RefTree(const int actual_number) {
        fActualNumber = actual_number;
        FillRefTree();
    }

    void FillRefTree() {
        std::div_t division = std::div(fActualNumber, 2);
        if (fActualNumber != 1) {
            fChildLeft = new RefTree(division.quot);
            fChildRight = new RefTree(division.quot + division.rem);
        }
    }
};

class TPZMultiscaleGridGen2D {

private:
    // Member variables
    const TPZManVector<int, 2> fNDivFineGrid{0};
    const int fNElemCoarseGrid{0};
    const TPZManVector<REAL, 3> fMinX{0};
    const TPZManVector<REAL, 3> fMaxX{0};

    RefTree fRefTree;
    std::map<int, std::map<int, TPZRefPattern>> fRefPatterns;

    TPZGeoMesh *fGeoMesh = nullptr;

    // Private member functions
    static TPZRefPattern CreateNonUniformLineRefPattern(int a, int b);

};

#endif // TPZMULTISCALEGRIDGEN2D_H
