//
// Created by Gustavo Batistela on 11/8/21.
//

#ifndef TPZMULTISCALEGRIDGEN2D_H
#define TPZMULTISCALEGRIDGEN2D_H

#include <pzgmesh.h>

struct RefTree {
    int fActualSize{0};
    RefTree *fChildLeft = nullptr;
    RefTree *fChildRight = nullptr;

    RefTree() = default;

    explicit RefTree(const int actual_number) {
        fActualSize = actual_number;
        FillRefTree();
    }

    void FillRefTree() {
        std::div_t division = std::div(fActualSize, 2);
        if (fActualSize != 1) {
            fChildLeft = new RefTree(division.quot);
            fChildRight = new RefTree(division.quot + division.rem);
        }
    }
};

class TPZMultiscaleGridGen2D {

public:
    // Constructors
    TPZMultiscaleGridGen2D() = delete;

    TPZMultiscaleGridGen2D(const TPZVec<REAL> &minX, const TPZVec<REAL> &maxX, const TPZVec<int> &NDivFineGrid,
                           int NElemCoarseGrid);

private:
    // Member variables
    const TPZManVector<int, 2> fNDivFineGrid{0};
    const int fNElemCoarseGrid{0};
    const TPZManVector<REAL, 3> fMinX{0};
    const TPZManVector<REAL, 3> fMaxX{0};

    RefTree *fRefTreeDesiredSize = nullptr;
    RefTree *fRefTreeRemainderX = nullptr;
    RefTree *fRefTreeRemainderY = nullptr;
    std::map<std::pair<int, int>, TPZRefPattern> fRefPatterns;
    std::map<int64_t, RefTree*> fSkelIdToRefTree;

    TPZGeoMesh *fGeoMesh = nullptr;

    // Private member functions
    static TPZRefPattern CreateNonUniformLineRefPattern(int a, int b);

    void GenerateRefPatterns();

    void CreateFineGridMesh();

    void CreateSkeletonElements();

    void RefineSkeletonElements();
};

#endif // TPZMULTISCALEGRIDGEN2D_H
