//
// Created by Gustavo Batistela on 11/8/21.
//

#include "TPZMultiscaleGridGen2D.h"
#include <tpzgeoelrefpattern.h>

TPZMultiscaleGridGen2D::TPZMultiscaleGridGen2D(const TPZVec<REAL> &minX, const TPZVec<REAL> &maxX,
                                               const TPZVec<int> &NDivFineGrid, int NElemCoarseGrid)
    : fMinX(minX), fMaxX(maxX), fNDivFineGrid(NDivFineGrid), fNElemCoarseGrid(NElemCoarseGrid) {

    if (NDivFineGrid[0] < NElemCoarseGrid || NDivFineGrid[1] < NElemCoarseGrid) DebugStop();

    fRefTreeDesiredSize = new RefTree(NElemCoarseGrid);

    std::div_t div_x = std::div(NDivFineGrid[0], NElemCoarseGrid);
    std::div_t div_y = std::div(NDivFineGrid[1], NElemCoarseGrid);
    if (div_x.rem != 0) fRefTreeRemainderX = new RefTree(div_x.rem);
    if (div_y.rem != 0) fRefTreeRemainderY = new RefTree(div_y.rem);

    GenerateRefPatterns();

}

TPZRefPattern TPZMultiscaleGridGen2D::CreateNonUniformLineRefPattern(const int a, const int b) {

    auto *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(1);

    TPZManVector<REAL, 3> coord(3, 0.);

    gmesh->NodeVec().Resize(3);
    gmesh->NodeVec()[0].Initialize(coord, *gmesh);

    coord[0] = a;
    gmesh->NodeVec()[1].Initialize(coord, *gmesh);
    coord[0] = a + b;
    gmesh->NodeVec()[2].Initialize(coord, *gmesh);

    constexpr int mat_id = 1;

    // Inserts line elements
    TPZManVector<int64_t, 4> nodes_id_vec(2, 0);
    nodes_id_vec[1] = 2;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodes_id_vec, mat_id, *gmesh);

    nodes_id_vec[0] = 0;
    nodes_id_vec[1] = 1;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodes_id_vec, mat_id, *gmesh);
    nodes_id_vec[0] = 1;
    nodes_id_vec[1] = 2;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodes_id_vec, mat_id, *gmesh);

    gmesh->BuildConnectivity();
    return TPZRefPattern(*gmesh);
}

void TPZMultiscaleGridGen2D::GenerateRefPatterns() {
    std::set<std::pair<int, int>> ref_levels;

    std::function<void(RefTree*)> VisitNode = [&](RefTree *node) -> void {
        if (node->fActualSize != 1) {
            const auto left_size = node->fChildLeft->fActualSize;
            const auto right_size = node->fChildRight->fActualSize;
            if (left_size != right_size) {
                ref_levels.insert({node->fChildLeft->fActualSize, node->fChildRight->fActualSize});
            } else {
                ref_levels.insert({1, 1});
            }
            VisitNode(node->fChildLeft);
            VisitNode(node->fChildRight);
        }
    };

    VisitNode(fRefTreeDesiredSize);
    if (fRefTreeRemainderX) VisitNode(fRefTreeRemainderX);
    if (fRefTreeRemainderY) VisitNode(fRefTreeRemainderY);

    for (auto it : ref_levels) {
        const auto a = it.first;
        const auto b = it.second;

        TPZRefPattern refPattern = CreateNonUniformLineRefPattern(a, b);
        fRefPatterns.insert({{it.first, it.second}, refPattern});
    }
}