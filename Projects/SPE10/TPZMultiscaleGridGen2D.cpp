//
// Created by Gustavo Batistela on 11/8/21.
//

#include "TPZMultiscaleGridGen2D.h"
#include <TPZGenGrid2D.h>
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
    CreateFineGridMesh();
    CreateSkeletonElements();
    RefineSkeletonElements();

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

void TPZMultiscaleGridGen2D::CreateFineGridMesh() {

    const TPZManVector<int, 4> bcIDs = {-2, -1, -2, -2};

    TPZGenGrid2D gen(fNDivFineGrid, fMinX, fMaxX, 1, 0);
    gen.SetRefpatternElements(true);

    fGeoMesh = new TPZGeoMesh;
    gen.Read(fGeoMesh);

    gen.SetBC(fGeoMesh, 4, bcIDs[0]);
    gen.SetBC(fGeoMesh, 5, bcIDs[1]);
    gen.SetBC(fGeoMesh, 6, bcIDs[2]);
    gen.SetBC(fGeoMesh, 7, bcIDs[3]);

    fGeoMesh->SetDimension(2);
}

void TPZMultiscaleGridGen2D::CreateSkeletonElements() {

    const auto node_vec = fGeoMesh->NodeVec();

    const std::div_t div_x = std::div(fNDivFineGrid[0], fNElemCoarseGrid);
    const std::div_t div_y = std::div(fNDivFineGrid[1], fNElemCoarseGrid);

    const int n_full_size_x = div_x.quot;
    const int n_full_size_y = div_y.quot;

    const int ny_correction = div_y.rem == 0 ? 1 : 0;
    const int nx_correction = div_x.rem == 0 ? 1 : 0;

    const int n_nodes_x = fNDivFineGrid[0] + 1;

    TPZManVector<int64_t, 2> node_id_vec(2, 0);
    int node0, node1;
    for (auto iy = 0; iy < n_full_size_y - ny_correction; iy++) {
        for (auto ix = 0; ix < n_full_size_x; ix++) {
            node0 =
                fNElemCoarseGrid * fNDivFineGrid[0] + (ix + 1) * fNElemCoarseGrid + fNElemCoarseGrid * n_nodes_x * iy;
            node1 = node0 + fNElemCoarseGrid;
            node_id_vec[0] = node0;
            node_id_vec[1] = node1;
            auto * gel = new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(node_id_vec, 666, *fGeoMesh);
            fSkelIdToRefTree.insert({gel->Id(), fRefTreeDesiredSize});
        }
        if (div_x.rem != 0) {
            node0 = fNElemCoarseGrid * fNDivFineGrid[0] + (n_full_size_x + 1) * fNElemCoarseGrid +
                    fNElemCoarseGrid * n_nodes_x * iy;
            node1 = node0 + div_x.rem;
            node_id_vec[0] = node0;
            node_id_vec[1] = node1;
            auto * gel = new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(node_id_vec, 666, *fGeoMesh);
            fSkelIdToRefTree.insert({gel->Id(), fRefTreeRemainderX});
        }
    }

    for (auto ix = 0; ix < n_full_size_x - nx_correction; ix++) {
        for (auto iy = 0; iy < n_full_size_y; iy++) {
            node0 = (fNElemCoarseGrid) * (ix + 1) + n_nodes_x * iy * fNElemCoarseGrid;
            node1 = node0 + n_nodes_x * fNElemCoarseGrid;
            node_id_vec[0] = node0;
            node_id_vec[1] = node1;
            auto * gel = new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(node_id_vec, 666, *fGeoMesh);
            fSkelIdToRefTree.insert({gel->Id(), fRefTreeDesiredSize});
        }
        if (div_y.rem != 0) {
            node0 = (fNElemCoarseGrid) * (ix + 1) + n_nodes_x * n_full_size_y * fNElemCoarseGrid;
            node1 = node0 + n_nodes_x * div_y.rem;
            node_id_vec[0] = node0;
            node_id_vec[1] = node1;
            auto * gel = new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(node_id_vec, 666, *fGeoMesh);
            fSkelIdToRefTree.insert({gel->Id(), fRefTreeRemainderY});
        }
    }
}

void TPZMultiscaleGridGen2D::RefineSkeletonElements() {

    std::function<void(const std::pair<int64_t, RefTree *>)> RefineSkeleton = [&](std::pair<int64_t, RefTree*> data) -> void {
        auto * node = data.second;
        if (node->fActualSize != 1) {

            auto right_size = node->fChildRight->fActualSize;
            auto left_size = node->fChildLeft->fActualSize;

            if (right_size == left_size) {
                right_size = 1;
                left_size = 1;
            }
            const auto ref_pattern = fRefPatterns.find({left_size, right_size});
            if (ref_pattern == fRefPatterns.end()) DebugStop();
            TPZAutoPointer<TPZRefPattern> ref_pattern_ptr(&ref_pattern->second);
            auto * skel = fGeoMesh->Element(data.first);
            skel->SetRefPattern(ref_pattern_ptr);
            TPZManVector<TPZGeoEl *, 2> sons(2, nullptr);
            skel->Divide(sons);

            RefineSkeleton({sons[0]->Id(), node->fChildLeft});
            RefineSkeleton({sons[1]->Id(), node->fChildRight});
        }
    };

    std::for_each(fSkelIdToRefTree.begin(), fSkelIdToRefTree.end(), RefineSkeleton);
}
