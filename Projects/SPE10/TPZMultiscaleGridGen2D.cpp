//
// Created by Gustavo Batistela on 11/8/21.
//

#include "TPZMultiscaleGridGen2D.h"
#include <tpzgeoelrefpattern.h>

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