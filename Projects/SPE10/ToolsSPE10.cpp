//
// Created by Gustavo Batistela on 10/26/21.
//

#include "ToolsSPE10.h"
#include <TPZGenGrid2D.h>
#include <pzgeoquad.h>
#include <tpzgeoelrefpattern.h>

[[maybe_unused]] TPZGeoMesh *SPE10::CreateFineGridGeoMesh() {

    TPZManVector<int, 4> bcIDs = {-2, -1, -2, -2};
    TPZManVector<int, 2> n_elements = {220, 60};
    TPZManVector<REAL, 3> x0 = {0, 0, 0};
    TPZManVector<REAL, 3> x1 = {220, 60, 0};

    TPZGenGrid2D gen(n_elements, x0, x1, 1, 0);
    gen.SetRefpatternElements(true);

    auto *gmesh = new TPZGeoMesh;
    gen.Read(gmesh);

    gen.SetBC(gmesh, 4, bcIDs[0]);
    gen.SetBC(gmesh, 5, bcIDs[1]);
    gen.SetBC(gmesh, 6, bcIDs[2]);
    gen.SetBC(gmesh, 7, bcIDs[3]);

    gmesh->SetDimension(2);

    return gmesh;
}

[[maybe_unused]] TPZGeoMesh *SPE10::CreateMHMGeoMesh() {

    auto* gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);

    TPZManVector<REAL, 3> coord(3, 0.);

    for (int y = 0; y <= 8; y++) {
        for (int x = 0; x <= 28; x++) {
            coord = {8.0 * x, 8.0 * y, 0};
            if (x == 28) coord[0] = 220.;
            if (y == 8) coord[1] = 60.;

            // Create new node
            const auto newID = gmesh->NodeVec().AllocateNewElement();
            std::cout << newID << ": " << coord[0] << ", " << coord[1] << ", " << coord[2] << '\n';
            gmesh->NodeVec()[newID].Initialize(coord, *gmesh);
        }
    }

    constexpr int porousMediaMatId = 1;
    constexpr int bcMatId1 = -1;
    constexpr int bcMatId2 = -2;

    // Inserts quad elements
    TPZManVector<int64_t, 4> nodesIdVec(4);
    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 28; x++) {
            nodesIdVec[0] = 29 * y + x;
            nodesIdVec[1] = 29 * y + x + 1;
            nodesIdVec[2] = 29 * (y + 1) + x + 1;
            nodesIdVec[3] = 29 * (y + 1) + x + 0;
            new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec, porousMediaMatId, *gmesh);
        }
    }

    nodesIdVec.resize(2);
    for (int x = 0; x < 28; x++) {
        nodesIdVec[0] = x;
        nodesIdVec[1] = x + 1;
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, bcMatId1, *gmesh);
        nodesIdVec[0] = x + 8 * 29;
        nodesIdVec[1] = x + 8 * 29 + 1;
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, bcMatId1, *gmesh);
    }

    for (int y = 0; y < 8; y++) {
        nodesIdVec[0] = y * 29;
        nodesIdVec[1] = (y + 1) * 29;
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, bcMatId1, *gmesh);
        nodesIdVec[0] = y * 29 + 28;
        nodesIdVec[1] = (y + 1) * 29 + 28;
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, bcMatId2, *gmesh);
    }

    //for (int64_t i = 0; i < 6; i++) {
    //    nodesIdVec[0] = 0;
    //    nodesIdVec[1] = 1 + 2 * i;
    //    nodesIdVec[2] = 3 + 2 * i;
    //    new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle>>(nodesIdVec, matIdTriangle, *gmesh);
    //}
    //// Inserts arc elements
    //for (int64_t i = 0; i < 6; i++) {
    //    nodesIdVec[0] = 1 + 2 * i;
    //    nodesIdVec[1] = 3 + 2 * i;
    //    nodesIdVec[2] = 2 + 2 * i;
    //    new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, matIdArc, *gmesh);
    //}
    //// Finally, inserts line elements to complete boundary
    //nodesIdVec.Resize(2);
    //nodesIdVec[0] = 0;
    //nodesIdVec[1] = 1;
    //new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, matIdArc, *gmesh);

    //nodesIdVec[0] = 0;
    //nodesIdVec[1] = 13;
    //new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, matIdArc, *gmesh);

    gmesh->BuildConnectivity();
    return gmesh;
}
