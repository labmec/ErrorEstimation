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

    auto * gmesh2x1 = CreateRefinementGeoMesh(2, 1);
    TPZRefPattern ref_pat2x1(*gmesh2x1);
    TPZAutoPointer<TPZRefPattern> ref2x1(&ref_pat2x1);

    auto * gmesh1x2 = CreateRefinementGeoMesh(1, 2);
    TPZRefPattern ref_pat1x2(*gmesh1x2);
    TPZAutoPointer<TPZRefPattern> ref1x2(&ref_pat1x2);

    auto * gmesh1x1 = CreateRefinementGeoMesh(1, 1);
    TPZRefPattern ref_pat1x1(*gmesh1x1);
    TPZAutoPointer<TPZRefPattern> ref1x1(&ref_pat1x1);

    TPZManVector<REAL, 3> coord(3, 0.);
    for (int y = 0; y <= 8; y++) {
        for (int x = 0; x <= 28; x++) {
            coord = {8.0 * x, 8.0 * y, 0};
            if (x == 28) coord[0] = 220.;
            if (y == 8) coord[1] = 60.;

            // Create new node
            const auto newID = gmesh->NodeVec().AllocateNewElement();
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
            auto * gel = new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec, porousMediaMatId, *gmesh);
            if (x == 27 && y != 7) gel->SetRefPattern(ref1x2);
            if (x != 27 && y == 7) gel->SetRefPattern(ref2x1);
            if (x == 27 && y == 7) gel->SetRefPattern(ref1x1);
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

    gmesh->BuildConnectivity();

    TPZManVector<TPZGeoEl*, 4> sons;
    for (int div = 0; div < 3; div++) {
        auto nelem = gmesh->NElements();
        for (int64_t i = 0; i < nelem; i++) {
            auto * gel = gmesh->ElementVec()[i];
            const int has_sub = gel->HasSubElement();
            if (has_sub == 0) {
                gel->Divide(sons);
            }
        }
    }
    return gmesh;
}

[[maybe_unused]] TPZGeoMesh *SPE10::CreateRefinementGeoMesh(const int nx, const int ny) {
    auto* gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    TPZManVector<REAL, 3> coord(3, 0.);

    for (int y = 0; y <= ny; y++) {
        for (int x = 0; x <= nx; x++) {
            coord = {static_cast<REAL>(x), static_cast<REAL>(y), 0.};
            // Create new node
            const auto newID = gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[newID].Initialize(coord, *gmesh);
        }
    }

    constexpr int matId = 1;

    // Inserts quad elements
    TPZManVector<int64_t, 4> nodesIdVec(4, 0);
    nodesIdVec[1] = nx;
    nodesIdVec[2] = (ny + 1) * (nx + 1) - 1;
    nodesIdVec[3] = ny * (nx + 1);
    auto * father_gel = new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec, matId, *gmesh);

    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            nodesIdVec[0] = (nx + 1) * y + x;
            nodesIdVec[1] = (nx + 1) * y + x + 1;
            nodesIdVec[2] = (nx + 1) * (y + 1) + x + 1;
            nodesIdVec[3] = (nx + 1) * (y + 1) + x + 0;
            auto * gel = new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec, matId, *gmesh);
            gel->SetFather(father_gel);
        }
    }

    gmesh->BuildConnectivity();
    return gmesh;
}