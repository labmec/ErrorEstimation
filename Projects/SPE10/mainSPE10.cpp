//
// Created by Gustavo Batistela on 4/5/21.
//

#include <TPZGenGrid3D.h>
#include <iostream>
#include <memory>
#include <pzgmesh.h>
#include "Tools.h"

TPZGeoMesh *CreateSPE10GeoMesh();
void ReadSPE10CellPermeabilities(TPZFMatrix<REAL>* perm_mat);

int main() {
    constexpr int n_cells = 60 * 220 * 85;
    auto perm_mat = std::make_unique<TPZFNMatrix<n_cells * 3, REAL>>(n_cells, 3);
    ReadSPE10CellPermeabilities(perm_mat.get());

    REAL max_perm = 0.;
    REAL min_perm = 1e15;
    for (int i = 0; i < perm_mat->Rows(); i++) {
        for (int j = 0; j < perm_mat->Cols(); j++) {
            const REAL perm = perm_mat->operator()(i, j);
            if (perm > max_perm) max_perm = perm;
            if (perm < min_perm) min_perm = perm;
        }
    }

    TPZGeoMesh *gmesh = CreateSPE10GeoMesh();

    int n_cells2 = 0;

    constexpr int c1 = 85;
    constexpr int c2 = 220;
    constexpr int c3 = 60;
    for (int i = 0; i < c1; i++) {
        for (int j = 0; j < c2; j++) {
            for (int k = 0; k < c3; k++) {
                const int cell_id = c2 * c3 * i + c3 * j + k;
                const double perm = perm_mat->operator()(cell_id, 0);
                const double relative_perm = (perm - min_perm) / (max_perm - min_perm);
                const int mat_id = static_cast<int>(std::round(255 * relative_perm));
                gmesh->Element(cell_id)->SetMaterialId(mat_id);
                n_cells2++;
            }
        }
    }
    Tools::PrintGeometry(gmesh, "SPE10Grid", false, true);
    std::cout << "SPE10 initial grid created. NElem: " << gmesh->NElements() << "\n";



    return 0;
}

TPZGeoMesh *CreateSPE10GeoMesh() {
    std::cout << "Creating SPE10 initial grid...\n";

    const TPZManVector<REAL, 3> x0 = {0, 0, 0};
    const TPZManVector<REAL, 3> x1 = {1, 11./3, 85.};
    const TPZManVector<int, 3> ndiv = {60, 220, 85};

    constexpr std::array bcIDs = {-1, -1, -1, -1, -1, -1};
    constexpr int domainID = 1;

    TPZGenGrid3D gen(x0, x1, ndiv, MMeshType::EHexahedral);
    TPZGeoMesh * gmesh = gen.BuildVolumetricElements(domainID);
    //TPZGeoMesh *gmesh = gen.BuildBoundaryElements(bcIDs[0], bcIDs[1], bcIDs[2], bcIDs[3], bcIDs[4], bcIDs[5]);

    std::cout << "SPE10 initial grid created. NElem: " << gmesh->NElements() << "\n";

    return gmesh;
}

void ReadSPE10CellPermeabilities(TPZFMatrix<REAL>* perm_mat) {

    std::cout << "Reading permeability data...\n";

    std::ifstream perm_file("InputData/spe_perm.dat", std::ios::in);
    if (!perm_file) {
        std::cerr << "Unable to open input file\n";
        DebugStop();
    }

    const int n_cells = perm_mat->Rows();
    int cell_id = 0;
    int coord_id = 0;
    int line_num = 0;
    while (perm_file) {
        line_num++;
        std::string line;
        std::getline(perm_file, line, '\n');
        if (line.length() <= 1) continue;
        std::stringstream stream(line);
        if (cell_id == n_cells) {
            cell_id = 0;
            coord_id++;
            if (coord_id > 2) DebugStop();
        }
        for (int i = 0; i < 6; i++) {
            stream >> perm_mat->operator()(cell_id, coord_id);
            cell_id++;
        }
    }
    std::cout << "Finished reading permeability data from input file!\n";
}
