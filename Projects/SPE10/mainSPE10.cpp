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

    //std::cout << "Creating SPE10 initial grid...\n";
    //TPZGeoMesh *gmesh = CreateSPE10GeoMesh();
    //Tools::PrintGeometry(gmesh, "SPE10Grid", false, true);
    //std::cout << "SPE10 initial grid created. NElem: " << gmesh->NElements() << "\n";



    return 0;
}

TPZGeoMesh *CreateSPE10GeoMesh() {
    const TPZManVector<REAL, 3> x0 = {0, 0, 0};
    const TPZManVector<REAL, 3> x1 = {1200, 2200, 170};
    const TPZManVector<int, 3> ndiv = {60, 220, 85};

    constexpr std::array bcIDs = {-1, -1, -1, -1, -1, -1};
    constexpr int domainID = 1;

    TPZGenGrid3D gen(x0, x1, ndiv, MMeshType::EHexahedral);
    TPZGeoMesh * gmesh = gen.BuildVolumetricElements(domainID);
    //TPZGeoMesh *gmesh = gen.BuildBoundaryElements(bcIDs[0], bcIDs[1], bcIDs[2], bcIDs[3], bcIDs[4], bcIDs[5]);

    return gmesh;
}

void ReadSPE10CellPermeabilities(TPZFMatrix<REAL>* perm_mat) {

    std::ifstream perm_file("InputData/spe_perm.dat", std::ios::in);
    if (!perm_file) {
        std::cerr << "Unable to open input file\n";
        DebugStop();
    }

    int n_cells = perm_mat->Rows();
    int cell_id = 0;
    int coord_id = 0;
    int line_num = 0;
    while (perm_file) {
        line_num++;
        std::string line;
        std::getline(perm_file, line, '\n');
        if (line.length() <= 1) continue;
        std::stringstream stream(line);
        //std::cout << stream.str() << std::endl;
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
