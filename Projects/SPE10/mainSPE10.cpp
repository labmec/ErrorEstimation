//
// Created by Gustavo Batistela on 4/5/21.
//

#include "Tools.h"
#include <TPZGenGrid2D.h>
#include <TPZGenGrid3D.h>
#include <iostream>
#include <memory>
#include <pzgmesh.h>
#include <libInterpolate/Interpolate.hpp>

typedef _2D::BicubicInterpolator<REAL> Interpolator;

// Global variables
Interpolator interpolator;

// Function declarations
TPZGeoMesh *CreateSPE10GeoMesh();
TPZCompMesh *CreateSPE10CompMesh(TPZGeoMesh *gmesh);
void ReadSPE10CellPermeabilities(TPZVec<REAL>*perm_vec, int layer);
void PermeabilityFunction(const TPZVec<REAL> &x, TPZVec<REAL> &res, TPZFMatrix<REAL> &res_mat);
void InsertMaterials(TPZCompMesh *cmesh);

int main() {

    constexpr int layer = 35;
    constexpr int nx = 220;
    constexpr int ny = 60;
    constexpr int n_cells = nx * ny;

    auto perm_vec = TPZManVector<REAL, n_cells>(n_cells, 1);
    ReadSPE10CellPermeabilities(&perm_vec, layer);

    TPZGeoMesh *gmesh = CreateSPE10GeoMesh();
    std::cout << "SPE10 initial grid created. NElem: " << gmesh->NElements() << "\n";

    std::vector<REAL> x, y, perm;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            const int cell_id = ny * i + j;
            const double cell_perm = perm_vec[cell_id + nx * ny * layer];
            x.push_back(0.5 + nx);
            y.push_back(0.5 + ny);
            perm.push_back(cell_perm);
        }
    }

    interpolator.setData(x.size(), x.data(), y.data(), perm.data());

    TPZCompMesh *cmesh = CreateSPE10CompMesh(gmesh);

    return 0;
}

TPZGeoMesh *CreateSPE10GeoMesh() {
    std::cout << "Creating SPE10 initial grid...\n";

    const TPZManVector<REAL, 3> x0 = {0, 0, 0};
    const TPZManVector<REAL, 3> x1 = {60., 220., 0.};
    const TPZManVector<int, 3> ndiv = {60, 220, 0};

    TPZGenGrid2D gen(ndiv, x0, x1);

    gen.SetRefpatternElements(true);
    auto gmesh = new TPZGeoMesh;
    gen.Read(gmesh);

    gen.SetBC(gmesh, 4, -2);
    gen.SetBC(gmesh, 5, -2);
    gen.SetBC(gmesh, 6, -1);
    gen.SetBC(gmesh, 7, -2);

    std::cout << "SPE10 initial grid created. NElem: " << gmesh->NElements() << "\n";

    return gmesh;
}

void ReadSPE10CellPermeabilities(TPZVec<REAL> *perm_vec, const int layer) {

    std::cout << "Reading permeability data...\n";

    std::ifstream perm_file("InputData/spe_perm.dat", std::ios::in);
    if (!perm_file) {
        std::cerr << "Unable to open input file\n";
        DebugStop();
    }

    int cell_id = 0;
    const int n_cells = perm_vec->size();
    const int start_line = 1 + n_cells * (layer - 1) / 6;

    int line_num = 0;
    int line_num2 = 0;
    while (perm_file) {
        line_num++;
        line_num2++;
        std::string line;
        std::getline(perm_file, line, '\n');

        if (line_num < start_line) continue;

        std::stringstream stream(line);
        for (int i = 0; i < 6; i++) {
            stream >> perm_vec->operator[](cell_id);
            cell_id++;
        }
        if (cell_id == n_cells) break;
    }
    std::cout << "Finished reading permeability data from input file!\n";
}

void PermeabilityFunction(const TPZVec<REAL> &x, TPZVec<REAL> &res, TPZFMatrix<REAL> &res_mat) {
    res.resize(1);
    res[0] = 1.5;
}

TPZCompMesh *CreateSPE10CompMesh(TPZGeoMesh *gmesh) {
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    int dim = gmesh->Dimension();
    cmesh->SetDimModel(dim);

    InsertMaterials(cmesh);

}

void InsertMaterials(TPZCompMesh *cmesh) {

    auto *mix = new TPZMixedPoisson(1, cmesh->Dimension());

    TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
    constexpr int dirichlet_bc = 0;
    constexpr int neumann_bc = 1;

    // Zero flux (reservoir boundary)
    TPZBndCond *zero_flux = mix->CreateBC(mix, -1, neumann_bc, val1, val2);
    // Unit pressure at left reservoir boundary
    val2(0, 0) = 1.;
    TPZBndCond *pressure_left = mix->CreateBC(mix, -2, dirichlet_bc, val1, val2);

    cmesh->InsertMaterialObject(mix);
    cmesh->InsertMaterialObject(zero_flux);
    cmesh->InsertMaterialObject(pressure_left);

}