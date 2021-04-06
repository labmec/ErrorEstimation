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
constexpr int n_cells = 60 * 220 * 85;
auto perm_mat = TPZFNMatrix<n_cells * 3, REAL>(n_cells, 3);
Interpolator interpolator;

// Function declarations
TPZGeoMesh *CreateSPE10GeoMesh();
TPZCompMesh *CreateSPE10CompMesh(TPZGeoMesh *gmesh);
void ReadSPE10CellPermeabilities();
void PermeabilityFunction(const TPZVec<REAL> &x, TPZVec<REAL> &res, TPZFMatrix<REAL> &res_mat);
void InsertMaterials(TPZCompMesh *cmesh);

int main() {
    ReadSPE10CellPermeabilities();

    REAL max_perm = 0.;
    REAL min_perm = 1e15;
    for (int i = 0; i < perm_mat.Rows(); i++) {
        for (int j = 0; j < perm_mat.Cols(); j++) {
            const REAL perm = perm_mat.operator()(i, j);
            if (perm > max_perm) max_perm = perm;
            if (perm < min_perm) min_perm = perm;
        }
    }

    TPZGeoMesh *gmesh = CreateSPE10GeoMesh();
    Tools::PrintGeometry(gmesh, "SPE10Grid", false, true);
    std::cout << "SPE10 initial grid created. NElem: " << gmesh->NElements() << "\n";

    //constexpr int layer = 35;
    //constexpr int nx = 220;
    //constexpr int ny = 60;
    //for (int i = 0; i < nx; i++) {
    //    for (int j = 0; j < ny; j++) {
    //        const int cell_id = ny * i + j;
    //        const double perm = perm_mat(cell_id + nx * ny * layer, 0);
    //        const double relative_perm = (perm - min_perm) / (max_perm - min_perm);
    //        const int mat_id = static_cast<int>(std::round(255 * relative_perm));
    //        gmesh->Element(cell_id)->SetMaterialId(mat_id);
    //    }
    //}

    TPZCompMesh *cmesh = CreateSPE10CompMesh(gmesh);

    return 0;
}

TPZGeoMesh *CreateSPE10GeoMesh() {
    std::cout << "Creating SPE10 initial grid...\n";

    const TPZManVector<REAL, 3> x0 = {0, 0, 0};
    const TPZManVector<REAL, 3> x1 = {1, 11./3, 0.};
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

void ReadSPE10CellPermeabilities() {

    std::cout << "Reading permeability data...\n";

    std::ifstream perm_file("InputData/spe_perm.dat", std::ios::in);
    if (!perm_file) {
        std::cerr << "Unable to open input file\n";
        DebugStop();
    }

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
            stream >> perm_mat(cell_id, coord_id);
            cell_id++;
        }
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