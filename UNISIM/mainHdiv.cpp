#include <Material/TPZVecL2.h>
#include <Material/pzbndcond.h>
#include <Matrix/pzstepsolver.h>
#include <Mesh/TPZCompMeshTools.h>
#include <Post/TPZVTKGeoMesh.h>
#include <Pre/TPZGmshReader.h>
#include <Pre/TPZHybridizeHDiv.h>
#include <StrMatrix/pzskylstrmatrix.h>
#include <libInterpolate/Interpolate.hpp>

#include <iostream>
#include <map>
#include <vector>

TPZGeoMesh *CreateFlatGeoMesh(std::string &geometry_file2D);

void ModifyZCoordinates(TPZGeoMesh *gmesh, std::string &filename);

void ReadReservoirGeometryData(const std::string &name, std::vector<double> &x,
                               std::vector<double> &y, std::vector<double> &z);

void PrintGeometry(TPZGeoMesh *gmesh, const std::string &file_name);

void UNISIMHDiv();

int main() {

    InitializePZLOG();
    UNISIMHDiv();

    return 0;
}

void UNISIMHDiv() {

    std::string geometry_file2D = "InputData/UNISIMFlatMesh.msh";
    TPZGeoMesh *gmesh = CreateFlatGeoMesh(geometry_file2D);

    std::string name = "InitialGeoMesh";
    PrintGeometry(gmesh, name);
}

TPZGeoMesh *CreateFlatGeoMesh(std::string &geometry_file2D) {

    TPZGmshReader Geometry;
    TPZGeoMesh *gmesh2d;
    REAL l = 1.0;
    Geometry.SetCharacteristiclength(l);
    Geometry.SetFormatVersion("4.1");
    gmesh2d = Geometry.GeometricGmshMesh(geometry_file2D);
    Geometry.PrintPartitionSummary(std::cout);

    std::string filename1 = "InputData/UNISIMPointCloud.txt";
    ModifyZCoordinates(gmesh2d, filename1);
    return gmesh2d;
}

void ModifyZCoordinates(TPZGeoMesh *gmesh, std::string &filename) {

    std::vector<double> x, y, z;
    ReadReservoirGeometryData(filename, x, y, z);

    _2D::ThinPlateSplineInterpolator<double> interp;
    interp.setData(x, y, z);

    int64_t nCoordinates = gmesh->NodeVec().NElements();
    double sum = 0.0;
    for (auto val : z) {
        sum += val;
    }
    double val_storage = sum / z.size();

    for (int icoord = 0; icoord < nCoordinates; icoord++) {
        TPZGeoNode node = gmesh->NodeVec()[icoord];
        TPZVec<REAL> co(3);
        node.GetCoordinates(co);
        double val_interp = interp(co[0], co[1]);

        if (val_interp == 0.0) {
            co[2] = val_storage;
        }
        if (val_interp > 1.0) {
            co[2] = val_interp;
        }

        gmesh->NodeVec()[icoord].SetCoord(co);
    }
}

void ReadReservoirGeometryData(const std::string &name, std::vector<double> &x,
                               std::vector<double> &y, std::vector<double> &z) {

    std::ifstream file;
    file.open(name);

    std::string line;
    int i = 1;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        char l = line[0];
        if (l != '/') {
            i = i + 1;
            int val = i % 15;
            if (val == 0) {
                double a, b, c;
                if (iss >> a >> b >> c)
                    ;
                x.push_back(a);
                y.push_back(b);
                z.push_back(c);
            }
        }
    }

    if (x.empty()) {
        std::cout << "No data read.\n";
        DebugStop();
    }

    std::cout << "File successfully read!\n";
}

void PrintGeometry(TPZGeoMesh *gmesh, const std::string &file_name) {
    std::stringstream txt_name;
    std::stringstream vtk_name;
    txt_name << file_name << ".txt";
    vtk_name << file_name << ".vtk";
    std::ofstream textfile(txt_name.str().c_str());
    gmesh->Print(textfile);
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
}
