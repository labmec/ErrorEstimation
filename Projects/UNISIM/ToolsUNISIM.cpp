//
// Created by Gustavo A. Batistela on 25/07/2020.
//

#include "ToolsUNISIM.h"
#include "Tools.h"

#include <fstream>
#include <string>
#include <vector>

#include <Post/TPZVTKGeoMesh.h>
#include <Pre/TPZGmshReader.h>
#include <Refine/TPZRefPatternTools.h>

#include <libInterpolate/Interpolate.hpp>

void ReadReservoirGeometryData(const std::string &name, std::vector<double> &x, std::vector<double> &y,
                               std::vector<double> &z) {
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

TPZGeoMesh *CreateUNISIMSurfaceGeoMesh(bool modifyZCoordinates) {

    std::string gmshFile = "InputData/UNISIMFlatMesh.msh";
#ifdef MACOSX
    gmshFile = "../" + gmshFile;
#endif

    TPZGmshReader gmeshReader;
    TPZGeoMesh *gmesh;

    TPZManVector<std::map<std::string, int>, 4> meshMaterialData(3);
    // BC materials
    meshMaterialData[1].insert({"ZeroFlux", -1});
    meshMaterialData[1].insert({"Productors", -2});
    meshMaterialData[1].insert({"Injectors", -3});
    meshMaterialData[1].insert({"Faults", 99}); // Won't be used
    // Domain materials
    meshMaterialData[2].insert({"RockMatrix", 1});
    meshMaterialData[2].insert({"RockMatrix2", 1});
    gmeshReader.SetDimNamePhysical(meshMaterialData);

    gmeshReader.SetFormatVersion("4.1");
    gmesh = gmeshReader.GeometricGmshMesh(gmshFile);

    std::string filename = "InputData/UNISIMPointCloud.txt";

    if (modifyZCoordinates) {
        ModifyZCoordinates(gmesh, filename);
    }

    MoveMeshToOrigin(gmesh);

    return gmesh;
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
    gmesh->SetDimension(3);
}

void MoveMeshToOrigin(TPZGeoMesh *gmesh) {
    TPZManVector<REAL, 3> nod0(3, 0.);
    gmesh->NodeVec()[0].GetCoordinates(nod0);
    int64_t nnodes = gmesh->NNodes();
    for (int64_t no = 0; no < nnodes; no++) {
        TPZManVector<REAL, 3> co(3);
        gmesh->NodeVec()[no].GetCoordinates(co);
        for (int ic = 0; ic < 3; ic++) {
            co[ic] -= nod0[ic];
        }
        gmesh->NodeVec()[no].SetCoord(co);
    }
}

void PrintGeometry(TPZGeoMesh *gmesh, const std::string &file_name, bool printTXT, bool printVTK) {
    if (printTXT) {
        std::stringstream txt_name;
        txt_name << file_name << ".txt";
        std::ofstream textfile(txt_name.str().c_str());
        gmesh->Print(textfile);
    }
    if (printVTK) {
        std::stringstream vtk_name;
        vtk_name << file_name << ".vtk";
        std::ofstream vtkfile(vtk_name.str().c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
    }
}

void hAdaptivity(TPZGeoMesh *gmesh, TPZVec<REAL> &elementErrors, REAL thresholdRatio) {

    REAL maxError = 0.;
    std::cout << "Starting h-adaptivity procedure...\n"
              << "Number of elements before refinement: " << gmesh->NElements() << '\n';
    int64_t nelem = gmesh->NElements();
    for (int64_t iel = 0; iel < nelem; iel++) {
        if (elementErrors[iel] > maxError) {
            maxError = elementErrors[iel];
        }
    }

    REAL threshold = thresholdRatio * maxError;

    for (int64_t iel = 0; iel < nelem; iel++) {
        REAL elementError = elementErrors[iel];
        if (elementError > threshold) {
            TPZGeoEl *gel = gmesh->Element(iel);
            if (!gel) DebugStop();
            if (gel->Dimension() != gmesh->Dimension()) DebugStop();
            TPZVec<TPZGeoEl *> sons;
            if (!gel->HasSubElement()) {
                gel->Divide(sons);
            }
        }
    }

    PrintGeometry(gmesh, "gmeshBeforeSpread");
    SpreadMeshRefinement(gmesh);
    DivideLowerDimensionalElements(gmesh);

    std::cout << "Number of elements after refinement: " << gmesh->NElements() << '\n';
}

// Refines elements that has two (or more) refined neighbours or a neighbour that has been refined twice
void SpreadMeshRefinement(TPZGeoMesh *gmesh) {
    bool hasChanged = true;
    int dim = gmesh->Dimension();
    while (hasChanged) {
        hasChanged = false;

        TPZStack<TPZGeoEl *> gelsToRefine;
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (!gel || gel->Dimension() != dim || gel->HasSubElement()) {
                continue;
            }

            int nsides = gel->NSides();
            int ncorner = gel->NCornerNodes();
            int nRefinedNeighbours = 0;
            bool needsRefinement = false;
            for (int side = ncorner; side < nsides; side++) {
                TPZGeoElSide gelside(gel, side);
                TPZGeoElSide neighbour(gelside.Neighbour());
                while (neighbour != gelside) {
                    if (neighbour.HasSubElement() && neighbour.NSubElements() >= 1) {
                        nRefinedNeighbours++;
                        // Check if two neighbours have been refined
                        if (nRefinedNeighbours >= 2) {
                            needsRefinement = true;
                        }

                        TPZStack<TPZGeoElSide> subNeighs;
                        neighbour.GetSubElements2(subNeighs);
                        for (int i = 0; i < subNeighs.size(); i++) {
                            // Check if an neighbour has been refined twice
                            if (subNeighs[i].NSubElements() > 1) {
                                needsRefinement = true;
                                break;
                            }
                            if (needsRefinement) break;
                        }
                    }
                    if (needsRefinement) break;
                    neighbour = neighbour.Neighbour();
                }
                if (needsRefinement) {
                    gelsToRefine.Push(gel);
                    break;
                }
            }
        }
        if (gelsToRefine.size()) {
            hasChanged = true;
            for (int64_t i = 0; i < gelsToRefine.size(); i++) {
                TPZManVector<TPZGeoEl *> subEls;
                gelsToRefine[i]->Divide(subEls);
            }
        }
    }
}

void ApplyDirectionalRefinement(TPZGeoMesh *gmesh, int nRef) {
    // Mat IDs of productors and injectors BCs
    set<int> matids{-2, -3};

    for (auto i = 0; i < nRef; i++) {
        int nelements = gmesh->NElements();
        for (auto el = 0; el < nelements; el++) {
            TPZGeoEl *element = gmesh->ElementVec()[el];
            if (!element) continue;
            TPZRefPatternTools::RefineDirectional(element, matids);
        }
        cout << "Refinement step: " << i << "\nNumber of elements = " << gmesh->NElements() << '\n';
        stringstream meshfilename;
        meshfilename << "DirectionalRefinementGMesh" << i;
        PrintGeometry(gmesh, meshfilename.str(), false, true);
    }
}
