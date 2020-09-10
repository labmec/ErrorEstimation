//
// Created by Gustavo A. Batistela on 25/07/2020.
//

#ifndef ERRORESTIMATION_TOOLSUNISIM_H
#define ERRORESTIMATION_TOOLSUNISIM_H

#include <Mesh/pzgmesh.h>

void ReadReservoirGeometryData(const std::string &name, std::vector<double> &x, std::vector<double> &y,
                               std::vector<double> &z);

TPZGeoMesh *CreateUNISIMSurfaceGeoMesh(bool modifyZCoordinates);

void ModifyZCoordinates(TPZGeoMesh *gmesh, std::string &filename);

void MoveMeshToOrigin(TPZGeoMesh *gmesh);

void PrintGeometry(TPZGeoMesh *gmesh, const std::string &file_name, bool printTXT = true, bool printVTK = true);

void hAdaptivity(TPZGeoMesh *gmesh, TPZVec<REAL> &elementErrors, REAL thresholdRatio);

void SpreadMeshRefinement(TPZGeoMesh *gmesh);

void ApplyDirectionalRefinement(TPZGeoMesh *gmesh, int nRef);

#endif // ERRORESTIMATION_TOOLSUNISIM_H
