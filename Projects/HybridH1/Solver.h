//
// Created by victor on 16/03/2021.
//

#ifndef ERRORESTIMATION_SOLVER_H
#define ERRORESTIMATION_SOLVER_H


#include <TPZMultiphysicsCompMesh.h>
#include "pzanalysis.h"
#include "DataStructure.h"
#include "ProblemConfig.h"

//// Call required methods to build a computational mesh for an Pryymal Hybrid approximation
void CreateHybridH1ComputationalMesh(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, int &InterfaceMatId, int &fluxMatID,PreConfig &eData, ProblemConfig &config,int hybridLevel);

//// Call required methods to build a computational mesh for a Mixed approximation
void CreateMixedComputationalMesh(TPZMultiphysicsCompMesh *cmesh_H1Mixed,PreConfig &eData, ProblemConfig &config);
void CreateCondensedMixedElements(TPZMultiphysicsCompMesh *cmesh_Mixed);

//// Solve classical H1 problem
void SolveH1Problem(TPZCompMesh *cmeshH1,struct ProblemConfig &config, struct PreConfig &eData);

//// Solve Primal Hybrid problem
void SolveHybridH1Problem(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, int InterfaceMatId, struct ProblemConfig config, struct PreConfig &eData,int hybridLevel);

//// Solve Mixed problem
void SolveMixedProblem(TPZMultiphysicsCompMesh *cmesh_Mixed,struct ProblemConfig config,struct PreConfig &eData);

//// Error Management
void StockErrorsH1(TPZAnalysis &an,TPZCompMesh *cmesh,std::ofstream &Erro, TPZVec<REAL> *Log, PreConfig &eData);

//// Error Management
void StockErrors(TPZAnalysis &an,TPZMultiphysicsCompMesh *cmesh,std::ofstream &Erro, TPZVec<REAL> *Log, PreConfig &eData);

//// Solve desired problem
void Solve(ProblemConfig &config, PreConfig &preConfig);

void EstimateError(ProblemConfig &config, PreConfig &preConfig, int fluxMatID, TPZMultiphysicsCompMesh *multiCmesh);

//// Draw geometric and computational mesh
void DrawMesh(ProblemConfig &config, PreConfig &preConfig, TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *multiCmesh);


#endif //ERRORESTIMATION_SOLVER_H
