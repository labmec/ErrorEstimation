//
// Created by victor on 16/03/2021.
//

#ifndef ERRORESTIMATION_SOLVER_H
#define ERRORESTIMATION_SOLVER_H


#include <TPZMultiphysicsCompMesh.h>
#include "TPZLinearAnalysis.h"
#include "DataStructure.h"
#include "ProblemConfig.h"
#include "TPZMixedErrorEstimate.h"

//// Call required methods to build a computational mesh for an Pryymal Hybrid approximation
void CreateHybridH1ComputationalMesh(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, int &InterfaceMatId, int &fluxMatID,PreConfig &eData, ProblemConfig &config,int hybridLevel);

//// Call required methods to build a computational mesh for a Mixed approximation
void CreateMixedComputationalMesh(TPZMultiphysicsCompMesh *cmesh_H1Mixed,PreConfig &eData, ProblemConfig &config);
void CreateCondensedMixedElements(TPZMultiphysicsCompMesh *cmesh_Mixed);
// Cretes multiphysics computational mesh to compute the solution difference between Hyb and Mix approximations.
void CreateHybMixCompMesh(TPZMultiphysicsCompMesh *multiHyb, TPZMultiphysicsCompMesh *multiMix, TPZMultiphysicsCompMesh *multiHybMix,PreConfig &hybConfig, ProblemConfig &ConfHyb);

//// Solve classical H1 problem
void SolveH1Problem(TPZCompMesh *cmeshH1,struct ProblemConfig &config, struct PreConfig &eData);

//// Solve Primal Hybrid problem
void SolveHybridH1Problem(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, int InterfaceMatId, struct ProblemConfig &config, struct PreConfig &eData,int hybridLevel);

//// Solve Mixed problem
void SolveMixedProblem(TPZMultiphysicsCompMesh *cmesh_Mixed,struct ProblemConfig config,struct PreConfig &eData);

//// Error Management
void StockErrorsH1(TPZAnalysis &an,TPZCompMesh *cmesh,std::ofstream &Erro, TPZVec<REAL> *Log, PreConfig &eData);

//// Error Management
void StockErrors(TPZAnalysis &an,TPZMultiphysicsCompMesh *cmesh,std::ofstream &Erro, TPZVec<REAL> *Log, PreConfig &eData);

//// Solve desired problem
void Solve(ProblemConfig &config, PreConfig &preConfig);
//// Solve hybrid and mixed problem and compute its difference
void SolveDiff(PreConfig &hybConfig, PreConfig &mixConfig,char *argv[]);
//// Print errors for simutaneously Hyb and Mix simulations
void PrintErrorsDiff(TPZVec<REAL> errorVec, ProblemConfig &config);

void EstimateError(ProblemConfig &config, PreConfig &preConfig, int fluxMatID, TPZMultiphysicsCompMesh *multiCmesh);

//// Draw geometric and computational mesh
void DrawMesh(ProblemConfig &config, PreConfig &preConfig, TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *multiCmesh);

void PostProcessHybMix(TPZMultiphysicsCompMesh *multiHybMix,PreConfig &hybConfig,ProblemConfig &config);

void FluxErrorCreateCompMesh(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, int &InterfaceMatId, int &fluxMatID,PreConfig &eData, ProblemConfig &config);

bool PostProcessing(TPZCompMesh * pressuremesh, TPZFMatrix<STATE> true_elerror, TPZFMatrix<STATE> estimate_elerror, ProblemConfig &config);

#endif //ERRORESTIMATION_SOLVER_H
