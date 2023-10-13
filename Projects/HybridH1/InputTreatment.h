//
// Created by victor on 16/03/2021.
//

#ifndef ERRORESTIMATION_INPUTTREATMENT_H
#define ERRORESTIMATION_INPUTTREATMENT_H

#include "DataStructure.h"
#include "ProblemConfig.h"

void DataInitialization(int argc, char *argv[],PreConfig &hybConfig,PreConfig &mixConfig);
void EvaluateEntry(int argc, char *argv[],PreConfig &pConfig);
void InitializeOutstream(PreConfig &pConfig,char *argv[]);
void IsInteger(char *argv);
void Configure(ProblemConfig &config,int ndiv,PreConfig &pConfig,char *argv[]);
void ReadEntry(ProblemConfig &config, PreConfig &pConfig);
TLaplaceExample1::EExactSol ChooseAnaliticSolution(PreConfig &preConfig);
void FluxErrorConfigure(ProblemConfig &config,PreConfig &pConfig);
void CopyHybSetup(PreConfig &hybConfig, PreConfig &mixConfig);

#endif //ERRORESTIMATION_INPUTTREATMENT_H
