//
// Created by victor on 16/03/2021.
//

#ifndef ERRORESTIMATION_INPUTTREATMENT_H
#define ERRORESTIMATION_INPUTTREATMENT_H

#include "DataStructure.h"
#include "ProblemConfig.h"

void EvaluateEntry(int argc, char *argv[],PreConfig &pConfig);
void InitializeOutstream(PreConfig &pConfig,char *argv[]);
void IsInteger(char *argv);
void Configure(ProblemConfig &config,int ndiv,PreConfig &pConfig,char *argv[]);
void ReadEntry(ProblemConfig &config, PreConfig &pConfig);
void ConfigureNFconvergence(ProblemConfig &config,PreConfig &pConfig);

#endif //ERRORESTIMATION_INPUTTREATMENT_H
