//
// Created by victor on 16/03/2021.
//



#ifndef ERRORESTIMATION_FORCINGFUNCTION_H
#define ERRORESTIMATION_FORCINGFUNCTION_H

#include "pzvec.h"
#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzfunction.h"

#include <string>

void LinearFunc( const TPZVec<REAL> &x, TPZVec<REAL> &g, TPZFMatrix<REAL> &dg);

void SingularityExact( const TPZVec<REAL> &x, TPZVec<REAL> &g, TPZFMatrix<REAL> &dg);

void SingularityForcingFunction( const TPZVec<REAL> &x, TPZVec<REAL> &f);

#endif //ERRORESTIMATION_FORCINGFUNCTION_H
