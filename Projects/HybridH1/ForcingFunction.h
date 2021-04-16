//
// Created by victor on 16/03/2021.
//

#include "pzvec.h"
#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzfunction.h"

#include <string>

#ifndef ERRORESTIMATION_FORCINGFUNCTION_H
#define ERRORESTIMATION_FORCINGFUNCTION_H

void LinearFunc( const TPZVec<REAL> &x, TPZVec<REAL> &g, TPZFMatrix<REAL> &dg)
{
    dg.Resize(3,1);
    dg(0,0) = 0;
    dg(1,0) = 0.5-x[0];
    dg(2,0) = 0;
}

#endif //ERRORESTIMATION_FORCINGFUNCTION_H
