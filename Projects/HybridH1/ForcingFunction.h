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

void ExactSolution( const TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix<REAL> &deriv)
{
    u[0] = 2*x[1]-1;
}

TPZAutoPointer<TPZFunction<STATE> > ExactSol = new TPZDummyFunction<STATE> (ExactSolution);

class BCfunction : public virtual TPZFunction {
public:
     void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df) override
    {

    }
}

struct LinearBCquadrilateral : public TPZAnalyticSolution{
    void LinearBCquadrilateral::Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &sigma) const
    {
        sigma.Resize(2,1);
        sigma(0,0) = 2*x[1]-1 // sigma_x = 2y -1
        sigma(0,1) = 0;
    }
};

#endif //ERRORESTIMATION_FORCINGFUNCTION_H
