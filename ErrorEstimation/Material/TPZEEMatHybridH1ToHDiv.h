//
// Created by victor on 07/10/2020.
//

#ifndef ERRORESTIMATION_TPZEEMATHYBRIDH1TOHIV_H
#define ERRORESTIMATION_TPZEEMATHYBRIDH1TOHIV_H

#include <stdio.h>
#include "mixedpoisson.h"
#include "TPZMatLaplacianHybrid.h"

class TPZEEMatHybridH1ToHDiv: public TPZMixedPoisson
{

public:


    TPZEEMatHybridH1ToHDiv(int matid, int dim);

    TPZEEMatHybridH1ToHDiv();

    TPZEEMatHybridH1ToHDiv(const TPZEEMatHybridH1ToHDiv &copy);

    TPZEEMatHybridH1ToHDiv( TPZMixedPoisson &copy);

    TPZEEMatHybridH1ToHDiv( TPZMatLaplacianHybrid matlaplacian);

    virtual ~TPZEEMatHybridH1ToHDiv();

    TPZEEMatHybridH1ToHDiv &operator=(const TPZEEMatHybridH1ToHDiv &copy);


    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;

    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override{
        DebugStop();
    }
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override{
        DebugStop();
    }


    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override;
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec) override;

    virtual void UpdateBCValues(TPZVec<TPZMaterialData> &datavec);




    bool fNeumannLocalProblem = true;

    virtual int NEvalErrors() override {return 5;}//erro de oscilacao de dados tbem

    /// Compute the error and error estimate
    // error[0] - error computed with exact pressure
    // error[1] - error computed with reconstructed pressure
    // error[2] - energy error computed with exact solution
    // error[3] - energy error computed with reconstructed solution
    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors) override;

    void ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors,TPZBndCond &bc)override;


    virtual int VariableIndex(const std::string &name) override;
    virtual int NSolutionVariables(int var) override;
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var,
                          TPZVec<STATE> &Solout) override;

    // Returns the first non-null approximation space index, which will be the
    // H1 reconstruction space
    int FirstNonNullApproxSpaceIndex(TPZVec<TPZMaterialData> &datavec);
};

#endif //ERRORESTIMATION_TPZEEMATHYBRIDH1TOHIV_H
