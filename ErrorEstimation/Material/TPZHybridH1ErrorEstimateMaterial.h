//
//  TPZHybridH1ErrorEstimateMaterial.hpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 11/06/19.
//


#ifndef ERRORESTIMATION_TPZEEMatHybridH1ToH1_H
#define ERRORESTIMATION_TPZEEMatHybridH1ToH1_H

#include <stdio.h>
#include "TPZMatLaplacianHybrid.h"
#include "mixedpoisson.h"

class TPZHybridH1ErrorEstimateMaterial: public TPZMixedPoisson
{
public:

    TPZHybridH1ErrorEstimateMaterial(int matid, int dim);

    TPZHybridH1ErrorEstimateMaterial();

    TPZHybridH1ErrorEstimateMaterial(TPZMatLaplacianHybrid &matlaplacian);

    TPZHybridH1ErrorEstimateMaterial(const TPZHybridH1ErrorEstimateMaterial &copy);

    TPZHybridH1ErrorEstimateMaterial(const TPZMixedPoisson &copy);

    virtual ~TPZHybridH1ErrorEstimateMaterial();

    TPZHybridH1ErrorEstimateMaterial &operator=(const TPZHybridH1ErrorEstimateMaterial &copy);

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

    virtual TPZMaterial * NewMaterial() override{
        return new TPZHybridH1ErrorEstimateMaterial(*this);
    }

    bool fNeumannLocalProblem = true;

    virtual int NEvalErrors() override {return 7;}

    /// Compute the error and error estimate
    // error[0] - error computed with exact pressure
    // error[1] - error computed with reconstructed pressure
    // error[2] - energy error computed with exact solution
    // error[3] - energy error computed with reconstructed solution
    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors) override;

    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<REAL> &errors) override;

    void ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors,TPZBndCond &bc)override;


    virtual int VariableIndex(const std::string &name) override;
    virtual int NSolutionVariables(int var) override;
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var,
                          TPZVec<STATE> &Solout) override;
};
#endif //ERRORESTIMATION_TPZEEMatHybridH1ToH1_H
