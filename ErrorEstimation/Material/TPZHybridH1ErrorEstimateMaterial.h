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
#include "DarcyFlow/TPZMixedDarcyFlow.h"

typedef TPZMixedDarcyFlow TPZMixedPoisson;

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

    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc) override;


    virtual void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;
    virtual void FillBoundaryConditionDataRequirements(int type,TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    virtual TPZMaterial * NewMaterial() const override{
        return new TPZHybridH1ErrorEstimateMaterial(*this);
    }

    bool fNeumannLocalProblem = true;

    virtual int NEvalErrors() const override {return 7;}

    /// Compute the error and error estimate
    // error[0] - error computed with exact pressure
    // error[1] - error computed with reconstructed pressure
    // error[2] - energy error computed with exact solution
    // error[3] - energy error computed with reconstructed solution
    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data,  TPZVec<REAL> &errors) override;

    void ErrorsBC(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors,TPZBndCondT<STATE> &bc);


    virtual int VariableIndex(const std::string &name) const override;
    virtual int NSolutionVariables(int var) const override;
    virtual void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var,
                          TPZVec<STATE> &Solout) override;
};
#endif //ERRORESTIMATION_TPZEEMatHybridH1ToH1_H
