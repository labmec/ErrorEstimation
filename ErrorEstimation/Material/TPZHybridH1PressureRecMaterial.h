//
// Created by victor on 14/02/23.
//

#ifndef ERRORESTIMATION_TPZHYBRIDH1PRESSURERECMATERIAL_H
#define ERRORESTIMATION_TPZHYBRIDH1PRESSURERECMATERIAL_H

#include <stdio.h>
#include "TPZMatLaplacianHybrid.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"

typedef TPZMixedDarcyFlow TPZMixedPoisson;

class TPZHybridH1PressureRecMaterial : public TPZMixedPoisson
{
public:
    TPZHybridH1PressureRecMaterial(int matid, int dim);

    TPZHybridH1PressureRecMaterial();

    TPZHybridH1PressureRecMaterial(TPZMatLaplacianHybrid &matlaplacian);

    TPZHybridH1PressureRecMaterial(const TPZHybridH1PressureRecMaterial &copy);

    TPZHybridH1PressureRecMaterial(const TPZMixedPoisson &copy);

    virtual ~TPZHybridH1PressureRecMaterial();

    TPZHybridH1PressureRecMaterial &operator=(const TPZHybridH1PressureRecMaterial &copy);

    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc) override;


    virtual void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    virtual void FillBoundaryConditionDataRequirements(int type,TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    virtual TPZMaterial * NewMaterial() const override{
        return new TPZHybridH1PressureRecMaterial(*this);
    }

    bool fNeumannLocalProblem = true;

    virtual int NEvalErrors() const override {return 7;}

    /// Compute the error and error estimate
    // error[0] - error computed with exact pressure
    // error[1] - error computed with reconstructed pressure
    // error[2] - energy error computed with exact solution
    // error[3] - energy error computed with reconstructed solution
    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data,  TPZVec<REAL> &errors) override;

    virtual int VariableIndex(const std::string &name) const override;

    virtual int NSolutionVariables(int var) const override;

    virtual void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var,
                          TPZVec<STATE> &Solout) override;
};

#endif // ERRORESTIMATION_TPZHYBRIDH1PRESSURERECMATERIAL_H
