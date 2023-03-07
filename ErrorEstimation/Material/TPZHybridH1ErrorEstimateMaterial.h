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

    [[nodiscard]] std::string Name() const override { return "TPZHybridH1ErrorEstimateMaterial"; }

    TPZHybridH1ErrorEstimateMaterial &operator=(const TPZHybridH1ErrorEstimateMaterial &copy);

    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc) override;


    virtual void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    virtual void FillBoundaryConditionDataRequirements(int type,TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    virtual TPZMaterial * NewMaterial() const override{
        return new TPZHybridH1ErrorEstimateMaterial(*this);
    }

    virtual int NEvalErrors() const override {return 7;}

    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data,  TPZVec<REAL> &errors) override;

    void ErrorsBC(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors,TPZBndCondT<STATE> &bc);


    virtual int VariableIndex(const std::string &name) const override;

    virtual int NSolutionVariables(int var) const override;

    virtual void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var,
                          TPZVec<STATE> &Solout) override;

    private:

    int fHDivConformPosition      =    0;
    int fH1conformPosition        =    1;
    int fLagrangeCoeffPosition    =    2;
    int fFEMbrokenH1Position      =    3;
    int fSourceProjectionPosition =    4;
};
#endif //ERRORESTIMATION_TPZEEMatHybridH1ToH1_H
