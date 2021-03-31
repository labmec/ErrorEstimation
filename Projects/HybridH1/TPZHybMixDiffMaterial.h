//
// Created by victor on 29/03/2021.
//

#ifndef ERRORESTIMATION_TPZHYBMIXDIFFMATERIAL_H
#define ERRORESTIMATION_TPZHYBMIXDIFFMATERIAL_H

#include "TPZMatLaplacianHybrid.h"
#include "mixedpoisson.h"


class TPZHybMixDiffMaterial: public TPZMixedPoisson {
    
public:

    TPZHybMixDiffMaterial(int matid, int dim);

    TPZHybMixDiffMaterial();

    TPZHybMixDiffMaterial(TPZMatLaplacianHybrid &matlaplacian);

    TPZHybMixDiffMaterial(const TPZHybMixDiffMaterial &copy);

    TPZHybMixDiffMaterial(const TPZMixedPoisson &copy);

    virtual ~TPZHybMixDiffMaterial();

    TPZHybMixDiffMaterial &operator=(const TPZHybMixDiffMaterial &copy);

    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override;
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec) override;

    virtual TPZMaterial * NewMaterial() override{
        return new TPZHybMixDiffMaterial(*this);
    }

    virtual int NEvalErrors() override {return 7;}

    /// Compute the error and error estimate
    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors) override;

    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<REAL> &errors) override;
    
    virtual int VariableIndex(const std::string &name) override;
    virtual int NSolutionVariables(int var) override;
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var,
                          TPZVec<STATE> &Solout) override;
};


#endif //ERRORESTIMATION_TPZHYBMIXDIFFMATERIAL_H
