//
// Created by victor on 29/03/2021.
//

#ifndef ERRORESTIMATION_TPZHYBMIXDIFFMATERIAL_H
#define ERRORESTIMATION_TPZHYBMIXDIFFMATERIAL_H

#include "TPZMatLaplacianHybrid.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"

typedef TPZMixedDarcyFlow TPZMixedPoisson;

class TPZHybMixDiffMaterial: public TPZMixedPoisson {
    
public:

    TPZHybMixDiffMaterial(int matid, int dim);

    TPZHybMixDiffMaterial();

    TPZHybMixDiffMaterial(TPZMatLaplacianHybrid &matlaplacian);

    TPZHybMixDiffMaterial(const TPZHybMixDiffMaterial &copy);

    TPZHybMixDiffMaterial(const TPZMixedPoisson &copy);

    virtual ~TPZHybMixDiffMaterial();

    TPZHybMixDiffMaterial &operator=(const TPZHybMixDiffMaterial &copy);

    virtual void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    virtual void FillBoundaryConditionDataRequirements(int type,TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    virtual TPZMaterial * NewMaterial() const override{
        return new TPZHybMixDiffMaterial(*this);
    }

    virtual int NEvalErrors() const override {return 6;}

    /// Compute the error and error estimate
//    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors) override;

    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;
    
    virtual int VariableIndex(const std::string &name) const override;
    virtual int NSolutionVariables(int var) const override;
    virtual void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var,
                          TPZVec<STATE> &Solout) override;
};


#endif //ERRORESTIMATION_TPZHYBMIXDIFFMATERIAL_H
