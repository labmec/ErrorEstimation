//
//  TPZMixedErrorEstimate.hpp
//  ErrorEstimate
//
//  Created by Philippe Devloo on 20/04/18.
//

#ifndef TPZMixedErrorEstimate_hpp
#define TPZMixedErrorEstimate_hpp

#include <cstdio>
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <TPZMaterialDataT.h>
#include "pzreal.h"


class TPZMixedHDivErrorEstimate : public virtual TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
        TPZMatErrorCombinedSpaces<STATE>, TPZDarcyFlowInterface>, public TPZMixedDarcyFlow {

public:
    
    TPZMixedHDivErrorEstimate();

    TPZMixedHDivErrorEstimate(int matid, int dim);
    
    ~TPZMixedHDivErrorEstimate() override;
    
    explicit TPZMixedHDivErrorEstimate(const TPZMixedDarcyFlow &cp);
    
    TPZMixedHDivErrorEstimate(const TPZMixedHDivErrorEstimate &cp);
    
    TPZMixedHDivErrorEstimate &operator=(const TPZMixedHDivErrorEstimate &copy);
    
    [[nodiscard]] TPZMaterial * NewMaterial() const override {
        return new TPZMixedHDivErrorEstimate(*this);
    }
    
    
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    
    /// make a contribution to the error computation
    void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;

    int NEvalErrors() override { return 5; }
    
    [[nodiscard]] int VariableIndex(const std::string &name) const override;
    
    [[nodiscard]] int NSolutionVariables(int var) const override;
    
    /**
     * @brief It return a solution to multiphysics simulation.
     * @param datavec [in] Data material vector
     * @param var [in] number of solution variables. See  NSolutionVariables() method
     * @param Solout [out] is the solution vector
     */
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) override;

    /// whether the post processing mesh will be H(div) or H1
    bool fPostProcesswithHDiv = true;
};


#endif /* TPZMixedErrorEstimate_hpp */
