//
//  TPZMixedErrorEstimate.hpp
//  ErrorEstimate
//
//  Created by Philippe Devloo on 20/04/18.
//

#ifndef TPZMixedErrorEstimate_hpp
#define TPZMixedErrorEstimate_hpp

#include <stdio.h>
#include "pzreal.h"

template<class TVar>
class TPZFMatrix;

class TPZMaterial;
class TPZMaterialData;

template<class TVar>
class TPZMaterialDataT;

template<class TVar>
class TPZVec;

template<class MixedMat>
class TPZMixedHDivErrorEstimate : public MixedMat
{
    
    
public:
    
    TPZMixedHDivErrorEstimate();

    TPZMixedHDivErrorEstimate(int matid, int dim);
    
    virtual ~TPZMixedHDivErrorEstimate();
    
    TPZMixedHDivErrorEstimate(const MixedMat &cp);
    
    TPZMixedHDivErrorEstimate(const TPZMixedHDivErrorEstimate<MixedMat> &cp);
    
    TPZMixedHDivErrorEstimate &operator=(const TPZMixedHDivErrorEstimate &copy);
    
    virtual TPZMaterial * NewMaterial() const override {
        return new TPZMixedHDivErrorEstimate(*this);
    }
    
    
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;
    
    /// make a contribution to the error computation
    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;
    
    virtual int NEvalErrors() const override {
        return 5;
        
    }
    
    virtual int VariableIndex(const std::string &name) const override;
    
    virtual int NSolutionVariables(int var) const override;
    
    /**
     * @brief It return a solution to multiphysics simulation.
     * @param datavec [in] Data material vector
     * @param var [in] number of solution variables. See  NSolutionVariables() method
     * @param Solout [out] is the solution vector
     */
    virtual void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) override;

    /// whether the post processing mesh will be H(div) or H1
    bool fPostProcesswithHDiv = true;
};


#endif /* TPZMixedErrorEstimate_hpp */
