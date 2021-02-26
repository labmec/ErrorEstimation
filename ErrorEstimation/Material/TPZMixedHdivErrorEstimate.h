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
class TPZVec;

template<class MixedMat>
class TPZMixedHDivErrorEstimate : public MixedMat
{
    
    
public:
    
    TPZMixedHDivErrorEstimate();

    TPZMixedHDivErrorEstimate(int matid, int dim);
    
    virtual ~TPZMixedHDivErrorEstimate();
    
    TPZMixedHDivErrorEstimate(const MixedMat &cp);
    
    TPZMixedHDivErrorEstimate(const TPZMixedHDivErrorEstimate &cp);
    
    TPZMixedHDivErrorEstimate &operator=(const TPZMixedHDivErrorEstimate &copy);
    
    virtual TPZMaterial * NewMaterial() override {
        return new TPZMixedHDivErrorEstimate(*this);
    }
    
    
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override;
    
    /// make a contribution to the error computation
    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors) override;
    
    virtual int NEvalErrors() override {
        return 5;
        
    }
    
    virtual int VariableIndex(const std::string &name) override;
    
    virtual int NSolutionVariables(int var) override;
    
    /**
     * @brief It return a solution to multiphysics simulation.
     * @param datavec [in] Data material vector
     * @param var [in] number of solution variables. See  NSolutionVariables() method
     * @param Solout [out] is the solution vector
     */
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override;

    /// whether the post processing mesh will be H(div) or H1
    bool fPostProcesswithHDiv = true;
};


#endif /* TPZMixedErrorEstimate_hpp */
