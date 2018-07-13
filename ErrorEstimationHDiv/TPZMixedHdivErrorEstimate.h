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
    
    /// sign convention adopted by MixedMat
    int fSignConvention = 1;
    
public:
    
    TPZMixedHDivErrorEstimate();
    
    TPZMixedHDivErrorEstimate(int matid, int dim);
    
    virtual ~TPZMixedHDivErrorEstimate();
    
    TPZMixedHDivErrorEstimate(const TPZMixedHDivErrorEstimate &cp);
    
    TPZMixedHDivErrorEstimate &operator=(const TPZMixedHDivErrorEstimate &copy);
    
    virtual TPZMaterial * NewMaterial(){
        return new TPZMixedHDivErrorEstimate(*this);
    }
    
    int SignConvention() const
    {
        return fSignConvention;
    }
    
    void SetSignConvention(int sign)
    {
        fSignConvention = sign;
    }
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    /// make a contribution to the error computation
    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors);
    
    virtual int NEvalErrors() {
        return 4;
        
    }
    
    virtual int VariableIndex(const std::string &name);
    
    virtual int NSolutionVariables(int var);
    
    /**
     * @brief It return a solution to multiphysics simulation.
     * @param datavec [in] Data material vector
     * @param var [in] number of solution variables. See  NSolutionVariables() method
     * @param Solout [out] is the solution vector
     */
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    

};


#endif /* TPZMixedErrorEstimate_hpp */
