//
//  TPZMixedErrorEstimate.hpp
//  ErrorEstimate
//
//  Created by Philippe Devloo on 20/04/18.
//

#ifndef TPZMixedErrorEstimate2_hpp
#define TPZMixedErrorEstimate2_hpp

#include <stdio.h>
#include "pzreal.h"
#include "Material/DarcyFlow/TPZMixedDarcyFlow.h"

template<class TVar>
class TPZFMatrix;

class TPZMaterial;
class TPZMaterialData;

template<class TVar>
class TPZVec;

/// class to compute the contribution of an HDiv reconstruction
template<class MixedMat>
class TPZMixedErrorEstimate2 : public MixedMat
{
 
    enum MMeshPositions {Eflux = 0, Epressure = 1, Epatch = 2, Eorigin = 3};
    /// sign convention adopted by MixedMat
    int fSignConvention = 1;
    
public:
    
    TPZMixedErrorEstimate2 ();
    
    TPZMixedErrorEstimate2 (int matid, int dim);
    
    virtual ~TPZMixedErrorEstimate2 ();
    
    TPZMixedErrorEstimate2 (const TPZMixedErrorEstimate2  &cp);
    
    TPZMixedErrorEstimate2  &operator=(const TPZMixedErrorEstimate2  &copy);
    
    virtual TPZMaterial * NewMaterial(){
        return new TPZMixedErrorEstimate2 (*this);
    }
    
    int SignConvention() const
    {
        return fSignConvention;
    }
    
    void SetSignConvention(int sign)
    {
        fSignConvention = sign;
    }
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec)const override;

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE> > &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    /// make a contribution to the error computation
    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE> > &data, TPZVec<REAL> &errors) override;
    

};


#endif /* TPZMixedErrorEstimate2 _hpp */
