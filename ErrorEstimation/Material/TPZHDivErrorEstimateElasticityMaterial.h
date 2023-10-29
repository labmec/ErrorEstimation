/* 
 * File:   TPZHDivErrorEstimateElasticityMaterial.h
 * Author: quinelato
 *
 * Created on September 6, 2023, 10:50 PM
 */

#ifndef TPZHDIVERRORESTIMATEELASTICITYMATERIAL_H
#define TPZHDIVERRORESTIMATEELASTICITYMATERIAL_H

#include "TPZHDivErrorEstimateMaterial.h"
#include "Elasticity/TPZMixedElasticityND.h"


class TPZHDivErrorEstimateElasticityMaterial : public TPZHDivErrorEstimateMaterial<TPZMixedElasticityND>{
public:
    TPZHDivErrorEstimateElasticityMaterial();
    virtual ~TPZHDivErrorEstimateElasticityMaterial();
    
    TPZHDivErrorEstimateElasticityMaterial(const TPZHDivErrorEstimateElasticityMaterial &copy): TPZHDivErrorEstimateMaterial<TPZMixedElasticityND>(copy){
    }
    
    TPZHDivErrorEstimateElasticityMaterial(const TPZHDivErrorEstimateMaterial<TPZMixedElasticityND> &copy): TPZHDivErrorEstimateMaterial<TPZMixedElasticityND>(copy){
        
    }
    
    TPZHDivErrorEstimateElasticityMaterial(const TPZMixedElasticityND &copy): TPZHDivErrorEstimateMaterial<TPZMixedElasticityND>(copy){
        
    }
        
    TPZHDivErrorEstimateElasticityMaterial &operator=(const TPZHDivErrorEstimateElasticityMaterial &copy) {
        TPZHDivErrorEstimateMaterial<TPZMixedElasticityND>::operator=(copy); 
        return *this;
    }



    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef) override;

    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
    
    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;
    
    virtual void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) override;
    virtual int VariableIndex(const std::string &name) const override;

private:


};

#endif /* TPZHDIVERRORESTIMATEELASTICITYMATERIAL_H */

