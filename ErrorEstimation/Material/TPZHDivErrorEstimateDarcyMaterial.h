/* 
 * File:   TPZHDivErrorEstimateDarcyMaterial.h
 * Author: quinelato
 *
 * Created on September 6, 2023, 10:50 PM
 */

#ifndef TPZHDIVERRORESTIMATEDARCYMATERIAL_H
#define TPZHDIVERRORESTIMATEDARCYMATERIAL_H

#include "TPZHDivErrorEstimateMaterial.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"

class TPZHDivErrorEstimateDarcyMaterial : public TPZHDivErrorEstimateMaterial<TPZMixedDarcyFlow> {
public:
    
    TPZHDivErrorEstimateDarcyMaterial(int matid, int dim) : TPZHDivErrorEstimateMaterial<TPZMixedDarcyFlow>(matid, dim) {
        
    }

    TPZHDivErrorEstimateDarcyMaterial() : TPZHDivErrorEstimateMaterial<TPZMixedDarcyFlow>(){
        
    }
    
    TPZHDivErrorEstimateDarcyMaterial(const TPZHDivErrorEstimateDarcyMaterial &copy): TPZHDivErrorEstimateMaterial<TPZMixedDarcyFlow>(copy){
        
    }
    
    TPZHDivErrorEstimateDarcyMaterial(const TPZHDivErrorEstimateMaterial<TPZMixedDarcyFlow> &copy): TPZHDivErrorEstimateMaterial<TPZMixedDarcyFlow>(copy){
        
    }
    
    TPZHDivErrorEstimateDarcyMaterial(const TPZMixedDarcyFlow &copy): TPZHDivErrorEstimateMaterial<TPZMixedDarcyFlow>(copy){
        
    }
        
    TPZHDivErrorEstimateDarcyMaterial &operator=(const TPZHDivErrorEstimateDarcyMaterial &copy) {
        TPZHDivErrorEstimateMaterial<TPZMixedDarcyFlow>::operator=(copy); 
        return *this;
    }
    
    virtual ~TPZHDivErrorEstimateDarcyMaterial() {
        
    }
    
    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef) override;

    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
    
    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;
    
    virtual void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) override;
private:

};

#endif /* TPZHDIVERRORESTIMATEDARCYMATERIAL_H */

