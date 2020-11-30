//
//  TPZSteklovMaterial.hpp
//  SteklovProblem
//
//  Created by Denise De Siqueira on 27/11/20.
//
/**
 This material implement the variational form for Steklov problem: find (sigma, u, lambda) such that
 sigma = -K grad u
 div(sigma)=0
 sigma.n = lamda u
 
 |A B^T||sigma|=Inverse(lambda)|C 0||sigma|
 |B 0  ||u    |                |0 0||u    |
 
 
 **/

#ifndef TPZSteklovMaterial_hpp
#define TPZSteklovMaterial_hpp

#include <stdio.h>
#include "mixedpoisson.h"

class TPZSteklovMaterial : public TPZMixedPoisson {
    
public:
    TPZSteklovMaterial(int matid, int dim);
    
    TPZSteklovMaterial();
    
    TPZSteklovMaterial(const TPZSteklovMaterial &copy);
    
    TPZSteklovMaterial(const TPZMixedPoisson &copy);
    
    virtual ~TPZSteklovMaterial();
    
    TPZSteklovMaterial &operator=(const TPZSteklovMaterial &copy);
    
 
    
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef) override;
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
//
//    virtual void FillDataRequirements(TPZVec<TPZMaterialData> &datavec) override;
//    virtual void FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec) override;
//
//    virtual void UpdateBCValues(TPZVec<TPZMaterialData> &datavec);
//
//
//    virtual int NEvalErrors() override { return 5; } // erro de oscilacao de dados tbem
//
//    /// Compute the error and error estimate
//    // error[0] - error computed with exact pressure
//    // error[1] - error computed with reconstructed pressure
//    // error[2] - energy error computed with exact solution
//    // error[3] - energy error computed with reconstructed solution
//    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact,
//                        TPZVec<REAL> &errors) override;
//
//    void ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact,
//                  TPZVec<REAL> &errors, TPZBndCond &bc) override;
//
//    virtual int VariableIndex(const std::string &name) override;
//    virtual int NSolutionVariables(int var) override;
//    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override;
    
};

#endif /* TPZSteklovMaterial_hpp */
