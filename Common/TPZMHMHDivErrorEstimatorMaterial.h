//
//  TPZMHMHDivErrorEstimateMaterial.hpp
//  ErrorEstimation
//
//  Created by Denise Siqueira on 22/07/19.
//
/**
 This material implement the MHM error estimation
 **/

#ifndef TPZMHMHDivErrorEstimateMaterial_hpp
#define TPZMHMHDivErrorEstimateMaterial_hpp

#include <stdio.h>
#include "TPZMatLaplacian.h"
#include "mixedpoisson.h"


class TPZMHMHDivErrorEstimateMaterial : public TPZMixedPoisson
{

public:
    
    
    TPZMHMHDivErrorEstimateMaterial(int matid, int dim);
    
    TPZMHMHDivErrorEstimateMaterial();
    
    TPZMHMHDivErrorEstimateMaterial(const TPZMHMHDivErrorEstimateMaterial &copy);
    
    TPZMHMHDivErrorEstimateMaterial(const TPZMixedPoisson &copy);
    
    virtual ~TPZMHMHDivErrorEstimateMaterial();
    
    TPZMHMHDivErrorEstimateMaterial &operator=(const TPZMHMHDivErrorEstimateMaterial &copy);
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override;
    

    virtual int NEvalErrors() override {return 5;}//erro de oscilacao de dados tbem
    
    /// Compute the error and error estimate
    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors) override;
    
    int VariableIndex(const std::string &name)override;
    int NSolutionVariables(int var)override;
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override;
    virtual void  Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)override;
    virtual int IsH1Position(TPZVec<TPZMaterialData> &datavec);
 

    
    
    

};

#endif /* TPZMHMHDivErrorEstimateMaterial_hpp */
