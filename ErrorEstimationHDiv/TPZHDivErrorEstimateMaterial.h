//
//  TPZHDivErrorEstimateMaterial.hpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 11/06/19.
//
/**
 This material implement the Ainsworth proposal
 **/

#ifndef TPZHDivErrorEstimateMaterial_hpp
#define TPZHDivErrorEstimateMaterial_hpp

#include <stdio.h>
#include "TPZMatLaplacian.h"
#include "mixedpoisson.h"

class TPZHDivErrorEstimateMaterial : public TPZMixedPoisson
{

public:
    
    TPZHDivErrorEstimateMaterial(int matid, int dim);
    
    TPZHDivErrorEstimateMaterial();
    
    TPZHDivErrorEstimateMaterial(const TPZHDivErrorEstimateMaterial &copy);
    
    TPZHDivErrorEstimateMaterial(const TPZMixedPoisson &copy);
    
    virtual ~TPZHDivErrorEstimateMaterial();
    
    TPZHDivErrorEstimateMaterial &operator=(const TPZHDivErrorEstimateMaterial &copy);
    

    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override;
    
    bool fNeumannLocalProblem = true;

};

#endif /* TPZHDivErrorEstimateMaterial_hpp */
