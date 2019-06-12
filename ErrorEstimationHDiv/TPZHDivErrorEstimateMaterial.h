//
//  TPZHDivErrorEstimateMaterial.hpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 11/06/19.
//

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

};

#endif /* TPZHDivErrorEstimateMaterial_hpp */
