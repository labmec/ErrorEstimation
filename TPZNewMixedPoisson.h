//
//  TPZNewMixedPoisson.hpp
//  Tools
//
//  Created by Denise De Siqueira on 18/05/20.
//

#ifndef TPZNewMixedPoisson_hpp
#define TPZNewMixedPoisson_hpp

#include <stdio.h>
#include "mixedpoisson.h"


class TPZNewMixedPoisson : public TPZMixedPoisson
{

public:
    
    
    TPZNewMixedPoisson(int matid, int dim);
    
    TPZNewMixedPoisson();
    
    TPZNewMixedPoisson(const TPZNewMixedPoisson &copy);
    
    TPZNewMixedPoisson(const TPZMixedPoisson &copy);
    
    virtual ~TPZNewMixedPoisson();
    
    TPZNewMixedPoisson &operator=(const TPZNewMixedPoisson &copy);
    

    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
     virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
    
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override;
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec) override;
    
    
    
    
};



#endif /* TPZNewMixedPoisson_hpp */


