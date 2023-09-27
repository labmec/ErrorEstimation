//
//  TPZPressureProjection.hpp
//  ErrorEstimationLib
//
//  Created by Denise De Siqueira on 21/05/20.
//

#ifndef TPZPressureProjection_hpp
#define TPZPressureProjection_hpp

#include <stdio.h>

#include "TPZHDivErrorEstimateDarcyMaterial.h"

class TPZPressureProjection : public TPZHDivErrorEstimateDarcyMaterial
{

public:
    
    
    TPZPressureProjection(int matid, int dim);
    
    TPZPressureProjection();
    
    TPZPressureProjection(const TPZPressureProjection &copy);
    
    TPZPressureProjection(const TPZHDivErrorEstimateDarcyMaterial &copy);
    
    virtual ~TPZPressureProjection();
    
    TPZPressureProjection &operator=(const TPZPressureProjection &copy);
    

    
     virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc) override;
    
    
};





#endif /* TPZPressureProjection_hpp */
