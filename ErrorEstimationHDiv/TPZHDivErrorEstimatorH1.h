//
// Created by gustavo on 30/05/19.
//

#include "TPZHybridHDivErrorEstimator.h"

#ifndef ERRORESTIMATION_TPZHDIVERRORESTIMATORH1_H
#define ERRORESTIMATION_TPZHDIVERRORESTIMATORH1_H


class TPZHDivErrorEstimatorH1 : public TPZHybridHDivErrorEstimator {

    /// number of orders the pressure polynomial order is increase
    int fUpliftOrder = 0;
protected:
    
    /// create the post processed multiphysics mesh (which is necessarily hybridized)
    virtual void CreatePostProcessingMesh() override;
    
    /// compute a more precise approximation for the pressure
    void UpliftPressure();
};


#endif //ERRORESTIMATION_TPZHDIVERRORESTIMATORH1_H
