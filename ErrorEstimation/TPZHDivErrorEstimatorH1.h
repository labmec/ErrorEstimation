//
// Created by gustavo on 30/05/19.
//
#ifndef ERRORESTIMATION_TPZHDIVERRORESTIMATORH1_H
#define ERRORESTIMATION_TPZHDIVERRORESTIMATORH1_H


#include "TPZHybridHDivErrorEstimator.h"



class TPZHDivErrorEstimatorH1 : public TPZHybridHDivErrorEstimator {

public:
    /// number of orders the pressure polynomial order is increase
    int fUpliftOrder = 0;
    //for solve a neumann problem
    bool fperformUplift = true;
public:
    
    TPZHDivErrorEstimatorH1(TPZMultiphysicsCompMesh &InputMesh) : TPZHybridHDivErrorEstimator(InputMesh)
    {
        
    }
    /// switch material object from mixed poisson to TPZMixedHdivErrorEstimate
    virtual void SwitchMaterialObjects() override;    

    /// create the post processed multiphysics mesh (which is necessarily hybridized)
    virtual void CreatePostProcessingMesh() override;
    
    /// Increase the polynomial orders of the pressures
    void IncreasePressureOrders(TPZCompMesh *pressuremesh);
    
    /// Compute an uplifted solution for the pressure
    void UpliftPressure();
    
    /// return a pointer to the pressure mesh
    virtual TPZCompMesh *PressureMesh() override;
    

    /// create a constant pressure mesh used for uplifting the pressure
    TPZCompMesh *CreateDiscontinuousMesh(const TPZCompMesh *pressuremesh);
    
    /// prepare the elements of postprocmesh to compute the pressures with increased accuracy
    void PreparePostProcElements();
    
    /// copy the solution from the neighbouring skeleton elements to the H1 pressure elements
    void CopySolutionFromSkeleton() override;
    virtual void ComputePressureWeights() override;
    
};


#endif //ERRORESTIMATION_TPZHDIVERRORESTIMATORH1_H
