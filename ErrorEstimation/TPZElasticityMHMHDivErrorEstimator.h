/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZElasticityMHMHDivErrorEstimator.h
 * Author: quinelato
 *
 * Created on September 27, 2023, 12:00 PM
 */

#ifndef TPZELASTICITYMHMHDIVERRORESTIMATOR_H
#define TPZELASTICITYMHMHDIVERRORESTIMATOR_H

#include "TPZMHMHDivErrorEstimator.h"
#include "Elasticity/TPZMixedElasticityND.h"
#include "Mesh/TPZMultiphysicsCompMesh.h"
#include "Pre/TPZMHMixedMeshControl.h"

class TPZElasticityMHMHDivErrorEstimator : public TPZMHMHDivErrorEstimator<TPZMixedElasticityND>{
public:
    TPZElasticityMHMHDivErrorEstimator(TPZMultiphysicsCompMesh &originalMesh, TPZMHMixedMeshControl *mhm, bool postProcWithHDiv = false);
    TPZElasticityMHMHDivErrorEstimator(const TPZElasticityMHMHDivErrorEstimator& orig)=delete;
    virtual ~TPZElasticityMHMHDivErrorEstimator();
    
    virtual void ComputePrimalWeights() override;
    virtual  void PostProcessing(TPZAnalysis &an, std::string &out) override;
    virtual void PlotState(const std::string& filename, int targetDim, TPZCompMesh* cmesh, bool atomic);
    
private:

};

#endif /* TPZELASTICITYMHMHDIVERRORESTIMATOR_H */

