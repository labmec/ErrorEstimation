/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZElasticityMHMHDivErrorEstimator.cc
 * Author: quinelato
 * 
 * Created on September 27, 2023, 12:00 PM
 */

#include "TPZElasticityMHMHDivErrorEstimator.h"
#include "Common/pzerror.h"

TPZElasticityMHMHDivErrorEstimator::TPZElasticityMHMHDivErrorEstimator(TPZMultiphysicsCompMesh &originalMesh, TPZMHMixedMeshControl *mhm, bool postProcWithHDiv) : TPZMHMHDivErrorEstimator<TPZMixedElasticityND>(originalMesh, mhm, postProcWithHDiv) {
}

TPZElasticityMHMHDivErrorEstimator::~TPZElasticityMHMHDivErrorEstimator() {
}

void TPZElasticityMHMHDivErrorEstimator::ComputePrimalWeights() {
    DebugStop();
}

