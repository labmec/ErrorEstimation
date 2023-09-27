/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZDarcyMHMHDivErrorEstimator.h
 * Author: quinelato
 *
 * Created on September 27, 2023, 11:47 AM
 */

#ifndef TPZDARCYMHMHDIVERRORESTIMATOR_H
#define TPZDARCYMHMHDIVERRORESTIMATOR_H

#include "TPZMHMHDivErrorEstimator.h"
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include "Mesh/TPZMultiphysicsCompMesh.h"
#include "Pre/TPZMHMixedMeshControl.h"


class TPZDarcyMHMHDivErrorEstimator : public TPZMHMHDivErrorEstimator<TPZMixedDarcyFlow> {
public:
    TPZDarcyMHMHDivErrorEstimator(TPZMultiphysicsCompMesh &originalMesh, TPZMHMixedMeshControl *mhm, bool postProcWithHDiv = false);
    TPZDarcyMHMHDivErrorEstimator(const TPZDarcyMHMHDivErrorEstimator& orig)=delete;
    virtual void ComputePrimalWeights();
    virtual ~TPZDarcyMHMHDivErrorEstimator();
private:

};

#endif /* TPZDARCYMHMHDIVERRORESTIMATOR_H */

