//
// Created by Gustavo on 30/05/19.
//

#ifndef TPZAPOSTERIORIERRORESTIMATOR_H
#define TPZAPOSTERIORIERRORESTIMATOR_H

#include "pzanalysis.h"
#include "TPZAnalyticSolution.h"

class TPZAPosterioriErrorEstimator {

private:
    TPZVec<REAL> fElementErrors;
    TPZAnalysis *fAnalysis;
    TPZAnalyticSolution *fAnalyticSolution;
    
    TPZCompMesh *originalPressureMesh;
    TPZCompMesh *originalFluxMesh;
    TPZCompMesh *reconstructedPressureMesh;
    TPZCompMesh *reconstructedFluxMesh;

public:
    virtual void CreateAuxiliarySpaces() = 0;
    virtual void ComputeErrors() = 0;

    void SetAnalyticSolution(TPZAnalyticSolution *analyticSolution);
    void CalculateEffectivityIndex();
    void PostProcess();
};


#endif // TPZAPOSTERIORIERRORESTIMATOR_H
