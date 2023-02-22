//
// Created by victor on 16/02/23.
//

#ifndef ERRORESTIMATION_TPZHYBRIDH1CREATERECMESHES_H
#define ERRORESTIMATION_TPZHYBRIDH1CREATERECMESHES_H

#include "TPZMultiphysicsCompMesh.h"

class TPZHybridH1ErrorEstimator;
class TPZAnalysis;

class TPZHybridH1CreateRecMeshes {
public:

    TPZHybridH1CreateRecMeshes() {
        DebugStop();
    }

    TPZHybridH1CreateRecMeshes(TPZHybridH1ErrorEstimator *pEstimator) : fHybridH1EE(pEstimator){};

    ~TPZHybridH1CreateRecMeshes();

    TPZMultiphysicsCompMesh *CreateFluxReconstructionMesh();

    TPZCompMesh *CreateFluxReconstructionL2Mesh();

    TPZCompMesh *CreateFluxReconstructionConstantMesh();

    // a method for generating the HDiv mesh
    TPZCompMesh *CreateFluxReconstructionHDivMesh();

public: // redundant description separates a subset of methods by functionality

    // Checks if lagrange coefficients are continuous
    virtual void VerifyBoundaryFluxConsistency(TPZCompMesh* cmesh);

    /// computing the element stifnesses will "automatically" compute the condensed form of the matrices
    void ComputeElementStiffnesses();

    /// computing the element stifnesses will "automatically" compute the condensed form of the matrices
    void ComputeElementStiffnesses(TPZCompMesh &cmesh);

    void PrintSolutionVTK(TPZAnalysis &an);

    void PostProcess(TPZMultiphysicsCompMesh *postProcMesh);

private:
    TPZHybridH1ErrorEstimator *fHybridH1EE;
};

#endif // ERRORESTIMATION_TPZHYBRIDH1CREATERECMESHES_H
