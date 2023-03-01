//
// Created by victor on 16/02/23.
//

#ifndef ERRORESTIMATION_TPZHYBRIDH1CREATEHDIVRECONSTRUCTION_H
#define ERRORESTIMATION_TPZHYBRIDH1CREATEHDIVRECONSTRUCTION_H

#include "TPZMultiphysicsCompMesh.h"
#include "TPZHybridH1ErrorEstimator.h"

class TPZHybridH1ErrorEstimator;
class TPZAnalysis;

class TPZHybridH1CreateHDivReconstruction {
public:
    TPZHybridH1CreateHDivReconstruction() {
        DebugStop();
    }

    TPZHybridH1CreateHDivReconstruction(TPZHybridH1ErrorEstimator *pEstimator) : fHybridH1EE(pEstimator){
         fOriginal = pEstimator->fOriginal;
     };

    ~TPZHybridH1CreateHDivReconstruction();

    TPZMultiphysicsCompMesh *CreateFluxReconstructionMesh();

    TPZCompMesh *CreateFluxReconstructionL2Mesh();

    TPZCompMesh *CreateFluxReconstructionConstantMesh();

    // a method for generating the HDiv mesh
    TPZCompMesh *CreateFluxReconstructionHDivMesh();

public: // redundant description separates a subset of methods by functionality

    // Checks if lagrange coefficients are continuous
    virtual void VerifyBoundaryFluxConsistency(TPZCompMesh* cmesh);

    /// computing the element stifnesses will "automatically" compute the condensed form of the matrices
    void ComputeElementStiffnesses(TPZCompMesh &cmesh);

    void PrintSolutionVTK(TPZAnalysis &an);

    void PostProcess(TPZMultiphysicsCompMesh *postProcMesh);

private:
    TPZHybridH1ErrorEstimator *fHybridH1EE;

    TPZMultiphysicsCompMesh *fOriginal;
};

#endif // ERRORESTIMATION_TPZHYBRIDH1CREATEHDIVRECONSTRUCTION_H
