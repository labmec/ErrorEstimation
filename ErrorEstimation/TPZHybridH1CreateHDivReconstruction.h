//
// Created by victor on 16/02/23.
//

#ifndef ERRORESTIMATION_TPZHYBRIDH1CREATEHDIVRECONSTRUCTION_H
#define ERRORESTIMATION_TPZHYBRIDH1CREATEHDIVRECONSTRUCTION_H

#include "TPZMultiphysicsCompMesh.h"
#include "ProblemConfig.h"
#include "TPZHybridH1ReconstructionBase.h"

class TPZAnalysis;

class TPZHybridH1CreateHDivReconstruction : public TPZHybridH1ReconstructionBase {
public:
    TPZHybridH1CreateHDivReconstruction() {
        DebugStop();
    }

    TPZHybridH1CreateHDivReconstruction(EstimatorConfig *pEstimator) : TPZHybridH1ReconstructionBase(pEstimator){
       
       fFolderOutput = "HDivReconstruction/";

       InitializeFolderOutput();

       fLagrangeMatId = pEstimator->fLagrangeMatId;
     };

    ~TPZHybridH1CreateHDivReconstruction(){}

    TPZMultiphysicsCompMesh *CreateFluxReconstructionMesh();

    TPZCompMesh *CreateFluxReconstructionL2Mesh();

    TPZCompMesh *CreateFluxReconstructionConstantMesh();

    // a method for generating the HDiv mesh
    TPZCompMesh *CreateFluxReconstructionHDivMesh();

public: // redundant description separates a subset of methods by functionality

    // Checks if lagrange coefficients are continuous
    virtual void VerifyBoundaryFluxConsistency(TPZCompMesh* cmesh);

    void PostProcess() override;

    inline TPZCompMesh* GetReconstructionMesh(){
        return fHDivReconstructionAtomicMesh;
    }

protected:

    void FillVTKoutputVariables(TPZStack<std::string> &scalnames,TPZStack<std::string> &vecnames) override;

private:

    int fLagrangeMatId;

protected:
 
    TPZCompMesh* fHDivReconstructionAtomicMesh = NULL;
};

#endif // ERRORESTIMATION_TPZHYBRIDH1CREATEHDIVRECONSTRUCTION_H
