//
// Created by victor on 16/02/23.
//

#ifndef ERRORESTIMATION_TPZHYBRIDH1CREATEHDIVRECONSTRUCTION_H
#define ERRORESTIMATION_TPZHYBRIDH1CREATEHDIVRECONSTRUCTION_H

#include "TPZMultiphysicsCompMesh.h"
#include "TPZHybridH1ErrorEstimator.h"
#include "ProblemConfig.h"
#include "TPZHybridH1ReconstructionBase.h"

class TPZHybridH1ErrorEstimator;
class TPZAnalysis;

class TPZHybridH1CreateHDivReconstruction : public TPZHybridH1ReconstructionBase {
public:
    TPZHybridH1CreateHDivReconstruction() {
        DebugStop();
    }

    TPZHybridH1CreateHDivReconstruction(EstimatorConfig *pEstimator) : TPZHybridH1ReconstructionBase(pEstimator){
       
       fFolderOutput = "HDivRecDebug/";

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

    void PostProcess(TPZMultiphysicsCompMesh *postProcMesh);

protected:

    void FillVTKoutputVariables(TPZStack<std::string> &scalnames,TPZStack<std::string> &vecnames) override;

private:

    int fLagrangeMatId;
};

#endif // ERRORESTIMATION_TPZHYBRIDH1CREATEHDIVRECONSTRUCTION_H
