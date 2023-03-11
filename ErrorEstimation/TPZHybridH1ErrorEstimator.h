//
//  TPZHDivErrorEstimator.hpp
//  ErrorEstimateHDiv
//
//  Created by Philippe Devloo on 10/06/18.
//

#ifndef TPZHybridH1ErrorEstimator_hpp
#define TPZHybridH1ErrorEstimator_hpp

#include "ProblemConfig.h"
#include "TPZAnalysis.h"
#include "TPZAnalyticSolution.h"
#include "TPZCreateMultiphysicsSpace.h"
#include "TPZHybridizeHDiv.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzmanvector.h"
#include "TPZHybridH1ReconstructionBase.h"
#include <stdio.h>


class TPZCompMesh;
class TPZSubCompMesh;

/// this class will compute the estimated value of the energy error of the input mesh
//  The error is estimated from a H1-conform and HDiv conform meshes
//  These meshes must be provided.
class TPZHybridH1ErrorEstimator : public TPZHybridH1ReconstructionBase {

public:

    TPZHybridH1ErrorEstimator() {
        DebugStop();
    }

    TPZHybridH1ErrorEstimator(EstimatorConfig *pEstimator) : TPZHybridH1ReconstructionBase(pEstimator){
       
        fFolderOutput = "ErrorEstimate/"; 

        InitializeFolderOutput();
     };

    
    TPZHybridH1ErrorEstimator(const TPZHybridH1ErrorEstimator &copy) :  TPZHybridH1ReconstructionBase(copy) 
    {
        // this method wont work because multiphysics meshes have no copy constructor (yet)
        DebugStop();
    }
    
    TPZHybridH1ErrorEstimator &operator=(const TPZHybridH1ErrorEstimator &cp)
    {
        // this method wont work because multiphysics meshes have no operator= (yet)
        DebugStop();
        return *this;
    }
    
    ~TPZHybridH1ErrorEstimator();

    /// create graphical output of estimated and true errors using the analysis
    void PostProcessing(TPZAnalysis &an);

    // Compute approximation error and generate VTK outputs
    void PostProcess() override;

    // Plots State solution of elements of target dimension
    void PlotState(const std::string& filename, int targetDim, TPZCompMesh* cmesh);

    /// create the post processed multiphysics mesh (which is necessarily hybridized)
    virtual void CreatePostProcessingMesh();

    inline void SetH1conformMesh(TPZCompMesh*cmesh){
        fH1conformMesh = cmesh;
    }

    inline void SetHDivConformMesh(TPZCompMesh* cmesh){
        fHDivconformMesh = cmesh;
    }

protected:

    // a method for generating the HDiv mesh
    virtual TPZCompMesh *ForceProjectionMesh();

    /// return a pointer to the pressure mesh
    virtual TPZCompMesh *PressureMesh();
    
    /// increase the side orders of the post processing mesh
    static void IncreaseSideOrders(TPZCompMesh *fluxmesh);

    /// Insert material for HDiv reconstruction
    /// Switch H1 material for H1 reconstruction material
    virtual void SwitchMaterialObjects();
    
    /// compute the effectivity indices of the pressure error and flux error and store in the element solution
    void ComputeEffectivityIndices(double &globalIndex);

    /// compute the effectivity indices of the pressure error and flux error and store in the element solution
    void ComputeEffectivityIndices(TPZSubCompMesh *cmesh);

    void FillVTKoutputVariables(TPZStack<std::string> &scalnames,TPZStack<std::string> &vecnames) override;

private:

    TPZCompMesh* fH1conformMesh = NULL;

    TPZCompMesh* fHDivconformMesh = NULL;

};

#endif /* TPZHybridH1ErrorEstimator_hpp */
