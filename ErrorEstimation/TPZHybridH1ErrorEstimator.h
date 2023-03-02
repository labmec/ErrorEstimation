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
#include <stdio.h>
// #include "../Projects/Tools/LoadCases.h"


class TPZCompMesh;
class TPZSubCompMesh;

/// this class will compute the estimated value of the energy error of the input mesh
// the class should work for any HybridSquared mesh
// first create the post processing mesh
// then compute the errors
struct TPZHybridH1ErrorEstimator
{
    /// The HybridSquared approximation mesh for which we will compute the error
    TPZMultiphysicsCompMesh *fOriginal;

    /// Computational mesh with pressure and flux reconstructions
    TPZMultiphysicsCompMesh fPostProcMesh;

    // material id of the dim-1 skeleton elements
    int fPressureSkeletonMatId;

    int fLagrangeMatId = -999;
    
    TPZAnalyticSolution *fExact;
    
    ProblemConfig fProblemConfig;

    std::string fDebugDirName = "HybridH1_ReconstructionDebug";

    TPZHybridH1ErrorEstimator(TPZMultiphysicsCompMesh &InputMesh) : fOriginal(&InputMesh),
    fPostProcMesh(0),fExact(NULL)
    {
        FindFreeMatID(fPressureSkeletonMatId);
    }

    TPZHybridH1ErrorEstimator(TPZMultiphysicsCompMesh &InputMesh, int skeletonMatId, int HDivMatId) : fOriginal(&InputMesh),
                                                                    fPostProcMesh(0),fExact(NULL),
                                                                    fPressureSkeletonMatId(fPressureSkeletonMatId)
    {

    }
    
    TPZHybridH1ErrorEstimator(const TPZHybridH1ErrorEstimator &copy) : fOriginal(copy.fOriginal),
        fPostProcMesh(copy.fPostProcMesh), fExact(copy.fExact), fProblemConfig(copy.fProblemConfig),fPressureSkeletonMatId(copy.fPressureSkeletonMatId)
    {
        // this method wont work because multiphysics meshes have no copy constructor (yet)
        DebugStop();
    }
    
    TPZHybridH1ErrorEstimator &operator=(const TPZHybridH1ErrorEstimator &cp)
    {
        fOriginal = cp.fOriginal;
        // this method wont work because multiphysics meshes have no operator= (yet)
        DebugStop();

        fPostProcMesh = cp.fPostProcMesh;
        fExact = cp.fExact;
        fProblemConfig = cp.fProblemConfig;
        fPressureSkeletonMatId = cp.fPressureSkeletonMatId;
        return *this;
    }
    
    ~TPZHybridH1ErrorEstimator();

    void SetLagrangeMatID(int lagrangeID){
        fLagrangeMatId = lagrangeID;
    }
    
    /// Set the analytic solution object
    void SetAnalyticSolution(TPZAnalyticSolution &exact)
    {
        fExact = &exact;
    }

    /// compute the element errors comparing the reconstructed solution based on average pressures
    /// with the original solution
    virtual void ComputeErrors(TPZVec<REAL> &errorVec, TPZVec<REAL> &elementerrors, bool store);

    /// Create pos processing mesh, call H1 and HDiv reconstructions
    void CreateReconstructionSpaces();

    /// create graphical output of estimated and true errors using the analysis
    void PostProcessing(TPZAnalysis &an);

    // Plots State solution of elements of target dimension
    void PlotState(const std::string& filename, int targetDim, TPZCompMesh* cmesh);

protected:

    /// create the post processed multiphysics mesh (which is necessarily hybridized)
    virtual void CreatePostProcessingMesh();

    // a method for generating the HDiv mesh
    virtual TPZCompMesh *ForceProjectionMesh();

    /// return a pointer to the pressure mesh
    virtual TPZCompMesh *PressureMesh();
    
    /// increase the side orders of the post processing mesh
    static void IncreaseSideOrders(TPZCompMesh *fluxmesh);

    /// Insert material for HDiv reconstruction
    /// Switch H1 material for H1 reconstruction material
    virtual void SwitchMaterialObjects();

    /// Find free matID number
    void FindFreeMatID(int &matID);
    
    /// compute the effectivity indices of the pressure error and flux error and store in the element solution
    void ComputeEffectivityIndices(double &globalIndex);

    /// compute the effectivity indices of the pressure error and flux error and store in the element solution
    void ComputeEffectivityIndices(TPZSubCompMesh *cmesh);


    /// identify the peripheral material objects and store the information in fHybridizer
    void IdentifyPeripheralMaterialIds();



    friend class TPZHybridH1CreateH1Reconstruction;
    friend class TPZHybridH1CreateHDivReconstruction;

};

#endif /* TPZHybridH1ErrorEstimator_hpp */
