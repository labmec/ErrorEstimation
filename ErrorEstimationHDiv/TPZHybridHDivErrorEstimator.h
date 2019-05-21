//
//  TPZHybridHDivErrorEstimator.hpp
//  ErrorEstimateHDiv
//
//  Created by Philippe Devloo on 10/06/18.
//

#ifndef TPZHybridHDivErrorEstimator_hpp
#define TPZHybridHDivErrorEstimator_hpp

#include <stdio.h>
#include "pzmanvector.h"
#include "TPZHybridizeHDiv.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZAnalyticSolution.h"
#include "ProblemConfig.h"
#include "pzanalysis.h"


class TPZCompMesh;

/// this class will compute the estimated value of the energy error of the input mesh
// it is expected that the estimated value of the error is improved if the input mesh is of H(div)++ type
// the class should work for any H(div) mesh
// first create the post processing mesh
// then compute the errors
struct TPZHybridHDivErrorEstimator
{
    
    /// The H(Div) approximation mesh for which we will compute the error
    TPZMultiphysicsCompMesh *fOriginal;
//    TPZManVector<TPZCompMesh *,3> fOriginal;
    
    bool fOriginalIsHybridized = true;
    
    /// order increase of the boundary flux (depending on the original mesh)
    int fUpliftPostProcessMesh = 0;
    
    /// Locally created computational mesh to compute the error
    TPZMultiphysicsCompMesh fPostProcMesh;
//    TPZManVector<TPZCompMesh *,7> fPostProcMesh;
    
    // object to assist in creating a hybridized version of the computational mesh
    TPZHybridizeHDiv fHybridizer;
    
    TPZAnalyticSolution *fExact;
    
    ProblemConfig fProblemConfig;
    
    TPZHybridHDivErrorEstimator(TPZMultiphysicsCompMesh &InputMesh, bool InputisHybridized = true) : fOriginal(&InputMesh),
    fOriginalIsHybridized(InputisHybridized), fPostProcMesh(0),fExact(NULL)
    {
        
    }
    
    TPZHybridHDivErrorEstimator(const TPZHybridHDivErrorEstimator &copy) : fOriginal(copy.fOriginal),
        fOriginalIsHybridized(copy.fOriginalIsHybridized),fUpliftPostProcessMesh(copy.fUpliftPostProcessMesh),
        fPostProcMesh(copy.fPostProcMesh), fExact(copy.fExact), fProblemConfig(copy.fProblemConfig)
    {
        
    }
    
    TPZHybridHDivErrorEstimator &operator=(const TPZHybridHDivErrorEstimator &cp)
    {
        fOriginal = cp.fOriginal;
        fOriginalIsHybridized = cp.fOriginalIsHybridized;
        fUpliftPostProcessMesh = cp.fUpliftPostProcessMesh;
        fPostProcMesh = cp.fPostProcMesh;
        fExact = cp.fExact;
        fProblemConfig = cp.fProblemConfig;
        return *this;
    }
    
    ~TPZHybridHDivErrorEstimator();
    
    /// Set the analytic solution object
    void SetAnalyticSolution(TPZAnalyticSolution &exact)
    {
        fExact = &exact;
    }
    
    /// compute the element errors comparing the reconstructed solution based on average pressures
    /// with the original solution
    void ComputeErrors(TPZVec<REAL> &elementerrors, bool store = true);
    
    //reconstruction of potential using hybrid solution on enrichement space
    void PotentialReconstruction();
    
    
    void PostProcessingHybridMesh();
    void CreateMultiphysicsHybridMesh();
    void PostProcessing(TPZAnalysis &an);
    
    
    void PlotLagrangeMultiplier(const std::string &filename, bool reconstructed = true);
private:
    
    /// create the post processed multiphysics mesh (which is necessarily hybridized)
    void CreatePostProcessingMesh();
    
    /// computing the element stifnesses will "automatically" compute the condensed form of the matrices
    void ComputeElementStiffnesses();
    
    /// increase the side orders of the post processing mesh
    static void IncreaseSideOrders(TPZCompMesh *fluxmesh);
    
    
    /// increase the order of the lower dimensional pressure elements
    void IncreasePressureSideOrders(TPZCompMesh *pressure_mesh);
    
    /// compute the average pressures of the hybridized form of the H(div) mesh
    void ComputeAveragePressures();
    
    /// set the cornernode values equal to the averages
    void ComputeNodalAverages();
    
    /// switch material object from mixed poisson to TPZMixedHdivErrorEstimate
    void SwitchMaterialObjects();
    
    /// clone the meshes into the post processing mesh
    void CloneMeshVec();
    
    /// create the multiphysics mesh using the TPZMixedErrorEstimate material
    void CreateMultiphysicsMesh();
    
    /// clone the fluxmesh but using TPZCompElReferred objects
    void CreateFluxMesh();
    
    /// clone the pressure mesh using TPZCompElReferred objects
    void CreatePressureMesh();
    
    /// compute the effectivity indices of the pressure error and flux error and store in the element solution
    void ComputeEffectivityIndices();
    
    /// returns true if the material associated with the element is a boundary condition
    /// and if the boundary condition is dirichlet type
    bool IsDirichletCondition(TPZGeoElSide gelside);
    
    /// return the value of the Dirichlet condition
    void GetDirichletValue(TPZGeoElSide gelside, TPZVec<STATE> &vals);

    /// identify the peripheral material objects and store the information in fHybridizer
    void IdentifyPeripheralMaterialIds();
};

#endif /* TPZHybridHDivErrorEstimator_hpp */
