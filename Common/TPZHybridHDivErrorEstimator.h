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
    
    // material id of the dim-2 skeleton elements
    int fSkeletonMatId = 6;
    
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
        // this method wont work because multiphysics meshes have no copy constructor (yet)
        DebugStop();
    }
    
    TPZHybridHDivErrorEstimator &operator=(const TPZHybridHDivErrorEstimator &cp)
    {
        fOriginal = cp.fOriginal;
        fOriginalIsHybridized = cp.fOriginalIsHybridized;
        fUpliftPostProcessMesh = cp.fUpliftPostProcessMesh;
        // this method wont work because multiphysics meshes have no operator= (yet)
        DebugStop();

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
    virtual void ComputeErrors(TPZVec<REAL> &elementerrors, bool store = true);
    
    //reconstruction of potential using hybrid solution on enrichement space
    virtual void PotentialReconstruction();
    
    /// create graphical output of estimated and true errors using the analysis
    void PostProcessing(TPZAnalysis &an);
    
    
    void PlotLagrangeMultiplier(const std::string &filename, bool reconstructed = true);
protected:
    
    /// create the post processed multiphysics mesh (which is necessarily hybridized)
    virtual void CreatePostProcessingMesh();
    
    /// computing the element stifnesses will "automatically" compute the condensed form of the matrices
    void ComputeElementStiffnesses();
    
    // a method for generating the HDiv mesh
    virtual TPZCompMesh *CreateFluxMesh();
    // a method for creating the pressure mesh
    virtual TPZCompMesh *CreatePressureMesh();

    /// return a pointer to the pressure mesh
    virtual TPZCompMesh *PressureMesh();
    
    /// increase the side orders of the post processing mesh
    static void IncreaseSideOrders(TPZCompMesh *fluxmesh);
    
    
    /// increase the order of the lower dimensional pressure elements
    void IncreasePressureSideOrders(TPZCompMesh *pressure_mesh);
    
    /// compute the average pressures of across faces of the H(div) mesh
    void ComputeAverageFacePressures();
    
    /// compute the average pressures of across edges of the H(div) mesh
    virtual void ComputeAveragePressures(int target_dim);
    
    /// transfer the solution of the edge functions to the face functions
    void TransferEdgeSolution();
    
    /// create dim-2 skeleton mesh based on the dim-1 faces
    // will do nothing if the dimension of the mesh == 2
    void CreateEdgeSkeletonMesh(TPZCompMesh *pressure_mesh);
    
    /// adjust the interpolation orders so as to create an H1/2 boundary mesh
    // this method is called by the CreateEdgeSkeletonMesh method
    void AdjustNeighbourPolynomialOrders(TPZCompMesh *pressure_mesh);
    
    /// restrain the edge elements that have larger elements as neighbours
    void RestrainSmallEdges(TPZCompMesh *pressuremesh);
    
    /// set the cornernode values equal to the averages
    virtual void ComputeNodalAverages();
    
    /// compute the nodal average of all elements that share a point
    void ComputeNodalAverage(TPZCompElSide &celside);
    
    /// copy the solution from the neighbouring skeleton elements
    // this is a placeholder for the derived class TPZHDivErrorEstimatorH1
    virtual void CopySolutionFromSkeleton()
    {
        
    }
    
    /// switch material object from mixed poisson to TPZMixedHdivErrorEstimate
    virtual void SwitchMaterialObjects();
    
    /// clone the meshes into the post processing mesh
    void CloneMeshVec();
    
    /// compute the effectivity indices of the pressure error and flux error and store in the element solution
    void ComputeEffectivityIndices();
    
    /// returns true if the material associated with the element is a boundary condition
    /// and if the boundary condition is dirichlet type
    bool IsDirichletCondition(TPZGeoElSide gelside);
    
    /// return the value of the Dirichlet condition
    void GetDirichletValue(TPZGeoElSide gelside, TPZVec<STATE> &vals);

    /// identify the peripheral material objects and store the information in fHybridizer
    void IdentifyPeripheralMaterialIds();
    

    // Checks if the solution is in fact continuous
    void VerifySolutionConsistency(TPZCompMesh* cmesh);
    
protected:
    
    // compute the average of an element iel in the pressure mesh looking at its neighbours
    void ComputeAverage(TPZCompMesh *pressuremesh, int64_t iel);
};

#endif /* TPZHybridHDivErrorEstimator_hpp */