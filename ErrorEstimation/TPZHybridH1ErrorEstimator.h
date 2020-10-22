//
//  TPZHybridHDivErrorEstimator.hpp
//  ErrorEstimateHDiv
//
//  Created by Philippe Devloo on 10/06/18.
//

#ifndef TPZHybridH1ErrorEstimator_hpp
#define TPZHybridH1ErrorEstimator_hpp

#include <stdio.h>
#include "pzmanvector.h"
#include "TPZHybridizeHDiv.h"
#include "TPZCreateMultiphysicsSpace.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZAnalyticSolution.h"
#include "ProblemConfig.h"
#include "pzanalysis.h"


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

    /// This method is contructed for HybridSquared meshes. Therefore this variable shall be redundant for now.
    bool fOriginalIsHybridized = true;
    
    /// order increase of the boundary flux (depending on the original mesh)
    int fUpliftPostProcessMesh = 0;
    
    /// whether the post processing mesh will be H(div) or H1
    bool fPostProcesswithHDiv = true;
    
    /// Locally created computational mesh to compute the error
    TPZMultiphysicsCompMesh fPostProcMesh;
//    TPZManVector<TPZCompMesh *,7> fPostProcMesh;
    
    /// weights associated with the dim-1 pressure elements to compute the averages
    TPZVec<REAL> fPressureweights;
    /// weights associated with material ids
    std::map<int,REAL> fMatid_weights;
    /// compute the pressure weights and material weights
    // fills in the data structure pressureweights and matid_weights
   virtual void ComputePressureWeights();
    // object to assist in creating a hybridized version of the computational mesh
    TPZHybridizeHDiv fHybridizer;
    //TPZCreateMultiphysicsSpace fHybrid;

    // material id of the dim-1 skeleton elements
    int fPressureSkeletonMatId;

    // material id of the HDiv reconstruction
    int fHDivResconstructionMatId;
    
    TPZAnalyticSolution *fExact;
    
    ProblemConfig fProblemConfig;

    std::string fDebugDirName = "HybridH1_ReconstructionDebug";

    TPZHybridH1ErrorEstimator(TPZMultiphysicsCompMesh &InputMesh) : fOriginal(&InputMesh),
    fPostProcMesh(0),fExact(NULL)
    {
        FindFreeMatID(fPressureSkeletonMatId);
        FindFreeMatID(fHDivResconstructionMatId);
    }

    TPZHybridH1ErrorEstimator(TPZMultiphysicsCompMesh &InputMesh, int skeletonMatId, int HDivMatId) : fOriginal(&InputMesh),
                                                                    fPostProcMesh(0),fExact(NULL),
                                                                    fPressureSkeletonMatId(fPressureSkeletonMatId),fHDivResconstructionMatId(HDivMatId)
    {

    }
    
    TPZHybridH1ErrorEstimator(const TPZHybridH1ErrorEstimator &copy) : fOriginal(copy.fOriginal),
        fOriginalIsHybridized(copy.fOriginalIsHybridized),fUpliftPostProcessMesh(copy.fUpliftPostProcessMesh),
        fPostProcMesh(copy.fPostProcMesh), fExact(copy.fExact), fProblemConfig(copy.fProblemConfig),fPressureSkeletonMatId(copy.fPressureSkeletonMatId)
    {
        // this method wont work because multiphysics meshes have no copy constructor (yet)
        DebugStop();
    }
    
    TPZHybridH1ErrorEstimator &operator=(const TPZHybridH1ErrorEstimator &cp)
    {
        fOriginal = cp.fOriginal;
        fOriginalIsHybridized = cp.fOriginalIsHybridized;
        fUpliftPostProcessMesh = cp.fUpliftPostProcessMesh;
        // this method wont work because multiphysics meshes have no operator= (yet)
        DebugStop();

        fPostProcMesh = cp.fPostProcMesh;
        fExact = cp.fExact;
        fProblemConfig = cp.fProblemConfig;
        fPressureSkeletonMatId = cp.fPressureSkeletonMatId;
        return *this;
    }
    
    ~TPZHybridH1ErrorEstimator();
    
    /// Set the analytic solution object
    void SetAnalyticSolution(TPZAnalyticSolution &exact)
    {
        fExact = &exact;
    }

    void SetHybridizer(TPZHybridizeHDiv &hybridizer) {
        fHybridizer = hybridizer;
    }
    /// compute the element errors comparing the reconstructed solution based on average pressures
    /// with the original solution
    virtual void ComputeErrors(TPZVec<REAL> &errorVec, TPZVec<REAL> &elementerrors, bool store);

    //reconstruction of potential using hybrid solution on enrichement space
    virtual void PotentialReconstruction();
    
    /// create graphical output of estimated and true errors using the analysis
    void PostProcessing(TPZAnalysis &an);
    
    void PlotLagrangeMultiplier(const std::string &filename, bool reconstructed = true);

    // Plots State solution of elements of target dimension
    void PlotState(const std::string& filename, int targetDim, TPZCompMesh* cmesh);

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

    // Creates skeleton geometric elements on which the average pressure will be calculated
    virtual void CreateSkeletonElements(TPZCompMesh *pressure_mesh);
    
    /// increase the side orders of the post processing mesh
    static void IncreaseSideOrders(TPZCompMesh *fluxmesh);
    
    /// increase the order of the lower dimensional pressure elements
    void IncreasePressureSideOrders(TPZCompMesh *pressure_mesh);
    
    /// compute the average pressures of across faces of the H(div) mesh
    void ComputeAverageFacePressures();
    
    void ComputeBoundaryL2Projection(TPZCompMesh *pressuremesh,int target_dim);
    void NewComputeBoundaryL2Projection(TPZCompMesh *pressuremesh,int target_dim);
    void BoundaryPressureProjection(TPZCompMesh *pressuremesh, int target_dim);
    
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

    /// Restrain a side of a small element to a side of a large element.
    void RestrainSkeletonSides(TPZCompMesh *pressure_mesh);
    
    /// restrain the edge elements that have larger elements as neighbours
    void RestrainSmallEdges(TPZCompMesh *pressuremesh);
    
    /// set the cornernode values equal to the averages
    virtual void ComputeNodalAverages();
    
    /// compute the nodal average of all elements that share a point
    void ComputeNodalAverage(TPZCompElSide &celside);
    //compute the global efectivity index using the CharacteristicSize() of element
    void GlobalEffectivityIndex();
    
    /// copy the solution from the neighbouring skeleton elements
    // this is a placeholder for the derived class TPZHDivErrorEstimatorH1
    virtual void CopySolutionFromSkeleton();

    /// Insert material for HDiv reconstruction
    /// Switch H1 material for H1 reconstruction material
    virtual void InsertEEMaterial();

    /// Find free matID number
    void FindFreeMatID(int &matID);

    ///  Insert BC material into the pressure mesh material vector,
    ///  Create computational element on BC elements
    void AddBC2PressureMesh(TPZCompMesh *pressureMesh);
    
    /// clone the meshes into the post processing mesh
    void CloneMeshVec();
    
    /// compute the effectivity indices of the pressure error and flux error and store in the element solution
    void ComputeEffectivityIndices();

    /// compute the effectivity indices of the pressure error and flux error and store in the element solution
    void ComputeEffectivityIndices(TPZSubCompMesh *cmesh);

    /// Compute skeleton averages;
    void MakeSkeletonContinuous();
    
    /// returns true if the material associated with the element is a boundary condition
    /// and if the boundary condition is dirichlet type
    bool IsDirichletCondition(TPZGeoElSide gelside);
    
    /// return the value of the Dirichlet condition
    void GetDirichletValue(TPZGeoElSide gelside, TPZVec<STATE> &vals);

    /// identify the peripheral material objects and store the information in fHybridizer
    void IdentifyPeripheralMaterialIds();
    

    // Checks if the solution is in fact continuous
    virtual void VerifySolutionConsistency(TPZCompMesh* cmesh);
    
  //  int FirstNonNullApproxSpaceIndex(TPZVec<TPZMaterialData> &datavec);

protected:
    
    // compute the average of an element iel in the pressure mesh looking at its neighbours
    void ComputeAverage(TPZCompMesh *pressuremesh, int64_t iel);


    void PrepareElementsForH1Reconstruction();

};

#endif /* TPZHybridH1ErrorEstimator_hpp */