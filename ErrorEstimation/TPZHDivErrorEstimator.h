//
//  TPZHDivErrorEstimator.hpp
//  ErrorEstimateHDiv
//
//  Created by Philippe Devloo on 10/06/18.
//

#ifndef TPZHybridHDivErrorEstimator_hpp
#define TPZHybridHDivErrorEstimator_hpp

#include "pzmanvector.h"
#include "TPZHybridizeHDiv.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZAnalyticSolution.h"
#include "ProblemConfig.h"
#include "TPZAnalysis.h"
#include "pzsubcmesh.h"

#include "Elasticity/TPZMixedElasticityND.h"

/// this class will compute the estimated value of the energy error of the input mesh
// it is expected that the estimated value of the error is improved if the input mesh is of H(div)++ type
// the class should work for any H(div) mesh
// first create the post processing mesh
// then compute the errors

template <typename MixedMaterial>
class TPZHDivErrorEstimator {
protected:
    /// The H(Div) approximation mesh for which we will compute the error
    TPZMultiphysicsCompMesh *fOriginal;

    /// Locally created computational mesh to compute the error
    TPZMultiphysicsCompMesh fPostProcMesh;

    /// weights associated with the dim-1 pressure/displacement elements to compute the averages
    TPZVec<REAL> fPrimalWeights;
    /// weights associated with material ids
    std::map<int,REAL> fMatid_weights;
    // object to assist in creating a hybridized version of the computational mesh
    TPZHybridizeHDiv fHybridizer;

    // material id of the dim-1 skeleton elements
    int fPrimalSkeletonMatId = 0;

    TPZAnalyticSolution *fExact;
    /// whether the post processing mesh will be H(div) or H1
    bool fPostProcesswithHDiv = false;

public:
    explicit TPZHDivErrorEstimator(TPZMultiphysicsCompMesh &originalMesh, bool postProcWithHDiv = false)
        : fOriginal(&originalMesh), fPostProcMesh(nullptr), fExact(nullptr) {
        fPostProcesswithHDiv = postProcWithHDiv;
    }

    // this method wont work because multiphysics meshes have no operator= (yet)
    TPZHDivErrorEstimator(const TPZHDivErrorEstimator<MixedMaterial> &copy) = delete;

    // this method wont work because multiphysics meshes have no operator= (yet)
    TPZHDivErrorEstimator &operator=(const TPZHDivErrorEstimator<MixedMaterial> &cp) = delete;

    ~TPZHDivErrorEstimator();

    //void SetProblemConfig(const ProblemConfig &cfg) { fProblemConfig = cfg; }

    /// Set the analytic solution object
    void SetAnalyticSolution(TPZAnalyticSolution &exact) { fExact = &exact; }

    void SetHybridizer(TPZHybridizeHDiv &hybridizer) { fHybridizer = hybridizer; }

    /**
     * Compute the element errors comparing the reconstructed solution based on 
     * average pressures/displacements
     * with the original solution
     */
    virtual void ComputeErrors(TPZVec<REAL> &error_vec, TPZVec<REAL> &element_errors, std::string& vtkPath);

    // reconstruction of pressure/displacement using hybrid solution on enrichment space
    virtual void PrimalReconstruction();

    /// create graphical output of estimated and true errors using the analysis
  virtual  void PostProcessing(TPZAnalysis &an, std::string &out);

    void PlotPrimalSkeleton(const std::string &filename, bool reconstructed = true);
    void PlotInterfaceFluxes(const std::string &filename, bool reconstructed = true);

    // Plots State solution of elements of target dimension
    static void PlotState(const std::string& filename, int targetDim, TPZCompMesh* cmesh, bool atomic = true);

    int PrimalSkeletonMatId() const { return fPrimalSkeletonMatId; }

    TPZMultiphysicsCompMesh *PostProcMesh() { return &fPostProcMesh; }
    TPZGeoMesh *GMesh() { return fOriginal->Reference(); }

protected:

    /// compute the pressure/displacement weights and material weights
    // fills in the data structure primalWeights and matid_weights
    virtual void ComputePrimalWeights() = 0;

    /// create the post processed multiphysics mesh (which is necessarily hybridized)
    virtual void CreatePostProcessingMesh();

    /// computing the element stifnesses will "automatically" compute the condensed form of the matrices
    void ComputeElementStiffnesses();

    // a method for generating the HDiv mesh
    virtual TPZCompMesh *CreateHDivMesh();

    // this method clones the original pressure/displacement mesh and performs the following:
    // 1. delete dim - 1 elements, where dim is the mesh dimension
    // 2. build BC elements // TODO improve comments
    // 3. create skeleton geometric elements and comp elements
    virtual TPZCompMesh *CreatePrimalMesh();

    /// return a pointer to the pressure/displacement mesh
    virtual TPZCompMesh *PrimalMesh();

    // Finds a material ID that has not been used yet
    static int FindFreeMatId(TPZGeoMesh *gmesh);

    // Creates skeleton geometric elements on which the average pressure/displacement will be calculated
    virtual void CreateSkeletonElements(TPZCompMesh *primal_mesh);

    /// increase the side orders of the post processing mesh
    static void IncreaseSideOrders(TPZCompMesh *fluxmesh);


    /// increase the order of the lower dimensional pressure/displacement elements
    static void IncreasePrimalSideOrders(TPZCompMesh *primal_mesh);

    void ComputeBoundaryL2Projection(int target_dim);

    /// compute the average pressures/displacements of across edges of the H(div) mesh
    virtual void ComputeAveragePrimal(int target_dim);

    /// transfer the solution of the edge functions to the face functions
    void TransferEdgeSolution();

    /// create dim-2 skeleton mesh based on the dim-1 faces
    // will do nothing if the dimension of the mesh == 2
    void CreateEdgeSkeletonMesh(TPZCompMesh *primal_mesh);

    /// adjusts the interpolation orders so as to create an H1/2 boundary mesh
    // this method is called by the CreateEdgeSkeletonMesh method
    void AdjustNeighbourPolynomialOrders(TPZCompMesh *primal_mesh);

    /// restrains the edge elements that have larger elements as neighbours
    void RestrainSmallEdges(TPZCompMesh *primalmesh);

    /// sets the cornernode values equal to the averages
    virtual void ComputeNodalAverages();

    /// computes the nodal average of all elements that share a point
    void ComputeNodalAverage(TPZCompElSide &node_celside);
    
    //computes the global effectivity index using the CharacteristicSize() of element
    void GlobalEffectivityIndex();

    // copies the solution from the neighbouring skeleton elements
    // this is a placeholder for the derived class TPZHDivErrorEstimatorH1
    virtual void CopySolutionFromSkeleton();

    /// switches material object from MixedMaterial to TPZMixedHdivErrorEstimate
    virtual void SwitchMaterialObjects();

    /// compute the effectivity indices of the pressure error and flux error and store in the element solution
    void ComputeEffectivityIndices();

    /// compute the effectivity indices of the pressure error and flux error and store in the element solution
    void ComputeEffectivityIndices(TPZSubCompMesh *subcmesh);

    /// returns true if the material associated with the element is a boundary condition
    /// and if the boundary condition is dirichlet type
    bool IsDirichletCondition(TPZGeoElSide gelside);

    void RestrainSkeletonSides(TPZCompMesh *pressure_mesh);

    // Checks if the solution is in fact continuous
    virtual void VerifySolutionConsistency(TPZCompMesh* cmesh);

    // computes the average of the element iel in the pressure/displacement mesh looking at its neighbours
    void ComputeAverage(TPZCompMesh *primalmesh, int64_t iel);

    void PrepareElementsForH1Reconstruction();

    bool IsAdjacentToHangingNode(const TPZCompElSide &celside);

    static std::set<int> GetBCMatIDs(const TPZCompMesh* cmesh);
};

#endif /* TPZHybridHDivErrorEstimator_hpp */
