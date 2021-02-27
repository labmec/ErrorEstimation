//
//  TPZHybridHDivErrorEstimator.hpp
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
#include "pzanalysis.h"

class TPZCompMesh;
class TPZSubCompMesh;

/// this class will compute the estimated value of the energy error of the input mesh
// it is expected that the estimated value of the error is improved if the input mesh is of H(div)++ type
// the class should work for any H(div) mesh
// first create the post processing mesh
// then compute the errors
class TPZHybridHDivErrorEstimator
{
protected:
    /// The H(Div) approximation mesh for which we will compute the error
    TPZMultiphysicsCompMesh *fOriginal;

    /// order increase of the boundary flux (depending on the original mesh)
    int fUpliftPostProcessOrder = 0;

    /// Locally created computational mesh to compute the error
    TPZMultiphysicsCompMesh fPostProcMesh;

    /// weights associated with the dim-1 pressure elements to compute the averages
    TPZVec<REAL> fPressureweights;
    /// weights associated with material ids
    std::map<int,REAL> fMatid_weights;
    // object to assist in creating a hybridized version of the computational mesh
    TPZHybridizeHDiv fHybridizer;

    // material id of the dim-1 skeleton elements
    int fPressureSkeletonMatId = 0;

    TPZAnalyticSolution *fExact;
    ProblemConfig fProblemConfig;
    /// whether the post processing mesh will be H(div) or H1
    bool fPostProcesswithHDiv = false;

public:
    explicit TPZHybridHDivErrorEstimator(TPZMultiphysicsCompMesh &originalMesh, bool postProcWithHDiv = false)
        : fOriginal(&originalMesh), fPostProcMesh(nullptr), fExact(nullptr) {
        fPostProcesswithHDiv = postProcWithHDiv;
    }

    // this method wont work because multiphysics meshes have no operator= (yet)
    TPZHybridHDivErrorEstimator(const TPZHybridHDivErrorEstimator &copy) = delete;

    // this method wont work because multiphysics meshes have no operator= (yet)
    TPZHybridHDivErrorEstimator &operator=(const TPZHybridHDivErrorEstimator &cp) = delete;

    ~TPZHybridHDivErrorEstimator();

    void SetProblemConfig(const ProblemConfig &cfg) { fProblemConfig = cfg; }

    void SetPostProcUpliftOrder(const int order) { fUpliftPostProcessOrder = order; }

    /// Set the analytic solution object
    void SetAnalyticSolution(TPZAnalyticSolution &exact) { fExact = &exact; }

    void SetHybridizer(TPZHybridizeHDiv &hybridizer) { fHybridizer = hybridizer; }

    /// compute the element errors comparing the reconstructed solution based on average pressures
    /// with the original solution
    virtual void ComputeErrors(TPZVec<REAL> &error_vec, TPZVec<REAL> &element_errors, std::string&vtkPath);

    // reconstruction of potential using hybrid solution on enrichment space
    virtual void PotentialReconstruction();

    /// create graphical output of estimated and true errors using the analysis
    void PostProcessing(TPZAnalysis &an, std::string &out);

    void PlotPressureSkeleton(const std::string &filename, bool reconstructed = true);
    void PlotInterfaceFluxes(const std::string &filename, bool reconstructed = true);

    // Plots State solution of elements of target dimension
    static void PlotState(const std::string& filename, int targetDim, TPZCompMesh* cmesh, bool atomic = true);

    int PressureSkeletonMatId() const { return fPressureSkeletonMatId; }

    TPZMultiphysicsCompMesh *PostProcMesh() { return &fPostProcMesh; }
    TPZGeoMesh *GMesh() { return fOriginal->Reference(); }

protected:

    /// compute the pressure weights and material weights
    // fills in the data structure pressureweights and matid_weights
    virtual void ComputePressureWeights();

    /// create the post processed multiphysics mesh (which is necessarily hybridized)
    virtual void CreatePostProcessingMesh();

    /// computing the element stifnesses will "automatically" compute the condensed form of the matrices
    void ComputeElementStiffnesses();

    // a method for generating the HDiv mesh
    virtual TPZCompMesh *CreateFluxMesh();

    // this method clones the original pressure mesh and perform the following:
    // 1. delete dim - 1 elements, where dim is the mesh dimension
    // 2. build BC elements // TODO improve comments
    // 3. create skeleton geometric elements and comp elements
    virtual TPZCompMesh *CreatePressureMesh();

    /// return a pointer to the pressure mesh
    virtual TPZCompMesh *PressureMesh();

    // Finds a material ID that has not been used yet
    static int FindFreeMatId(TPZGeoMesh *gmesh);

    // Creates skeleton geometric elements on which the average pressure will be calculated
    virtual void CreateSkeletonElements(TPZCompMesh *pressure_mesh);

    /// increase the side orders of the post processing mesh
    static void IncreaseSideOrders(TPZCompMesh *fluxmesh);


    /// increase the order of the lower dimensional pressure elements
    static void IncreasePressureSideOrders(TPZCompMesh *pressure_mesh);

    /// compute the average pressures of across faces of the H(div) mesh
    void ComputeAverageFacePressures();

    void ComputeBoundaryL2Projection(int target_dim);
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

    /// restrain the edge elements that have larger elements as neighbours
    void RestrainSmallEdges(TPZCompMesh *pressuremesh);

    /// set the cornernode values equal to the averages
    virtual void ComputeNodalAverages();

    /// compute the nodal average of all elements that share a point
    void ComputeNodalAverage(TPZCompElSide &node_celside);
    //compute the global efectivity index using the CharacteristicSize() of element
    void GlobalEffectivityIndex();

    /// copy the solution from the neighbouring skeleton elements
    // this is a placeholder for the derived class TPZHDivErrorEstimatorH1
    virtual void CopySolutionFromSkeleton();

    /// switch material object from mixed poisson to TPZMixedHdivErrorEstimate
    virtual void SwitchMaterialObjects();

    /// clone the meshes into the post processing mesh
    void CloneMeshVec();

    /// compute the effectivity indices of the pressure error and flux error and store in the element solution
    void ComputeEffectivityIndices();

    /// compute the effectivity indices of the pressure error and flux error and store in the element solution
    void ComputeEffectivityIndices(TPZSubCompMesh *cmesh);

    /// returns true if the material associated with the element is a boundary condition
    /// and if the boundary condition is dirichlet type
    bool IsDirichletCondition(TPZGeoElSide gelside);

    /// return the value of the Dirichlet condition
    void GetDirichletValue(TPZGeoElSide gelside, TPZVec<STATE> &vals);

    /// identify the peripheral material objects and store the information in fHybridizer
    void IdentifyPeripheralMaterialIds();

    void RestrainSkeletonSides(TPZCompMesh *pressure_mesh);

    // Checks if the solution is in fact continuous
    virtual void VerifySolutionConsistency(TPZCompMesh* cmesh);

  //  int FirstNonNullApproxSpaceIndex(TPZVec<TPZMaterialData> &datavec);

    // compute the average of an element iel in the pressure mesh looking at its neighbours
    void ComputeAverage(TPZCompMesh *pressuremesh, int64_t iel);

    void PrepareElementsForH1Reconstruction();

    bool IsAdjacentToHangingNode(const TPZCompElSide &celside);
};

#endif /* TPZHybridHDivErrorEstimator_hpp */
