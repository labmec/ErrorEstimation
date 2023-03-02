//
// Created by victor on 22/02/23.
//

#ifndef ERRORESTIMATION_TPZHYBRIDH1CREATEH1RECONSTRUCTION_H
#define ERRORESTIMATION_TPZHYBRIDH1CREATEH1RECONSTRUCTION_H

#include "pzerror.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZHybridH1ErrorEstimator.h"

class TPZHybridH1ErrorEstimator;
class TPZAnalysis;


class TPZHybridH1CreateH1Reconstruction{
public:
    TPZHybridH1CreateH1Reconstruction() {
        DebugStop();
    }

    TPZHybridH1CreateH1Reconstruction(TPZHybridH1ErrorEstimator *pEstimator);

    ~TPZHybridH1CreateH1Reconstruction();

    TPZCompMesh *CreateH1ReconstructionMesh();

    ///  Insert BC material into the pressure mesh material vector,
    ///  Create computational element on BC elements
    void AddBC2PressureMesh();

    // Creates skeleton geometric elements on which the average pressure will be calculated
    virtual void CreateSkeletonElements();

    /// increase the order of the lower dimensional pressure elements
    void IncreasePressureSideOrders();

    /// Restrain a side of a small element to a side of a large element.
    void RestrainSkeletonSides();

    /// create dim-2 skeleton mesh based on the dim-1 faces
    // will do nothing if the dimension of the mesh == 2
    void CreateEdgeSkeletonMesh();

    /// adjust the interpolation orders so as to create an H1/2 boundary mesh
    // this method is called by the CreateEdgeSkeletonMesh method
    void AdjustNeighbourPolynomialOrders();

    void BuildMultiphysicsSpace();

    /// restrain the edge elements that have larger elements as neighbours
    void RestrainSmallEdges();

    /// Create and delete geometric elements
    void PrepareGeometricElements();

    /// Compute skeleton averages;
    void MakeSkeletonContinuous();

    void CreateGroupedAndCondensedElements();

    // Checks if the skeleton is continuous
    void VerifySkeletonContinuity();

    /// computing the element stifnesses will "automatically" compute the condensed form of the matrices
    void ComputeElementStiffnesses();

    /// compute the average pressures of across edges of the H(div) mesh
    virtual void ComputeAveragePressures(int target_dim);

    /// transfer the solution of the edge functions to the face functions
    void TransferEdgeSolution();

    /// set the cornernode values equal to the averages
    virtual void ComputeNodalAverages();

    // If all nodal (0-dimensional) connects linked to this node have dependency, it lies* on a hanging side.
    bool LiesOnHangingSide(TPZCompElSide &node_celside);

    // Impose the solution on a hanging side
    void ImposeHangingNodeSolution(TPZCompElSide &node_celside);

    /// compute the nodal average of all elements that share a point
    void ComputeNodalAverage(TPZCompElSide &node_celside);

    /// compute the average pressures of across faces of the H(div) mesh
    void ComputeAverageFacePressures();

    ///Legacy code. Is this required?
    void NewComputeBoundaryL2Projection(TPZCompMesh *pressuremesh,int target_dim);
    void BoundaryPressureProjection(TPZCompMesh *pressuremesh, int target_dim);

    // compute the average of an element iel in the pressure mesh looking at its neighbours
    void ComputeAverage(TPZCompMesh *pressuremesh, int64_t iel);

    void ComputeBoundaryL2Projection(TPZCompMesh *pressuremesh,int target_dim);


    // Verify if average were well executed
    void VerifyAverage(int target_dim);

    /// copy the solution from the neighbouring skeleton elements
    //     this is a placeholder for the derived class TPZHDivErrorEstimatorH1
    virtual void CopySolutionFromSkeleton();

    /// returns true if the material associated with the element is a boundary condition
    /// and if the boundary condition is dirichlet type
    bool IsDirichletCondition(TPZGeoElSide gelside);

    /// compute the pressure weights and material weights
    // fills in the data structure pressureweights and matid_weights
    virtual void ComputePressureWeights();

    void AdditionalAverageWeights(TPZGeoEl *large, TPZGeoEl *small, REAL &large_weight, REAL &small_weight, REAL sum);

    static TPZGeoElSide HasNeighbour(const TPZGeoElSide &gelside, int matid);

    void PlotLagrangeMultiplier(const std::string &filename, bool reconstructed = true);

    void GetPressureMatIDs(std::set<int> &matIDs);

private:
    TPZHybridH1ErrorEstimator *fHybridH1EE;

    TPZCompMesh* fPressureMesh = NULL;

    TPZMultiphysicsCompMesh *fMultiphysicsH1reconstructionMesh = NULL;

    TPZMultiphysicsCompMesh *fOriginal;

    /// weights associated with the dim-1 pressure elements to compute the averages
    TPZVec<REAL> fPressureweights;

    /// weights associated with material ids
    std::map<int,REAL> fMatid_weights;

    std::string fPressureReconstructionFolderOutput = "H1RecDebug/";

    /// Compute average pressure between skeletons on hanging nodes mode;
    /// 1 : average between both skeletons; 2 : Just small skeleton; 3 : weighted average w~1/h^(p+1)
    int fAverageMode = 0;
};

#endif // ERRORESTIMATION_TPZHYBRIDH1CREATEH1RECONSTRUCTION_H