//
//  TPZHDivErrorEstimator.hpp
//  ErrorEstimateHDiv
//
//  Created by Philippe Devloo on 10/06/18.
//

#ifndef TPZMHMHDIVERRORESTIMATOR_H
#define TPZMHMHDIVERRORESTIMATOR_H

#include "TPZHDivErrorEstimator.h"
#include "TPZMHMixedMeshControl.h"

/// this class will compute the estimated value of the energy error of an MHM mesh
// the class should work for any MHM-H(div) mesh
// first create the post processing mesh
// then compute the errors
// compute a reference solution without MHM
class TPZMHMHDivErrorEstimator : public TPZHDivErrorEstimator {

public:

    TPZMHMHDivErrorEstimator(TPZMultiphysicsCompMesh &originalMesh, TPZMHMixedMeshControl *mhm, bool postProcWithHDiv = false)
        : TPZHDivErrorEstimator(originalMesh, postProcWithHDiv), fMHM(mhm) {
    }

    // this method wont work because multiphysics meshes have no copy constructor (yet)
    TPZMHMHDivErrorEstimator(const TPZMHMHDivErrorEstimator &copy) = delete;

    // this method wont work because multiphysics meshes have no copy constructor (yet)
    TPZMHMHDivErrorEstimator &operator=(const TPZMHMHDivErrorEstimator &cp) = delete;

    virtual ~TPZMHMHDivErrorEstimator() = default;

private:
    // material id of the dim-1 multiphysic inerface and wrap elements
    // (only used for HDiv reconstruction)
    int fMultiPhysicsInterfaceMatId = 0;
    int fHDivWrapMatId = 0;

    /// a pointer to the datastructure used to generate the MHM mesh
    TPZMHMixedMeshControl *fMHM = nullptr;
    // a method for generating the HDiv mesh
    TPZCompMesh *CreateFluxMesh() override;
    // a method for creating the pressure mesh
    TPZCompMesh *CreatePressureMesh() override;
    // method fro creating a discontinuous pressure mesh
    TPZCompMesh *CreateDiscontinuousPressureMesh();
    // method for creating a pressure mesh that is continuous in each MHM domain, but not globally
    TPZCompMesh *CreateInternallyContinuousPressureMesh();

    // Creates skeleton geometric elements on which the average pressure will be calculated
    void CreateSkeletonElements(TPZCompMesh *pressure_mesh) override;
    // Creates H1 discontinuous space on skeleton elements
    void CreateSkeletonApproximationSpace(TPZCompMesh *pressure_mesh);
    // Creates TPZMultiphysicsInterface elements between subdomains if reconstructing with H(div)
    void CreateFluxSkeletonElements(TPZCompMesh *flux_mesh);
    void CreateMultiphysicsInterfaces();

    // a method for generating the hybridized multiphysics post processing mesh
    void CreatePostProcessingMesh() override;
    // a method for transferring the multiphysics elements in submeshes
    void SubStructurePostProcessingMesh();
    // transfer embedded elements to submeshes
    // the substructuring method will only transfer the H(div) and surrounding elements to the submesh
    // it does not detect that the boundary pressure elements belong to the submesh. This is done in the
    // following method
    void TransferEmbeddedElements();
    // remove the materials that are not listed in MHM
    void RemoveMaterialObjects(std::map<int,TPZMaterial *> &matvec);
    // a method for computing the pressures between subdomains as average pressures
    /// compute the average pressures of across edges of the H(div) mesh
    void ComputeAveragePressures(int target_dim) override;

    /// compute the average pressure over corners
    /// set the cornernode values equal to the averages
    void ComputeNodalAverages() override;

    void CopySolutionFromSkeleton() override;

    void VerifySolutionConsistency(TPZCompMesh* cmesh) override;

    // Fill a list with the connect indexes of volumetric elements sides
    // in the neighbourhood of the skeleton
    void ComputeConnectsNextToSkeleton(std::set<int64_t>& connectList);
};

#endif