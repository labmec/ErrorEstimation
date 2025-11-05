#ifndef TPZPOSTPROCESSERRORSBFEM_H
#define TPZPOSTPROCESSERRORSBFEM_H

#include "TPZPostProcessError.h"
#include "TPZBuildSBFemHybrid.h"
#include "TPZSBFemElementGroup.h"
#include "TPZGeoPatch.h"
#include "TPZMaterialT.h"
class TPZMultiPhysicsMeshWindow;

/// @brief class organizing computing the post processing error for SBFem simulations
class TPZPostProcessErrorSBFem : public TPZPostProcessError {
public:
    enum MeshPosition {EFlux = 0, EPressure = 1, EAverage = 2, EOrigin = 3};
    // Constructor
    TPZPostProcessErrorSBFem(TPZCompMesh *origin, TPZCompMesh *postprocess, TPZBuildSBFemHybrid &buildSBFemHybrid);


    TPZPostProcessErrorSBFem(TPZCompMesh * origin) : TPZPostProcessError(origin) {
        DebugStop();
    }

    TPZPostProcessErrorSBFem(TPZCompMesh * origin,ProblemConfig &config, bool useHDiv)  : TPZPostProcessError(origin,config,useHDiv) {
        DebugStop();
    }

    TPZPostProcessErrorSBFem(TPZVec<TPZCompMesh *> &meshvec) : TPZPostProcessError(meshvec) {
        DebugStop();
    }

    // Destructor
    ~TPZPostProcessErrorSBFem();

    virtual void BuildPatchStructures2() override; //one patch by color
    
    // create the mesh that allow us to compute the error estimate
    void CreateMultiphysicsMesh();
    
    // create a multiphysics window mesh
    void CreateMultiphysicsWindowMesh();
    
    /// @brief Compute a locally conservative hybrid H1 approximation using reconstruction
    void ReconstructHybridH1();

    // compute the estimated H1 seminorm errors
    virtual void ComputeElementErrors(TPZVec<STATE> &elementerrors) override;
    
    /// Add the Interface elements to the multiphysics mesh
    virtual void AddInterfaceElements(TPZMultiphysicsCompMesh *mfmesh) override;

    

    
    // print the relevant information of the patches
    virtual void PrintPatchInformation(std::ostream &out = std::cout) override;
    
    // plot the relevant patch information to a plot file
    void PlotPatches(const std::string rootname);
    


private:
    // Add any additional members or methods specific to this class
    TPZBuildSBFemHybrid fBuildSBFemHybrid;
    
    /// add coarse scale skeleton elements so that the intermediate nodes will be restrained
    void AddLevel0SkeletonElements(TPZCompMesh *partition_unitiy_mesh);
    
    /// restrain the central node of each patch
    void RestrainCentralConnectofPatches();
    
    /// Expand the geometric patch data to include the skeleton elements, the interface elements and flux elements
    /// These will be used to create a TPZMultiphysicsWindow object
    void ExpandGeoPatchMeshes();
    
    /// insert post processing materials in the multiphysics mesh
    void InsertPostProcessingMaterials(TPZMultiphysicsCompMesh *mphys);
    
    /// Creating an average pressure mesh corresponding to elements of mesh dimension
    virtual void CreateAveragePressureMesh() override;
    
    /// Identify the element groups in the H1 and Hybrid H1 mesh
    void IdentifyElementGroups();
    
    
    /// Expose the SBFemVolume elements and hide the element groups
    /// this method is necessary to build the multiphysics mesh
    void ExposeSBFemVolume();
    
    /// Expose the SBFemVolume elements and hide the element groups
    /// this method is necessary to build the multiphysics mesh
    void HideSBFemVolume();
    
    /// create and initialize the interface material objects
    void InitializeInterfaceMaterialObjects();
    
    /// HideInterfaceMaterial objects
    void HideInterfaceMaterial(TPZCompMesh *cmesh);
    
    /// Insert the interface material objects
    void InsertInterfaceMaterial(TPZCompMesh *cmesh);
    
    /// put the multiphysics elements corresponding to TPZSBFemVolume into postprocess groups
    void GroupSBFemMultiphysics(TPZMultiphysicsCompMesh *mphys, bool hasboundary);
    
    /// Element groups present in the original mesh
    TPZStack<TPZSBFemElementGroup *> fElementGroupsH1;

    /// Element groups present in the hybrid mesh
    TPZStack<TPZSBFemElementGroup *> fElementGroupsHybrid;
    
    /// vector of TPZGeoPatch objects, defining the patch for each flux reconstruction
    TPZVec<TPZGeoPatch> fElementPatches;
    
    /// vector of interface material objects. these pointers are inserted and hidden when creating a multiphysics mesh
    TPZManVector<TPZMaterialT<STATE> *> fInterfaceMaterials;

    /// @brief identify connect multiplicator values
    /// if the connect has no dependency, the multiplicator value is 1
    /// other wise it will be a list of connects and scalar values
    ///  the sum of the multiplicator values should always be one
    void ComputeConnectMultiplicators(TPZCompMesh &cmesh, int64_t cindex, REAL depval, std::map<int64_t,REAL> &mult);
    
    /// Compute the element graph as a function of connect indices for the partition of unity mesh
    void ComputeElementGraph(TPZStack<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex);
    
    /// compute Connect to element graph
    void ComputeConnectToElementGraph(TPZVec<int64_t> &connectgraph, TPZVec<int64_t> &connectgraphindex);
    
    /// set the center nodes of the discontinuous elements to the center of the patch
    void LoadElementCenters(TPZMultiPhysicsMeshWindow *cmesh, TPZGeoPatch &patch);
    
};

#endif // TPZPOSTPROCESSERRORSBFEM_H
