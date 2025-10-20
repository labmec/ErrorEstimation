
/**
 * @file
 * @brief Contains the TPZMultiPhysicsMeshWindow class which is a specialized TPZCompMesh for multiphysics simulations.
 */
#ifndef TPZMULTIPHYSICSMESHWINDOW_H
#define TPZMULTIPHYSICSMESHWINDOW_H

#include "TPZMultiphysicsCompMesh.h"

class TPZMultiPhysicsMeshWindow : public TPZMultiphysicsCompMesh {

    /// @brief Correspondence between connect indices in the multiphysics mesh and the individual physics meshes
    // Each entry in the outer vector corresponds to a physics space
    // Each inner vector contains pairs of (multiphysics connect index, physics mesh connect index)
    // This structure helps in mapping solutions and boundary conditions across different physics
    TPZManVector<TPZVec<int64_t>, 7> m_connect_correspondence;

public:
    TPZMultiPhysicsMeshWindow();
    
    /// Constructor based on TPZGeoMesh pointer and vector of meshes
    TPZMultiPhysicsMeshWindow(TPZGeoMesh * gmesh, bool isComplex=false) : TPZMultiphysicsCompMesh(gmesh,isComplex) {
        
    }
    
    /// Constructor based on TPZGeoMesh pointer and vector of meshes
    TPZMultiPhysicsMeshWindow(TPZAutoPointer<TPZGeoMesh>  gmesh, bool isComplex=false) : TPZMultiphysicsCompMesh(gmesh,isComplex) {
        
    }
    


    TPZMultiPhysicsMeshWindow(const TPZMultiPhysicsMeshWindow &other);
    TPZMultiPhysicsMeshWindow& operator=(const TPZMultiPhysicsMeshWindow &other);
    virtual ~TPZMultiPhysicsMeshWindow();

    /// Automatic builder for the computational mesh structure
    /// Build the multiphysics space using the previously provided mesh vector and active spaces
    virtual void AutoBuild() override {
        std::cout << __PRETTY_FUNCTION__ << " should not be called\n";
    }

        /// Build the multiphysics space where the elements are associated with a material with memory
    /// This signature assumes the mesh vector and active spaces have been already set
    virtual void BuildMultiphysicsSpaceWithMemory() override{
        std::cout << __PRETTY_FUNCTION__ << " has not been implemented. Use BuildMultiphysicsSpace instead\n";
        std::cout << "Please implement this method if you intend to use materials with memory in TPZMultiPhysicsMeshWindow.\n";
        DebugStop();
    }

    /// Build the multiphysics space where the elements are associated with a material with memory
    /// This signature assumes the mesh vector and active spaces have been already set
    virtual void BuildMultiphysicsSpaceWithMemory(std::set<int> matsIdWithMem, std::set<int> matsIdNoMem) override {
        std::cout << __PRETTY_FUNCTION__ << " has not been implemented. Use BuildMultiphysicsSpace instead\n";
        std::cout << "Please implement this method if you intend to use materials with memory in TPZMultiPhysicsMeshWindow.\n";
        DebugStop();
    }

    /// @brief Build the multiphysics space for the corresponding geometric elements
    virtual void BuildMultiphysicsSpace(const TPZVec<int64_t> &gelindexes) override;

    /// @brief build the map between connect indices of the atomic mesh and the multiphysics mesh
    void BuildConnectMap(int imesh, std::map<int64_t, int64_t> &connectmap);

        /// @brief Load the solution from the meshes in the multiphysics mesh
    virtual void LoadSolutionFromMeshes() override {
        DebugStop();
    }
    
    /// @brief Load the solution from the multiphysics mesh to the meshes in the multiphysics mesh
    virtual void LoadSolutionFromMultiPhysics() override {
        DebugStop();
    }

    int64_t AtomicConnectIndex(int imesh, int64_t connectindex) {
        return m_connect_correspondence[imesh][connectindex];
    }

protected:
    /// @brief Add the connects from the atomic meshes
    virtual void AddConnects() override;

    /// @brief  Add the atomic elements to the multiphysics elements
    virtual void AddElements() override;
};

#endif // TPZMULTIPHYSICSMESHWINDOW_H
