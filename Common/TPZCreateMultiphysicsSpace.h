//
//  TPZCreateMultiPhysicsSpace.hpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 13/07/19.
//

#ifndef TPZCreateMultiPhysicsSpace_hpp
#define TPZCreateMultiPhysicsSpace_hpp

#include <stdio.h>
#include <set>
#include "pzmanvector.h"
class TPZCompMesh;
class TPZGeoMesh;
class TPZMultiphysicsCompMesh;
class TPZCompEl;


class TPZCreateMultiPhysicsSpace
{
public:
    /// types of spaces this class can create
    enum MSpaceType {Enone, EH1Hybrid, EH1HybridHybrid};
    
private:
    
    /// the type of space this object will generate
    MSpaceType fSpaceType = Enone;
    
    /// the materialids which will be used to create the atomic meshes
    std::set<int> fMaterialIds;
    
    /// the boundary condition material ids
    std::set<int> fBCMaterialIds;
    
    /// default internal order for the H1 elements
    int fDefaultPOrder = 2;
    
    /// default Lagrange multiplier order
    int fDefaultLagrangeOrder = 1;
    
    /// the dimension of the geometric elements that will be used to generate computational elements
    int fDimension = -1;
    
    /// the geometric mesh which will generate the computational mesh
    TPZGeoMesh *fGeoMesh = 0;
    
public:
    
    /// All parameters needed for creating a hybrid H1 space
    struct TConfigH1Hybrid
    {
        /// the material id of the boundary elements that serve as a placeholder for the interface elements
        std::pair<int,int> fMatWrapId = {-1,-1};
        
        /// material id of the dim-1 flux elements
        int fFluxMatId = -1;
        
        /// material id of left and right lagrange multipliers
        std::pair<int, int> fLagrangeMatid = {-1,-1};
        
        /// indicated whether the boundary conditions should be hybridized as well
        bool fHybridizeBC = false;
        
        /// default constructor
        TConfigH1Hybrid(){}
        
        /// copy constructor
        TConfigH1Hybrid(const TConfigH1Hybrid &copy);
        
        /// copy operator
        TConfigH1Hybrid &operator=(const TConfigH1Hybrid &copy);
    };
    
    /// object which contains the relevant information for create a hybrid H1 mesh
    TConfigH1Hybrid fH1Hybrid;
    
    /// default constructor
    TPZCreateMultiPhysicsSpace(TPZGeoMesh *gmesh);
    
    /// copy constructor
    TPZCreateMultiPhysicsSpace(const TPZCreateMultiPhysicsSpace &copy);
    
    /// = operator
    TPZCreateMultiPhysicsSpace &operator=(const TPZCreateMultiPhysicsSpace &copy);
    
    /// Indicate to create Hybridized H1 meshes
    void SetH1Hybridized(const TConfigH1Hybrid &config);
    
    /// create meshes and elements for all geometric elements
    void CreateAtomicMeshes(TPZVec<TPZCompMesh *> &meshvec);
    
    /// add interface elements to the multiphysics space
    void AddInterfaceElements(TPZMultiphysicsCompMesh *mphys);
    
    /// group and condense the elements
    void GroupandCondenseElements(TPZMultiphysicsCompMesh *mphys);
    
private:
    
    /// Create the pressure mesh
    TPZCompMesh *CreatePressureMesh();
    
    /// Create the flux mesh
    TPZCompMesh *CreateFluxMesh();

    /// create the geometric elements for the lagrange multipliers
    // these elements will go with the largest H1 element
    void CreateLagrangeGeometricElements(TPZCompMesh *H1mesh);
    
    /// create the pressure boundary elements if the boundary is not hybridized
    void CreatePressureBoundaryElements(TPZCompMesh *pressure);
    
    /// insert the pressure material ids
    void InsertPressureMaterialIds(TPZCompMesh *pressure);
    
    /// insert flux material ids
    void InsertFluxMaterialIds(TPZCompMesh *fluxmesh);
    
    /// insert materialids for the null space
    void InsertNullSpaceMaterialIds(TPZCompMesh *nullspace);
    
    /// Find the neighbouring flux element
    TPZCompEl *FindFluxElement(TPZCompEl *wrapelement);
};


#endif /* TPZCreateMultiPhysicsSpace_hpp */
