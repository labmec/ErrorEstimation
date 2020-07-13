//
//  TPZCreateMultiphysicsSpace.hpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 13/07/19.
//

#ifndef TPZCreateMultiphysicsSpace_hpp
#define TPZCreateMultiphysicsSpace_hpp

#include <stdio.h>
#include <set>
#include "pzmanvector.h"
class TPZCompMesh;
class TPZGeoMesh;
class TPZMultiphysicsCompMesh;
class TPZCompEl;
class TPZGeoElSide;


class TPZCreateMultiphysicsSpace
{
public:
    /// types of spaces this class can create
    enum MSpaceType {Enone, EH1Hybrid, EH1HybridSquared};
    
private:
    
    /// the type of space this object will generate
    MSpaceType fSpaceType = Enone;
    
    /// the materialids which will be used to create the atomic meshes
    std::set<int> fMaterialIds;
    
    /// the boundary condition material ids
    std::set<int> fBCMaterialIds;
    
    /// default internal order for the H1 elements
    int fDefaultPOrder = 3;
    
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
        /// the material id of  elements created on the sides that serve as a placeholder for the interface elements
        int fMatWrapId = -1;
        
        /// material id of the dim-1 flux elements
        int fFluxMatId = -1;

        /// material id of left and right lagrange multipliers
        std::pair<int, int> fLagrangeMatid = {-1,-1};
        
        /// material for the lagrange multiplier between flux and interface pressure
        int fSecondLagrangeMatid = {-1};
        
        /// material id of the pressure interface element
        int fInterfacePressure = -1;
        
        /// indicated whether the boundary conditions should be hybridized as well
        int fHybridizeBCLevel = 0;
        
        bool fHybridSquared = false;
        /// indicates whether a second hybridizations will be applied
        
        /// default constructor
        TConfigH1Hybrid(){}
        
        /// copy constructor
        TConfigH1Hybrid(const TConfigH1Hybrid &copy);
        
        /// copy operator
        TConfigH1Hybrid &operator=(const TConfigH1Hybrid &copy);
    };
    
    void SetPOrder(int order){
        fDefaultPOrder = order;
    }
    
    void SetLagrangeOrder(int order)
    {
        fDefaultLagrangeOrder = order;
    }
    
    
    /// object which contains the relevant information for create a hybrid H1 mesh
    TConfigH1Hybrid fH1Hybrid;
    
    /// default constructor
    TPZCreateMultiphysicsSpace(TPZGeoMesh *gmesh, MSpaceType spacetype = EH1Hybrid);
    
    /// copy constructor
    TPZCreateMultiphysicsSpace(const TPZCreateMultiphysicsSpace &copy);
    
    /// = operator
    TPZCreateMultiphysicsSpace &operator=(const TPZCreateMultiphysicsSpace &copy);
    
    /// Configure the Hybridized H1 meshes
    void SetH1Hybridized(const TConfigH1Hybrid &config);
    
    /// Compute Periferal Material ids
    // the material ids will be computed from a number whose modulus by base is zero
    void ComputePeriferalMaterialIds(int base = 10);
    
    /// Insert the periferal material objects (for wrapmatid, fluxmatid and lagrange matid
    void InsertPeriferalMaterialObjects(TPZMultiphysicsCompMesh *mphys);
    
    /// insert the lagrange material objects
    void InsertLagranceMaterialObjects(TPZMultiphysicsCompMesh *mphys);
    
    /// Initialize the material ids and bc material ids
    void SetMaterialIds(const std::set<int> &matids, const std::set<int> &bc_matids)
    {
        fMaterialIds = matids;
        fBCMaterialIds = bc_matids;
    }
    /// create meshes and elements for all geometric elements
    void CreateAtomicMeshes(TPZVec<TPZCompMesh *> &meshvec,int pressureOrder, int lagrangeorder);
    
    /// add interface elements to the multiphysics space
    void AddInterfaceElements(TPZMultiphysicsCompMesh *mphys);
    
    /// group and condense the elements
    void GroupandCondenseElements(TPZMultiphysicsCompMesh *mphys);
    
private:
    
    /// Create geometric elements needed for the computational elements
    void AddGeometricWrapElements();
    
    /// Create the pressure mesh with order fDefaultOrder
    TPZCompMesh *CreatePressureMesh();
    
    /// Create the flux mesh with order fDefaultLagrangeOrder
    TPZCompMesh *CreateBoundaryFluxMesh();

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
    
    /// if there a neighbouring element with matid == lagrangematid -> return true
    bool ShouldCreateFluxElement(TPZGeoElSide &gelside, int lagrangematid);
    
    /// associate an element group index with the computational elements
    // the grouping of elements depends on the type of mesh created
    // in all cases the volumetric elements nucleate groups
    // interface elements that shape connects with a group are incorporated to the group
    // if(hybridizeSquare) then the elements that have flux degrees of freedom will be incorporated
    void AssociateElements(TPZCompMesh *cmesh, TPZVec<int64_t> &elementgroup);

};


#endif /* TPZCreateMultiphysicsSpace_hpp */
