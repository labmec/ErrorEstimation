//
//  TPZGeoPatch.hpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 22/09/25.
//

#ifndef TPZGeoPatch_hpp
#define TPZGeoPatch_hpp

#include <stdio.h>
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzgmesh.h"
#include <set>

class TPZMultiphysicsCompMesh;

/// @brief class destined to define a patch of elements to be used in generating a multi physics window object
/// an additional data item define the connect indexes and values for the value of the hat functions in the partition of unit mesh
class TPZGeoPatch {
    
    std::set<int64_t> fElementIndexes; // list of element indexes in the patch
    std::map<int64_t,STATE> fConnectIndexValues; // list of connect indexes in the patch mesh
    /// @brief geometric node indexes corresponding to the connect indexes of the patch
    /// This datastructure has the same size of fConnectIndexes
    /// If a flux element (internal or boundary) does not contain a fGeoNodes element, then it will not contribute to the patch
    /// This information will be used in order not to expand the patch
    std::set<int64_t> fGeoNodes;
    
    /// coordinate of the center node of the patch
    TPZManVector<REAL,3> fCenter = TPZManVector<REAL,3>(3,0.);
    
    /// area of the patch
    REAL fArea = 0.;
public:
    // Default constructor
    TPZGeoPatch();
    
    // Copy constructor
    TPZGeoPatch(const TPZGeoPatch &other);
    // Move constructor
    TPZGeoPatch(TPZGeoPatch &&other) noexcept;
    
    // Destructor
    ~TPZGeoPatch();
    
    /// empty the datastructure of a geopatch
    void Empty() {
        fElementIndexes.clear();
        fConnectIndexValues.clear();
        fGeoNodes.clear();
    }
    // Assignment operator
    TPZGeoPatch &operator=(const TPZGeoPatch &other);
    
    // Move assignment operator
    TPZGeoPatch &operator=(TPZGeoPatch &&other) noexcept;
    
    /// Sum two patches
    TPZGeoPatch operator+(const TPZGeoPatch &other) const;
    
    TPZGeoPatch &operator+=(const TPZGeoPatch &other);
    
    // Access method for element indexes
    const std::set<int64_t> &ElementIndexes() const {
        return fElementIndexes;
    }
    
    // Access method for connect indexes
    const std::map<int64_t,REAL> &ConnectIndexValues() const {
        return fConnectIndexValues;
    }
    
    const TPZVec<REAL> &Center() const {
        return fCenter;
    }
    
    // Mutable access method for element indexes
    void SetElementIndexes(TPZGeoMesh *gmesh, const std::set<int64_t> &elindices) {
        fElementIndexes = elindices;
        // compute the area and centernode
        ComputeAreaCenter(gmesh, elindices, fCenter, fArea);
    }
    
    // Mutable access method for element indexes
    void AddElementIndexes(TPZGeoMesh *gmesh, const std::set<int64_t> &elindices) {
        // compute the center node and area of the added elements
        REAL addarea = 0.;
        TPZManVector<REAL,3> addcenter(3,0.);
        ComputeAreaCenter(gmesh, elindices, addcenter, addarea);
        fElementIndexes.insert(elindices.begin(),elindices.end());
        // adjust the area and centernode
        if(fArea > 0. || addarea > 0.) {
            REAL sumarea = fArea+addarea;
            for(int i=0; i<3; i++) {
                fCenter[i] = (fCenter[i]*fArea + addcenter[i]*addarea)/sumarea;
            }
            fArea += addarea;
        }
    }
    
    // Access method for connect indexes
    void SetConnectIndexValues(std::map<int64_t,REAL> &connectvalues)  {
        fConnectIndexValues = connectvalues;
    }
    
    
    /// @brief this method will expand the patch datastructure to include the data
    void AddConnect(int64_t connectindex, STATE value, int64_t geonodeindex);
    
    // load the patch values in the patch mesh
    void LoadPatchVelues(TPZCompMesh *cmesh);
    
    // zero the patch values in the patch mesh
    void ZeroPatchValues(TPZCompMesh *cmesh);
    
    /// @brief verify is the connectindex is internal to the patch
    bool isInternalNode(int64_t nodeindex) {
        return (fGeoNodes.find(nodeindex) != fGeoNodes.end());
    }
    
    bool HasBoundary(TPZGeoMesh &gmesh, std::set<int> &matids) {
        for(auto elindex : fElementIndexes) {
            auto el = gmesh.Element(elindex);
            if(!el) continue;
            int elmatid = el->MaterialId();
            if(matids.find(elmatid) != matids.end()) return true;
        }
        return false;
    }
    
    void Print(std::ostream &out = std::cout) const;
    
    void PlotPatch(const std::string &plotname, int step, TPZMultiphysicsCompMesh &cmesh);
    
    void ComputeAreaCenter(TPZGeoMesh *gmesh, const std::set<int64_t> &elindices, TPZVec<REAL> &center, REAL &area);
};


#endif /* TPZGeoPatch_hpp */
