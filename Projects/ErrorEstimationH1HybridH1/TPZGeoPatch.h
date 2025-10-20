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
    
    TPZVec<int64_t> fElementIndexes; // list of element indexes in the patch
    TPZManVector<int64_t> fConnectIndexes; // list of connect indexes in the patch mesh
    TPZManVector<STATE> fHatFunctionValues; // list of values of the hat functions in the connect indexes
    /// @brief geometric node indexes corresponding to the connect indexes of the patch
    /// This datastructure has the same size of fConnectIndexes
    /// If a flux element (internal or boundary) does not contain a fGeoNodes element, then it will not contribute to the patch
    /// This information will be used in order not to expand the patch
    std::set<int64_t> fGeoNodes;
public:
    // Default constructor
    TPZGeoPatch();
    
    // Copy constructor
    TPZGeoPatch(const TPZGeoPatch &other);
    // Move constructor
    TPZGeoPatch(TPZGeoPatch &&other) noexcept;
    
    // Destructor
    ~TPZGeoPatch();
    
    // Assignment operator
    TPZGeoPatch &operator=(const TPZGeoPatch &other);
    
    // Move assignment operator
    TPZGeoPatch &operator=(TPZGeoPatch &&other) noexcept;
    
    // Access method for element indexes
    const TPZVec<int64_t> &ElementIndexes() const {
        return fElementIndexes;
    }
    
    // Access method for connect indexes
    const TPZVec<int64_t> &ConnectIndexes() const {
        return fConnectIndexes;
    }
    
    // Access method for hat function values
    const TPZVec<STATE> &HatFunctionValues() const {
        return fHatFunctionValues;
    }
    
    // Mutable access method for element indexes
    TPZVec<int64_t> &ElementIndexes() {
        return fElementIndexes;
    }
    
    // Mutable access method for connect indexes
    TPZVec<int64_t> &ConnectIndexes() {
        return fConnectIndexes;
    }
    
    // Mutable access method for hat function values
    TPZVec<STATE> &HatFunctionValues() {
        return fHatFunctionValues;
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
};


#endif /* TPZGeoPatch_hpp */
