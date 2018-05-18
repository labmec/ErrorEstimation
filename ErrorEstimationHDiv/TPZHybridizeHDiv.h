//
//  TPZHybridizeHDiv.hpp
//  ErrorEstimate
//
//  Created by Philippe Devloo on 16/05/18.
//

#ifndef TPZHybridizeHDiv_hpp
#define TPZHybridizeHDiv_hpp

#include <stdio.h>
#include <map>

template<class T>
class TPZVec;

class TPZCompMesh;
class TPZMaterial;
class TPZCompElSide;
class TPZInterpolatedElement;

struct TPZHybridizeHDiv {
    int HDivWrapMatid = 0;
    int LagrangeInterface = 0;
    int InterfaceMatid = 0;
    int NState = 1;
    
    TPZHybridizeHDiv() = default;
    
    TPZHybridizeHDiv(TPZVec<TPZCompMesh *> &meshvec);
    
    void ComputeMatIds(TPZVec<TPZCompMesh*>& meshvec);

    void ComputeNState(TPZVec<TPZCompMesh*>& meshvec);


    /// split the connects between flux elements and create a dim-1 pressure element
    void HybridizeInternalSides(TPZVec<TPZCompMesh *> &meshvec);

    /// Create interface elements with material id InterfaceMatid
    void CreateInterfaceElements(TPZCompMesh *cmesh, TPZVec<TPZCompMesh *> &meshvec);

    /// create a multiphysics mesh using the materials pointed to in the vector
    void CreateMultiphysicsMesh(TPZCompMesh *cmesh, TPZVec<TPZCompMesh *> &meshvec);
    
    /// group and condense the elements
    static void GroupElements(TPZCompMesh *cmesh);
    
    /// insert the material objects for HDivWrap, LagrangeInterface and InterfaceMatid in the atomic meshes
    void InsertPeriferalMaterialObjects(TPZVec<TPZCompMesh *> &meshvec);
    
    /// insert the material objects for HDivWrap, LagrangeInterface and InterfaceMatid in the multiphysics mesh
    void InsertPeriferalMaterialObjects(TPZCompMesh *cmesh);
    
private:
    
    std::tuple<int64_t,int> SplitConnects(const TPZCompElSide &left, const TPZCompElSide &right, TPZVec<TPZCompMesh *> &meshvec);

    static TPZCompElSide RightElement(TPZInterpolatedElement *intel, int side);

};

#endif /* TPZHybridizeHDiv_hpp */
