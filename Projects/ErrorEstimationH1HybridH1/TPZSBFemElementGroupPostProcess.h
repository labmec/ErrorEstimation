//
//  TPZSBFemElementGroupPostProcess.hpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 19/09/25.
//

#ifndef TPZSBFemElementGroupPostProcess_hpp
#define TPZSBFemElementGroupPostProcess_hpp

#include <stdio.h>
#include "pzelementgroup.h"
#include "TPZSBFemElementGroup.h"

class TPZSBFemElementGroupPostProcess : public TPZElementGroup {

    
    TPZSBFemElementGroup *fReferred;
    
    bool fHasBoundary = false;
    
    /// right hand side of the bubble function with respect to patch reconstruction
    TPZFMatrix<CSTATE> fRhsBubble;
    
public:
    TPZSBFemElementGroupPostProcess();
    
    TPZSBFemElementGroupPostProcess(TPZCompMesh &mesh) : TPZElementGroup(mesh)
    {
        
    }
    
    /** @brief create a copy of the element group in the other mesh */
    TPZSBFemElementGroupPostProcess(TPZCompMesh &mesh, const TPZElementGroup &copy) : TPZElementGroup(mesh,copy) {
        std::cout << __PRETTY_FUNCTION__ << " please implement me\n";
        DebugStop();
    }
    
    TPZSBFemElementGroupPostProcess(TPZCompMesh &mesh, TPZSBFemElementGroup *elgr) : TPZElementGroup(mesh),
        fReferred(elgr) {
        
    }
    
    virtual ~TPZSBFemElementGroupPostProcess();
    
    bool HasBoundary() const{
        return fHasBoundary;
    }
    void SetHasBoundary(bool hasbdry){
        fHasBoundary = hasbdry;
    }
    /**
     * @brief Computes the element stifness matrix and right hand side
     * @param ek element stiffness matrix
     * @param ef element load vector
     */
    virtual void CalcStiff(TPZElementMatrixT<CSTATE> &ek,TPZElementMatrixT<CSTATE> &ef) override{
        DebugStop();
    }
    /**
     * @brief Computes the element stifness matrix and right hand side
     * @param ek element stiffness matrix
     * @param ef element load vector
     */
    virtual void CalcStiff(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef) override;

    /** @brief Loads the solution within the internal data structure of the element */
    /**
     * Is used to initialize the solution of connect objects with dependency \n
     * Is also used to load the solution within SuperElements
     */
    virtual void LoadSolution() override;
private:
    
    /// reorganize the connect indexes to correspond to the original sbfem group
    void ReorganizeConnectOrder();
    
    /// cpmpute the right hand side contribution of the hybrid h1 reconstruction for sbfem volume elements
    void ComputeRhs(TPZFMatrix<CSTATE> &rhssbfem, TPZFMatrix<CSTATE> &rhsbubble);
};

#endif /* TPZSBFemElementGroupPostProcess_hpp */
