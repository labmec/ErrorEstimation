//
//  TPZHybridHDivErrorEstimator.hpp
//  ErrorEstimateHDiv
//
//  Created by Philippe Devloo on 10/06/18.
//

#ifndef TPZMHMHDivErrorEstimator_hpp
#define TPZMHMHDivErrorEstimator_hpp

#include <stdio.h>
#include "TPZHybridHDivErrorEstimator.h"
#include "TPZMHMixedMeshControl.h"

class TPZCompMesh;

/// this class will compute the estimated value of the energy error of an MHM mesh
// the class should work for any MHM-H(div) mesh
// first create the post processing mesh
// then compute the errors
// compute a reference solution without MHM
struct TPZMHMHDivErrorEstimator : public TPZHybridHDivErrorEstimator
{
    
    /// a pointer to the datastructure used to generate the MHM mesh
    TPZMHMixedMeshControl *fMHM = 0;
    
    TPZMHMHDivErrorEstimator(TPZMultiphysicsCompMesh &InputMesh, TPZMHMixedMeshControl *mhm) : TPZHybridHDivErrorEstimator(InputMesh,true),
        fMHM(mhm)
    {
        fProblemConfig.makepressurecontinuous = true;
    }
    
    TPZMHMHDivErrorEstimator(const TPZMHMHDivErrorEstimator &copy) : TPZHybridHDivErrorEstimator(copy), fMHM(copy.fMHM)
    {
        // this method wont work because multiphysics meshes have no copy constructor (yet)
        DebugStop();
    }
    
    TPZMHMHDivErrorEstimator &operator=(const TPZMHMHDivErrorEstimator &cp)
    {
        TPZHybridHDivErrorEstimator::operator=(cp);
        fMHM = cp.fMHM;
        return *this;
    }
    
    virtual ~TPZMHMHDivErrorEstimator()
    {
        
    }

    /// compute the element errors comparing the reconstructed solution based on average pressures
    /// with the original solution
    virtual void ComputeErrors(TPZVec<REAL> &elementerrors, bool store = true) override;
    

    // a method for generating the HDiv mesh
    virtual TPZCompMesh *CreateFluxMesh() override;
    // a method for creating the pressure mesh
    virtual TPZCompMesh *CreatePressureMesh() override;
    // a method for generating the hybridized multiphysics post processing mesh
    virtual void CreatePostProcessingMesh() override;
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
    virtual void ComputeAveragePressures(int target_dim) override;
    /// compute the average pressure over corners
    /// set the cornernode values equal to the averages
    virtual void ComputeNodalAverages() override;

    // a method for computing a reference solution
};

#endif /* TPZHybridHDivErrorEstimator_hpp */
