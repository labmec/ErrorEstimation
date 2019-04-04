//
//  TPZPostProcessError.hpp
//  PZ
//
//  Created by Philippe Devloo on 6/30/16.
//
//

#ifndef TPZPostProcessError_hpp
#define TPZPostProcessError_hpp

#include <stdio.h>
#include <iterator>

#include "pzmanvector.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzblock.h"
#include "pzfunction.h"

struct TPZPatch
{
    // connect index of the partition of unity mesh
    int64_t fPartitionConnectIndex;
    // location of the partition connect
    TPZManVector<REAL,3> fCo;
    // vector of element indices of HDiv elements
    TPZManVector<int64_t,20> fElIndices;
    // vector of open set of connect indices that will be used for flux and pressure computations
    TPZManVector<int64_t,25> fConnectIndices;
    
    // vector of closed set of connect indexes included in the elements
    TPZManVector<int64_t,30> fBoundaryConnectIndices;
    
    void ClosedSet(std::set<int64_t> &closed)
    {
//        std::copy (bar.begin(),bar.end(),std::inserter(foo,it));
        std::copy(&(fConnectIndices[0]),(&(fConnectIndices[0])+fConnectIndices.size()),std::inserter(closed,closed.begin()));
    }
    
    TPZPatch() : fPartitionConnectIndex(-1), fCo(3,-1.)
    {
        
    }
    
    TPZPatch(const TPZPatch &copy) : fPartitionConnectIndex(copy.fPartitionConnectIndex), fCo(copy.fCo), fElIndices(copy.fElIndices),
    fConnectIndices(copy.fConnectIndices), fBoundaryConnectIndices(copy.fBoundaryConnectIndices)
    {
        
    }
    TPZPatch &operator=(const TPZPatch &copy)
    {
        fPartitionConnectIndex = copy.fPartitionConnectIndex;
        fCo = copy.fCo;
        fElIndices = copy.fElIndices;
        fConnectIndices = copy.fConnectIndices;
        fBoundaryConnectIndices = copy.fBoundaryConnectIndices;
        return *this;
    }
    
    void Print(std::ostream &out)
    {
        out << "The generating partitionindex = " << fPartitionConnectIndex << std::endl;
        out << "Coordinate of the partition node " << fCo << std::endl;
        out << "Element indices " << fElIndices << std::endl;
        out << "Open set connect indices " << fConnectIndices << std::endl;
        out << "Boundary set connect indices " << fBoundaryConnectIndices << std::endl;
    }
    
    // return the first equation associated with a lagrange multiplier
    int64_t FirstLagrangeEquation(TPZCompMesh *cmesh) const;
    

};

class TPZPostProcessError
{
public:
    
    TPZPostProcessError(TPZCompMesh * origin);
    
    TPZPostProcessError(TPZVec<TPZCompMesh *> &meshvec);
    
private:
    
    // mesh vector
    TPZManVector<TPZCompMesh *,6> fMeshVector;
    
    // vector of vector of patches
    // each vector of patches corresponds to one color
    TPZManVector<TPZStack<TPZPatch>, 10> fVecVecPatches;
    
    // build vector of patches of a same color
    void BuildPatchStructures();
    
    // print the relevant information of the patches
    void PrintPatchInformation(std::ostream &out);
    
    // original connect sequence numbers
    TPZVec<int64_t> fConnectSeqNumbers;
    
    // multiplying coefficients of the reconstructed fluxes and pressures
    TPZFMatrix<STATE> fSolution;
    
    // block corresponding to the original connect sequence numbers
    TPZBlock<STATE> fBlock;
    
    /// size of the connects in the multiphysics mesh
    TPZVec<int64_t> fConnectSizes;
    
    // plot the reconstructed fluxes
    void PlotFluxes(const std::string &filename);
    
    // solve for the reconstructed fluxes of a given color. Add the flux coefficients
    void ComputePatchFluxes();
    
    // determine if a given patch is boundary or not
    bool PatchHasBoundary(TPZPatch &patch) const;
    
    // Sum the solution stored in fSolution of the second mesh to the fSolution vector
    void TransferAndSumSolution(TPZCompMesh *cmesh);
    
    // Reset the state of the HDiv mesh to its original structure
    void ResetState();
    
    // check whether the connectsizes have changed
    void CheckConnectSizes();
    
    // create the meshes that allow us to compute the error estimate
    void CreateAuxiliaryMeshes();

    /// create a fluxmesh based on the original H1 mesh
    // the flux mesh will be put in the second position of the mesh vector
    void CreateFluxMesh();
    
    /// create the lagrange mesh corresponding to the flux mesh
    void CreatePressureMesh();

    /// create the multiphysics mesh that will compute the projection matrix
    void CreateMixedMesh();

    /// create the partition of unity mesh
    void CreatePartitionofUnityMesh();
    
public:
    
    // print partition diagnostics
    void PrintPartitionDiagnostics(int color, std::ostream &out) const ;
    
    // Collect the connect indices and elements which will contribute to the patch caracterized by the set of nodes
    // generally each node will form a patch
    TPZPatch BuildPatch(TPZCompElSide &seed);

    // compute the estimated H1 seminorm errors
    void ComputeHDivSolution();
    
    // compute the estimated H1 seminorm errors
    void ComputeElementErrors(TPZVec<STATE> &elementerrors);
    
    // compute the exact element errors
    void ComputeExactH1SemiNormErrors(TPZFunction<STATE> &exact, TPZVec<STATE> &exacterror)
    {
        DebugStop();
    }
    
    TPZCompMesh *MultiPhysicsMesh()
    {
        return fMeshVector[1];
    }
    
};

#endif /* TPZPostProcessError_hpp */
