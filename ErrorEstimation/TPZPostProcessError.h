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

#include <TPZMultiphysicsCompMesh.h>
#include "pzmanvector.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzblock.h"
#include "pzfunction.h"
#include "TPZAnalyticSolution.h"
#include "ProblemConfig.h"

struct TPZPatch
{
    // connect index of the partition of unity mesh
    TPZManVector<std::pair<int64_t,REAL> > fPartitionConnectIndices;
    // location of the partition connect
    TPZManVector<REAL,3> fCo;
    // vector of element indices of multiphysics elements
    TPZManVector<int64_t,20> fElIndices;
    // vector of open set of connect indices that will be used for flux and pressure computations
    TPZManVector<int64_t,25> fConnectIndices;
    
    // vector of closed set of connect indexes included in the elements
    TPZManVector<int64_t,30> fBoundaryConnectIndices;
    
    bool fPatchIsBoundary;
    
    void ClosedSet(std::set<int64_t> &closed)
    {
//        std::copy (bar.begin(),bar.end(),std::inserter(foo,it));
        std::copy(&(fConnectIndices[0]),(&(fConnectIndices[0])+fConnectIndices.size()),std::inserter(closed,closed.begin()));
    }

    TPZPatch() : fPartitionConnectIndices(1, std::make_pair(-1, 1.)), fCo(3,-1.)
    {
        
    }
    
    TPZPatch(const TPZPatch &copy) : fPartitionConnectIndices(copy.fPartitionConnectIndices), fCo(copy.fCo), fElIndices(copy.fElIndices),
    fConnectIndices(copy.fConnectIndices), fBoundaryConnectIndices(copy.fBoundaryConnectIndices), fPatchIsBoundary(copy.fPatchIsBoundary)
    {
        
    }
    
    TPZPatch &operator=(const TPZPatch &copy)
    {
        fPartitionConnectIndices = copy.fPartitionConnectIndices;
        fCo = copy.fCo;
        fElIndices = copy.fElIndices;
        fConnectIndices = copy.fConnectIndices;
        fBoundaryConnectIndices = copy.fBoundaryConnectIndices;
        fPatchIsBoundary = copy.fPatchIsBoundary;
        return *this;
    }
    
    void Print(std::ostream &out)
    {
        out << "The generating partitionindex = " << fPartitionConnectIndices<< std::endl;
        out << "Coordinate of the partition node " << fCo << std::endl;
        out << "Patch is boundary " << fPatchIsBoundary << std::endl;
        out << "Element indices " << fElIndices << std::endl;
        out << "Open set connect indices " << fConnectIndices << std::endl;
        out << "Boundary set connect indices " << fBoundaryConnectIndices << std::endl;
    }
    
    // return the first equation associated with a lagrange multiplier
    int64_t FirstLagrangeEquation(TPZCompMesh *cmesh) const;
    
};

enum MMeshPositions {Emulti = 0, Eflux = 1, Epressure = 2, Epatch = 3, Eorigin = 4, Epressureaverage = 5};

class TPZPostProcessError
{
public:

    TPZPostProcessError(TPZCompMesh * origin);

    TPZPostProcessError(TPZCompMesh * origin,ProblemConfig &config, bool useHDiv);
    
    TPZPostProcessError(TPZVec<TPZCompMesh *> &meshvec);
    
    

protected:
    // mesh vector
    TPZManVector<TPZCompMesh *,6> fMeshVector;
    
    // vector of vector of patches
    // each vector of patches corresponds to one color
    TPZManVector<TPZStack<TPZPatch>, 10> fVecVecPatches;
    
    /// use HDiv or hybrid H1 to construct a conservative approximation
    bool fuseHDiv = true;
    
    int fNState = 1;
    
    /// @brief material ids associated with error computation
    std::set<int> fMaterialIds;
    /// material ids for building the hybrid H1 mesh
    int fMatWrap = 10;
    int fInterfacePositive = 11;
    int fInterfaceNegative = 12;
    int fMatFlux = 15;
    
    // build vector of patches of a same color
    void BuildPatchStructures();
    virtual void BuildPatchStructures2();//one patch by color

    
    // print the relevant information of the patches
    virtual void PrintPatchInformation(std::ostream &out);
    
    // original connect sequence numbers
    TPZVec<int64_t> fConnectSeqNumbers;
    
    // multiplying coefficients of the reconstructed fluxes and pressures
    TPZFMatrix<STATE> fSolution;
    
    // block corresponding to the original connect sequence numbers
    TPZBlock fBlock;
    
    /// size of the connects in the multiphysics mesh
    TPZVec<int64_t> fConnectSizes;
    
    // plot the reconstructed fluxes
    void PlotFluxes(const std::string &filename);
    
    // solve for the reconstructed fluxes of a given color. Add the flux coefficients
    void ComputePatchFluxes();// not implemented
    
    // determine if a given patch is boundary or not
    bool PatchHasBoundary(TPZPatch &patch, const std::set<int64_t> &internalconnects) const;
    
    // Sum the solution stored in fSolution of the multiphysics mesh to the fSolution vector
    void TransferAndSumSolution(TPZCompMesh *cmesh); // what is second mesh?

    // Reset the state of the HDiv mesh to its original structure
    void ResetState();

    // check whether the connectsizes have changed
    void CheckConnectSizes();
    
    /// identify the material ids of the boundary conditions in the root mesh
    std::set<int> BCMaterialIds() const;

    // create the meshes that allow us to compute the error estimate
    void CreateMultiphysicsMesh();
    
    /// add geometric wrappers, interface and interface flux elements
    void AddWrapperElements();

    /// create a fluxmesh based on the original H1 mesh
    // the flux mesh will be put in position EFlux of the mesh vector
    void CreateFluxMesh();
    
    /// create a boundary flux mesh based on the original H1 mesh
    /// the boundary flux mesh will be put in position EFlux
    void CreateBoundaryFluxMesh();
    
    /// create the lagrange mesh corresponding to the flux mesh
    void CreatePressureMesh();
    
    /// create the hybrid H1 mesh corresponding to the H1 mesh
    void CreateDiscontinuousPressureMesh();

    virtual void CreateAveragePressureMesh();
    
    /// create the multiphysics mesh combining hdiv elements that will compute the projection matrix
    void CreateMultiphysicsHdivMesh();

    /// create the multiphysics mesh using hybrid H1 mesh for reconstruction
    void CreateMultiphysicsHybridH1Mesh();

    /// Add the Interface elements to the multiphysics mesh
    virtual void AddInterfaceElements(TPZMultiphysicsCompMesh *mfmesh);

    /// create the partition of unity mesh
    void CreatePartitionofUnityMesh();
    
    TPZAnalyticSolution *fExact;
    //TPZAutoPointer<TLaplaceExample1> fExact;

public:
    
    // print partition diagnostics
    void PrintPartitionDiagnostics(int64_t color, std::ostream &out) const ;
    
    // include the wrap, interface and flux element in the connected element list
    void IncludeDim1Neighbours(int64_t seednodeindex, TPZGeoEl *gel, std::set<TPZCompEl *> &patchwrappers);
    // Collect the connect indices and elements which will contribute to the patch caracterized by the set of nodes
    // generally each node will form a patch
    TPZPatch BuildPatch(TPZCompElSide &seed);

    // compute the estimated H1 seminorm errors
    void ComputeHDivSolution();//Not used
    
    // compute the estimated H1 seminorm errors
    virtual void ComputeElementErrors(TPZVec<STATE> &elementerrors);
    
    // compute the exact element errors
    void ComputeExactH1SemiNormErrors(TPZFunction<STATE> &exact, TPZVec<STATE> &exacterror)
    {
        DebugStop();
    }
    
    TPZMultiphysicsCompMesh *MultiPhysicsMesh()
    {
        return dynamic_cast<TPZMultiphysicsCompMesh*>(fMeshVector[Emulti]) ;
    }
    
    void SetAnalyticSolution(TPZAnalyticSolution &exact)
    {
        fExact = &exact;
    }
    
};

#endif /* TPZPostProcessError_hpp */
