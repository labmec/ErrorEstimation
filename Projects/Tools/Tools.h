//
//  Tools.h
//  ErrorEstimation
//
//  Created by Denise De Siqueira on 28/03/19.
//

#ifndef Tools_h
#define Tools_h

#include "pzlog.h"
#include "TPZRefPatternDataBase.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"
#include "ProblemConfig.h"

#include "mixedpoisson.h"
#include "TPZVecL2.h"
#include "pzbndcond.h"

#include "pzintel.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZHybridizeHDiv.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZMultiphysicsCompMesh.h"

#include "TPZHybridHDivErrorEstimator.h"

class TPZMultiphysicsCompMesh;

#include <tuple>
#include <memory>

#include <cstdio>

namespace Tools {

    // Create a geometric mesh on a unit square domain with boundary condition ids defined by bcids
    TPZGeoMesh *CreateGeoMesh(int nelem, TPZVec<int> &bcids);

    TPZGeoMesh *CreateLCircleGeoMesh();

    TPZGeoMesh *CreateLShapeMesh(TPZVec<int> &bcids);

    TPZGeoMesh *CreateQuadLShapeMesh(TPZVec<int> &bcids);

    TPZGeoMesh *CreateSingleTriangleMesh(TPZVec<int> &bcids);

    TPZGeoMesh *CreateSingleQuadMesh(TPZVec<int> &bcids);

    TPZGeoMesh *CreateQuadMeshRefTriang(TPZVec<int> &bcids);

    TPZGeoMesh *ReadGeometricMesh(struct ProblemConfig &config, bool IsgmeshReader);


    TPZCompMesh *CreateFluxHDivMesh(const ProblemConfig &problem);

    TPZCompMesh *CreatePressureMesh(const ProblemConfig &problem);

    TPZMultiphysicsCompMesh *CreateHDivMesh(const ProblemConfig &problem);

    void CloneMeshVec(TPZVec<TPZCompMesh *> &meshvec, TPZVec<TPZCompMesh *> &meshvec_clone);

/// Increase the approximation orders of the sides of the flux elements
    void IncreaseSideOrders(TPZCompMesh *fluxmesh);

/// Set the interface pressure to the average pressure
    void ComputeAveragePressure(TPZCompMesh *pressure, TPZCompMesh *pressureHybrid, int InterfaceMatid);

    void UniformRefinement(int nDiv, TPZGeoMesh *gmesh);

    void UniformRefinement(int nDiv, int dim, TPZGeoMesh *gmesh);

    TPZGeoMesh *CreateTrapezoidalMesh(int nelx, int nely, REAL Lx, REAL Ly, TPZVec<int> &bcids);


/// Divide lower dimensional elements
    void DivideLowerDimensionalElements(TPZGeoMesh *gmesh);

    void MultiPhysicsCompel(const ProblemConfig &config);

/// numelrefine : number of elements to refine
/// depth : refinement depth of the mesh
    void RandomRefinement(TPZGeoMesh *gmesh, int64_t numelrefine, int depth);

    // Refine elements given a set of indexes
    void RefineElements(TPZGeoMesh *gmesh, std::set<int64_t> elsToRefine);

    std::tuple<TPZCompMesh *, TPZVec<TPZCompMesh *> >
    CreatePostProcessingMesh(TPZCompMesh *cmesh_HDiv, TPZVec<TPZCompMesh *> &meshvec_HDiv, TPZHybridizeHDiv &hybridize);

    void PrintSolAndDerivate(const ProblemConfig& config);

    void FunctionTest();

    void MultiPhysicsHybrid(const ProblemConfig &config);

    void Prefinamento(TPZCompMesh *cmesh, int ndiv, int porder);


    void SolveHybridProblem(TPZCompMesh *Hybridmesh, std::pair<int, int> InterfaceMatId, const ProblemConfig &problem,
                            bool PostProcessingFEM);

    void SolveMixedProblem(TPZCompMesh *cmesh_HDiv, const ProblemConfig &config);

    void PlotLagrangeMultiplier(TPZCompMesh *cmesh, const ProblemConfig &problem);


    TPZMultiphysicsCompMesh *HybridSolveProblem(TPZMultiphysicsCompMesh *cmesh_HDiv, struct ProblemConfig &config);

    TPZCompMesh *CMeshH1(ProblemConfig problem);

    void hAdaptivity(TPZCompMesh *postProcessMesh, TPZGeoMesh *gmeshToRefine, ProblemConfig &config);

    void ComputeError(TPZCompMesh *Hybridmesh, std::ofstream &out, const ProblemConfig &config);

    void VectorEnergyNorm(TPZCompMesh *hdivmesh, std::ostream &out, const ProblemConfig &problem);

    void Print(const FADREAL& a, std::ostream& out);

    void Print(const FADFADREAL& a, std::ostream& out);


}

#endif
