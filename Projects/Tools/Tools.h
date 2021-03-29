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

#include "TPZHDivErrorEstimator.h"

#include <tuple>
#include <memory>
#include "DataStructure.h"

#include <cstdio>

namespace Tools {

    // Create a geometric mesh on a unit square domain with boundary condition ids defined by bcids
    TPZGeoMesh *CreateGeoMesh(int nelem, TPZVec<int> &bcids);
    TPZGeoMesh *CreateNewGeoMesh(int nelem, TPZVec<int> &bcids);

    TPZGeoMesh *CreateCubeGeoMesh(const TPZVec<int> &nelDiv, const TPZVec<int> &bcids);

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

    void UniformRefinement(int nDiv, TPZGeoMesh *gmesh);

    void UniformRefinement(int nDiv, int dim, TPZGeoMesh *gmesh);

    TPZGeoMesh *CreateTrapezoidalMesh(int nelx, int nely, REAL Lx, REAL Ly, TPZVec<int> &bcids);


/// Divide lower dimensional elements
    void DivideLowerDimensionalElements(TPZGeoMesh *gmesh);

/// numelrefine : number of elements to refine
/// depth : refinement depth of the mesh
    void RandomRefinement(TPZGeoMesh *gmesh, int64_t numelrefine, int depth);

    // Refine elements given a set of indexes
    void RefineElements(TPZGeoMesh *gmesh, const std::set<int64_t>& elsToRefine);

    std::tuple<TPZCompMesh *, TPZVec<TPZCompMesh *> >
    CreatePostProcessingMesh(TPZCompMesh *cmesh_HDiv, TPZVec<TPZCompMesh *> &meshvec_HDiv, TPZHybridizeHDiv &hybridize);


    void FunctionTest();

    void Prefinamento(TPZCompMesh *cmesh, int ndiv, int porder);


    void SolveHybridProblem(TPZCompMesh *Hybridmesh, std::pair<int, int> InterfaceMatId, const ProblemConfig &problem,
                            bool PostProcessingFEM);

    void SolveMixedProblem(TPZCompMesh *cmesh_HDiv, const ProblemConfig &config);

    TPZCompMesh *CMeshH1(ProblemConfig problem);

    void hAdaptivity(TPZCompMesh *postProcessMesh, TPZGeoMesh *gmeshToRefine, ProblemConfig &config);

    void Print(const FADREAL& a, std::ostream& out);

    void Print(const FADFADREAL& a, std::ostream& out);

    void DrawGeoMesh(ProblemConfig &config, PreConfig &preConfig);

    void DrawCompMesh(ProblemConfig &config, PreConfig &preConfig, TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *multiCmesh);

    TPZGeoMesh* CreateGeoMesh(int nelem, TPZVec<int>& bcids,int dim, bool isOriginCentered, int topologyMode);
}

#endif
