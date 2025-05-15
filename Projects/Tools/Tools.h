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

#include "DarcyFlow/TPZDarcyFlow.h"

#include "pzintel.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZHybridizeHDiv.h"

#include "TPZAnalysis.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZMultiphysicsCompMesh.h"

#include "TPZHDivErrorEstimator.h"

#include <tuple>
#include <memory>
#include "DataStructure.h"
#include "TPZHDivApproxCreator.h"

#include <cstdio>

namespace Tools {

    void PrintGeometry(TPZGeoMesh *gmesh, const std::string &file_name, bool printTXT, bool printVTK);

// Create a geometric mesh on a unit square domain with boundary condition ids defined by bcids
    TPZGeoMesh *CreateGeoMesh(int nelem, TPZVec<int> &bcids, REAL distortion=0);
    TPZGeoMesh *CreateNewGeoMesh(int nelem, TPZVec<int> &bcids);

    TPZGeoMesh *CreateCubeGeoMesh(const TPZVec<int> &nelDiv, const TPZVec<int> &bcids);

    TPZGeoMesh *CreateLCircleGeoMesh();

    TPZGeoMesh *CreateLShapeMesh(TPZVec<int> &bcids);

    TPZGeoMesh *CreateQuadLShapeMesh(TPZVec<int> &bcids);

    TPZCompMesh *CreateFluxHDivMesh(ProblemConfig &problem);

    TPZCompMesh *CreatePressureMesh(ProblemConfig &problem);

    TPZMultiphysicsCompMesh *CreateHDivMesh(ProblemConfig &problem);

    void UniformRefinement(int nDiv, TPZGeoMesh *gmesh);

    void UniformRefinement(int nDiv, int dim, TPZGeoMesh *gmesh);

/// Divide lower dimensional elements
    void DivideLowerDimensionalElements(TPZGeoMesh *gmesh);

/// numelrefine : number of elements to refine
/// depth : refinement depth of the mesh
    void RandomRefinement(TPZGeoMesh *gmesh, int64_t numelrefine, int depth);

    // Refine elements given a set of indexes
    void RefineElements(TPZGeoMesh *gmesh, const std::set<int64_t>& elsToRefine);

    void Prefinamento(TPZCompMesh *cmesh, int ndiv, int porder);
    void PRefinementNew(TPZMultiphysicsCompMesh *&cmesh, ProblemConfig &config, TPZHDivApproxCreator &hdivCreator);


    void SolveHybridProblem(TPZCompMesh *Hybridmesh, std::pair<int, int> InterfaceMatId, ProblemConfig &problem,
                            bool PostProcessingFEM);

    void SolveMixedProblem(TPZCompMesh *cmesh_HDiv, ProblemConfig &config);

    TPZCompMesh *CMeshH1(ProblemConfig problem);

    void hAdaptivity(TPZCompMesh *postProcessMesh, TPZGeoMesh *gmeshToRefine, TPZCompMesh* cmeshToRefine, ProblemConfig &config);

    void Print(const FADREAL& a, std::ostream& out);

    void Print(const FADFADREAL& a, std::ostream& out);

    void DrawGeoMesh(ProblemConfig &config, PreConfig &preConfig);

    void DrawCompMesh(ProblemConfig &config, PreConfig &preConfig, TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *multiCmesh);

    TPZGeoMesh* CreateGeoMesh(int nelem, TPZVec<int>& bcids, int dim, bool isOriginCentered, int topologyMode);

    void PrintErrors(std::ofstream& out, ProblemConfig& config, const TPZVec<REAL>& error_vec);
void PrintElasticityErrors(std::ofstream& out, ProblemConfig& config, const TPZVec<REAL>& error_vec);
void PrintElasticityErrorsFEM(std::ofstream& out, ProblemConfig& config, const TPZVec<REAL>& error_vec);
void EstimateErrorElasticity(ProblemConfig &config, TPZMultiphysicsCompMesh *originalMesh);
}

#endif
