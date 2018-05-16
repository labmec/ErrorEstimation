/**
 * @file Poisson 3D in hexahedra with shock problem
 */
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

#include <tuple>

TPZCompMesh *CreateFluxHDivMesh(const ProblemConfig &problem);
TPZCompMesh *CreatePressureMesh(const ProblemConfig &problem);
TPZCompMesh *CreateHDivMesh(const ProblemConfig &problem, TPZVec<TPZCompMesh *> &meshvec);
void CloneMeshVec(TPZVec<TPZCompMesh *> &meshvec, TPZVec<TPZCompMesh *> &meshvec_clone);
/// Increase the approximation orders of the sides of the flux elements
void IncreaseSideOrders(TPZCompMesh *fluxmesh);

std::tuple<TPZCompMesh *, TPZVec<TPZCompMesh *> > CreatePostProcessingMesh(TPZCompMesh *cmesh_orig, TPZVec<TPZCompMesh *> &meshvec);

int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    ProblemConfig config;
    {
        TPZGmshReader gmsh;
        gmsh.fPZMaterialId[1]["dirichlet"] = -1;
        gmsh.fPZMaterialId[1]["neuman"] = -2;
        gmsh.fPZMaterialId[2]["domain"] = 1;
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        config.bcmaterialids.insert(-2);
        TPZGeoMesh *gmesh = 0;
        gmesh = gmsh.GeometricGmshMesh("BasicMesh.msh");
        config.gmesh = gmesh;
        {
            std::ofstream out("gmesh.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        }
    }
    TPZManVector<TPZCompMesh*, 2> meshvec_orig(2, 0), meshvec(2, 0);
    TPZCompMesh *cmesh_orig = CreateHDivMesh(config, meshvec_orig);
//    {
//        std::ofstream out("cmesh_orig.txt");
//        cmesh_orig->Print(out);
//    }
//    {
//        std::ofstream out("meshvec_orig_flux.txt");
//        meshvec_orig[0]->Print(out);
//    }
//    {
//        std::ofstream out("meshvec_orig_pres.txt");
//        meshvec_orig[1]->Print(out);
//    }

    TPZCompMesh *cmesh;
    tie(cmesh, meshvec) = CreatePostProcessingMesh(cmesh_orig, meshvec_orig);
//    {
//        std::ofstream out("cmesh.txt");
//        cmesh->Print(out);
//    }
//    {
//        std::ofstream out("meshvec_flux.txt");
//        meshvec[0]->Print(out);
//    }
//    {
//        std::ofstream out("meshvec_pres.txt");
//        meshvec[1]->Print(out);
//    }
    return 0;
}

TPZCompMesh *CreatePressureMesh(const ProblemConfig &problem) {
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial *mat = 0;
    for (auto matid : problem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    cmesh->SetDefaultOrder(problem.porder);
    if (problem.hdivmais) {
        cmesh->SetDefaultOrder(problem.porder + 1);
    }
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->AutoBuild();
    int64_t n_connects = cmesh->NConnects();
    for (int64_t i = 0; i < n_connects; ++i) {
        cmesh->ConnectVec()[i].SetLagrangeMultiplier(1);
    }
    return cmesh;
}

TPZCompMesh *CreateFluxHDivMesh(const ProblemConfig &problem) {
    int dim = problem.gmesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial *mat = NULL;
    problem.gmesh->ResetReference();
    for (auto matid : problem.materialids) {
        TPZVecL2 *mix = new TPZVecL2(matid);
        mix->SetDimension(dim);
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        TPZBndCond *bc = mat->CreateBC(mat, matid, 0, val1, val2);
        cmesh->InsertMaterialObject(bc);
    }
    cmesh->SetDefaultOrder(problem.porder);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmesh->AutoBuild();
    if (problem.hdivmais) {
        int64_t nel = cmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if (gel->Dimension() == dim) {
                int side = gel->NSides() - 1;
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
                intel->SetSideOrder(side, problem.porder + 1);
            }
        }
    }
    cmesh->InitializeBlock();
    return cmesh;

}

TPZCompMesh *CreateHDivMesh(const ProblemConfig &problem, TPZVec<TPZCompMesh *> &meshvector) {
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial *mat = NULL;
    for (auto matid : problem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        mix->SetForcingFunctionExact(problem.exact.ForcingFunction());
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        if (matid == -2) {
            bctype = 1;
        }
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunctionExact(problem.exact.Exact());
        cmesh->InsertMaterialObject(bc);
    }
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    cmesh->AutoBuild();

    meshvector[0] = CreateFluxHDivMesh(problem);
    meshvector[1] = CreatePressureMesh(problem);
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh);
    cmesh->LoadReferences();
    bool keepmatrix = false;
    //TPZCompMeshTools::CreatedCondensedElements(cmesh, true, keepmatrix);

    return cmesh;
}

void CloneMeshVec(TPZVec<TPZCompMesh *> &meshvec, TPZVec<TPZCompMesh *> &meshvec_clone) {
    for (int i = 0; i < meshvec.size(); i++) {
        meshvec_clone[i] = meshvec[i]->Clone();
    }
    IncreaseSideOrders(meshvec_clone[0]);
}

/// Increase the approximation orders of the sides of the flux elements

void IncreaseSideOrders(TPZCompMesh *fluxmesh) {
    int64_t nel = fluxmesh->NElements();
    int dim = fluxmesh->Dimension();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fluxmesh->Element(el);
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != dim) {
            continue;
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        int nc = cel->NConnects();
        int order = cel->Connect(nc - 1).Order();
        int nsides = gel->NSides();
        int ncorner = gel->NCornerNodes();
        for (int side = ncorner; side < nsides - 1; side++) {
            if (intel->NSideConnects(side)) {
                intel->SetSideOrder(side, order);
            }
        }
    }
    fluxmesh->InitializeBlock();
}

std::tuple<TPZCompMesh *, TPZVec<TPZCompMesh *> > CreatePostProcessingMesh(TPZCompMesh *cmesh_orig, TPZVec<TPZCompMesh *> &meshvec_orig) {
    TPZGeoMesh *gmesh = cmesh_orig->Reference();
    TPZManVector<TPZCompMesh *, 2> meshvec(2, 0);
    CloneMeshVec(meshvec_orig, meshvec);
    //   IncreaseSideOrders(meshvec[0]);
    TPZHybridizeHDiv hybridizer;
    /// insert the material objects for HDivWrap, LagrangeInterface and InterfaceMatid
    hybridizer.InsertPeriferalMaterialObjects(meshvec);


    hybridizer.HybridizeInternalSides(meshvec);
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh_orig->CopyMaterials(*cmesh);
    hybridizer.CreateMultiphysicsMesh(cmesh, meshvec);
    hybridizer.CreateInterfaceElements(cmesh, meshvec);
    hybridizer.GroupElements(cmesh);
    return std::make_tuple(cmesh, meshvec);
}

