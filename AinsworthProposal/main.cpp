// This program implements Ainsworth proposal for an a posteriori error estimator as seen in:
// M. Ainsworth, X. Ma, Non-uniform order mixed FEM approximation: Implementation, post-processing,
// computable error bound and adaptivity. J. Comput. Physc. 231(2) (2012) 436-453.


#include <iostream>


#include "TPZMultiphysicsCompMesh.h"
#include "TPZGmshReader.h"
#include "mixedpoisson.h"

TPZGeoMesh* createGeoMesh();

TPZCompMesh* createSolutionMesh(TPZGeoMesh* gmesh);

int main() {

    TPZGeoMesh* gmesh = createGeoMesh();
    TPZMultiphysicsCompMesh mesh_test(gmesh);

    std::cout << "blob";
}

TPZGeoMesh* createGeoMesh() {

    TPZGeoMesh *gmesh = nullptr;

    std::string meshfilename = "LMesh.msh";

    TPZGmshReader gmsh;

    gmsh.GetDimNamePhysical()[1]["dirichlet"] = -1;
    gmsh.GetDimNamePhysical()[1]["neumann"] = -2;
    gmsh.GetDimNamePhysical()[2]["domain"] = 1;

    gmsh.SetFormatVersion("4.1");

#ifdef MACOSX
    meshfilename = "../" + meshfilename;
    gmesh = gmsh.GeometricGmshMesh(meshfilename);
    gmsh.PrintPartitionSummary(std::cout);
#else
    gmesh = gmsh.GeometricGmshMesh(meshfilename);
#endif
    gmesh->SetDimension(2); // TODO Is this really needed?
    return gmesh;
}

TPZCompMesh *CreateSolutionMesh(TPZGeoMesh* gmesh) {

    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    TPZMaterial *mat = nullptr;

    //for (auto matid : problem.materialids) { // TODO multiple mat ids
    int matid = 1;
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        TPZAutoPointer<TPZFunction<STATE>> force1 = new TPZDummyFunction<STATE>(Forcing, 0);
        mix->SetForcingFunction(force1);
        mix->SetInternalFlux(1);
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    //}
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        if (matid == -2) {
            bctype = 0;
            val2.Zero();
        }
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        //   bc->TPZMaterial::SetForcingFunction(problem.exact.Exact());
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
    TPZCompMeshTools::CreatedCondensedElements(cmesh, true, keepmatrix);

    return cmesh;
}
