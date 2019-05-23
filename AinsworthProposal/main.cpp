//
// Created by Gustavo on 02/04/19.
//
// This program implements Ainsworth proposal for an a posteriori error estimator as seen in:
// M. Ainsworth, X. Ma, Non-uniform order mixed FEM approximation: Implementation, post-processing,
// computable error bound and adaptivity. J. Comput. Physc. 231(2) (2012) 436-453.
//

#include <iostream>

#include "TPZGmshReader.h"
#include "mixedpoisson.h"
#include "TPZCompMeshTools.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzbndcond.h"

#include "ProblemConfig.h"

// Creates TPZGeoMesh mesh from a Gmsh .msh file
TPZGeoMesh *createGeoMesh(const std::string &meshFileName);

// Creates mixed mesh to calculate the solution
TPZMultiphysicsCompMesh *createSolutionMixedMesh(ProblemConfig &config);

int main() {

    ProblemConfig config;
    config.gmesh = createGeoMesh("LMesh.msh");

    TPZMultiphysicsCompMesh *solutionMixedMesh = createSolutionMixedMesh(config);

    std::cout << "blob";
}

TPZGeoMesh *createGeoMesh(const std::string &meshFileName) {

    TPZGeoMesh *gmesh = nullptr;

    TPZGmshReader reader;

    reader.GetDimNamePhysical()[1]["dirichlet"] = -1;
    reader.GetDimNamePhysical()[1]["neumann"] = -2;
    reader.GetDimNamePhysical()[2]["domain"] = 1;

    reader.SetFormatVersion("4.1");

#ifdef MACOSX
    meshfilename = "../" + meshfilename;
    gmesh = reader.GeometricGmshMesh(meshfilename);
    reader.PrintPartitionSummary(std::cout);
#else
    gmesh = reader.GeometricGmshMesh(meshFileName);
#endif
    gmesh->SetDimension(2); // TODO Is this really needed?
    return gmesh;
}

TPZMultiphysicsCompMesh *createSolutionMixedMesh(ProblemConfig &config) {

    TPZMultiphysicsCompMesh *mixedMesh = new TPZMultiphysicsCompMesh(config.gmesh);

    TPZMaterial *mat = nullptr;
    for (auto matID : config.materialIDs) {

        TPZMixedPoisson *mixedMat = new TPZMixedPoisson(matID, mixedMesh->Dimension());

        mixedMat->SetForcingFunction(config.exactSolution.ForcingFunction());
        mixedMat->SetForcingFunctionExact(config.exactSolution.Exact());
        mixedMesh->InsertMaterialObject(mixedMat);

        if (!mat) mat = mixedMat;
    }

    for (auto matID : config.bcMaterialIDs) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);

        int bctype = 0;
        if (matID == -2) {
            val2.Zero();
        }

        TPZBndCond *bc = mat->CreateBC(mat, matID, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(config.exactSolution.Exact());

        mixedMesh->InsertMaterialObject(bc);
    }

    mixedMesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

    std::set<int> matid;
    matid.insert(1);
    matid.insert(-1);

    TPZManVector<TPZCompMesh *> meshvector(2,0);
    meshvector[0] = CreateFluxHDivMesh(config);
    meshvector[1] = CreatePressureMesh(config);

    TPZManVector<int> active(2,1);
    mixedMesh->BuildMultiphysicsSpace(active, meshvector);

    mixedMesh->LoadReferences();

    bool keepMatrix = false;
    TPZCompMeshTools::CreatedCondensedElements(mixedMesh, true, keepMatrix);

    return mixedMesh;
}
