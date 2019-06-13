//
// Created by Gustavo on 02/04/19.
//
// This program implements Ainsworth proposal for an a posteriori error estimator as seen in:
// M. Ainsworth, X. Ma, Non-uniform order mixed FEM approximation: Implementation, post-processing,
// computable error bound and adaptivity. J. Comput. Physc. 231(2) (2012) 436-453.
//

#include <iostream>

#include "TPZGmshReader.h"
#include "TPZAnalyticSolution.h"
#include "TPZMultiphysicsCompMesh.h"
#include "mixedpoisson.h"
#include "TPZVecL2.h"
#include "TPZCompMeshTools.h"
#include "pzbndcond.h"
#include "pzintel.h"

#include "pzanalysis.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "ProblemConfig.h"

// Creates TPZGeoMesh mesh from a Gmsh .msh file
TPZGeoMesh *CreateGeoMesh(const std::string &meshFileName);

// Creates mixed mesh to calculate the solution
TPZMultiphysicsCompMesh *CreateSolutionMixedMesh(ProblemConfig &config);
TPZCompMesh *CreateSolutionFluxMesh(ProblemConfig &config);
TPZCompMesh *CreateSolutionPressureMesh(ProblemConfig &config);

void Solve(TPZAnalysis &analysis);
void PostProcessSolution(TPZAnalysis &analysis);

void UniformRefinement(int nDiv, TPZGeoMesh *gmesh);

int main() {

    ProblemConfig config;

    config.meshFileName = "BasicMesh.msh";
    config.InsertDomainMat("domain", 1, 2);
    config.InsertBCMat("dirichlet", -1, 1);
    config.InsertBCMat("neumann", -2, 1);
    config.CreateGeoMesh();

    config.exactSolution.fExact = TLaplaceExample1::ESinSin;
    config.nDivisions = 2;
    config.pOrder = 2;

    UniformRefinement(config.nDivisions, config.gmesh);

    TPZMultiphysicsCompMesh *solutionMixedMesh = CreateSolutionMixedMesh(config);

    TPZAnalysis an(solutionMixedMesh);
    Solve(an);
    PostProcessSolution(an);

    
    
    
    std::cout << "hell yeah!";
}

TPZMultiphysicsCompMesh *CreateSolutionMixedMesh(ProblemConfig &config) {

    TPZMultiphysicsCompMesh *mixedMesh = new TPZMultiphysicsCompMesh(config.gmesh);

    TPZMaterial *mat = nullptr;
    for (const auto & domainMat: config.DomainMats) {
        int matID = std::get<1>(domainMat);
        TPZMixedPoisson *mixedMat = new TPZMixedPoisson(matID, config.gmesh->Dimension());

        mixedMat->SetForcingFunction(config.exactSolution.ForcingFunction());
        mixedMat->SetForcingFunctionExact(config.exactSolution.Exact());
        mixedMesh->InsertMaterialObject(mixedMat);

        if (!mat) mat = mixedMat;
    }

    for (const auto & bcMat: config.BCMats) {
        int matID = std::get<1>(bcMat);

        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);

        int bctype = 0; // TODO hardcoded bc type here
        if (matID == -2) {
            val2.Zero();
        }

        TPZBndCond *bc = mat->CreateBC(mat, matID, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(config.exactSolution.Exact());

        mixedMesh->InsertMaterialObject(bc);
    }

    mixedMesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

    TPZManVector<TPZCompMesh *, 2> meshvector(2, 0);
    meshvector[0] = CreateSolutionFluxMesh(config);
    meshvector[1] = CreateSolutionPressureMesh(config);

    TPZManVector<int> active(2,1);
    mixedMesh->BuildMultiphysicsSpace(active, meshvector);

    mixedMesh->LoadReferences();

    bool keepMatrix = false;
    TPZCompMeshTools::CreatedCondensedElements(mixedMesh, true, keepMatrix);

    mixedMesh->InitializeBlock();

    return mixedMesh;
}

TPZCompMesh *CreateSolutionFluxMesh(ProblemConfig &config) {

    TPZCompMesh *cmesh = new TPZCompMesh(config.gmesh);

    config.gmesh->ResetReference();

    TPZMaterial *mat = nullptr;

    int dim = config.gmesh->Dimension();

    for (const auto & domainMat: config.DomainMats) {
        int matID = std::get<1>(domainMat);
        TPZVecL2 *mix = new TPZVecL2(matID);
        mix->SetDimension(dim);
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }

    for (const auto & bcMat: config.BCMats) {
        int matID = std::get<1>(bcMat);
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 1.);
        if (!mat) DebugStop();
        TPZBndCond *bc = mat->CreateBC(mat, matID, 0, val1, val2);
        bc->TPZMaterial::SetForcingFunction(config.exactSolution.Exact());
        cmesh->InsertMaterialObject(bc);
    }

    cmesh->SetDefaultOrder(config.pOrder);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmesh->AutoBuild();

    if (config.hdivplus) {

        int64_t nel = cmesh->NElements();

        for (int64_t el = 0; el < nel; el++) {

            TPZCompEl *cel = cmesh->Element(el);
            TPZGeoEl *gel = cel->Reference();

            if (gel->Dimension() == dim) {
                // Increases intern side (of highest dimension)
                int side = gel->NSides() - 1;
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
                intel->SetSideOrder(side, config.pOrder + config.hdivplus);
                intel->SetPreferredOrder(config.pOrder + config.hdivplus); // TODO perguntar pro phil o que isso faz
            }
        }
    }

    //  if(config.pRefine){
    //      Prefinamento(cmesh, config.nDivisions, config.pOrder);
    //  } TODO: do something about this.

    cmesh->InitializeBlock();
    std::ofstream out("fluxsolmesh.txt");
    cmesh->Print(out);

    return cmesh;
}

TPZCompMesh *CreateSolutionPressureMesh(ProblemConfig &config) {

    TPZCompMesh *cmesh = new TPZCompMesh(config.gmesh);

    TPZMaterial *mat = nullptr;

    for (const auto & domainMat: config.DomainMats) {
        int matID = std::get<1>(domainMat);
        TPZMixedPoisson *mix = new TPZMixedPoisson(matID, cmesh->Dimension());
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    // TODO doesnt this need a conditional like if(problem.hdivplus) before increasing pol. order?
    cmesh->SetDefaultOrder(config.pOrder + config.hdivplus);

    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    cmesh->AutoBuild();

    int64_t n_connects = cmesh->NConnects();
    // TODO what is this
    for (int64_t i = 0; i < n_connects; ++i) {
        cmesh->ConnectVec()[i].SetLagrangeMultiplier(1);
    }

//    if(config.pRefine){
//        Prefinamento(cmesh, config.nDivisions, config.porder);
//    } TODO handle this

    return cmesh;
}

void UniformRefinement(int nDiv, TPZGeoMesh *gmesh) {

    TPZManVector<TPZGeoEl*> children;
    for(int division = 0; division < nDiv; division++) {

        int64_t nels = gmesh->NElements();

        for(int64_t elem = 0; elem < nels; elem++) {

            TPZGeoEl * gel = gmesh->ElementVec()[elem];

            if(!gel || gel->HasSubElement()) continue;
            if(gel->Dimension() == 0) continue;
            gel->Divide(children);
        }
    }
}

void Solve(TPZAnalysis &analysis) {

#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(analysis.Mesh());
    strmat.SetNumThreads(0);
    analysis.SetStructuralMatrix(strmat);
#else
    TPZSkylineStructMatrix strmat(analysis.Mesh());
    strmat.SetNumThreads(0);
#endif

    // Solves mixed problem
    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    analysis.SetSolver(*direct);
    delete direct;
    direct = nullptr;
    analysis.Assemble();
    analysis.Solve();
}

void PostProcessSolution(TPZAnalysis &analysis) {
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    vecnames.Push("Flux");
    int dim = analysis.Mesh()->Reference()->Dimension();
    analysis.DefineGraphMesh(dim, scalnames, vecnames, "Original.vtk");
    analysis.PostProcess(0, dim);
}

