//
// Created by gus on 14/03/19.
//

#ifndef TPZMARKERRORESTIMATION_H
#define TPZMARKERRORESTIMATION_H

#include "ProblemConfig.h"

#include "mixedpoisson.h"
#include "TPZVecL2.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZCompMeshTools.h"
#include "pzintel.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzskylstrmatrix.h"

class TPZMarkErrorEstimation {

private:
    ProblemConfig fProblem;

    TPZManVector<TPZCompMesh *,6> fMeshVector;

    void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

public:
    TPZMarkErrorEstimation(ProblemConfig problem);

    // Creates Hdiv x L2 meshes to solve the mixed problem.
    void CreateMultiphysicsMesh();
    void CreateSolutionFluxMesh();
    void CreateSolutionPressureMesh();

    // Solves mixed problem
    void SolveMixedProblem();

    //
    void CreatePostProcessMeshes();
};


#endif //TPZMARKERRORESTIMATION_H
