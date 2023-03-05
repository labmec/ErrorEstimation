//
// Created by victor on 16/02/23.
//

#ifndef ERRORESTIMATION_TPZHYBRIDH1RECONSTRUCTIONBASE_H
#define ERRORESTIMATION_TPZHYBRIDH1RECONSTRUCTIONBASE_H

#include "TPZMultiphysicsCompMesh.h"
#include "ProblemConfig.h"
#include "TPZLinearAnalysis.h"

class TPZAnalysis;

class TPZHybridH1ReconstructionBase {
public:
    TPZHybridH1ReconstructionBase() {
        DebugStop();
    }

    TPZHybridH1ReconstructionBase(EstimatorConfig *pEstimator);

    ~TPZHybridH1ReconstructionBase(){
    }

protected:

    TPZMultiphysicsCompMesh *fOriginal;

    TPZMultiphysicsCompMesh *fMultiphysicsReconstructionMesh = NULL;

    std::string fFolderOutput;

    int forderFEM_k;

    int forderFEM_n;

    int fnDivisions;

    int fAdaptivityStep;

    /// set of materialids in the mesh
   std::set<int> fmaterialids;
   /// set of boundary condition material ids
   std::set<int> fbcmaterialids;

   TPZAutoPointer<TLaplaceExample1> fExact;

   /// name identifying the problem
   std::string *fproblemname;

public:

    virtual void PostProcess() =0;

protected:

   TPZHybridH1ReconstructionBase(const TPZHybridH1ReconstructionBase &copy) 
    {
        // this method wont work because multiphysics meshes have no copy constructor (yet)
        DebugStop();
    }

    /// computing the element stifnesses will "automatically" compute the condensed form of the matrices
    void ComputeElementStiffnesses();

    void InitializeFolderOutput();

    virtual void FillVTKoutputVariables(TPZStack<std::string> &scalnames,TPZStack<std::string> &vecnames) = 0;

    void PrintSolutionVTK(TPZAnalysis &an);

    TPZVec<REAL>* ComputeErrors(TPZLinearAnalysis *an, int numberErrors);


};

#endif // ERRORESTIMATION_TPZHYBRIDH1RECONSTRUCTIONBASE_H
