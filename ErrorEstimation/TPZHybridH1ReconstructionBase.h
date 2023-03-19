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

    ~TPZHybridH1ReconstructionBase();

protected:

    TPZMultiphysicsCompMesh *fOriginal;

    TPZMultiphysicsCompMesh *fMultiphysicsReconstructionMesh = NULL;

    std::string fFolderOutput;

    int forderFEM_k;

    int forderFEM_n;

    int fnDivisions;

    int fAdaptivityStep;

    int fvtkResolution = 0;

    /// set of materialids in the mesh
   std::set<int> fmaterialids;
   /// set of boundary condition material ids
   std::set<int> fbcmaterialids;
    
    /// material id of the skeleton elements
    int fSkeletonMatId;

   TPZAutoPointer<TLaplaceExample1> fExact;

   /// name identifying the problem
   std::string *fproblemname;

public:

    [[maybe_unused]] virtual TPZVec<REAL> PostProcess() = 0;

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

    TPZVec<REAL> ComputeErrors(TPZLinearAnalysis *an, int numberErrors);

    /**
     * File management: Introduce the approximation errors corresponding to the current adaptive step
     * besides the errors from previous steps
     * @param errorVec : Approximation errors corresponding to the current adaptivity step;
     * @param errorVec : Complementary data that should be flushed (e.g. DOF, adapativity step,...)  ;
     * @param filepath : File path leading to "filename". Should end with a "/"  character;
     * @param filename : Name of the file containing errors. 
     **/ 
    static void FlushErrorDataIntoFile(const TPZVec<REAL> &errorVec,const TPZVec<int> &complementaryVec,
                                        const std::string &filePath,const std::string &filename);

};

#endif // ERRORESTIMATION_TPZHYBRIDH1RECONSTRUCTIONBASE_H
