#include "TPZHybridH1ReconstructionBase.h"
#include "ProblemConfig.h"
#include "pzcondensedcompel.h"
#include "pzsubcmesh.h"
#include "TPZMaterial.h"


TPZHybridH1ReconstructionBase::TPZHybridH1ReconstructionBase(EstimatorConfig *pEstimator){
       fOriginal = pEstimator->fOriginal;
       forderFEM_k = pEstimator->fk;
       forderFEM_n = pEstimator->fn;
       fmaterialids = pEstimator->fmaterialids;
       fbcmaterialids = pEstimator->fbcmaterialids;
       fExact = pEstimator->fExact;
       fproblemname = pEstimator->fproblemname;
       fnDivisions = pEstimator->fnDivisions;
       fAdaptivityStep = pEstimator->fAdaptivityStep;

       fMultiphysicsReconstructionMesh = new TPZMultiphysicsCompMesh(fOriginal->Reference());
};

void TPZHybridH1ReconstructionBase::InitializeFolderOutput(){
       std::string foldername = fFolderOutput;
       foldername.pop_back();
       std::string command = "mkdir -p " + foldername;
       system(command.c_str());
}

/// computing the element stifnesses will "automatically" compute the condensed form of the matrices
void TPZHybridH1ReconstructionBase::ComputeElementStiffnesses() {
    std::cout << "Solving local Dirichlet problem " << std::endl;
#ifdef ERRORESTIMATION_DEBUG2

    {
        std::ofstream out("MeshToComputeStiff.txt");
        fMultiphysicsH1reconstructionMesh->Print(out);
    }
#endif
    for (auto cel: fMultiphysicsReconstructionMesh->ElementVec()) {
        if (!cel) continue;
        TPZCondensedCompEl *condense = dynamic_cast<TPZCondensedCompEl *>(cel);
        if (condense) {
            // for Mark proposal ek correspond to local dirichlet problem
            condense->Assemble();
        }
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (subcmesh) {
            subcmesh->Assemble();
        }
#ifdef ERRORESTIMATION_DEBUG
        if(subcmesh && condense)
        {
            DebugStop();
        }
#endif
    }
}

void TPZHybridH1ReconstructionBase::PrintSolutionVTK(TPZAnalysis &an){

    TPZMaterial *mat = fMultiphysicsReconstructionMesh->FindMaterial(*fmaterialids.begin());
    int varindex = -1;
    if (mat) varindex = mat->VariableIndex("PressureFEM");
    if (varindex != -1) {
        TPZStack<std::string> scalnames, vecnames;
        
        FillVTKoutputVariables(scalnames,vecnames);

        int dim = fMultiphysicsReconstructionMesh->Reference()->Dimension();

        std::stringstream out;
        out << fFolderOutput << *fproblemname
            << "_k_" << forderFEM_k << "_n_"
            << forderFEM_n;
        if (fnDivisions != -1) {
            out << "_Ndiv_" << fnDivisions;
        }
        if (fAdaptivityStep != -1) {
            out << "_AdaptivityStep_" << fAdaptivityStep;
        }
        out << ".vtk";

        int res =2;

        if(fOriginal->NEquations()<100){
            res=4;
        }
        an.DefineGraphMesh(dim, scalnames, vecnames, out.str());
        an.PostProcess(res, dim);
    }
    else {
        std::cout << __PRETTY_FUNCTION__ << "\nPost Processing variable not found!\n";
    }
}

TPZVec<REAL>* TPZHybridH1ReconstructionBase::ComputeErrors(TPZLinearAnalysis *an, int numberErrors){
    
    if (fExact) {
        an->SetExact(fExact->ExactSolution());
    }


    auto errorVec = new TPZVec<REAL>;
    int64_t nErrorCols = numberErrors;
    errorVec->resize(nErrorCols);
    errorVec->Fill(0);
    for (int64_t i = 0; i < nErrorCols; i++) {
        (*errorVec)[i] = 0;
    }

    int64_t nelem = fMultiphysicsReconstructionMesh->NElements();
    fMultiphysicsReconstructionMesh->LoadSolution(fMultiphysicsReconstructionMesh->Solution());
    fMultiphysicsReconstructionMesh->ExpandSolution();
    fMultiphysicsReconstructionMesh->ElementSolution().Redim(nelem, nErrorCols-1);
    
    for(int64_t el = 0; el<nelem; el++)
    {
        TPZCompEl *cel = fMultiphysicsReconstructionMesh->Element(el);
        TPZSubCompMesh *subc = dynamic_cast<TPZSubCompMesh *>(cel);
        if(subc)
        {
            int64_t nelsub = subc->NElements();
            subc->ElementSolution().Redim(nelsub, 6);
        }
    }

    bool store=true;
    std::ofstream myDummyOfs;
    an->PostProcessError(*errorVec, store, myDummyOfs);//calculo do erro com sol exata e aprox e armazena no elementsolution

    return errorVec;
}

