#include "TPZHybridH1ReconstructionBase.h"
#include "ProblemConfig.h"
#include "pzcondensedcompel.h"
#include "pzsubcmesh.h"
#include "TPZMaterial.h"

#include <iostream>
#include <filesystem>


TPZHybridH1ReconstructionBase::TPZHybridH1ReconstructionBase(EstimatorConfig *pEstimator){
       fOriginal = pEstimator->fOriginal;
       forderFEM_k = pEstimator->fk;
       forderFEM_n = pEstimator->fn;
       fmaterialids = pEstimator->fmaterialids;
       fbcmaterialids = pEstimator->fbcmaterialids;
       fSkeletonMatId = pEstimator->fSkeletonMatId;
       fExact = pEstimator->fExact;
       fproblemname = pEstimator->fproblemname;
       fnDivisions = pEstimator->fnDivisions;
       fAdaptivityStep = pEstimator->fAdaptivityStep;
       fvtkResolution = pEstimator->fvtkResolution;
       fProblemFolderOutput = *pEstimator->fproblemname;

       fMultiphysicsReconstructionMesh = new TPZMultiphysicsCompMesh(fOriginal->Reference());
}

TPZHybridH1ReconstructionBase::~TPZHybridH1ReconstructionBase(){
    auto meshvec = fMultiphysicsReconstructionMesh->MeshVector();
    int64_t nmeshes = meshvec.size();
    TPZGeoMesh *gmesh = fMultiphysicsReconstructionMesh->Reference();
    gmesh->ResetReference();
    for(int64_t im=0; im<nmeshes; im++) {
        if(meshvec[im]) {
            int64_t nc = meshvec[im]->ConnectVec().NElements();
            for (int64_t ic = 0; ic<nc; ic++) {
                meshvec[im]->ConnectVec()[ic].RemoveDepend();
            }
            delete meshvec[im];
        }
    }
    delete fMultiphysicsReconstructionMesh;
}


void TPZHybridH1ReconstructionBase::InitializeFolderOutput(){
       std::string foldername = fFolderOutput;
       foldername.pop_back();
       std::string command = "mkdir -p " + foldername;
       system(command.c_str());
}

/// computing the element stifnesses will "automatically" compute the condensed form of the matrices
void TPZHybridH1ReconstructionBase::ComputeElementStiffnesses() {

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
    if (mat) varindex = mat->VariableIndex("uh");
    if (varindex != -1) {
        TPZStack<std::string> scalnames, vecnames;
        
        FillVTKoutputVariables(scalnames,vecnames);

        int dim = fMultiphysicsReconstructionMesh->Reference()->Dimension();

        std::stringstream out;
        out << fProblemFolderOutput << *fproblemname
            << "_k_" << forderFEM_k << "_n_"
            << forderFEM_n;
        if (fnDivisions != -1) {
            out << "_Ndiv_" << fnDivisions;
        }
        if (fAdaptivityStep != -1) {
            out << "_AdaptivityStep_" << fAdaptivityStep;
        }
        out << ".vtk";

        int res = fvtkResolution;
        if(res<0){
            DebugStop();
        }

        an.DefineGraphMesh(dim, scalnames, vecnames, out.str());
        an.PostProcess(res, dim);
    }
    else {
        std::cout << __PRETTY_FUNCTION__ << "\nPost Processing variable not found!\n";
    }
}

TPZVec<REAL> TPZHybridH1ReconstructionBase::ComputeErrors(TPZLinearAnalysis *an, int numberErrors){
    
    if (fExact) {
        an->SetExact(fExact->ExactSolution());
    }


    TPZVec<REAL> errorVec;
    int64_t nErrorCols = numberErrors;
    errorVec.resize(nErrorCols);
    errorVec.Fill(0);
    for (int64_t i = 0; i < nErrorCols; i++) {
        errorVec[i] = 0;
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
    an->PostProcessError(errorVec, store, myDummyOfs);//calculo do erro com sol exata e aprox e armazena no elementsolution

    return errorVec;
}

void TPZHybridH1ReconstructionBase::FlushErrorDataIntoFile(const TPZVec<REAL> &errorVec,const TPZVec<std::string> &complementaryVec,
                                    const std::string &filePath,const std::string &filename){

    int numErrors = errorVec.size();
    int numComplementary = complementaryVec.size();
    int numTargetLines = numErrors+numComplementary;
    std::string fullPath = filePath+filename;

//    if(!std::filesystem::exists(fullPath)){
//        DebugStop();
//    }

    std::ifstream fin(fullPath); 
    const int numFileLines = std::count(std::istreambuf_iterator<char>(fin), 
             std::istreambuf_iterator<char>(), '\n');
    fin.close();

    // Filling a list with the string content of every line;
    std::string outFileName = filePath + "ErrorReconstructionTemp.txt";
    std::ofstream fout(outFileName); 
    fin.open(fullPath);
    if(fin.is_open()) {
        int iLine;
        std::string line;
        for(iLine = 0; iLine < numFileLines - numTargetLines; iLine++){
            std::getline(fin,line);                      // Read the current line
            line += "\n";
            fout << line;
        }
        //auto read_pos = fin.tellg();
        int maxLineSize = 0;
        std::string extraChars = "\t&\t";
        std::list<std::pair<std::string,std::string>> lineList;
        std::pair<std::string,std::string> linePair;
        for(iLine =0; iLine < numTargetLines; iLine++){
            std::stringstream ss;
            std::getline(fin,line);                      // Read the current line
            linePair.first = line;
            if(iLine < numComplementary){
                ss << extraChars <<complementaryVec[iLine];
            }else{
                ss << extraChars << errorVec[iLine-numComplementary];
            }
            linePair.second = ss.str();
            lineList.push_back(linePair);
            int lineLength = ss.str().length();
            if(maxLineSize < lineLength){
                maxLineSize = lineLength;
            }
        }
         for(iLine =0; iLine < numTargetLines; iLine++){
            linePair = lineList.front();
            lineList.pop_front();
            int lineSize = linePair.second.length();
            line = linePair.first;
            if(lineSize < maxLineSize){
                int sizediff = maxLineSize-lineSize;
                for(int iWhiteSpace = 0; iWhiteSpace < sizediff;iWhiteSpace++){
                    line += " ";
                }
            }
            line += linePair.second+"\n";
            fout << line;
         }
    }
    fin.close();
    if(std::rename(outFileName.c_str(),fullPath.c_str())){
        std::cout << "outFileName: " << outFileName.c_str()<<std::endl;
        std::cout << "fullPath: " << fullPath.c_str() <<std::endl;
        std::cout << "Adaptivity step: " << complementaryVec[1]<<std::endl;
        DebugStop();
    }
}

void TPZHybridH1ReconstructionBase::InitializeProblemFolderOutput(std::string &problemName, std::string &folderOutput,const int &k, const int &n, const REAL threshold){
       std::stringstream ss;
       ss << "__k-" << k << "__n-" << n;

       int th = (int)(100.*threshold);
       if(threshold!=-1)
            ss << "__tal-"<< th;
       problemName += ss.str();
       std::string command = "mkdir -p " + folderOutput + problemName;
       system(command.c_str());
       problemName = folderOutput + problemName + "/";
}
