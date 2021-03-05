//
// Created by Gustavo A. Batistela on 25/07/2020.
//

#include "Tools.h"
#include "ToolsMHM.h"
#include "ToolsUNISIM.h"

#include <TPZMFSolutionTransfer.h>
#include <TPZMHMHDivErrorEstimator.h>

#include <Analysis/pzanalysis.h>
#include <Geom/pzgeoquad.h>
#include <Material/REAL/mixedpoisson.h>
#include <Material/TPZVecL2.h>
#include <Material/pzbndcond.h>
#include <Matrix/pzstepsolver.h>
#include <Mesh/TPZMultiphysicsCompMesh.h>
#include <Mesh/tpzgeoelrefpattern.h>
#include <ProblemConfig.h>
#include <StrMatrix/TPZSSpStructMatrix.h>
#include <StrMatrix/pzstrmatrix.h>
#include <Pre/TPZMHMixedMeshControl.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#define NOPE_FORCING_FUNCTION_DEBUG


void UNISIMMHM(TPZGeoMesh *gmesh);

void CreateMHMCompMesh(TPZMHMixedMeshControl *mhm, const ProblemConfig &config, int nInternalRef,
                       bool definePartitionByCoarseIndex, TPZManVector<int64_t> mhmIndexes);

void SolveMHMProblem(TPZMHMixedMeshControl *mhm, const ProblemConfig &config);

void InsertMaterialsInMHMMesh(TPZMHMixedMeshControl &control, const ProblemConfig &config);

void EstimateError(ProblemConfig &config, TPZMHMixedMeshControl *mhm);

void RemoveGelsOfGivenMaterial(TPZGeoMesh *gmesh, int matId);

int main() {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    gRefDBase.InitializeRefPatterns(2);
    gRefDBase.InitializeAllUniformRefPatterns();

    bool modifyZCoords = false;
    TPZGeoMesh *gmesh = CreateUNISIMSurfaceGeoMesh(modifyZCoords);

    std::string meshFileName{"UNISIMMesh"};
    PrintGeometry(gmesh, meshFileName, false, true);

    RemoveGelsOfGivenMaterial(gmesh, 99);
    meshFileName = "UNISIMMAfterRemovingMaterial";
    PrintGeometry(gmesh, meshFileName, false, true);

    UNISIMMHM(gmesh);

    delete gmesh;
    return 0;
}

void UNISIMMHM(TPZGeoMesh *gmesh) {

    ProblemConfig config;

    config.porder = 1;
    config.hdivmais = 1;
    config.dimension = 2;
    config.makepressurecontinuous = true;

    config.dir_name = "CILAMCE_UNISIM_MHM";
    config.problemname = "UNISIM_MHM";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

#ifdef FORCING_FUNCTION_DEBUG
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::EX;
    config.dir_name = "CILAMCE_UNISIM_MHM_DBEUG";
#endif

    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.bcmaterialids.insert(-2);
    config.bcmaterialids.insert(-3);

    config.gmesh = gmesh;
    TPZVec<int64_t> coarseIndexes;
    ComputeCoarseIndices(gmesh, coarseIndexes);

    int nInternalRef = 0;
    Tools::UniformRefinement(nInternalRef, 2, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    std::stringstream gmeshFileName;
    gmeshFileName << config.dir_name << "/GeoMesh";
    PrintGeometry(gmesh, gmeshFileName.str(), false, true);

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    bool definePartitionByCoarseIndexes = true;
    CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, coarseIndexes);
    {
        std::string fileName = config.dir_name + "/" + config.problemname + "GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    SolveMHMProblem(mhm, config);
    EstimateError(config, mhm);
}

void SolveMHMProblem(TPZMHMixedMeshControl *mhm, const ProblemConfig &config) {

    TPZAutoPointer<TPZCompMesh> cmesh = mhm->CMesh();

    bool shouldrenumber = true;
    TPZAnalysis an(cmesh, shouldrenumber);

#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(0 /*config.n_threads*/);
#else
    TPZSkylineStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(config.n_threads);
#endif

    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "Assembling\n";
    an.Assemble();

    std::cout << "Solving\n";
    an.Solve();
    std::cout << "Finished\n";
    an.LoadSolution(); // compute internal dofs

    TPZMFSolutionTransfer transfer;
    transfer.BuildTransferData(cmesh.operator->());
    transfer.TransferFromMultiphysics();

    TPZStack<std::string> scalnames, vecnames;
    TPZMaterial *mat = cmesh->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }

    scalnames.Push("Pressure");
    vecnames.Push("Flux");

    int resolution = 0;
    std::string plotname = config.dir_name + "/" + config.problemname + "Results.vtk";
    an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotname);
    an.PostProcess(resolution, cmesh->Dimension());

    //TPZManVector<REAL> errors(4, 0.);
    //an.SetThreadsForError(2);
    //an.SetExact(analytic->ExactSolution());
    //an.PostProcessError(errors, false);
}

void InsertMaterialsInMHMMesh(TPZMHMixedMeshControl &control, const ProblemConfig &config) {
    TPZCompMesh &cmesh = control.CMesh();

    int dim = control.GMesh()->Dimension();
    cmesh.SetDimModel(dim);

    auto *mix = new TPZMixedPoisson(1, dim);

    TPZFMatrix<REAL> K(3, 3, 0), invK(3, 3, 0);
    K.Identity();
    invK.Identity();

    mix->SetPermeabilityTensor(K, invK);


    TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
    int dirichlet = 0;

    // Zero flux (reservoir boundary)
    TPZBndCond *zeroFlux = mix->CreateBC(mix, -1, dirichlet, val1, val2);
    // Productors
    val2(0, 0) = -10.;
    TPZBndCond *productors = mix->CreateBC(mix, -2, dirichlet, val1, val2);
    // Injectors
    val2(0, 0) = 20.;
    TPZBndCond *injectors = mix->CreateBC(mix, -3, dirichlet, val1, val2);

#ifdef FORCING_FUNCTION_DEBUG
    mix->SetForcingFunctionExact(config.exact.operator*().Exact());
    mix->SetForcingFunction(config.exact.operator*().ForcingFunction());
    zeroFlux->TPZMaterial::SetForcingFunction(config.exact.operator*().Exact());
    productors->TPZMaterial::SetForcingFunction(config.exact.operator*().Exact());
    injectors->TPZMaterial::SetForcingFunction(config.exact.operator*().Exact());
#endif

    cmesh.InsertMaterialObject(mix);
    cmesh.InsertMaterialObject(zeroFlux);
    cmesh.InsertMaterialObject(productors);
    cmesh.InsertMaterialObject(injectors);
}

void EstimateError(ProblemConfig &config, TPZMHMixedMeshControl *mhm) {
    cout << "Error Estimation processing for MHM-Hdiv problem " << endl;

    // Error estimation
    TPZMultiphysicsCompMesh *InputMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm->CMesh().operator->());
    if (!InputMesh) DebugStop();

    bool postProcWithHDiv = false;
    TPZMHMHDivErrorEstimator ErrorEstimator(*InputMesh, mhm, postProcWithHDiv);
    ErrorEstimator.SetAnalyticSolution(config.exact);
    ErrorEstimator.SetProblemConfig(config);
    ErrorEstimator.PotentialReconstruction();

    {
        string command = "mkdir " + config.dir_name;
        system(command.c_str());

        TPZManVector<REAL, 6> errors;
        TPZManVector<REAL> elementerrors;
        std::string vtkPath = config.dir_name + "/" + config.problemname + "Errors.vtk";
        ErrorEstimator.ComputeErrors(errors, elementerrors, vtkPath);
    }
}

void CreateMHMCompMesh(TPZMHMixedMeshControl *mhm, const ProblemConfig &config, int nInternalRef,
                       bool definePartitionByCoarseIndex, TPZManVector<int64_t> mhmIndexes) {

    if (definePartitionByCoarseIndex) {
        mhm->DefinePartitionbyCoarseIndices(mhmIndexes);
    } else {
        mhm->DefinePartition(mhmIndexes);
    }

    // Indicate material indices to the MHM control structure
    mhm->fMaterialIds = config.materialids;
    mhm->fMaterialBCIds = config.bcmaterialids;

    // Insert the material objects in the multiphysics mesh
    InsertMaterialsInMHMMesh(*mhm, config);

    // General approximation order settings
    mhm->SetInternalPOrder(config.porder);
    mhm->SetSkeletonPOrder(config.porder);
    mhm->SetHdivmaismaisPOrder(config.hdivmais);

    // Refine skeleton elements
    mhm->DivideSkeletonElements(nInternalRef);
    mhm->DivideBoundarySkeletonElements();

    // Creates MHM mesh
    bool substructure = true;
    mhm->BuildComputationalMesh(substructure);
    {
        std::string fileName = "CompMesh.txt";
        std::ofstream file(fileName);
        mhm->CMesh()->Print(file);
    }
}

void RemoveGelsOfGivenMaterial(TPZGeoMesh *gmesh, int matId) {
    int64_t nElem = gmesh->NElements();

    for (int64_t i = 0; i < nElem; i++) {
        TPZGeoEl* gel = gmesh->Element(i);
        if (gel->MaterialId() == matId) {
            gmesh->DeleteElement(gel, i);
            std::cout << "i: " << i << '\n';
        }
    }
}
