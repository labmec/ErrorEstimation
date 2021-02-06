//
// Created by gus on 1/19/21.
//

#include <Mesh/pzgmesh.h>
#include <Pre/TPZGenGrid3D.h>
#include <Pre/TPZMHMixedMeshControl.h>
#include <TPZMFSolutionTransfer.h>
#include <ToolsMHM.h>
#include <Util/pzlog.h>

TPZGeoMesh *CreateQuadGeoMesh(int nCoarseDiv, int nInternalRef);
void CreateMHMCompMesh(TPZMHMixedMeshControl *mhm, const ProblemConfig &config, int nInternalRef,
                       bool definePartitionByCoarseIndex, TPZManVector<int64_t>& mhmIndexes);
void InsertMaterialsInMHMMesh(TPZMHMixedMeshControl &control, const ProblemConfig &config);
void SolveMHMProblem(TPZMHMixedMeshControl *mhm, const ProblemConfig &config);
void SolveH1Problem(TPZCompMesh *cmesh, const ProblemConfig &config);

void RunOscillatoryProblemMHM();
void RunOscillatoryProblemHDiv();
void RunOscillatoryProblemH1();

int main() {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    gRefDBase.InitializeAllUniformRefPatterns();

    RunOscillatoryProblemMHM();
    RunOscillatoryProblemHDiv();
    RunOscillatoryProblemH1();

    return 0;
}

void RunOscillatoryProblemHDiv() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::EArcTan;
    config.problemname = "HDivOscillatory";
    config.dir_name = "DifferentMethods";
    config.porder = 2;
    config.hdivmais = 2;
    config.ndivisions = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nCoarseDiv = 4;
    int nInternalRef = 0;
    config.gmesh = CreateQuadGeoMesh(nCoarseDiv, nInternalRef);

    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    TPZManVector<int64_t> coarseIndexes;
    ComputeCoarseIndices(config.gmesh, coarseIndexes);
    bool definePartitionByCoarseIndexes = true;
    CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, coarseIndexes);

    SolveMHMProblem(mhm, config);
}

void RunOscillatoryProblemMHM() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::EArcTan;
    config.problemname = "MHMOscillatory";
    config.dir_name = "DifferentMethods";
    config.porder = 2;
    config.hdivmais = 2;
    config.ndivisions = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nCoarseDiv = 3;
    int nInternalRef = 3;
    config.gmesh = CreateQuadGeoMesh(nCoarseDiv, nInternalRef);

    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    TPZManVector<int64_t> coarseIndexes;
    ComputeCoarseIndices(config.gmesh, coarseIndexes);
    bool definePartitionByCoarseIndexes = true;
    CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, coarseIndexes);

    SolveMHMProblem(mhm, config);
}

void RunOscillatoryProblemH1() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::EArcTan;
    config.problemname = "H1Oscillatory";
    config.dir_name = "DifferentMethods";
    config.porder = 3;
    config.hdivmais = 2;
    config.ndivisions = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nCoarseDiv = 12;
    int nInternalRef = 0;
    config.gmesh = CreateQuadGeoMesh(nCoarseDiv, nInternalRef);

    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    TPZCompMesh *cmesh = new TPZCompMesh(config.gmesh);

    cmesh->SetDimModel(config.dimension);
    cmesh->SetDefaultOrder(config.porder);

    TPZMatPoisson3d *mat = new TPZMatPoisson3d(1, config.dimension);

    TPZFMatrix<REAL> K(3, 3, 0), invK(3, 3, 0);
    K.Identity();
    invK.Identity();

    mat->SetForcingFunctionExact(config.exact.operator*().Exact());
    mat->SetForcingFunction(config.exact.operator*().ForcingFunction());
    mat->SetPermeabilityTensor(K, invK);

    cmesh->InsertMaterialObject(mat);

    for (auto matid : config.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(config.exact.operator*().Exact());
        cmesh->InsertMaterialObject(bc);
    }
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    SolveH1Problem(cmesh, config);
}

TPZGeoMesh *CreateQuadGeoMesh(int nCoarseDiv, int nInternalRef) {

    TPZManVector<int, 4> bcIDs(4, -1);
    TPZGeoMesh *gmesh = Tools::CreateGeoMesh(nCoarseDiv, bcIDs);
    gmesh->SetDimension(2);

    {
        std::ofstream file("GMeshCoarse.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    }

    Tools::UniformRefinement(nInternalRef, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    {
        std::ofstream file("GMeshFine.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    }

    return gmesh;
}

void CreateMHMCompMesh(TPZMHMixedMeshControl *mhm, const ProblemConfig &config, int nInternalRef,
                       bool definePartitionByCoarseIndex, TPZManVector<int64_t>& mhmIndexes) {

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
    mhm->DivideSkeletonElements(0);
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

void InsertMaterialsInMHMMesh(TPZMHMixedMeshControl &control, const ProblemConfig &config) {
    TPZCompMesh &cmesh = control.CMesh();

    int dim = control.GMesh()->Dimension();
    cmesh.SetDimModel(dim);

    auto *mat = new TPZMixedPoisson(1, dim);

    TPZFMatrix<REAL> K(3, 3, 0), invK(3, 3, 0);
    K.Identity();
    invK.Identity();

    mat->SetForcingFunctionExact(config.exact.operator*().Exact());
    mat->SetForcingFunction(config.exact.operator*().ForcingFunction());
    mat->SetPermeabilityTensor(K, invK);

    cmesh.InsertMaterialObject(mat);

    for (auto matid : config.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(config.exact.operator*().Exact());
        cmesh.InsertMaterialObject(bc);
    }
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

    TLaplaceExample1 *analytic = &config.exact.operator*();
    an.SetExact(analytic->ExactSolution());

    scalnames.Push("Pressure");
    scalnames.Push("ExactPressure");
    vecnames.Push("Flux");
    vecnames.Push("ExactFlux");

    int resolution = 5;
    std::string plotname = config.dir_name + "/" + config.problemname + "Results.vtk";
    an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotname);
    an.PostProcess(resolution, cmesh->Dimension());
}

void SolveH1Problem(TPZCompMesh *cmesh, const ProblemConfig &config) {

    bool shouldrenumber = true;
    TPZAnalysis an(cmesh, shouldrenumber);

#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh);
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

    TPZStack<std::string> scalnames, vecnames;
    TPZMaterial *mat = cmesh->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }

    TLaplaceExample1 *analytic = &config.exact.operator*();
    an.SetExact(analytic->ExactSolution());

    scalnames.Push("Pressure");
    scalnames.Push("ExactPressure");
    vecnames.Push("Flux");
    vecnames.Push("MinusKGradU");
    vecnames.Push("ExactFlux");

    int resolution = 3;
    std::string plotname = config.dir_name + "/" + config.problemname + "Results.vtk";
    an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotname);
    an.PostProcess(resolution, cmesh->Dimension());
}

#include "TPZRefPatternDataBase.h"
#include "TPZVTKGeoMesh.h"
#include "tpzgeoelrefpattern.h"

#include "ProblemConfig.h"
#include "mixedpoisson.h"
#include "TPZMatLaplacianHybrid.h"
#include "TPZVecL2.h"
#include "pzbndcond.h"

#include "TPZHybridizeHDiv.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZMultiphysicsCompMesh.h"
#include "Tools.h"
#include "TPZBFileStream.h"
#include "TPZCreateMultiphysicsSpace.h"
#include <memory>

bool neumann = true;
bool h1solution = false;
bool hybridh1 = true;

void InsertMaterialObjectsH1Hybrid(TPZMultiphysicsCompMesh *cmesh, ProblemConfig &config);
void SolveHybridH1Problem(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,int InterfaceMatId,struct ProblemConfig config);

int main2(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    ProblemConfig config;

    config.porder = 3;
    config.ndivisions = 3;
    config.dimension = 2;
    config.prefine = false;

    int orderlagrange = 1;

    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::EConst;

    config.problemname = "ESinSin k=1 e lagrange order 2";

    config.dir_name = "HybridH1_ESinSin";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    // geometric mesh
    TPZManVector<int, 4> bcids(4, -1);
    TPZGeoMesh *gmesh = Tools::CreateGeoMesh(2, bcids);

    config.gmesh = gmesh;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);



    Tools::UniformRefinement(config.ndivisions, gmesh);
    int refinement_depth = 3;
    Tools::RandomRefinement(config.gmesh, 5, refinement_depth);

#ifdef PZDEBUG
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmeshInitial.txt");
        gmesh->Print(out2);

    }
#endif
    if(hybridh1){
        config.exact.operator*().fSignConvention = 1;
        TPZCreateMultiphysicsSpace createspace(gmesh,TPZCreateMultiphysicsSpace::EH1Hybrid);

        createspace.SetMaterialIds({1}, {-2,-1});
        createspace.fH1Hybrid.fHybridizeBCLevel = 1;//opcao de hibridizar o contorno
        createspace.ComputePeriferalMaterialIds();


        std::cout<<"---Original PerifericalMaterialId --- "<<std::endl;
        std::cout <<" fMatWrapId + = "<<createspace.fH1Hybrid.fMatWrapId<<std::endl;
        std::cout <<" fLagrangeMatid + = "<<createspace.fH1Hybrid.fLagrangeMatid.first<<std::endl;
        std::cout <<" fLagrangeMatid - = "<<createspace.fH1Hybrid.fLagrangeMatid.second<<std::endl;
        std::cout <<" fFluxMatId = "<<createspace.fH1Hybrid.fFluxMatId<<std::endl;
        std::cout << "fSecond Lagrange MatID = " <<createspace.fH1Hybrid.fSecondLagrangeMatid<<std::endl;
        std::cout << "fInterfacePressure = " <<createspace.fH1Hybrid.fInterfacePressure<<std::endl;
        std::cout << "fIBC hybridization level = " <<createspace.fH1Hybrid.fHybridizeBCLevel<<std::endl;

        TPZManVector<TPZCompMesh *> meshvec;

        createspace.CreateAtomicMeshes(meshvec,config.porder,orderlagrange);

        TPZMultiphysicsCompMesh *cmesh_H1Hybrid = new TPZMultiphysicsCompMesh(gmesh);
        InsertMaterialObjectsH1Hybrid(cmesh_H1Hybrid, config);
        createspace.InsertPeriferalMaterialObjects(cmesh_H1Hybrid);

#ifdef PZDEBUG
        {
            std::ofstream out3("pressure.txt");
            meshvec[0]->Print(out3);
            std::ofstream out4("flux.txt");
            meshvec[1]->Print(out4);
            std::ofstream out("gmeshIncremented.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
            std::ofstream out2("gmeshIncremented.txt");
            gmesh->Print(out2);

        }
#endif
        cmesh_H1Hybrid->BuildMultiphysicsSpace(meshvec);

#ifdef PZDEBUG
        {
            std::map<int,int> matelem;
            int64_t nel = cmesh_H1Hybrid->NElements();
            for (int64_t el = 0; el<nel; el++) {
                TPZCompEl *cel = cmesh_H1Hybrid->Element(el);
                TPZGeoEl *gel = cel->Reference();
//                TPZManVector<REAL,3> center(3);
//                TPZGeoElSide gelside(gel);
//                gelside.CenterX(center);
//                std::cout << "Matid " << gel->MaterialId() << " center " << center << std::endl;
                matelem[gel->MaterialId()]++;
            }
            std::cout << __PRETTY_FUNCTION__ << " number of computational elements by material \n";
            for (auto it : matelem) {
                std::cout << "Material id " << it.first << " number of elements " << it.second << std::endl;
            }
        }
#endif


        createspace.InsertLagranceMaterialObjects(cmesh_H1Hybrid);
        createspace.AddInterfaceElements(cmesh_H1Hybrid);
#ifdef PZDEBUG
        {
            std::map<int,int> matelem;
            int64_t nel = cmesh_H1Hybrid->NElements();
            for (int64_t el = 0; el<nel; el++) {
                TPZCompEl *cel = cmesh_H1Hybrid->Element(el);
                TPZGeoEl *gel = cel->Reference();
                matelem[gel->MaterialId()]++;
            }
            std::cout << __PRETTY_FUNCTION__ << " number of computational elements by material \n";
            for (auto it : matelem) {
                std::cout << "Material id " << it.first << " number of elements " << it.second << std::endl;
            }
        }
#endif
        cmesh_H1Hybrid->ComputeNodElCon();
#ifdef PZDEBUG
        {
            std::ofstream out("mphysicsmeshBeforeCondense.txt");
            cmesh_H1Hybrid->Print(out);
        }
#endif
        createspace.GroupandCondenseElements(cmesh_H1Hybrid);

        cmesh_H1Hybrid->InitializeBlock();
        cmesh_H1Hybrid->ComputeNodElCon();

#ifdef PZDEBUG
        {
            std::ofstream out("mphysicsmesh.txt");
            cmesh_H1Hybrid->Print(out);
        }
#endif
        //Solve Hybrid problem

        SolveHybridH1Problem(cmesh_H1Hybrid,createspace.fH1Hybrid.fLagrangeMatid.first,config);

        // Post Processing for Lagrange Multiplier
        {
            TPZAnalysis an(meshvec[1], false);

            TPZStack<std::string> scalnames, vecnames;
            scalnames.Push("State");

            int dim = 1;
            std::string plotname("LagrangeMultiplier.vtk");
            an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
            an.PostProcess(0, dim);
        }
    }
    return 0.;
}

void InsertMaterialObjectsH1Hybrid(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, ProblemConfig &config) {
    TPZGeoMesh *gmesh = cmesh_H1Hybrid->Reference();
    int dim = gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    // Creates Poisson material
    TPZMatLaplacianHybrid *material = new TPZMatLaplacianHybrid(matID, dim);

    cmesh_H1Hybrid->InsertMaterialObject(material);
    if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
        material->SetForcingFunction(config.exact.operator*().ForcingFunction());
        material->SetForcingFunctionExact(config.exact.operator*().Exact());
    }
    //    TPZMaterial * mat(material);
    //    cmesh->InsertMaterialObject(mat);

    // Inserts boundary conditions
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 1.);
    TPZMaterial *BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2);
    if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
        BCond0->SetForcingFunction(config.exact.operator*().Exact());
    }
    val2.Zero();
    TPZMaterial *BCond1 = material->CreateBC(material, -2, neumann, val1, val2);

    cmesh_H1Hybrid->InsertMaterialObject(BCond0);
    cmesh_H1Hybrid->InsertMaterialObject(BCond1);
}

void SolveHybridH1Problem(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, int InterfaceMatId, struct ProblemConfig config) {

    TPZAnalysis an(cmesh_H1Hybrid);

#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh_H1Hybrid);
    strmat.SetNumThreads(0);
    //        strmat.SetDecomposeType(ELDLt);
#else
    //    TPZFrontStructMatrix<TPZFrontSym<STATE> > strmat(Hybridmesh);
    //    strmat.SetNumThreads(2);
    //    strmat.SetDecomposeType(ELDLt);
    TPZSkylineStructMatrix strmat(cmesh_H1Hybrid);
    strmat.SetNumThreads(0);
#endif

    std::set<int> matIds;

    for (auto matid : config.materialids) {
        matIds.insert(matid);
    }

    for (auto matidbc : config.bcmaterialids) {
        matIds.insert(matidbc);
    }

    matIds.insert(InterfaceMatId);
    strmat.SetMaterialIds(matIds);
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();

    // Pos processamento
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    scalnames.Push("PressureExact");

    int dim = 2;
    std::string plotname;
    {
        std::stringstream out;
        out << config.dir_name << "/"
            << "HybridH1" << config.porder << "_" << dim << "D_" << config.problemname << "Ndiv_ " << config.ndivisions
            << "HdivMais" << config.hdivmais << ".vtk";
        plotname = out.str();
    }
    int resolution = 3;
    an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
    an.PostProcess(resolution, dim);

    // Calculo do Erro
    an.SetExact(config.exact.operator*().ExactSolution());

    TPZManVector<REAL> errorvec(5, 0.);
    int64_t nelem = cmesh_H1Hybrid->NElements();
    cmesh_H1Hybrid->LoadSolution(cmesh_H1Hybrid->Solution());
    cmesh_H1Hybrid->ExpandSolution();
    cmesh_H1Hybrid->ElementSolution().Redim(nelem, 5);

    an.PostProcessError(errorvec); // calculo do erro com sol exata e aprox

    std::cout << "Computed errors " << errorvec << std::endl;

    // Erro
    ofstream myfile;
    myfile.open("ArquivosErrosH1Hibrido.txt", ios::app);
    myfile << "\n\n Estimator errors for Problem " << config.problemname;
    myfile << "\n-------------------------------------------------- \n";
    myfile << "Ndiv = " << config.ndivisions << " Order = " << config.porder << "\n";
    myfile << "DOF Total = " << cmesh_H1Hybrid->NEquations() << "\n";
    myfile << "error norm L2 = " << errorvec[0] << "\n";
    myfile << "semi norm H1 = " << errorvec[1] << "\n";
    myfile << "H1 norm = " << errorvec[2] << "\n";
    myfile << "energy norm = " << errorvec[3] << "\n";
    myfile.close();
}
