//
// Created by Gustavo A. Batistela on 06/07/2020.
//

#include <Mesh/pzgmesh.h>
#include <Pre/TPZGenGrid3D.h>
#include <Pre/TPZMHMixedMeshControl.h>
#include <TPZMFSolutionTransfer.h>
//#include <Tools.h>
#include <ToolsMHM.h>
#include <Util/pzlog.h>
#include <Material/DarcyFlow/TPZMixedDarcyFlow.h>

void RunSmoothProblem();
void RunHighGradientProblem();
void RunOscillatoryProblem();
void RunNonConvexProblem();
void Run3DProblem();
void RunInnerSingularityProblem();
void RunAdaptivityProblem();

TPZGeoMesh *CreateQuadGeoMesh(int nCoarseDiv, int nInternalRef);
TPZGeoMesh *CreateCubeGeoMesh(int nCoarseRef, int nInternalRef, TPZStack<int64_t> &coarseIndexes);
TPZGeoMesh *CreateLShapeGeoMesh(int nCoarseRef, int nInternalRef, TPZStack<int64_t> &mhmIndexes);

void InsertMaterialsInMHMMesh(TPZMHMixedMeshControl &control, const ProblemConfig &config);
void CreateMHMCompMesh(TPZMHMixedMeshControl *mhm, const ProblemConfig &config, int nInternalRef,
                       bool definePartitionByCoarseIndex, TPZManVector<int64_t>& mhmIndexes);

void SolveMHMProblem(TPZMHMixedMeshControl *mhm, const ProblemConfig &config);

void EstimateError(ProblemConfig &config, TPZMHMixedMeshControl *mhm);

void MHMAdaptivity(TPZMHMixedMeshControl *mhm, TPZGeoMesh* gmeshToRefine, ProblemConfig& config);

int main() {
    TPZLogger::InitializePZLOG();
    gRefDBase.InitializeAllUniformRefPatterns();
    //RunSmoothProblem();
    //RunHighGradientProblem();
    RunOscillatoryProblem();
    //RunNonConvexProblem();
    //Run3DProblem();
    //RunInnerSingularityProblem();
    
   // RunAdaptivityProblem();

    return 0;
}

void RunSmoothProblem() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
    config.problemname = "SinSin";
    config.dir_name = "MHM";
    config.porder = 1;
    config.hdivmais = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nCoarseDiv = 2;
    int nInternalRef = 0;

    config.ndivisions = nCoarseDiv;
    config.gmesh = CreateQuadGeoMesh(nCoarseDiv, nInternalRef);

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    TPZManVector<int64_t> coarseIndexes;
    ComputeCoarseIndices(config.gmesh, coarseIndexes);
    bool definePartitionByCoarseIndexes = true;
    CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, coarseIndexes);

    SolveMHMProblem(mhm, config);
    EstimateError(config, mhm);
}

void RunHighGradientProblem() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::EConst;
    config.problemname = "Constant";
    config.dir_name = "MHM";
    config.porder = 1;
    config.hdivmais = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nCoarseDiv = 2;
    int nInternalRef = 0;

    config.ndivisions = nCoarseDiv;
    config.gmesh = CreateQuadGeoMesh(nCoarseDiv, nInternalRef);

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    TPZManVector<int64_t> coarseIndexes;
    ComputeCoarseIndices(config.gmesh, coarseIndexes);
    bool definePartitionByCoarseIndexes = true;
    CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, coarseIndexes);

    SolveMHMProblem(mhm, config);
    EstimateError(config, mhm);
}

void RunOscillatoryProblem() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::EArcTan;
    config.problemname = "EArcTan";
    config.dir_name = "HdivRecTestOscillatory";
    config.porder = 1;
    config.hdivmais = 3;
    config.ndivisions = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nCoarseDiv = 3;
    int nInternalRef = 1;
    config.gmesh = CreateQuadGeoMesh(nCoarseDiv, nInternalRef);

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    TPZManVector<int64_t> coarseIndexes;
    ComputeCoarseIndices(config.gmesh, coarseIndexes);
    bool definePartitionByCoarseIndexes = true;
    CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, coarseIndexes);

    SolveMHMProblem(mhm, config);

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    EstimateError(config, mhm);
}

void RunNonConvexProblem() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
    config.problemname = "NonConvex";
    config.dir_name = "MHM";
    config.porder = 1;
    config.hdivmais = 3;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nDiv = 1;
    int nInternalRef = 0;

    TPZManVector<int64_t> coarseIndexes;
    config.ndivisions = nDiv;
    config.gmesh = CreateLMHMMesh(nDiv, coarseIndexes);

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    bool definePartitionByCoarseIndexes = false;
    CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, coarseIndexes);

    SolveMHMProblem(mhm, config);
    EstimateError(config, mhm);
}

void Run3DProblem() {
    ProblemConfig config;
    config.dimension = 3;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
    config.problemname = "SinSinCube";
    config.dir_name = "MHM";
    config.porder = 1;
    config.hdivmais = 1;
    config.ndivisions = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nCoarseDiv = 3;
    int nInternalRef = 1;

    TPZStack<int64_t> mhmIndexes;
    config.gmesh = CreateCubeGeoMesh(nCoarseDiv, nInternalRef, mhmIndexes);

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    bool definePartitionByCoarseIndexes = true;
    CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, mhmIndexes);

    SolveMHMProblem(mhm, config);
    EstimateError(config, mhm);
}

void RunInnerSingularityProblem() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinMark;
    config.problemname = "SinMarkLShape";
    config.dir_name = "MHM";
    config.porder = 1;
    config.hdivmais = 3;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nCoarseRef = 0;
    int nInternalRef = 0;

    config.ndivisions = nCoarseRef;
    TPZStack<int64_t> mhmIndexes;
    config.gmesh = CreateLShapeGeoMesh(nCoarseRef, nInternalRef, mhmIndexes);

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    bool definePartitionByCoarseIndexes = true;
    CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, mhmIndexes);

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }

    SolveMHMProblem(mhm, config);


    EstimateError(config, mhm);
}

TPZGeoMesh *CreateQuadGeoMesh(int nCoarseDiv, int nInternalRef) {

    TPZManVector<int, 4> bcIDs(4, -1);
    TPZGeoMesh *gmesh = Tools::CreateGeoMesh(nCoarseDiv, bcIDs);
    gmesh->SetDimension(2);
    
    Tools::UniformRefinement(nInternalRef, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    return gmesh;
}

TPZGeoMesh *CreateCubeGeoMesh(int nCoarseRef, int nInternalRef, TPZStack<int64_t> &coarseIndexes) {

    auto *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(3);
    int matID = 1;

    // Creates matrix with node coordinates
    const int NodeNumber = 8;
    REAL coordinates[NodeNumber][3] = {
        {0., 0., 0.},
        {1., 0., 0.},
        {1., 1., 0.},
        {0., 1., 0.},
        {0., 0., 1.},
        {1., 0., 1.},
        {1., 1., 1.},
        {0., 1., 1.}
    };

    // Inserts coordinates in the TPZGeoMesh object
    for (int i = 0; i < NodeNumber; i++) {
        int64_t nodeID = gmesh->NodeVec().AllocateNewElement();

        TPZVec<REAL> nodeCoord(3);
        nodeCoord[0] = coordinates[i][0];
        nodeCoord[1] = coordinates[i][1];
        nodeCoord[2] = coordinates[i][2];

        gmesh->NodeVec()[nodeID] = TPZGeoNode(i, nodeCoord, *gmesh);
    }

    // Creates cube element
    TPZManVector<int64_t> nodeIDs(8);
    nodeIDs[0] = 0;
    nodeIDs[1] = 1;
    nodeIDs[2] = 2;
    nodeIDs[3] = 3;
    nodeIDs[4] = 4;
    nodeIDs[5] = 5;
    nodeIDs[6] = 6;
    nodeIDs[7] = 7;
    new TPZGeoElRefPattern<pzgeom::TPZGeoCube>(nodeIDs, matID, *gmesh);

    // Creates boundary faces
    nodeIDs.Resize(4);
    matID = -1;
    nodeIDs[0] = 0;
    nodeIDs[1] = 1;
    nodeIDs[2] = 2;
    nodeIDs[3] = 3;
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeIDs, matID, *gmesh);
    nodeIDs[0] = 4;
    nodeIDs[1] = 5;
    nodeIDs[2] = 6;
    nodeIDs[3] = 7;
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeIDs, matID, *gmesh);

    for (int i = 0; i < 4; i++) {
        nodeIDs[0] = (0 + i) % 4;
        nodeIDs[1] = (1 + i) % 4;
        nodeIDs[2] = nodeIDs[1] + 4;
        nodeIDs[3] = nodeIDs[0] + 4;
        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeIDs, matID, *gmesh);
    }

    gmesh->BuildConnectivity();

    Tools::UniformRefinement(nCoarseRef, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    int64_t nElem = gmesh->NElements();

    coarseIndexes.clear();
    for (int64_t i = 0; i < nElem; i++) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (gel->Dimension() != gmesh->Dimension() || gel->NSubElements() > 0) continue;
        coarseIndexes.Push(i);
    }

    Tools::UniformRefinement(nInternalRef, gmesh);

    return gmesh;
}

TPZGeoMesh *CreateLShapeGeoMesh(int nCoarseRef, int nInternalRef, TPZStack<int64_t> &mhmIndexes) {

    TPZVec<int> bcIDs(8, -1);
    TPZGeoMesh *gmesh = Tools::CreateQuadLShapeMesh(bcIDs);
    gmesh->SetDimension(2);
    gmesh->BuildConnectivity();

    Tools::UniformRefinement(nCoarseRef, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    int64_t nElem = gmesh->NElements();
    for (int64_t i = 0; i < nElem; i++) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (gel->Dimension() != gmesh->Dimension() || gel->NSubElements() > 0) continue;
        mhmIndexes.Push(i);
    }

    Tools::UniformRefinement(nInternalRef, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    for (int64_t i = 0; i < mhmIndexes.size(); i++) {
        std::cout << mhmIndexes[i] << '\n';
    }
    std::cout << '\n';

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
}

void SolveMHMProblem(TPZMHMixedMeshControl *mhm, const ProblemConfig &config) {

    TPZAutoPointer<TPZCompMesh> cmesh = mhm->CMesh();

    bool shouldrenumber = true;
    TPZLinearAnalysis an(cmesh, shouldrenumber);

#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> strmat(cmesh.operator->());
    strmat.SetNumThreads(0 /*config.n_threads*/);
#else
    TPZSkylineStructMatrix<STATE> strmat(cmesh.operator->());
    strmat.SetNumThreads(0);
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

    int resolution = 0;
    std::string plotname = config.dir_name + "/" + config.problemname + "Results.vtk";
    an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotname);
    an.PostProcess(resolution, cmesh->Dimension());

    {
        TPZCompMesh* cmesh =  mhm->CMesh().operator->();
        std::ofstream out("SolvedMHMMesh.txt");
        cmesh->Print(out);

        std::ofstream outFlux("SolvedFluxMHMMesh.txt");
        mhm->FluxMesh()->Print(outFlux);

        std::stringstream solByElement;
        TPZCompMeshTools::PrintSolutionByGeoElement(mhm->PressureMesh().operator->(), solByElement);
        std::ofstream solByElementFile("SolByElementPressureMHM.txt");
        solByElementFile << solByElement.str();
    }

    //TPZManVector<REAL> errors(4, 0.);
    //an.SetThreadsForError(2);
    //an.SetExact(analytic->ExactSolution());
    //an.PostProcessError(errors, false);
}

void EstimateError(ProblemConfig &config, TPZMHMixedMeshControl *mhm) {

    std::cout << "\nError Estimation processing for MHM-Hdiv problem " << std::endl;

    // Error estimation
    TPZMultiphysicsCompMesh *originalMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm->CMesh().operator->());
    if (!originalMesh) DebugStop();

    bool postProcWithHDiv = false;
    TPZMHMHDivErrorEstimator ErrorEstimator(*originalMesh, mhm, postProcWithHDiv);
    ErrorEstimator.SetAnalyticSolution(config.exact);
    ErrorEstimator.PotentialReconstruction();

    {
        std::string command = "mkdir -p " + config.dir_name;
        system(command.c_str());

        TPZManVector<REAL, 6> errors;
        TPZManVector<REAL> elementerrors;
        std::string outVTK = config.dir_name + "/out.vtk";
        ErrorEstimator.ComputeErrors(errors, elementerrors, outVTK);
    }
}

void InsertMaterialsInMHMMesh(TPZMHMixedMeshControl &control, const ProblemConfig &config) {
    TPZCompMesh &cmesh = control.CMesh();

    int dim = control.GMesh()->Dimension();
    cmesh.SetDimModel(dim);

    auto *mat = new TPZMixedDarcyFlow(1, dim);

    auto exact_lambda = [config](const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv) {
        config.exact.operator*().Exact()->Execute(loc, result, deriv);
    };

    auto ff_lambda = [config](const TPZVec<REAL> &loc, TPZVec<STATE> &result) {
        config.exact.operator*().ForcingFunction()->Execute(loc, result);
    };

    mat->SetPermeabilityFunction(1);
    mat->SetExactSol(exact_lambda, 8);
    mat->SetForcingFunction(ff_lambda, 8);

    cmesh.InsertMaterialObject(mat);

    for (auto matid : config.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.);
        TPZManVector<REAL, 1> val2(1, 0.);
        int bctype = 0;
        auto *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->SetForcingFunctionBC(exact_lambda);
        cmesh.InsertMaterialObject(bc);
    }
}
void RunAdaptivityProblem(){
    
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::EBoundaryLayer;
    config.problemname = "EBoundaryLayer";
    config.dir_name = "Adaptivity";
    config.porder = 1;
    config.hdivmais = 3;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;
    config.ndivisions = 3;
    int refinementSteps = 1;

    
    
    for(int hsk=2; hsk < config.ndivisions; hsk++){
    
        for(int hin=1 ; hin < 2 ; hin ++){
            
            //    int nCoarseRef = 2;
            //    int nInternalRef = 0;
            //
            //    config.ndivisions = nCoarseRef;
            config.gmesh = CreateQuadGeoMesh(hsk, hin);
            
            
            
            
            std::string command = "mkdir -p " + config.dir_name;
            system(command.c_str());
            
            {
                std::string fileName = config.dir_name + "/" + config.problemname + "GeoMesh.vtk";
                std::ofstream file(fileName);
                TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
            }
            
            auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
            TPZManVector<int64_t> coarseIndexes;
            ComputeCoarseIndices(config.gmesh, coarseIndexes);
            bool definePartitionByCoarseIndexes = true;
            CreateMHMCompMesh(mhm, config, hin, definePartitionByCoarseIndexes, coarseIndexes);
            
            {
                std::string fileName = config.dir_name + "/" + config.problemname + "GeoMeshAfterPartition.vtk";
                std::ofstream file(fileName);
                TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
            }
            //
            //    SolveMHMProblem(mhm, config);
            //    EstimateError(config, mhm);
            
          //  for (int iSteps = 0; iSteps < refinementSteps; iSteps++) {
                
//                {
//                    std::string fileName = config.dir_name + "/" + config.problemname + "GeoMeshAdapt.vtk";
//                    std::ofstream file(fileName);
//                    TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
//                }
                
 //               config.adaptivityStep = iSteps;
                
                //        auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
                //        TPZManVector<int64_t> coarseIndexes;
                //        ComputeCoarseIndices(config.gmesh, coarseIndexes);
                //        bool definePartitionByCoarseIndexes = true;
                //        CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, coarseIndexes);
                
                SolveMHMProblem(mhm, config);
                EstimateError(config, mhm);
                
 //               MHMAdaptivity(mhm,  config.gmesh, config);
//#ifdef PZDEBUG
//                {
//                    std::ofstream out("GmeshAfterAdapty.vtk");
//                    TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, out);
//                }
//#endif
           // }
        }
        
    }
            
    
}


    
void MHMAdaptivity(TPZMHMixedMeshControl *mhm, TPZGeoMesh* gmeshToRefine, ProblemConfig& config) {

    // Column of the flux error estimate on the element solution matrix
    const int fluxErrorEstimateCol = 3;
    
    
    
    TPZMultiphysicsCompMesh *cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm->CMesh().operator->());
       if (!cmesh) DebugStop();

    
   // TPZCompMesh &cmesh = mhm->CMesh();
  
    TPZFMatrix<STATE> &elsol = cmesh->ElementSolution();
    int64_t nelem = elsol.Rows();

    //postProcessMesh->ElementSolution().Print("ElSolutionForAdaptivity",std::cout);

    // Iterates through element errors to get the maximum value
    REAL maxError = 0.;
    for (int64_t iel = 0; iel < nelem; iel++) {
        TPZCompEl* cel = cmesh->ElementVec()[iel];
        if (!cel) continue;
        if (cel->Dimension() != cmesh->Dimension()) continue;
        REAL elementError = elsol(iel, fluxErrorEstimateCol);


        if (elementError > maxError) {
            maxError = elementError;
        }
    }

    std::cout << "max error " << maxError << "\n";

    // Refines elements which error are bigger than 30% of the maximum error
    REAL threshold = 0.2 * maxError;

    for (int64_t iel = 0; iel < nelem; iel++) {
        TPZCompEl* cel = cmesh->ElementVec()[iel];
        if (!cel) continue;
        if (cel->Dimension() != cmesh->Dimension()) continue;

        REAL elementError = elsol(iel, fluxErrorEstimateCol);
        //prefinement
        if (elementError > threshold) {

            std::cout << "element error " << elementError << "el " << iel << "\n";
            TPZGeoEl* gel = cel->Reference();
            int iel = gel->Id();

            TPZVec<TPZGeoEl*> sons;
            TPZGeoEl* gelToRefine = gmeshToRefine->ElementVec()[iel];
            if (gelToRefine && !gelToRefine->HasSubElement()) {
                gelToRefine->Divide(sons);
#ifdef LOG4CXX2
                int nsides = gelToRefine->NSides();
                TPZVec<REAL> loccenter(gelToRefine->Dimension());
                TPZVec<REAL> center(3);
                gelToRefine->CenterPoint(nsides - 1, loccenter);
                
                gelToRefine->X(loccenter, center);
                static LoggerPtr logger(Logger::getLogger("HDivErrorEstimator"));
                if (logger->isDebugEnabled()) {
                    std::stringstream sout;
                    sout << "\nCenter coord: = " << center[0] << " " << center[1] << "\n";
                    sout << "Error = " << elementError << "\n\n";
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
            }
        } else {
            std::cout << "como refinar em p? " << "\n";
//            TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
//            if(!sp) continue;
//            int level = sp->Reference()->Level();
//            int ordem = config.porder + (config.adaptivityStep -1 ) + (level);
//            std::cout<<"level "<< level<<" ordem "<<ordem<<std::endl;
//            sp->PRefine(ordem);

        }
    }
    Tools::DivideLowerDimensionalElements(gmeshToRefine);
}
