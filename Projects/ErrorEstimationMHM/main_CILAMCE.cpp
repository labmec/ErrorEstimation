//
// Created by Gustavo A. Batistela on 06/07/2020.
//
// This file contains the numerical tests to be shown in the CILAMCE 2020 article.
//

#include <Mesh/pzgmesh.h>
#include <Pre/TPZGenGrid3D.h>
#include <Pre/TPZMHMixedMeshControl.h>
#include <TPZMFSolutionTransfer.h>
//#include <Tools.h>
#include <ToolsMHM.h>
#include <Util/pzlog.h>

void RunSinSinProblem();
void RunOscillatoryProblem();
void RunNonConvexProblem();
void Run3DProblem();
void RunSingularProblem();
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
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    gRefDBase.InitializeAllUniformRefPatterns();

    RunSinSinProblem();
    //RunOscillatoryProblem();
    //RunNonConvexProblem();
    //Run3DProblem();
    //RunSingularProblem();
    
    //RunAdaptivityProblem();

    return 0;
}

void RunSinSinProblem() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
    config.problemname = "SinSin";
    config.dir_name = "CILAMCE";
    config.porder = 1;
    config.hdivmais = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nCoarseDiv = 2;
    int nInternalRef = 0;

    config.ndivisions = nCoarseDiv;
    config.gmesh = CreateQuadGeoMesh(nCoarseDiv, nInternalRef);

    std::string command = "mkdir " + config.dir_name;
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
    config.problemname = "Oscillatory";
    config.dir_name = "CILAMCE";
    config.porder = 1;
    config.hdivmais = 1;
    config.ndivisions = 2;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;

    int nCoarseDiv = 2;
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
    config.dir_name = "CILAMCE";
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

    std::string command = "mkdir " + config.dir_name;
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
    config.dir_name = "CILAMCE";
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

    std::string command = "mkdir " + config.dir_name;
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

void RunSingularProblem() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinMark;
    config.problemname = "SinMarkLShape";
    config.dir_name = "CILAMCE";
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

    std::string command = "mkdir " + config.dir_name;
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
    {
        std::string fileName = "CompMesh.txt";
        std::ofstream file(fileName);
        mhm->CMesh()->Print(file);
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


    {
        std::cout << "GMesh elements after building MHM structure:\n";
        for (int i = 0; i < config.gmesh->NElements();i++) {
            TPZGeoEl * gel = config.gmesh->Element(i);
            std::cout << "Gel: " << i;
            if (gel) {
                std::cout << ", matid: " << gel->MaterialId() << "\n";
            }
            else {
                std:: cout << " NULL\n";
            }
        }
    }

    cout << "\nError Estimation processing for MHM-Hdiv problem " << endl;

    // Error estimation
    TPZMultiphysicsCompMesh *originalMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm->CMesh().operator->());
    if (!originalMesh) DebugStop();

    bool postProcWithHDiv = true;
    TPZMHMHDivErrorEstimator ErrorEstimator(*originalMesh, mhm, postProcWithHDiv);
    ErrorEstimator.SetAnalyticSolution(config.exact);
    ErrorEstimator.SetProblemConfig(config);
    ErrorEstimator.PotentialReconstruction();

    {
        string command = "mkdir " + config.dir_name;
        system(command.c_str());

        TPZManVector<REAL, 6> errors;
        TPZManVector<REAL> elementerrors;
        bool store_errors = true;
        ErrorEstimator.ComputeErrors(errors, elementerrors, store_errors);
    }
}

void InsertMaterialsInMHMMesh(TPZMHMixedMeshControl &control, const ProblemConfig &config) {
    TPZCompMesh &cmesh = control.CMesh();

    int dim = control.GMesh()->Dimension();
    cmesh.SetDimModel(dim);

    TPZMixedPoisson *mat = new TPZMixedPoisson(1, dim);

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
void RunAdaptivityProblem(){
    
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
    config.problemname = "SinSin";
    config.dir_name = "Adaptivity";
    config.porder = 1;
    config.hdivmais = 3;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.makepressurecontinuous = true;
    
    int nCoarseRef = 2;
    int nInternalRef = 0;
    
    config.ndivisions = nCoarseRef;
    config.gmesh = CreateQuadGeoMesh(nCoarseRef, nInternalRef);
    
    int refinementSteps = 2;
    
    
    std::string command = "mkdir " + config.dir_name;
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
//
//    SolveMHMProblem(mhm, config);
//    EstimateError(config, mhm);
    
    for (int iSteps = 0; iSteps < refinementSteps; iSteps++) {
        
        {
            std::string fileName = config.dir_name + "/" + config.problemname + "GeoMeshAdapt.vtk";
            std::ofstream file(fileName);
            TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
        }
        
        config.adaptivityStep = iSteps;
        
//        auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
//        TPZManVector<int64_t> coarseIndexes;
//        ComputeCoarseIndices(config.gmesh, coarseIndexes);
//        bool definePartitionByCoarseIndexes = true;
//        CreateMHMCompMesh(mhm, config, nInternalRef, definePartitionByCoarseIndexes, coarseIndexes);
        
        SolveMHMProblem(mhm, config);
        EstimateError(config, mhm);
        
        MHMAdaptivity(mhm,  config.gmesh, config);
#ifdef PZDEBUG
        {
            std::ofstream out("GmeshAfterAdapty.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, out);
        }
#endif
    }
            
    
}


    
void MHMAdaptivity(TPZMHMixedMeshControl *mhm, TPZGeoMesh* gmeshToRefine, ProblemConfig& config) {

    // Column of the flux error estimate on the element solution matrix
    const int fluxErrorEstimateCol = 3;
    
    
    
    TPZMultiphysicsCompMesh *cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm->CMesh().operator->());
       if (!cmesh) DebugStop();

    
   // TPZCompMesh &cmesh = mhm->CMesh();
  

    int64_t nelem = cmesh->ElementSolution().Rows();

    //postProcessMesh->ElementSolution().Print("ElSolutionForAdaptivity",std::cout);

    // Iterates through element errors to get the maximum value
    REAL maxError = 0.;
    for (int64_t iel = 0; iel < nelem; iel++) {
        TPZCompEl* cel = cmesh->ElementVec()[iel];
        if (!cel) continue;
        if (cel->Dimension() != cmesh->Dimension()) continue;
        REAL elementError = cmesh->ElementSolution()(iel, fluxErrorEstimateCol);


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

        REAL elementError = cmesh->ElementSolution()(iel, fluxErrorEstimateCol);
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
