#include "Tools.h"

#include <Analysis/pzanalysis.h>
#include <Geom/pzgeoquad.h>
#include <Material/REAL/mixedpoisson.h>
#include <Material/TPZNullMaterial.h>
#include <Material/TPZVecL2.h>
#include <Material/pzbndcond.h>
#include <Matrix/pzstepsolver.h>
#include <Mesh/TPZCompMeshTools.h>
#include <Mesh/TPZMultiphysicsCompMesh.h>
#include <Mesh/tpzgeoelrefpattern.h>
#include <Post/TPZVTKGeoMesh.h>
#include <Pre/TPZGmshReader.h>
#include <Pre/TPZHybridizeHDiv.h>
#include <ProblemConfig.h>
#include <Refine/TPZRefPatternTools.h>
#include <StrMatrix/TPZSSpStructMatrix.h>
#include <StrMatrix/pzstrmatrix.h>

#include <libInterpolate/Interpolate.hpp>

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>

#define DEBUGTEST

TPZGeoMesh *CreateFlatGeoMesh();

TPZGeoMesh *CreateDebugGeoMesh();

void ModifyZCoordinates(TPZGeoMesh *gmesh, std::string &filename);

void ReadReservoirGeometryData(const std::string &name, std::vector<double> &x, std::vector<double> &y,
                               std::vector<double> &z);

void PrintGeometry(TPZGeoMesh *gmesh, const std::string &file_name, bool printTXT = true, bool printVTK = true);

void UNISIMHDiv(TPZGeoMesh *gmesh);

TPZMultiphysicsCompMesh *CreateMixedCMesh(const ProblemConfig &problem);

TPZCompMesh *CreateFluxCMesh(const ProblemConfig &problem);

TPZCompMesh *CreatePressureCMesh(const ProblemConfig &problem);

void SolveMixedHybridProblem(TPZCompMesh *Hybridmesh, const ProblemConfig &problem);

void hAdaptivity(TPZGeoMesh* gmesh, TPZVec<REAL>& elementErrors, REAL thresholdRatio);

void ApplyDirectionalRefinement(TPZGeoMesh *gmesh, int nRef);

void MoveMeshToOrigin(TPZGeoMesh *gmesh);

void SpreadMeshRefinement(TPZGeoMesh* gmesh);

void RotateGeoMesh(TPZGeoMesh* gmesh, std::array<int, 2> coords_ids_to_swap);


int main() {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    gRefDBase.InitializeRefPatterns(2);
    gRefDBase.InitializeAllUniformRefPatterns();
#ifdef DEBUGTEST
    TPZGeoMesh *gmesh = CreateDebugGeoMesh();
    std::string meshFileName{"DebugMesh"};

    // TODO: if nDirectionalRefinements is equal to zero, the code fails during run time.
    //  If not (i.e. equal to 1), the code does not break, but the results are wrong.
    RotateGeoMesh(gmesh, {0, 2});
    int nDirectionalRefinements = 0;
#else
    TPZGeoMesh *gmesh = CreateFlatGeoMesh();
    std::string meshFileName{"UNISIMMesh"};
    int nDirectionalRefinements = 3;
#endif
    PrintGeometry(gmesh, meshFileName, false, true);
    ApplyDirectionalRefinement(gmesh, nDirectionalRefinements);
    meshFileName.append("AfterDirectionalRef");
    PrintGeometry(gmesh, meshFileName, false, true);

#ifdef DEBUGTEST
    int nSteps = 1;
#else
    int nSteps = 5;
#endif
    for (int i = 0; i < nSteps; i++) {
        UNISIMHDiv(gmesh);
    }
    delete gmesh;
    return 0;
}

void MoveMeshToOrigin(TPZGeoMesh *gmesh) {
    TPZManVector<REAL, 3> nod0(3, 0.);
    gmesh->NodeVec()[0].GetCoordinates(nod0);
    int64_t nnodes = gmesh->NNodes();
    for (int64_t no = 0; no < nnodes; no++) {
        TPZManVector<REAL, 3> co(3);
        gmesh->NodeVec()[no].GetCoordinates(co);
        for (int ic = 0; ic < 3; ic++) {
            co[ic] -= nod0[ic];
        }
        gmesh->NodeVec()[no].SetCoord(co);
    }
}

void ApplyDirectionalRefinement(TPZGeoMesh *gmesh, int nRef) {
    // Mat IDs of productors and injectors BCs
#ifdef DEBUGTEST
    set<int> matids{-2};
#else
    set<int> matids{-2, -3};
#endif

    for (auto i = 0; i < nRef; i++) {
        int nelements = gmesh->NElements();
        for (auto el = 0; el < nelements; el++) {
            TPZGeoEl *element = gmesh->ElementVec()[el];
            if (!element) continue;
            TPZRefPatternTools::RefineDirectional(element, matids);
        }
        cout << "Refinement step: " << i << "\nNumber of elements = " << gmesh->NElements() << '\n';
        stringstream meshfilename;
        meshfilename << "DirectionalRefinementGMesh" << i;
        PrintGeometry(gmesh, meshfilename.str(), false, true);
    }
}

void UNISIMHDiv(TPZGeoMesh *gmesh) {

    static int adaptivityStep = 0;

    ProblemConfig config;
    config.porder = 1;
    config.hdivmais = 1;
    config.dimension = 2;
    config.adaptivityStep = adaptivityStep;
    config.makepressurecontinuous = true;

#ifdef DEBUGTEST
    config.dir_name = "DebugTest";
    {
        config.exact = new TLaplaceExample1;
        config.exact.operator*().fExact = TLaplaceExample1::EConst;
    }
#else
    config.dir_name = "UNISIMSpreadLevel2Refinement";
#endif
    config.problemname = "UNISIM_Errors";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    std::stringstream gmeshFileName;
    gmeshFileName << config.dir_name << "/GeoMesh" << adaptivityStep;
    PrintGeometry(gmesh, gmeshFileName.str(), false, true);

    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.bcmaterialids.insert(-2);
    config.bcmaterialids.insert(-3);

    config.gmesh = new TPZGeoMesh(*gmesh);

    TPZMultiphysicsCompMesh *cmesh_HDiv = CreateMixedCMesh(config);
    cmesh_HDiv->InitializeBlock();

    // Hybridizes mesh
    TPZHybridizeHDiv hybrid;
    auto HybridMesh = hybrid.Hybridize(cmesh_HDiv);
    delete cmesh_HDiv->MeshVector()[0];
    delete cmesh_HDiv->MeshVector()[1];
    delete cmesh_HDiv;
    HybridMesh->CleanUpUnconnectedNodes();
    HybridMesh->AdjustBoundaryElements();

    cmesh_HDiv = HybridMesh;

    // Solves FEM problem
    SolveMixedHybridProblem(cmesh_HDiv, config);
    std::cout << "Finished simulation!\n";
    std::cout << "Starting error estimation procedure...\n";

    TPZManVector<REAL> elementerrors;
    {
        // Estimates error
        TPZHybridHDivErrorEstimator HDivEstimate(*cmesh_HDiv);
        HDivEstimate.fProblemConfig = config;
        HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;

        HDivEstimate.fPostProcesswithHDiv = false;
#ifdef DEBUGTEST
        HDivEstimate.SetAnalyticSolution(config.exact);
#endif
        std::cout << "Reconstructing potential...\n";
        HDivEstimate.PotentialReconstruction();

        std::cout << "Computing errors...\n";
        HDivEstimate.ComputeErrors(elementerrors);
    }
    delete HybridMesh->MeshVector()[0];
    delete HybridMesh->MeshVector()[1];
    delete HybridMesh;
    delete config.gmesh;
    // h-refinement on elements with bigger errors
    hAdaptivity(gmesh, elementerrors, 0.3);
    adaptivityStep++;
}

TPZGeoMesh *CreateFlatGeoMesh() {

    std::string gmshFile = "InputData/UNISIMFlatMesh.msh";
#ifdef MACOSX
    gmshFile = "../" + gmshFile;
#endif

    TPZGmshReader gmeshReader;
    TPZGeoMesh *gmesh;

    TPZManVector<std::map<std::string, int>, 4> meshMaterialData(3);
    // BC materials
    meshMaterialData[1].insert({"ZeroFlux", -1});
    meshMaterialData[1].insert({"Productors", -2});
    meshMaterialData[1].insert({"Injectors", -3});
    meshMaterialData[1].insert({"Faults", 99}); // Won't be used
    // Domain materials
    meshMaterialData[2].insert({"RockMatrix", 1});
    meshMaterialData[2].insert({"RockMatrix2", 1});
    gmeshReader.SetDimNamePhysical(meshMaterialData);

    gmeshReader.SetFormatVersion("4.1");
    gmesh = gmeshReader.GeometricGmshMesh(gmshFile);

    // std::string filename = "InputData/UNISIMPointCloud.txt";
    // ModifyZCoordinates(gmesh, filename);

    MoveMeshToOrigin(gmesh);

    return gmesh;
}

void ModifyZCoordinates(TPZGeoMesh *gmesh, std::string &filename) {

    std::vector<double> x, y, z;
    ReadReservoirGeometryData(filename, x, y, z);

    _2D::ThinPlateSplineInterpolator<double> interp;
    interp.setData(x, y, z);

    int64_t nCoordinates = gmesh->NodeVec().NElements();
    double sum = 0.0;
    for (auto val : z) {
        sum += val;
    }
    double val_storage = sum / z.size();

    for (int icoord = 0; icoord < nCoordinates; icoord++) {
        TPZGeoNode node = gmesh->NodeVec()[icoord];
        TPZVec<REAL> co(3);
        node.GetCoordinates(co);
        double val_interp = interp(co[0], co[1]);

        if (val_interp == 0.0) {
            co[2] = val_storage;
        }
        if (val_interp > 1.0) {
            co[2] = val_interp;
        }

        gmesh->NodeVec()[icoord].SetCoord(co);
    }
}

void ReadReservoirGeometryData(const std::string &name, std::vector<double> &x, std::vector<double> &y,
                               std::vector<double> &z) {
    std::ifstream file;
    file.open(name);

    std::string line;
    int i = 1;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        char l = line[0];
        if (l != '/') {
            i = i + 1;
            int val = i % 15;
            if (val == 0) {
                double a, b, c;
                if (iss >> a >> b >> c)
                    ;
                x.push_back(a);
                y.push_back(b);
                z.push_back(c);
            }
        }
    }

    if (x.empty()) {
        std::cout << "No data read.\n";
        DebugStop();
    }

    std::cout << "File successfully read!\n";
}

void PrintGeometry(TPZGeoMesh *gmesh, const std::string &file_name, bool printTXT, bool printVTK) {
    if (printTXT) {
        std::stringstream txt_name;
        txt_name << file_name << ".txt";
        std::ofstream textfile(txt_name.str().c_str());
        gmesh->Print(textfile);
    }
    if (printVTK) {
        std::stringstream vtk_name;
        vtk_name << file_name << ".vtk";
        std::ofstream vtkfile(vtk_name.str().c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
    }
}

TPZCompMesh *CreatePressureCMesh(const ProblemConfig &problem) {
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);

    TPZNullMaterial *pressureMat = new TPZNullMaterial(1);
    cmesh->InsertMaterialObject(pressureMat);

    cmesh->SetDefaultOrder(problem.porder + problem.hdivmais);

    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    cmesh->AutoBuild();
    // This puts pressure equations after flux equations. It is necessary so we
    // don't end up with '0' values in the diagonal of the global matrix.
    int64_t n_connects = cmesh->NConnects();
    for (int64_t i = 0; i < n_connects; ++i) {
        cmesh->ConnectVec()[i].SetLagrangeMultiplier(1);
    }

    return cmesh;
}

TPZCompMesh *CreateFluxCMesh(const ProblemConfig &problem) {
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);

    TPZVecL2 *fluxMat = new TPZVecL2(1);
    fluxMat->SetDimension(cmesh->Dimension());
    cmesh->InsertMaterialObject(fluxMat);

    for (auto bcID : {-1, -2, -3}) {
        // The information here is not important.
        // The materials are needed only so the mesh creates BC elements.
        TPZFNMatrix<1, REAL> val1, val2;
        int bctype = 999;
        TPZBndCond *bc = fluxMat->CreateBC(fluxMat, bcID, bctype, val1, val2);

        cmesh->InsertMaterialObject(bc);
    }

    cmesh->SetDefaultOrder(problem.porder);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(cmesh->Dimension());
    cmesh->AutoBuild();

    cmesh->InitializeBlock();
    return cmesh;
}

TPZMultiphysicsCompMesh *CreateMixedCMesh(const ProblemConfig &problem) {
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(problem.gmesh);

    TPZMixedPoisson *mix = new TPZMixedPoisson(1, cmesh->Dimension());

#ifdef DEBUGTEST
    {
        mix->SetForcingFunction(problem.exact.operator*().ForcingFunction());
        mix->SetForcingFunctionExact(problem.exact.operator*().Exact());
    }
#endif
    TPZFMatrix<REAL> K(3, 3, 0), invK(3, 3, 0);
    K.Identity();
    invK.Identity();
    mix->SetPermeabilityTensor(K, invK);

    cmesh->InsertMaterialObject(mix);

    // Insert boundary conditions
    TPZFNMatrix<1, REAL> val1(1, 1, 1.e12); // Not used by the material
    TPZFNMatrix<1, REAL> val2(1, 1, 0.);
    const int dirichlet = 0;

    // Zero flux (reservoir boundary)
    TPZBndCond *zeroFlux = mix->CreateBC(mix, -1, dirichlet, val1, val2);
    // Productors
    val2(0, 0) = -10.;
    TPZBndCond *productors = mix->CreateBC(mix, -2, dirichlet, val1, val2);
    // Injectors
    val2(0, 0) = 20.;
    TPZBndCond *injectors = mix->CreateBC(mix, -3, dirichlet, val1, val2);
#ifdef DEBUGTEST
    productors->TPZMaterial::SetForcingFunction(problem.exact.operator*().Exact());
    injectors->TPZMaterial::SetForcingFunction(problem.exact.operator*().Exact());
#endif

    cmesh->InsertMaterialObject(zeroFlux);
    cmesh->InsertMaterialObject(productors);
    cmesh->InsertMaterialObject(injectors);

    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

    TPZManVector<int> active(2, 1);
    TPZManVector<TPZCompMesh *> meshVector(2, 0);

    meshVector[0] = CreateFluxCMesh(problem);
    meshVector[1] = CreatePressureCMesh(problem);

    TPZCompMeshTools::AdjustFluxPolynomialOrders(meshVector[0], problem.hdivmais);
    TPZCompMeshTools::SetPressureOrders(meshVector[0], meshVector[1]);

    cmesh->BuildMultiphysicsSpace(active, meshVector);
    cmesh->LoadReferences();
    bool keepmatrix = false;
    bool keeponelagrangian = true;
    TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian, keepmatrix);

    return cmesh;
}

void SolveMixedHybridProblem(TPZCompMesh *Hybridmesh, const ProblemConfig &problem) {

    TPZAnalysis an(Hybridmesh);

#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(Hybridmesh);
    strmat.SetNumThreads(0);
#else
    TPZSkylineStructMatrix strmat(Hybridmesh);
    strmat.SetNumThreads(4);
#endif

    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;

    an.Assemble();
    an.Solve();
}

void hAdaptivity(TPZGeoMesh* gmesh, TPZVec<REAL>& elementErrors, REAL thresholdRatio) {

    REAL maxError = 0.;
    std::cout << "Starting h-adaptivity procedure...\n"
              << "Number of elements before refinement: " << gmesh->NElements() << '\n';
    int64_t nelem = gmesh->NElements();
    for (int64_t iel = 0; iel < nelem; iel++) {
        if (elementErrors[iel] > maxError) {
            maxError = elementErrors[iel];
        }
    }

    REAL threshold = thresholdRatio * maxError;

    for (int64_t iel = 0; iel < nelem; iel++) {
        REAL elementError = elementErrors[iel];
        if (elementError > threshold) {
            TPZGeoEl *gel = gmesh->Element(iel);
            if (!gel) DebugStop();
            if (gel->Dimension() != gmesh->Dimension()) DebugStop();
            TPZVec<TPZGeoEl *> sons;
            if (!gel->HasSubElement()) {
                gel->Divide(sons);
            }
        }
    }

    PrintGeometry(gmesh, "gmeshBeforeSpread");
    SpreadMeshRefinement(gmesh);
    DivideLowerDimensionalElements(gmesh);

    std::cout << "Number of elements after refinement: " << gmesh->NElements() << '\n';
}

TPZGeoMesh *CreateDebugGeoMesh() {

    // I'll try making the orientation messy to test refinement
    // robustness, so the indexes are all out of place
    TPZGeoMesh* gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    int matID = 1;
    int innerBCMatId = -2;
    int outerBCMatId = -3;

    // Creates matrix with node coordinates
    const int NodeNumber = 8;
    constexpr REAL coordinates[NodeNumber][3] = {
        {-1.0, 1.0, 0.},
        { 1.0, 1.0, 0.},
        { 0.5, 0.5, 0.},
        {-1.0,-1.0, 0.},
        { 0.5,-0.5, 0.},
        {-0.5, 0.5, 0.},
        { 1.0,-1.0, 0.},
        {-0.5,-0.5, 0.}
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

    // Creates 2D elements
    TPZManVector<int64_t> nodeIDs(4);
    nodeIDs[0] = 5;
    nodeIDs[1] = 2;
    nodeIDs[2] = 1;
    nodeIDs[3] = 0;
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeIDs, matID, *gmesh);
    nodeIDs[0] = 6;
    nodeIDs[1] = 1;
    nodeIDs[2] = 2;
    nodeIDs[3] = 4;
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeIDs, matID, *gmesh);
    nodeIDs[0] = 4;
    nodeIDs[1] = 7;
    nodeIDs[2] = 3;
    nodeIDs[3] = 6;
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeIDs, matID, *gmesh);
    nodeIDs[0] = 5;
    nodeIDs[1] = 7;
    nodeIDs[2] = 3;
    nodeIDs[3] = 0;
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeIDs, matID, *gmesh);

    // Creates outer BC elements
    nodeIDs.Resize(2);
    nodeIDs[0] = 6;
    nodeIDs[1] = 1;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodeIDs, outerBCMatId, *gmesh);
    nodeIDs[0] = 0;
    nodeIDs[1] = 1;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodeIDs, outerBCMatId, *gmesh);
    nodeIDs[0] = 3;
    nodeIDs[1] = 0;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodeIDs, outerBCMatId, *gmesh);
    nodeIDs[0] = 6;
    nodeIDs[1] = 3;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodeIDs, outerBCMatId, *gmesh);

    // Creates inner BC elements
    nodeIDs[0] = 5;
    nodeIDs[1] = 2;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodeIDs, innerBCMatId, *gmesh);
    nodeIDs[0] = 7;
    nodeIDs[1] = 5;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodeIDs, innerBCMatId, *gmesh);
    nodeIDs[0] = 7;
    nodeIDs[1] = 4;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodeIDs, innerBCMatId, *gmesh);
    nodeIDs[0] = 2;
    nodeIDs[1] = 4;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodeIDs, innerBCMatId, *gmesh);

    gmesh->BuildConnectivity();

    return gmesh;
}

// Refines elements that has two (or more) refined neighbours or a neighbour that has been refined twice
void SpreadMeshRefinement(TPZGeoMesh *gmesh) {
    bool hasChanged = true;
    int dim = gmesh->Dimension();
    while (hasChanged) {
        hasChanged = false;

        TPZStack<TPZGeoEl *> gelsToRefine;
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (!gel || gel->Dimension() != dim || gel->HasSubElement()) {
                continue;
            }

            int nsides = gel->NSides();
            int ncorner = gel->NCornerNodes();
            int nRefinedNeighbours = 0;
            bool needsRefinement = false;
            for (int side = ncorner; side < nsides; side++) {
                TPZGeoElSide gelside(gel, side);
                TPZGeoElSide neighbour(gelside.Neighbour());
                while (neighbour != gelside) {
                    if (neighbour.HasSubElement() && neighbour.NSubElements() >= 1) {
                        nRefinedNeighbours++;
                        // Check if two neighbours have been refined
                        if (nRefinedNeighbours >= 2) {
                            needsRefinement = true;
                        }

                        TPZStack<TPZGeoElSide> subNeighs;
                        neighbour.GetSubElements2(subNeighs);
                        for (int i = 0; i < subNeighs.size(); i++) {
                            // Check if an neighbour has been refined twice
                            if (subNeighs[i].NSubElements() > 1) {
                                needsRefinement = true;
                                break;
                            }
                            if (needsRefinement) break;
                        }
                    }
                    if (needsRefinement) break;
                    neighbour = neighbour.Neighbour();
                }
                if (needsRefinement) {
                    gelsToRefine.Push(gel);
                    break;
                }
            }
        }
        if (gelsToRefine.size()) {
            hasChanged = true;
            for (int64_t i = 0; i < gelsToRefine.size(); i++) {
                TPZManVector<TPZGeoEl *> subEls;
                gelsToRefine[i]->Divide(subEls);
            }
        }
    }
}

void RotateGeoMesh(TPZGeoMesh *gmesh, const std::array<int, 2> swap) {
    // sanity checks
    //assert(swap[0] != swap[1]);
    assert(swap[0] >= 0 && swap[0] < 3);
    assert(swap[1] >= 0 && swap[1] < 3);

    int64_t nnodes = gmesh->NNodes();
    TPZManVector<REAL, 3> coord(3);

    for (int64_t inode = 0; inode < nnodes; inode++) {
        gmesh->NodeVec()[inode].GetCoordinates(coord);
        coord[swap[0]] = std::exchange(coord[swap[1]], coord[swap[0]]);
        gmesh->NodeVec()[inode].SetCoord(coord);
    }
}
