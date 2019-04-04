/**
 * @file Poisson 3D in hexahedra with shock problem
 */
#include "TPZMarkErrorEstimation.h"
#include "pzlog.h"
#include "TPZRefPatternDataBase.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"
#include "ProblemConfig.h"

#include "mixedpoisson.h"
#include "TPZVecL2.h"
#include "pzbndcond.h"

#include "pzintel.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZHybridizeHDiv.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzskylstrmatrix.h"
#include <tuple>
#include <memory>

TPZCompMesh *CreateFluxHDivMesh(const ProblemConfig &problem);

TPZCompMesh *CreatePressureMesh(const ProblemConfig &problem);

TPZCompMesh *CreateHDivMesh(const ProblemConfig &problem, TPZVec<TPZCompMesh *> &meshvec);

void CloneMeshVec(TPZVec<TPZCompMesh *> &meshvec, TPZVec<TPZCompMesh *> &meshvec_clone);

/// Increase the approximation orders of the sides of the flux elements
void IncreaseSideOrders(TPZCompMesh *fluxmesh);

/// Set the interface pressure to the average pressure
void ComputeAveragePressure(TPZCompMesh *pressure, TPZCompMesh *pressureHybrid, int InterfaceMatid);

void UniformRefinement(int nDiv, TPZGeoMesh *gmesh);

TPZGeoMesh *LMesh();

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

void LocalNeumannProblem(TPZAnalysis &an, TPZCompMesh *mesh, int degreepK, int degreeqK);

void CalcStiff(TPZElementMatrix ek, TPZElementMatrix ef);

TPZCompMesh *EnrichedMesh(const ProblemConfig &problem,int degreeqk);
void InterElementSmoothing(TPZCompMesh *pressure, int InterfaceMatid);
void LocalDirichletProblem(TPZCompMesh *enrmesh);

std::tuple<TPZCompMesh *, TPZVec<TPZCompMesh *>>
CreatePostProcessingMesh(TPZCompMesh *cmesh_HDiv, TPZVec<TPZCompMesh *> &meshvec_HDiv, TPZHybridizeHDiv &hybridize);

bool gmshreader = false;

int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif

// Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);

    TPZGeoMesh *gmesh = nullptr;

    ProblemConfig config;
    config.porder = 1;
    config.hdivmais = 0;
    config.makepressurecontinuous = true;
    // config.forcingCte;
    config.exact.fExact = TLaplaceExample1::EArcTanSingular;//EArcTan;//ESinSinDirNonHom;//
    config.problemname = "EArcTanSingular";//"ArcTang";//"SinSin";//"SinSinNonHom";//

    //TPZMarkErrorEstimation estimation(config);

    if(gmshreader){
    std::string meshfilename = "../LMesh.msh";
    
    
    
        TPZGmshReader gmsh;
        gmsh.GetDimNamePhysical()[1]["dirichlet"] = -1;
        gmsh.GetDimNamePhysical()[1]["neuman"] = -2;
        gmsh.GetDimNamePhysical()[2]["domain"] = 1;
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        config.bcmaterialids.insert(-2);

        gmsh.SetFormatVersion("4.0");

#ifdef MACOSX
        gmesh = gmsh.GeometricGmshMesh(meshfilename);
        gmsh.PrintPartitionSummary(std::cout);
        //gmesh = gmsh.GeometricGmshMesh("../BasicMesh.msh");
#else
        gmesh = gmsh.GeometricGmshMesh("LMesh.msh");
#endif
        gmesh->SetDimension(2);
        config.gmesh = gmesh;

    } else {
        gmesh = LMesh();
    }
    {
        std::ofstream out("Sgmesh_SemRef.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    }

    int nDiv = 0;

    UniformRefinement(nDiv, gmesh);
    {
        std::ofstream out("Sgmesh_Ref.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    }


    TPZManVector<TPZCompMesh *, 2> meshvec_HDiv(2, 0);
    TPZCompMesh *cmesh_HDiv = CreateHDivMesh(config, meshvec_HDiv);//Hdiv x L2
    // cmesh_HDiv->InitializeBlock();
//    {
//        std::ofstream out("meshvec_HDiv_flux.txt");
//        meshvec_HDiv[0]->Print(out);
//    }
//    {
//        std::ofstream out("meshvec_HDiv_pres.txt");
//        meshvec_HDiv[1]->Print(out);
//    }


    cmesh_HDiv->InitializeBlock(); //IM BATMAN
//        {
//            std::ofstream out("cmesh_HDiv.txt");
//            cmesh_HDiv->Print(out);
//        }
    TPZAnalysis an(cmesh_HDiv, false);
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh_HDiv);
    strmat.SetNumThreads(0);
//        strmat.SetDecomposeType(ELDLt);
    an.SetStructuralMatrix(strmat);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmesh_HDiv);
    strmat.SetNumThreads(0);
    //        TPZSkylineStructMatrix strmat3(cmesh_HDiv);
    //        strmat3.SetNumThreads(8);
#endif

    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();//resolve o problema misto ate aqui


    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    vecnames.Push("Flux");
    an.DefineGraphMesh(2, scalnames, vecnames, "Original.vtk");
    //        meshvec_Hybrid[1]->Solution().Print("Press");

    // Post processing
    an.PostProcess(0, 2);

    //Solve a local Neumann Problem
    int degreepK = config.porder;
    int degreeqK = config.porder + 1;
    TPZCompMesh *enrmesh = EnrichedMesh(config, degreeqK);

    LocalNeumannProblem(an, enrmesh, degreepK, degreeqK);


    return 0;
}

TPZCompMesh *CreatePressureMesh(const ProblemConfig &problem) {
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial *mat = 0;
    for (auto matid : problem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    cmesh->SetDefaultOrder(problem.porder + problem.hdivmais);//ordem + hdivmais
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->AutoBuild();
    int64_t n_connects = cmesh->NConnects();
    for (int64_t i = 0; i < n_connects; ++i) {
        cmesh->ConnectVec()[i].SetLagrangeMultiplier(1);
    }
    return cmesh;
}

TPZCompMesh *CreateFluxHDivMesh(const ProblemConfig &problem) {
    int dim = problem.gmesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial *mat = NULL;
    problem.gmesh->ResetReference();
    for (auto matid : problem.materialids) {
        TPZVecL2 *mix = new TPZVecL2(matid);
        mix->SetDimension(dim);
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 1.);
        TPZBndCond *bc = mat->CreateBC(mat, matid, 0, val1, val2);
        cmesh->InsertMaterialObject(bc);
    }
    cmesh->SetDefaultOrder(problem.porder);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmesh->AutoBuild();
    if (problem.hdivmais) {
        int64_t nel = cmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if (gel->Dimension() == dim) {
                int side = gel->NSides() - 1;
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
                intel->SetSideOrder(side, problem.porder + problem.hdivmais);//seta ordem +hdivmais
                intel->SetPreferredOrder(problem.porder + problem.hdivmais);
            }
        }
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if (gel->Dimension() == dim - 1) {
                int side = gel->NSides() - 1;
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
                intel->SetSideOrder(side, problem.porder + problem.hdivmais);
                intel->SetPreferredOrder(problem.porder + problem.hdivmais);
            }
        }
    }
    cmesh->InitializeBlock();
    return cmesh;

}

TPZCompMesh *CreateHDivMesh(const ProblemConfig &problem, TPZVec<TPZCompMesh *> &meshvector) {

    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial *mat = nullptr;

    for (auto matid : problem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        // mix->SetForcingFunction(problem.exact.ForcingFunction());
        TPZAutoPointer<TPZFunction<STATE>> force1 = new TPZDummyFunction<STATE>(Forcing, 0);
        mix->SetForcingFunction(force1);
        //mix->SetForcingFunction(problem.forcingCte.);
        mix->SetInternalFlux(1);
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        if (matid == -2) {
            bctype = 0;
            val2.Zero();
        }
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        //   bc->TPZMaterial::SetForcingFunction(problem.exact.Exact());
        cmesh->InsertMaterialObject(bc);
    }
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    cmesh->AutoBuild();

    meshvector[0] = CreateFluxHDivMesh(problem);
    meshvector[1] = CreatePressureMesh(problem);
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh);
    cmesh->LoadReferences();
    bool keepmatrix = false;
    TPZCompMeshTools::CreatedCondensedElements(cmesh, true, keepmatrix);

    return cmesh;
}

void CloneMeshVec(TPZVec<TPZCompMesh *> &meshvec, TPZVec<TPZCompMesh *> &meshvec_clone) {
    for (int i = 0; i < meshvec.size(); i++) {
        meshvec_clone[i] = meshvec[i]->Clone();
    }
}

/// Increase the approximation orders of the sides of the flux elements
void IncreaseSideOrders(TPZCompMesh *fluxmesh) {
    int64_t nel = fluxmesh->NElements();
    int dim = fluxmesh->Dimension();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fluxmesh->Element(el);
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != dim) {
            continue;
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        int nc = cel->NConnects();
        int order = cel->Connect(nc - 1).Order();
        int nsides = gel->NSides();
        int ncorner = gel->NCornerNodes();
        intel->SetPreferredOrder(order);
        for (int side = ncorner; side < nsides - 1; side++) {
            if (intel->NSideConnects(side)) {
                intel->SetSideOrder(side, order);
            }
        }
//        intel->Print();
    }
    fluxmesh->InitializeBlock();
}

std::tuple<TPZCompMesh *, TPZVec<TPZCompMesh *>>
CreatePostProcessingMesh(TPZCompMesh *cmesh_HDiv, TPZVec<TPZCompMesh *> &meshvec_HDiv, TPZHybridizeHDiv &hybridizer) {
    TPZManVector<TPZCompMesh *, 2> meshvec_Hybrid(2, 0);
    CloneMeshVec(meshvec_HDiv, meshvec_Hybrid);
    IncreaseSideOrders(meshvec_Hybrid[0]);
    hybridizer.ComputePeriferalMaterialIds(meshvec_Hybrid);
    hybridizer.ComputeNState(meshvec_Hybrid);
    /// insert the material objects for HDivWrap and LagrangeInterface
    hybridizer.InsertPeriferalMaterialObjects(meshvec_Hybrid);
    hybridizer.HybridizeInternalSides(meshvec_Hybrid);
    TPZCompMesh *cmesh_Hybrid = hybridizer.CreateMultiphysicsMesh(cmesh_HDiv, meshvec_Hybrid);
    hybridizer.CreateInterfaceElements(cmesh_Hybrid, meshvec_Hybrid);
    hybridizer.GroupElements(cmesh_Hybrid);
    return std::make_tuple(cmesh_Hybrid, meshvec_Hybrid);
}


void UniformRefinement(int nDiv, TPZGeoMesh *gmesh) {

    TPZManVector<TPZGeoEl *> children;
    for (int division = 0; division < nDiv; division++) {

        int64_t nels = gmesh->NElements();

        for (int64_t elem = 0; elem < nels; elem++) {

            TPZGeoEl *gel = gmesh->ElementVec()[elem];

            if (!gel || gel->HasSubElement()) continue;
            if (gel->Dimension() == 0) continue;
            gel->Divide(children);
        }
    }
}

TPZGeoMesh *LMesh() {

    int dim = 2;
    int nodes = 8;
    int nel = 6;
    //nos
    REAL coord[8][3] =
            {
                    {0,   0,  0},
                    {1.,  0,  0},
                    {1,   1., 0},
                    {0.,  1., 0},
                    {-1., 1., 0},
                    {-1,  0,  0},
                    {-1,  -1, 0},
                    {0,   -1, 0},

            };
    //elementos
    int elno[6][3] =
            {
                    {0, 1, 2},
                    {0, 2, 3},
                    {0, 3, 5},
                    {3, 4, 5},
                    {0, 5, 6},
                    {0, 6, 7},
            };


    TPZGeoMesh *geomesh = new TPZGeoMesh;


    geomesh->NodeVec().Resize(nodes);
    geomesh->SetDimension(dim);

    for (int no = 0; no < nodes; no++) {
        TPZManVector<REAL, 3> xco(3, 0.);
        xco[0] = coord[no][0];
        xco[1] = coord[no][1];
        xco[2] = coord[no][2];
        geomesh->NodeVec()[no].Initialize(xco, *geomesh);
    }
    int64_t index;
    for (int el = 0; el < nel; el++) {

        TPZVec<int64_t> elnodes(3);
        for (int i = 0; i < 3; i++) {
            elnodes[i] = elno[el][i];
            REAL valor = elno[el][i];

        }
        geomesh->CreateGeoElement(ETriangle, elnodes, 1, index);

    }
    //condicao de contorno
    TPZManVector<int64_t> nodeIDs(3);

    //el 0
    nodeIDs[0] = 0;
    nodeIDs[1] = 1;
    geomesh->CreateGeoElement(EOned, nodeIDs, -1, index);
    //el 0

    nodeIDs[0] = 1;
    nodeIDs[1] = 2;
    geomesh->CreateGeoElement(EOned, nodeIDs, -2, index);
    //el 01

    nodeIDs[0] = 2;
    nodeIDs[1] = 3;
    geomesh->CreateGeoElement(EOned, nodeIDs, -3, index);
    //el 03

    nodeIDs[0] = 3;
    nodeIDs[1] = 4;
    geomesh->CreateGeoElement(EOned, nodeIDs, -4, index);
    //el 03

    nodeIDs[0] = 4;
    nodeIDs[1] = 5;
    geomesh->CreateGeoElement(EOned, nodeIDs, -5, index);
    //el 04

    nodeIDs[0] = 5;
    nodeIDs[1] = 6;
    geomesh->CreateGeoElement(EOned, nodeIDs, -6, index);
    //el 05

    nodeIDs[0] = 6;
    nodeIDs[1] = 7;
    geomesh->CreateGeoElement(EOned, nodeIDs, -7, index);
    //el 05
    nodeIDs[0] = 7;
    nodeIDs[1] = 0;
    geomesh->CreateGeoElement(EOned, nodeIDs, -8, index);


    geomesh->BuildConnectivity();

    return geomesh;


}

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp) {
    disp[0] = 1.;
    return;
}

void LocalNeumannProblem(TPZAnalysis &an, TPZCompMesh *mesh, int degreepK, int degreeqK) {

    int nelcomp = mesh->NElements();

    TPZElementMatrix ek, ef;
    for (int iel = 0; iel < nelcomp; iel++) {
        TPZCompEl *compEl = mesh->ElementVec()[iel];
        if (!compEl) continue;
        TPZGeoEl *gel = compEl->Reference();

        if (gel->Dimension() != 2) continue;

        compEl->CalcStiff(ek, ef);
    }
}


TPZCompMesh *EnrichedMesh(const ProblemConfig &problem, int degreeqk) {
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);

    TPZMatPoisson3d *mix = new TPZMixedPoisson(1, 2);
    mix->SetNeumannProblem();
    cmesh->InsertMaterialObject(mix);
    cmesh->SetDefaultOrder(degreeqk);
    cmesh->SetAllCreateFunctionsContinuous();

    TPZMatPoisson3d *material = new TPZMixedPoisson(1, degreeqk);
    material->SetNeumannProblem();
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);

    cmesh-> SetDefaultOrder(degreeqk);
    cmesh-> SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();

    return cmesh;
}

void InterElementSmoothing(TPZCompMesh *pressure, int InterfaceMatid){
    //tomar a media da pressao calcula no NeumannProblem
    
    TPZGeoMesh *gmesh = pressure->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    pressure->LoadReferences();
    int64_t nel = pressure->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = pressure->Element(el);
        if(!cel || !cel->Reference() || cel->Reference()->Dimension() != dim-1)
        {
            continue;
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        TPZGeoEl *gel = cel->Reference();
        if (gel->MaterialId() != InterfaceMatid) {
            continue;
        }
        if (!intel || gel->Dimension() != dim-1) {
            DebugStop();
        }
        int nc = cel->NConnects();
        int order = cel->Connect(nc-1).Order();
    }

    
    
    
}
void LocalDirichletProblem(TPZCompMesh *enrmesh){
    //REsolver um problema de Dirichlet em que a condicao de contorno é a solucao suavizada
    
    
    
    
}
