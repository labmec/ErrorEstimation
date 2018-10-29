/**
 * @file Poisson 3D in hexahedra with shock problem
 */
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

#include "TPZHybridHDivErrorEstimator.h"

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

std::tuple<TPZCompMesh *, TPZVec<TPZCompMesh *> > CreatePostProcessingMesh(TPZCompMesh *cmesh_HDiv, TPZVec<TPZCompMesh *> &meshvec_HDiv, TPZHybridizeHDiv &hybridize);

void Print(const FADREAL &a, std::ostream &out)
{
    out << " val " << a.val() << std::endl;
    for (int i=0; i< a.dx().size(); i++) {
        out << a.d(i) << " ";
    }
    out << std::endl;
}
void Print(const FADFADREAL &a, std::ostream &out)
{
    out << "Value ";
    Print(a.val(),out);
    out << "Derivatives\n";
    for (int i=0; i< a.dx().size(); i++) {
        Print(a.d(i),out);
    }
    out << "End\n";

}

int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    std::string meshfilename = "../BasicMesh.msh";
    ProblemConfig config;
    config.porder = 1;
    config.hdivmais = 2;
    config.makepressurecontinuous = true;
    config.exact.fExact = TLaplaceExample1::ESinSin;
    config.problemname = "SinSin_9x9_HDiv1";
    bool random_refine = false;
    // printing the solution value and its derivatives at a certain points
    if(0)
    {
        TPZManVector<REAL,3> x(3,0.25);
        
        TPZManVector<Fad<REAL>,3> xfad(x.size()), graduxy(x.size());
        TPZManVector<FADFADREAL,3> xfadfad(x.size()), uxyfadfad(1);
        for(int i=0; i<3; i++)
        {
            xfad[i] = Fad<REAL>(3,i,x[i]);
            xfadfad[i] = FADFADREAL(3,i,xfad[i]);
            for(int j=0; j<3; j++)
            {
                xfadfad[i].fastAccessDx(j) = Fad<REAL>(3,xfadfad[i].val().dx(j));
            }
        }
        std::cout << "xfadfad = \n";
        for(int i=0; i<3; i++)
        {
            Print(xfadfad[i],std::cout);
        }
        std::cout << std::endl;
        config.exact.graduxy(xfad, graduxy);
        config.exact.uxy(xfadfad, uxyfadfad);
        for(int i=0; i<3; i++)
        {
            std::cout << "xfad = ";
            Print(xfad[i],std::cout);
            std::cout << std::endl;
        }
        std::cout << "graduxy = \n";
        for(int i=0; i<3; i++)
        {
            Print(graduxy[i],std::cout);
        }
        std::cout << std::endl;
        std::cout << "uxyfadfad = \n";
        for(int i=0; i<uxyfadfad.size(); i++)
        {
            Print(uxyfadfad[i],std::cout);
        }
        REAL laplace = uxyfadfad[0].dx(0).dx(0)+uxyfadfad[0].dx(1).dx(1)+uxyfadfad[0].dx(2).dx(2);
        std::cout << "Laplacian " << laplace << std::endl;
    }
    {
        TPZGmshReader gmsh;
        gmsh.fPZMaterialId[1]["dirichlet"] = -1;
        gmsh.fPZMaterialId[1]["neuman"] = -2;
        gmsh.fPZMaterialId[2]["domain"] = 1;
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        config.bcmaterialids.insert(-2);
        TPZGeoMesh *gmesh = 0;
#ifdef MACOSX
        gmesh = gmsh.GeometricGmshMesh(meshfilename);
//        gmesh = gmsh.GeometricGmshMesh("../BasicMesh.msh");
#else
        gmesh = gmsh.GeometricGmshMesh("BasicMesh.msh");
#endif
        gmesh->SetDimension(2);
        config.gmesh = gmesh;
        {
            std::ofstream out("gmesh.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        }
    }
    if(random_refine)
    {
        int numelrefine = 10;
        int64_t nel = config.gmesh->NElements();
        if (numelrefine > nel/2) {
            numelrefine = 1;
        }
        int count = 0;
        while(count < numelrefine) {
            int64_t elindex = random()%nel;
            TPZGeoEl *gel = config.gmesh->Element(elindex);
            if(gel && gel->Dimension() == config.gmesh->Dimension() && !gel->Father())
            {
                TPZStack<TPZGeoEl *> subels;
                gel->Divide(subels);
                count++;
            }
        }
        nel = config.gmesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = config.gmesh->Element(el);
            if(gel && gel->Dimension() < config.gmesh->Dimension())
            {
                TPZGeoElSide gelside(gel,gel->NSides()-1);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while (neighbour != gelside) {
                    if (neighbour.Element()->NSubElements() != 0) {
                        TPZStack<TPZGeoEl *> subels;
                        gel->Divide(subels);
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
        }

    }
    TPZManVector<TPZCompMesh*, 2> meshvec_HDiv(2, 0);
    TPZCompMesh *cmesh_HDiv = CreateHDivMesh(config, meshvec_HDiv);
    cmesh_HDiv->InitializeBlock();
//    {
//        std::ofstream out("meshvec_HDiv_flux.txt");
//        meshvec_HDiv[0]->Print(out);
//    }
//    {
//        std::ofstream out("meshvec_HDiv_pres.txt");
//        meshvec_HDiv[1]->Print(out);
//    }

    {
        cmesh_HDiv->InitializeBlock();
        {
            std::ofstream out("cmesh_HDiv.txt");
            cmesh_HDiv->Print(out);
        }
        TPZAnalysis an(cmesh_HDiv);
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
        an.Solve();
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Pressure");
        vecnames.Push("Flux");
        an.DefineGraphMesh(2, scalnames, vecnames, "Original.vtk");
        //        meshvec_Hybrid[1]->Solution().Print("Press");
        // Post processing
        an.PostProcess(1,2);

    }
    {
        TPZManVector<TPZCompMesh *> MeshesHDiv(3);
        MeshesHDiv[0] = cmesh_HDiv;
        MeshesHDiv[1] = meshvec_HDiv[0];
        MeshesHDiv[2] = meshvec_HDiv[1];
        TPZHybridHDivErrorEstimator HDivEstimate(MeshesHDiv);
        HDivEstimate.fProblemConfig = config;
        HDivEstimate.fUpliftPostProcessMesh = 0;
        HDivEstimate.SetAnalyticSolution(&config.exact);
        TPZManVector<REAL> elementerrors;
        HDivEstimate.ComputeErrors(elementerrors);
    }
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
    cmesh->SetDefaultOrder(problem.porder+problem.hdivmais);
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
                intel->SetSideOrder(side, problem.porder + problem.hdivmais);
                intel->SetPreferredOrder(problem.porder+problem.hdivmais);
            }
        }
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if (gel->Dimension() == dim-1) {
                int side = gel->NSides() - 1;
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
                intel->SetSideOrder(side, problem.porder + problem.hdivmais);
                intel->SetPreferredOrder(problem.porder+problem.hdivmais);
            }
        }
   }
    cmesh->InitializeBlock();
    return cmesh;

}

TPZCompMesh *CreateHDivMesh(const ProblemConfig &problem, TPZVec<TPZCompMesh *> &meshvector) {
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial *mat = NULL;
    for (auto matid : problem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        mix->SetForcingFunction(problem.exact.ForcingFunction());
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
        bc->TPZMaterial::SetForcingFunction(problem.exact.Exact());
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

std::tuple<TPZCompMesh *, TPZVec<TPZCompMesh *> > CreatePostProcessingMesh(TPZCompMesh *cmesh_HDiv, TPZVec<TPZCompMesh *> &meshvec_HDiv, TPZHybridizeHDiv &hybridizer) {
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

/// Set the interface pressure to the average pressure
void ComputeAveragePressure(TPZCompMesh *pressure, TPZCompMesh *pressureHybrid, int InterfaceMatid)
{
    TPZGeoMesh *gmesh = pressure->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    pressure->LoadReferences();
    int64_t nel = pressureHybrid->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = pressureHybrid->Element(el);
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
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 1, 0);
        TPZManVector<TPZTransform<REAL> ,2> tr(2);
        tr[0] = gelside.NeighbourSideTransform(celstack[0].Reference());
        {
            TPZGeoEl *right = celstack[0].Element()->Reference();
            TPZTransform<REAL> tmp = right->SideToSideTransform(celstack[0].Side(), right->NSides()-1);
            tr[0] = tmp.Multiply(tr[0]);
        }
        if (celstack.size() == 1) {
            TPZCompElSide lowlevel = gelside.LowerLevelCompElementList2(1);
            if (!lowlevel) {
                DebugStop();
            }
            celstack.Push(lowlevel);
            tr[1] = TPZTransform<REAL>(gelside.Dimension());
            gel->BuildTransform2(gelside.Side(), lowlevel.Reference().Element(), tr[1]);
        }
        else if(celstack.size() == 2)
        {
            tr[1] = gelside.NeighbourSideTransform(celstack[1].Reference());
        }
        else
        {
            DebugStop();
        }
        {
            TPZGeoEl *right = celstack[1].Element()->Reference();
            TPZTransform<REAL> tmp = right->SideToSideTransform(celstack[1].Side(), right->NSides()-1);
            tr[1] = tmp.Multiply(tr[1]);
        }

        std::unique_ptr<TPZIntPoints> intp( gel->CreateSideIntegrationRule(gel->NSides()-1, 2*order));
        int nshape = intel->NShapeF();
        TPZFNMatrix<20,REAL> L2Mat(nshape,nshape,0.), L2Rhs(nshape,1,0.);
        TPZFNMatrix<220,REAL> phi(nshape,1,0.), dshape(dim,nshape);
        int64_t npoints = intp->NPoints();
        for (int64_t ip=0; ip<npoints; ip++) {
            TPZManVector<REAL,3> pt(dim-1,0.),pt1(dim,0.), pt2(dim,0.),sol1(1),sol2(1);
            REAL weight;
            intp->Point(ip, pt, weight);
            intel->Shape(pt, phi, dshape);
            tr[0].Apply(pt, pt1);
            tr[1].Apply(pt, pt2);
            celstack[0].Element()->Solution(pt1, 0, sol1);
            celstack[1].Element()->Solution(pt2, 0, sol2);
//            std::cout << "Values " << sol1 << " " << sol2 << std::endl;
            for (int ishape=0; ishape<nshape; ishape++) {
                L2Rhs(ishape,0) += weight*phi(ishape,0)*(sol1[0]+sol2[0])/2.;
                for (int jshape = 0; jshape<nshape; jshape++) {
                    L2Mat(ishape,jshape) += weight*phi(ishape,0)*phi(jshape,0);
                }
            }
        }
        L2Mat.SolveDirect(L2Rhs, ECholesky);
//        L2Rhs.Print("Average pressure");
        int count = 0;
        for (int ic=0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            int64_t seqnum = c.SequenceNumber();
            int64_t pos = pressureHybrid->Block().Position(seqnum);
            int ndof = c.NShape()*c.NState();
            for (int idf = 0; idf<ndof; idf++) {
                pressureHybrid->Solution()(pos+idf,0) = L2Rhs(count++);
            }
        }
    }
}
