//
//  Tools.cpp
//  ErrorEstimation
//
//  Created by Denise De Siqueira on 28/03/19.
//

#include "Tools.h"
#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "TPZGeoLinear.h"
#include "TPZGenGrid2D.h"
#include <tuple>
#include <memory>
#include <Pre/TPZGenGrid3D.h>
#include "DataStructure.h"

#include "pzelementgroup.h"

void Tools::PrintGeometry(TPZGeoMesh *gmesh, const std::string &file_name, bool printTXT, bool printVTK) {
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

TPZCompMesh* Tools::CreatePressureMesh(const ProblemConfig& problem) {
    TPZCompMesh* cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial* mat = 0;
    for (auto matid : problem.materialids) {
        TPZMixedPoisson* mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    
    cmesh->SetDefaultOrder(problem.porder + problem.hdivmais);
    //cmesh->SetDefaultOrder(problem.porder);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->AutoBuild();
    int64_t n_connects = cmesh->NConnects();
    for (int64_t i = 0; i < n_connects; ++i) {
        cmesh->ConnectVec()[i].SetLagrangeMultiplier(1);
    }
    
    if (problem.prefine) {
        Prefinamento(cmesh, problem.ndivisions, problem.porder);
    }
    
    
    return cmesh;
}

TPZCompMesh* Tools::CreateFluxHDivMesh(const ProblemConfig& problem) {
    int dim = problem.gmesh->Dimension();
    TPZCompMesh* cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial* mat = NULL;
    problem.gmesh->ResetReference();
    for (auto matid : problem.materialids) {
        TPZVecL2* mix = new TPZVecL2(matid);
        mix->SetDimension(dim);
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 1.);
        int bctype;
        if (matid == -1 || matid == 2) {
            bctype = 0;
        } else {
            bctype = 1;
        }
        TPZBndCond* bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(problem.exact.operator*().Exact());
        cmesh->InsertMaterialObject(bc);
    }
    cmesh->SetDefaultOrder(problem.porder);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmesh->AutoBuild();
    
    if (problem.prefine) {
        Prefinamento(cmesh, problem.ndivisions, problem.porder);
    }
    
    
    cmesh->InitializeBlock();
    return cmesh;
    
}

TPZMultiphysicsCompMesh* Tools::CreateHDivMesh(const ProblemConfig& problem) {

    TPZMultiphysicsCompMesh* cmesh = new TPZMultiphysicsCompMesh(problem.gmesh);
    TPZMaterial* mat = NULL;
    TPZFMatrix<REAL> K(3, 3, 0), invK(3, 3, 0);
    K.Identity();
    invK.Identity();
    
    STATE Km = problem.Km;

    
    if (problem.TensorNonConst && problem.gmesh->Dimension() == 3) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (i == j) {
                    K(i, j) = 2.;
                    invK(i, j) = 3. / 4.;
                } else {
                    K(i, j) = 1.;
                    invK(i, j) = (-1.) / 4.;
                }
            }
        }
        
    }

//    K.Print(std::cout);
//    invK.Print(std::cout);
    
    for (auto matid : problem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        mix->SetForcingFunction(problem.exact.operator*().ForcingFunction());
        mix->SetExactSol(problem.exact.operator*().Exact());
        mix->SetPermeabilityTensor(K, invK);

        if (!mat) mat = mix;

        cmesh->InsertMaterialObject(mix);
    }
        
    
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype;
    
        switch (matid) {
            case -1 :{
            bctype = 0;
                break;
            }
                
                
            case -2:{
            bctype = 1;
    
            break;
            }
            case -3:{
            bctype = 4;// different from mixed (bctype 2) already implemented on TPZMixedPoisson3d
            val1(0,0) = Km ;

                
            break;
            }
        }
        TPZBndCond* bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(problem.exact.operator*().Exact());
        cmesh->InsertMaterialObject(bc);
    }
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    
    TPZManVector<int> active(2, 1);
    TPZManVector<TPZCompMesh*> meshvector(2, 0);
    
    meshvector[0] = CreateFluxHDivMesh(problem);
    meshvector[1] = CreatePressureMesh(problem);
    
    TPZCompMeshTools::AdjustFluxPolynomialOrders(meshvector[0], problem.hdivmais);
    TPZCompMeshTools::SetPressureOrders(meshvector[0], meshvector[1]);
    
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    bool keepmatrix = false;
    bool keeponelagrangian = true;
    TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian, keepmatrix);
    cmesh->InitializeBlock();
    return cmesh;
}

void Tools::UniformRefinement(int nDiv, TPZGeoMesh* gmesh) {
    
    TPZManVector<TPZGeoEl*> children;
    for (int division = 0; division < nDiv; division++) {
        
        int64_t nels = gmesh->NElements();
        
        for (int64_t elem = 0; elem < nels; elem++) {
            
            TPZGeoEl* gel = gmesh->ElementVec()[elem];
            
            if (!gel || gel->HasSubElement()) continue;
            if (gel->Dimension() == 0) continue;
            gel->Divide(children);
        }
    }
}

// This overload takes the dimension of the elements to be refined.
// The function DivideLowerDimensionalElements must be called afterwards to guarantee mesh consistency.
void Tools::UniformRefinement(int nDiv, int dim, TPZGeoMesh* gmesh) {

    TPZManVector<TPZGeoEl*> children;
    for (int division = 0; division < nDiv; division++) {

        int64_t nels = gmesh->NElements();

        for (int64_t elem = 0; elem < nels; elem++) {

            TPZGeoEl* gel = gmesh->Element(elem);

            if (!gel || gel->HasSubElement()) continue;
            if (gel->Dimension() != dim) continue;
            gel->Divide(children);
        }
    }
}

TPZGeoMesh* Tools::CreateNewGeoMesh(int nel, TPZVec<int>& bcids) {

    TPZManVector<int> nx(2, nel);
    TPZManVector<REAL> x0(3, -1.), x1(3, 1.);
    x1[2] = 0.;
    TPZGenGrid2D gen(nx, x0, x1, 1, 0);

    gen.SetRefpatternElements(true);
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    gen.Read(gmesh);
    gen.SetBC(gmesh, 4, bcids[0]);
    gen.SetBC(gmesh, 5, bcids[1]);
    gen.SetBC(gmesh, 6, bcids[2]);
    gen.SetBC(gmesh, 7, bcids[3]);

    gmesh->SetDimension(2);

    return gmesh;
}

TPZGeoMesh* Tools::CreateGeoMesh(int nel, TPZVec<int>& bcids) {
    
    TPZManVector<int> nx(2, nel);
    TPZManVector<REAL> x0(3, 0.), x1(3, 1.);
    x1[2] = 0.;
    TPZGenGrid2D gen(nx, x0, x1, 1, 0);
    
    gen.SetRefpatternElements(true);
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    gen.Read(gmesh);
    gen.SetBC(gmesh, 4, bcids[0]);
    gen.SetBC(gmesh, 5, bcids[1]);
    gen.SetBC(gmesh, 6, bcids[2]);
    gen.SetBC(gmesh, 7, bcids[3]);
    
    gmesh->SetDimension(2);
    
    return gmesh;
}

TPZGeoMesh *Tools::CreateCubeGeoMesh(const TPZVec<int> &nelDiv, const TPZVec<int> &bcids) {
    TPZManVector<REAL, 3> x0(3, 0.), x1(3, 1.);
    x1[2] = 0.5;
    TPZGenGrid3D gen(x0, x1, nelDiv, MMeshType::EHexahedral);

    TPZGeoMesh *gmesh{nullptr};
    gmesh = gen.BuildVolumetricElements(1);
    gmesh = gen.BuildBoundaryElements(bcids[0], bcids[1], bcids[2], bcids[3], bcids[4], bcids[5]);

    return gmesh;
}

void Tools::RandomRefinement(TPZGeoMesh *gmesh, int64_t numelrefine, int depth) {
    
    int64_t nel = gmesh->NElements();
    if (numelrefine > nel / 2) {
        numelrefine = nel/2;
    }
    for (int id=0; id<depth; id++)
    {
        nel = gmesh->NElements();
        int count = 0;
        while (count < numelrefine) {
            int64_t elindex = rand() % nel;
            TPZGeoEl* gel = gmesh->Element(elindex);
            int level = gel->Level();
            if (gel && gel->Dimension() == gmesh->Dimension() && gel->Level() == id) {
                TPZStack<TPZGeoEl*> subels;
                gel->Divide(subels);
                count++;
            }
        }
    }
    nel = gmesh->NElements();
    bool changed = true;
    while (changed) {
        changed = false;
        nel = gmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZGeoEl* gel = gmesh->Element(el);
            if (gel && gel->Dimension() < gmesh->Dimension()) {
                TPZGeoElSide gelside(gel, gel->NSides() - 1);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while (neighbour != gelside) {
                    if (neighbour.Element()->HasSubElement() != 0 && !gel->HasSubElement()) {
                        TPZStack<TPZGeoEl*> subels;
                        gel->Divide(subels);
                        changed = true;
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
        }
    }
}

void Tools::RefineElements(TPZGeoMesh *gmesh, const std::set<int64_t>& elsToRefine) {

    for (const auto it : elsToRefine) {
        TPZGeoEl *gel = gmesh->Element(it);
        if (!gel) DebugStop();
        if (gel->Dimension() == gmesh->Dimension()) {
            TPZStack<TPZGeoEl *> subels;
            gel->Divide(subels);
        }
    }
    int64_t nel;
    bool changed = true;
    while (changed) {
        changed = false;
        nel = gmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (gel && gel->Dimension() < gmesh->Dimension()) {
                TPZGeoElSide gelside(gel, gel->NSides() - 1);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while (neighbour != gelside) {
                    if (neighbour.Element()->HasSubElement() != 0 && !gel->HasSubElement()) {
                        TPZStack<TPZGeoEl *> subels;
                        gel->Divide(subels);
                        changed = true;
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
        }
    }
}

void Tools::Print(const FADREAL& a, std::ostream& out) {
    out << " val " << a.val() << std::endl;
    for (int i = 0; i < a.dx().size(); i++) {
        out << a.d(i) << " ";
    }
    out << std::endl;
}

void Tools::Print(const FADFADREAL& a, std::ostream& out) {
    out << "Value ";
    Print(a.val(), out);
    out << "Derivatives\n";
    for (int i = 0; i < a.dx().size(); i++) {
        Print(a.d(i), out);
    }
    out << "End\n";

}

void Tools::Prefinamento(TPZCompMesh* cmesh, int ndiv, int porder) {
    if (ndiv < 1) return;
    int nel = cmesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZCompEl* cel = cmesh->ElementVec()[iel];
        if (!cel) continue;

        TPZInterpolationSpace* sp = dynamic_cast<TPZInterpolationSpace*>(cel);
        if (!sp) continue;
        int level = sp->Reference()->Level();
        TPZGeoEl* gel = sp->Reference();
        if ((gel->Dimension() == 2) && (iel % 2 == 0)) {
            int ordem = 0;
            ordem = porder + (ndiv - 1) + (level);
            std::cout << "level " << level << " ordem " << ordem << std::endl;
            sp->PRefine(ordem);
        }
    }
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();

    std::stringstream sout;
    sout << "malha computacional apos pRefinamento\n";
    cmesh->Print(sout);

}

void Tools::SolveHybridProblem(TPZCompMesh *Hybridmesh, std::pair<int, int> InterfaceMatId, const ProblemConfig &problem,
                   bool PostProcessingFEM) {

    TPZAnalysis an(Hybridmesh);

#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(Hybridmesh);
    strmat.SetNumThreads(0);
    //        strmat.SetDecomposeType(ELDLt);
#else
    //    TPZFrontStructMatrix<TPZFrontSym<STATE> > strmat(Hybridmesh);
    //    strmat.SetNumThreads(2);
    //    strmat.SetDecomposeType(ELDLt);
    TPZSkylineStructMatrix strmat(Hybridmesh);
    strmat.SetNumThreads(0);
#endif

    std::set<int> matIds;

    for (auto matid : problem.materialids) {
        matIds.insert(matid);
    }

    for (auto matidbc : problem.bcmaterialids) {
        matIds.insert(matidbc);
    }

    matIds.insert(InterfaceMatId.first);
    matIds.insert(InterfaceMatId.second);

    strmat.SetMaterialIds(matIds);

    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();

    if (PostProcessingFEM) {
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("ExactPressure");
        scalnames.Push("Pressure");
        vecnames.Push("ExactFlux");
        vecnames.Push("Flux");

        std::stringstream sout;
        sout << problem.dir_name << "/" << "OriginalHybrid_Order_" << problem.porder << "Nref_" << problem.ndivisions
             << "NAdapStep_" << problem.adaptivityStep << ".vtk";
        an.DefineGraphMesh(2, scalnames, vecnames, sout.str());
        int resolution = 0;
        an.PostProcess(resolution, Hybridmesh->Dimension());

        if (problem.exact.operator*().Exact()) {
            TPZManVector<REAL> errors(5, 0.);
            an.SetThreadsForError(0);
            an.SetExact(problem.exact.operator*().ExactSolution());
            an.PostProcessError(errors, false);
            /*Error on MixedPoisson
             [0] L2 for pressure
             [1] L2 for flux
             [2] L2 for div(flux)
             [3] Grad pressure (Semi H1)
             [4] Hdiv norm
             */

            // Erro
            std::ofstream myfile;
            myfile.open("ErrorBCFemProblem.txt", std::ios::app);

            myfile << "\n\n Error for Mixed formulation ";
            myfile << "\n-------------------------------------------------- \n";
            myfile << "Ndiv = " << problem.ndivisions
                   << " Order k = " << problem.porder << " n "<<problem.hdivmais<< " K_R = "<<problem.Km<<" Ndofs = "<<Hybridmesh->NEquations() <<"\n";
            myfile << "L2 pressure = " << errors[0] << "\n";
            myfile << "L2 flux= " << errors[1] << "\n";
            myfile << "L2 div(flux) = " << errors[2] << "\n";
            myfile << "L2 flux +L2 pressure "<< errors[0] + errors[1]<<"\n";
          //  myfile << "Semi H1 = " << errors[3] << "\n";
           // myfile << "Hdiv norm = " << errors[4] << "\n";

            myfile.close();
        }

    }

    {
        std::ofstream out("SolvedHDivMesh.txt");
        Hybridmesh->Print(out);
        std::ofstream outFlux("SolvedFluxMesh.txt");
        TPZMultiphysicsCompMesh *multiphysicsMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Hybridmesh);
        if (!multiphysicsMesh) DebugStop();
        multiphysicsMesh->MeshVector()[0]->Print(outFlux);

        std::stringstream solByElement;
        TPZCompMeshTools::PrintSolutionByGeoElement(multiphysicsMesh->MeshVector()[1], solByElement);

        std::ofstream solByElementFile("SolByElementPressureHDiv.txt");
        solByElementFile << solByElement.str();
    }

}

void Tools::SolveMixedProblem(TPZCompMesh* cmesh_HDiv, const ProblemConfig& config) {
#ifdef PZDEBUG
    {
        std::ofstream out("gmeshSolve.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, out);

    }
#endif


    TPZAnalysis an(cmesh_HDiv, false);


    TPZSymetricSpStructMatrix strmat(cmesh_HDiv);
    strmat.SetNumThreads(0);
    an.SetStructuralMatrix(strmat);

    std::set<int> matids;

    for (auto mat:config.materialids) {
        matids.insert(mat);
    }

    for (auto mat:config.bcmaterialids) {
        matids.insert(mat);
    }

    strmat.SetMaterialIds(matids);
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE>* direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();

    an.Solve();//resolve o problema misto ate aqui

    TPZFNMatrix<20, REAL> pressure(7,1,0);

    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    scalnames.Push("ExactPressure");
    vecnames.Push("Flux");
    vecnames.Push("ExactFlux");

    int dim = config.gmesh->Dimension();

    std::stringstream sout;

    sout << config.dir_name
         << "/"
            "OriginalMixed_Order_"
         << config.problemname << "Order" << config.porder << "NAdapStep_"
         << config.adaptivityStep << ".vtk";

    an.DefineGraphMesh(dim, scalnames, vecnames, sout.str());
    int resolution = 2;
    an.PostProcess(resolution, dim);

    if (config.exact.operator*().Exact()) {
        TPZManVector<REAL> errors(4, 0.);
        an.SetThreadsForError(0);
        an.SetExact(config.exact.operator*().ExactSolution());
        an.PostProcessError(errors, false);

        // Erro
        std::ofstream myfile;
        /*Error on MixedPoisson
           [0] L2 for pressure
           [1] L2 for flux
           [2] L2 for div(flux)
           [3] Grad pressure (Semi H1)
           [4] Hdiv norm
           */

          // Erro
          myfile.open("ErrorMixed.txt", std::ios::app);
          myfile << "\n\n Error for Mixed formulation ";
          myfile << "\n-------------------------------------------------- \n";
          myfile << "Ndiv = " << config.ndivisions
                 << " Order k = " << config.porder << " n "<<config.hdivmais<< " K_R = "<<config.Km<<" Ndofs = "<<cmesh_HDiv->NEquations() <<"\n";
          myfile << "L2 pressure = " << errors[0] << "\n";
          myfile << "L2 flux= " << errors[1] << "\n";
          myfile << "L2 div(flux) = " << errors[2] << "\n";
        //  myfile << "Semi H1 = " << errors[3] << "\n";
         // myfile << "Hdiv norm = " << errors[4] << "\n";
        myfile.close();
    }
}

/// Divide lower dimensional elements
void Tools::DivideLowerDimensionalElements(TPZGeoMesh* gmesh) {
    bool haschanged = true;
    int dim = gmesh->Dimension();
    while (haschanged) {
        haschanged = false;
        int64_t nel = gmesh->NElements();
        TPZStack<TPZGeoEl*> geldivide;
        for (int64_t el = 0; el < nel; el++) {
            TPZGeoEl* gel = gmesh->Element(el);
            if (!gel || gel->Dimension() == dim || gel->Dimension() == 0) {
                continue;
            }
            if (gel->HasSubElement()) continue;
            int nsides = gel->NSides();
            int ncorner = gel->NCornerNodes();
            for (int side = ncorner; side < nsides; side++) {
                TPZGeoElSide gelside(gel, side);
                TPZGeoElSide neighbour(gelside.Neighbour());
                bool found = false;
                while (neighbour != gelside) {
                    if (neighbour.HasSubElement())
                    {
                        if(neighbour.NSubElements() > 1)
                        {
                            geldivide.Push(gel);
                            found = true;
                            break;
                        }
                    }
                    neighbour = neighbour.Neighbour();
                }
                if (found) break;
            }
        }
        if (geldivide.size()) {
            haschanged = true;
            for (int64_t i = 0; i < geldivide.size(); i++) {
                TPZManVector<TPZGeoEl*> sub;
                geldivide[i]->Divide(sub);
            }
        }
    }
}

TPZCompMesh* Tools::CMeshH1(ProblemConfig problem) {

    TPZCompMesh* cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial* mat = 0;


    for (auto matid : problem.materialids) {
        TPZMatPoisson3d *mix = new TPZMatPoisson3d(matid, cmesh->Dimension());
        mix->SetExactSol(problem.exact.operator*().Exact());
        mix->SetForcingFunction(problem.exact.operator*().ForcingFunction());

        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }

    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        val2.Zero();
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(problem.exact.operator*().Exact());

        cmesh->InsertMaterialObject(bc);
    }

    cmesh->SetDefaultOrder(problem.porder);//ordem
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();

    cmesh->AutoBuild();


    return cmesh;
}

void Tools::hAdaptivity(TPZCompMesh* postProcessMesh, TPZGeoMesh* gmeshToRefine, ProblemConfig& config) {

    // Column of the flux error estimate on the element solution matrix
    const int fluxErrorEstimateCol = 3;

    int64_t nelem = postProcessMesh->ElementSolution().Rows();

    //postProcessMesh->ElementSolution().Print("ElSolutionForAdaptivity",std::cout);

    // Iterates through element errors to get the maximum value
    REAL maxError = 0.;
    for (int64_t iel = 0; iel < nelem; iel++) {
        TPZCompEl* cel = postProcessMesh->ElementVec()[iel];
        if (!cel) continue;
        if (cel->Dimension() != postProcessMesh->Dimension()) continue;
        REAL elementError = postProcessMesh->ElementSolution()(iel, fluxErrorEstimateCol);


        if (elementError > maxError) {
            maxError = elementError;
        }
    }

    std::cout << "max error " << maxError << "\n";

    // Refines elements which error are bigger than 30% of the maximum error
    REAL threshold = 0.2 * maxError;

    for (int64_t iel = 0; iel < nelem; iel++) {
        TPZCompEl* cel = postProcessMesh->ElementVec()[iel];
        if (!cel) continue;
        if (cel->Dimension() != postProcessMesh->Dimension()) continue;

        REAL elementError = postProcessMesh->ElementSolution()(iel, fluxErrorEstimateCol);
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
    DivideLowerDimensionalElements(gmeshToRefine);
}

TPZGeoMesh* Tools::CreateLCircleGeoMesh() {

    TPZGeoMesh* gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);

    TPZVec<REAL> coord(3, 0.);

    // Inserts node at origin
    gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[0].Initialize(coord, *gmesh);

    // Inserts circumference nodes
    for (int64_t i = 0; i < 13; i++) {
        const REAL step = M_PI / 8;
        coord[0] = cos(i * step);
        coord[1] = sin(i * step);
        const int64_t newID = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[newID].Initialize(coord, *gmesh);
    }

    int matIdTriangle = 1, matIdArc = 2;

    // Inserts triangle elements
    TPZManVector<int64_t, 3> nodesIdVec(3);
    for (int64_t i = 0; i < 6; i++) {
        nodesIdVec[0] = 0;
        nodesIdVec[1] = 1 + 2 * i;
        nodesIdVec[2] = 3 + 2 * i;
        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle>>(nodesIdVec, matIdTriangle, *gmesh);
    }
    // Inserts arc elements
    for (int64_t i = 0; i < 6; i++) {
        nodesIdVec[0] = 1 + 2 * i;
        nodesIdVec[1] = 3 + 2 * i;
        nodesIdVec[2] = 2 + 2 * i;
        new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, matIdArc, *gmesh);
    }
    // Finally, inserts line elements to complete boundary
    nodesIdVec.Resize(2);
    nodesIdVec[0] = 0;
    nodesIdVec[1] = 1;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, matIdArc, *gmesh);

    nodesIdVec[0] = 0;
    nodesIdVec[1] = 13;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, matIdArc, *gmesh);

    gmesh->BuildConnectivity();
    return gmesh;
}

TPZGeoMesh* Tools::CreateLShapeMesh(TPZVec<int>& bcids) {

    TPZGeoMesh* gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    int matID = 1;

    // Creates matrix with node coordinates
    const int NodeNumber = 8;
    REAL coordinates[NodeNumber][3] = {
        {0.,  0., 0.},
        {1.,  0., 0.},
        {1.,  1., 0.},
        {0.,  1., 0.},
        {-1., 1., 0.},
        {-1., 0., 0.},
        {-1.,-1., 0.},
        {0., -1., 0.}
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
    bool blob = false;
    if (blob) {
        // Creates 2D elements
        TPZManVector<int64_t> nodeIDs(3);
        nodeIDs[0] = 0;
        nodeIDs[1] = 1;
        nodeIDs[2] = 3;
        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodeIDs, matID, *gmesh);
        nodeIDs[0] = 2;
        nodeIDs[1] = 3;
        nodeIDs[2] = 1;
        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodeIDs, matID, *gmesh);
        nodeIDs[0] = 3;
        nodeIDs[1] = 4;
        nodeIDs[2] = 0;
        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodeIDs, matID, *gmesh);
        nodeIDs[0] = 5;
        nodeIDs[1] = 0;
        nodeIDs[2] = 4;
        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodeIDs, matID, *gmesh);
        nodeIDs[0] = 0;
        nodeIDs[1] = 5;
        nodeIDs[2] = 7;
        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodeIDs, matID, *gmesh);
        nodeIDs[0] = 6;
        nodeIDs[1] = 7;
        nodeIDs[2] = 5;
        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodeIDs, matID, *gmesh);
    }
    else {
        // Creates 2D elements
        TPZManVector<int64_t> nodeIDs(4);
        nodeIDs[0] = 0;
        nodeIDs[1] = 1;
        nodeIDs[2] = 2;
        nodeIDs[3] = 3;
        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeIDs, matID, *gmesh);
        nodeIDs[0] = 0;
        nodeIDs[1] = 3;
        nodeIDs[2] = 4;
        nodeIDs[3] = 5;
        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeIDs, matID, *gmesh);
        nodeIDs[0] = 0;
        nodeIDs[1] = 5;
        nodeIDs[2] = 6;
        nodeIDs[3] = 7;
        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeIDs, matID, *gmesh);
    }

    // Creates line elements where boundary conditions will be inserted
    TPZManVector<int64_t> nodeIDs(2);
    for (int i = 0; i < NodeNumber; i++) {
        nodeIDs[0] = i % NodeNumber;
        nodeIDs[1] = (i + 1) % NodeNumber;
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodeIDs, bcids[i], *gmesh);
    }

    gmesh->BuildConnectivity();

    return gmesh;

}

TPZGeoMesh* Tools::CreateQuadLShapeMesh(TPZVec<int>& bcids) {

    TPZGeoMesh* gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    int matID = 1;

    // Creates matrix with node coordinates
    const int NodeNumber = 8;
    REAL coordinates[NodeNumber][3] = {
            {0.,  0., 0.},
            {1.,  0., 0.},
            {1.,  1., 0.},
            {0.,  1., 0.},
            {-1., 1., 0.},
            {-1., 0., 0.},
            {-1.,-1., 0.},
            {0., -1., 0.}
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
    for (int i = 0; i < 3; i++) {
        nodeIDs[0] = 0;
        nodeIDs[1] = (2 * i + 1) % NodeNumber;
        nodeIDs[2] = (2 * i + 2) % NodeNumber;
        nodeIDs[3] = (2 * i + 3) % NodeNumber;
        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeIDs, matID, *gmesh);
    }

    // Creates line elements where boundary conditions will be inserted
    nodeIDs.Resize(2);
    for (int i = 0; i < NodeNumber; i++) {
        nodeIDs[0] = i % NodeNumber;
        nodeIDs[1] = (i + 1) % NodeNumber;
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodeIDs, bcids[i], *gmesh);
    }

    gmesh->BuildConnectivity();

    return gmesh;

}

TPZGeoMesh* Tools::CreateGeoMesh(int nel, TPZVec<int>& bcids, int dim, bool isOriginCentered, int topologyMode) {

    if (dim == 2){
        TPZManVector<int> nx(2, nel);
        TPZManVector<REAL> x0(3, 0.), x1(3, 1.);
        if(isOriginCentered == 1){
            x0[0]= x0[1] = -1;
        }
        x1[2] = x0[2] = 0.;

        TPZGenGrid2D gen(nx, x0, x1, 1, 0);
        MMeshType eltype;
        if(topologyMode == 1) eltype = MMeshType::ETriangular;
        else if(topologyMode == 2) eltype = MMeshType::EQuadrilateral;
        else DebugStop();
        gen.SetElementType(eltype);

        //TPZGenGrid2D gen(nx, x0, x1);
        gen.SetRefpatternElements(true);
        TPZGeoMesh* gmesh = new TPZGeoMesh;
        gen.Read(gmesh);
        gen.SetBC(gmesh, 4, bcids[0]);
        gen.SetBC(gmesh, 5, bcids[1]);
        gen.SetBC(gmesh, 6, bcids[2]);
        gen.SetBC(gmesh, 7, bcids[3]);

        gmesh->SetDimension(2);

        return gmesh;
    }
    if (dim == 3){
        int volID = 1;
        int bcID = -1;
        TPZManVector<int> nx(2, nel);
        TPZManVector<REAL> x0(3, 0.), x1(3, 1.);
        if(isOriginCentered == 1){
            x0[0]= x0[1] = x0[2]  = -1;
        }

        TPZManVector<int> nelDiv(3, 1);
        MMeshType  eltype;
        if(topologyMode == 3) eltype = MMeshType::ETetrahedral;
        else if(topologyMode == 4) eltype = MMeshType::EHexahedral;
        else if(topologyMode == 5) eltype = MMeshType::EPrismatic;
        else DebugStop();
        TPZGenGrid3D *gen = new TPZGenGrid3D(x0, x1, nelDiv, eltype);

        TPZGeoMesh* gmesh = new TPZGeoMesh;
        gmesh = gen->BuildVolumetricElements(volID);
        gmesh = gen->BuildBoundaryElements(bcID,bcID,bcID,bcID,bcID,bcID);

        gmesh->SetDimension(3);
        return gmesh;
    }
    DebugStop(); // Dim should be 2 or 3
}

void Tools::DrawGeoMesh(ProblemConfig &config, PreConfig &preConfig) {

    std::stringstream ref;
    ref << "_ref-" << 1/preConfig.h <<" x " << 1/preConfig.h;
    std::string refinement =  ref.str();

    std::ofstream out(preConfig.plotfile + "/gmesh_"+ preConfig.topologyFileName + refinement + ".vtk");
    std::ofstream out2(preConfig.plotfile + "/gmesh"+ refinement + ".txt");

    TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, out);
    config.gmesh->Print(out2);
}

void Tools::DrawCompMesh(ProblemConfig &config, PreConfig &preConfig, TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *multiCmesh) {

    std::stringstream ref;
    ref << "_ref-" << 1/preConfig.h <<" x " << 1/preConfig.h;
    std::string refinement =  ref.str();

    std::ofstream out(preConfig.plotfile + "/cmesh" + refinement + ".txt");

    if (preConfig.mode == 0) cmesh->Print(out);
    else multiCmesh->Print(out);
}

void Tools::PrintErrors(std::ofstream& out, const ProblemConfig& config, const TPZVec<REAL>& error_vec) {

    std::stringstream ss;
    ss << "\nEstimator errors for Problem " << config.problemname;
    ss << "\n-------------------------------------------------- \n";
    ss << "Ndiv = " << config.ndivisions << ", AdaptivityStep = " << config.adaptivityStep
       << ", Order k = " << config.porder << ", Order n = " << config.hdivmais
       << ", K_R = " << config.Km << "\n";
    //ss << "DOF Total = " << config.fPostProcMesh.NEquations() << "\n";
    ss << "Global estimator = " << error_vec[3] << "\n";
    ss << "|ufem-urec| = " << error_vec[1] << "\n";
    ss << "Residual Error L2 = " << error_vec[4] << "\n";
    if (config.exact) {
        ss << "Global exact error = " << error_vec[2] << "\n";
        ss << "|uex-ufem| = " << error_vec[0] << "\n";
        REAL global_index = 1;
        if (!IsZero(error_vec[4] + error_vec[3]) && !IsZero(error_vec[2])) {
            global_index = sqrt(error_vec[4] + error_vec[3]) / sqrt(error_vec[2]);
        }
        ss << "Global Index = " << global_index;
    } else {
        ss << "[Unknown exact solution and errors]\n";
    }

    out << ss.str();
    std::cout << ss.str();
}
