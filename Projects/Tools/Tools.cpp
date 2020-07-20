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

#include "pzcondensedcompel.h"
#include "pzelementgroup.h"

TPZCompMesh* CreatePressureMesh(const ProblemConfig& problem) {
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

TPZCompMesh* CreateFluxHDivMesh(const ProblemConfig& problem) {
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

TPZMultiphysicsCompMesh* CreateHDivMesh(const ProblemConfig& problem) {
    
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
        mix->SetForcingFunctionExact(problem.exact.operator*().Exact());
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
    
    return cmesh;
}

void CloneMeshVec(TPZVec<TPZCompMesh*>& meshvec, TPZVec<TPZCompMesh*>& meshvec_clone) {
    for (int i = 0; i < meshvec.size(); i++) {
        meshvec_clone[i] = meshvec[i]->Clone();
    }
}

/// Increase the approximation orders of the sides of the flux elements


void UniformRefinement(int nDiv, TPZGeoMesh* gmesh) {
    
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


TPZGeoMesh* CreateGeoMesh(int nel, TPZVec<int>& bcids) {
    
    TPZManVector<int> nx(2, nel);
    TPZManVector<REAL> x0(3, 0.), x1(3, 1.);
    x1[2] = 0.;
    TPZGenGrid2D gen(nx, x0, x1, 1, 0);
    
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


void MultiPhysicsCompel(const ProblemConfig& config) {
    
    TPZManVector<TPZCompMesh*, 2> MeshesHDiv(2);
    TPZMultiphysicsCompMesh* mixed_cmesh = CreateHDivMesh(config);
    MeshesHDiv = mixed_cmesh->MeshVector();
    
    TPZMultiphysicsCompMesh* mphysicCompMesh = new TPZMultiphysicsCompMesh(config.gmesh);
    std::ofstream outgeo("geometria.txt");
    mphysicCompMesh->Reference()->Print(outgeo);
    
    
    //Have to include the materials. Here we just did a copy of previous materials
    TPZCompMesh* cmesh = dynamic_cast<TPZCompMesh*>(mphysicCompMesh);
    mixed_cmesh->CopyMaterials(*cmesh);
    
    TPZManVector<TPZCompMesh*, 3> mp_meshes_vec(3);
    mp_meshes_vec[0] = mixed_cmesh;
    mp_meshes_vec[1] = MeshesHDiv[0];
    mp_meshes_vec[2] = MeshesHDiv[1];
    
    mphysicCompMesh->SetDimModel(2);
    TPZManVector<int, 5> active_approx_spaces(3, 1);//teste usando todos os espaços
    mphysicCompMesh->BuildMultiphysicsSpace(active_approx_spaces, mp_meshes_vec);
    
    {
        std::ofstream out("mixed.txt");
        mphysicCompMesh->MeshVector()[0]->Print(out);
        
        std::ofstream out2("hdiv.txt");
        mphysicCompMesh->MeshVector()[1]->Print(out2);
        
        std::ofstream out3("L2.txt");
        mphysicCompMesh->MeshVector()[2]->Print(out3);
        
    }
    
    
}

void MultiPhysicsHybrid(const ProblemConfig& config) {
    
    
    TPZManVector<TPZCompMesh*, 2> MeshesHDiv(2, 0);
    TPZMultiphysicsCompMesh* mixed_cmesh = CreateHDivMesh(config);//Hdiv x L2
    MeshesHDiv = mixed_cmesh->MeshVector();
    mixed_cmesh->InitializeBlock();
    
    //cria malha hibrida
    TPZHybridizeHDiv hybrid;
    auto HybridMesh = hybrid.Hybridize(mixed_cmesh);
    (HybridMesh)->CleanUpUnconnectedNodes();//enumerar adequadamente os connects
    delete mixed_cmesh;
    delete MeshesHDiv[0];
    delete MeshesHDiv[1];
    
    mixed_cmesh = (HybridMesh);//malha hribrida
    MeshesHDiv[0] = (HybridMesh)->MeshVector()[0];//malha Hdiv
    MeshesHDiv[1] = (HybridMesh)->MeshVector()[1];//malha L2
    
    //////
    
    
    
    TPZMultiphysicsCompMesh* mphysicCompMesh = new TPZMultiphysicsCompMesh(config.gmesh);
    std::ofstream outgeo("geometria.txt");
    mphysicCompMesh->Reference()->Print(outgeo);
    
    
    //Have to include the materials. Here we just did a copy of previous materials
    TPZCompMesh* cmesh = dynamic_cast<TPZCompMesh*>(mphysicCompMesh);
    mixed_cmesh->CopyMaterials(*cmesh);
    
    TPZManVector<TPZCompMesh*, 3> mp_meshes_vec(3);
    mp_meshes_vec[0] = mixed_cmesh;
    mp_meshes_vec[1] = MeshesHDiv[0];
    mp_meshes_vec[2] = MeshesHDiv[1];
    
    mphysicCompMesh->SetDimModel(2);
    TPZManVector<int, 5> active_approx_spaces(3, 1);//teste usando todos os espaços
    mphysicCompMesh->BuildMultiphysicsSpace(active_approx_spaces, mp_meshes_vec);
    
    {
        std::ofstream out("mixed.txt");
        mphysicCompMesh->MeshVector()[0]->Print(out);
        
        std::ofstream out2("hdiv.txt");
        mphysicCompMesh->MeshVector()[1]->Print(out2);
        
        std::ofstream out3("L2.txt");
        mphysicCompMesh->MeshVector()[2]->Print(out3);
        
    }
    
    
}

void RandomRefine(ProblemConfig& config, int numelrefine, int depth) {
    
    int64_t nel = config.gmesh->NElements();
    if (numelrefine > nel / 2) {
        numelrefine = nel/2;
    }
    for (int id=0; id<depth; id++)
    {
        nel = config.gmesh->NElements();
        int count = 0;
        while (count < numelrefine) {
            int64_t elindex = rand() % nel;
            TPZGeoEl* gel = config.gmesh->Element(elindex);
            int level = gel->Level();
            if (gel && gel->Dimension() == config.gmesh->Dimension() && gel->Level() == id) {
                TPZStack<TPZGeoEl*> subels;
                gel->Divide(subels);
                count++;
            }
        }
    }
    nel = config.gmesh->NElements();
    bool changed = true;
    while (changed) {
        changed = false;
        nel = config.gmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZGeoEl* gel = config.gmesh->Element(el);
            if (gel && gel->Dimension() < config.gmesh->Dimension()) {
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


void Print(const FADREAL& a, std::ostream& out) {
    out << " val " << a.val() << std::endl;
    for (int i = 0; i < a.dx().size(); i++) {
        out << a.d(i) << " ";
    }
    out << std::endl;
}

void Print(const FADFADREAL& a, std::ostream& out) {
    out << "Value ";
    Print(a.val(), out);
    out << "Derivatives\n";
    for (int i = 0; i < a.dx().size(); i++) {
        Print(a.d(i), out);
    }
    out << "End\n";
    
}

void PrintSolAndDerivate(const ProblemConfig config) {
    
    TPZManVector<REAL, 3> x(3, 0.25);
    
    TPZManVector<Fad<REAL>, 3> xfad(x.size()), graduxy(x.size());
    TPZManVector<FADFADREAL, 3> xfadfad(x.size()), uxyfadfad(1);
    for (int i = 0; i < 3; i++) {
        xfad[i] = Fad<REAL>(3, i, x[i]);
        xfadfad[i] = FADFADREAL(3, i, xfad[i]);
        for (int j = 0; j < 3; j++) {
            xfadfad[i].fastAccessDx(j) = Fad<REAL>(3, xfadfad[i].val().dx(j));
        }
    }
    std::cout << "xfadfad = \n";
    for (int i = 0; i < 3; i++) {
        Print(xfadfad[i], std::cout);
    }
    std::cout << std::endl;
    config.exact.operator*().graduxy(xfad, graduxy);
    config.exact.operator*().uxy(xfadfad, uxyfadfad);
    for (int i = 0; i < 3; i++) {
        std::cout << "xfad = ";
        Print(xfad[i], std::cout);
        std::cout << std::endl;
    }
    std::cout << "graduxy = \n";
    for (int i = 0; i < 3; i++) {
        Print(graduxy[i], std::cout);
    }
    std::cout << std::endl;
    std::cout << "uxyfadfad = \n";
    for (int i = 0; i < uxyfadfad.size(); i++) {
        Print(uxyfadfad[i], std::cout);
    }
    REAL laplace = uxyfadfad[0].dx(0).dx(0) + uxyfadfad[0].dx(1).dx(1) + uxyfadfad[0].dx(2).dx(2);
    std::cout << "Laplacian " << laplace << std::endl;
}


void FunctionTest() {
    TLaplaceExample1 Denise;
    Denise.fExact = TLaplaceExample1::ESinMark;//ESinSinDirNonHom;//TLaplaceExample1::
    TPZVec<FADFADREAL> x(3);
    FADFADREAL x0 = (FADFADREAL) 0.013;
    FADFADREAL x1 = (FADFADREAL) 0.25;
    FADFADREAL x2 = (FADFADREAL) 0;
    x[0] = x0;
    x[1] = x1;
    x[2] = x2;
    TPZVec<FADFADREAL> disp(1);
    Denise.uxy(x, disp);
    std::cout << "Pto x[0] " << x[0] << std::endl;
    std::cout << "Pto x[1] " << x[1] << std::endl;
    std::cout << "Pto x[2] " << x[2] << std::endl;
    
    std::cout << "valor de ur0 " << disp[0] << std::endl;
    
    TPZVec<REAL> x_r(3);
    x_r[0] = x[0].val().val();
    x_r[1] = x[1].val().val();
    x_r[2] = x[2].val().val();
    TPZManVector<REAL, 3> grad(3);
    Denise.graduxy(x_r, grad);
    
    std::cout << "valor de grad " << grad[0] <<", "<<grad[1]<<","<< grad[2]<< std::endl;
    
    
    
    REAL force;
    Denise.DivSigma(x_r, force);
    
    std::cout << "valor de div " << force << std::endl;
    
}


void Prefinamento(TPZCompMesh* cmesh, int ndiv, int porder) {
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

void
SolveHybridProblem(TPZCompMesh* Hybridmesh, std::pair<int,int> InterfaceMatId, const ProblemConfig& problem, bool PostProcessingFEM) {
    
    
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
    
    
    TPZStepSolver<STATE>* direct = new TPZStepSolver<STATE>;
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
            
            //
            
            
            
            //
            
            ofstream myfile;
            myfile.open("ErrorBCFemProblem.txt", ios::app);
            
            
            
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
    
    
}

void ComputeError(TPZCompMesh *Hybridmesh, std::ofstream &out,const ProblemConfig &config)
{
    long nel = Hybridmesh->NElements();
    int dim = Hybridmesh->Dimension();
    TPZManVector<STATE,10> globerrors(10,0.);
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = Hybridmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        
        TPZManVector<REAL,10> elerror(5,0.);
        
        int matId = gel->MaterialId();
        

        TPZMaterial *mat = Hybridmesh->FindMaterial(matId);
        
        if(matId != 1 || matId != 4) continue;
        
        
        TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
        

        
        cel->EvaluateError(config.exact->ExactSolution(), elerror, NULL);

        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    // out << "Errors associated with HDiv space\n";
    out << "L2 Norm for flux = "    << sqrt(globerrors[1]) << std::endl;
    out << "L2 Norm for divergence = "    << sqrt(globerrors[2]) << std::endl;
    out << "Hdiv Norm for flux = "    << sqrt(globerrors[3])  <<std::endl;
    
}


void PlotLagrangeMultiplier(TPZCompMesh* cmesh, const ProblemConfig& problem) {
    
    TPZAnalysis an(cmesh, false);
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("State");
    
    int dim = cmesh->Reference()->Dimension() - 1;
    std::string plotname;
    {
        std::stringstream out;
        out << problem.dir_name << "/" << "OriginalLagrangeMultiplier" << ".vtk";
        plotname = out.str();
    }
    an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
    an.PostProcess(2, dim);
    
}

void SolveMixedProblem(TPZCompMesh* cmesh_HDiv, const ProblemConfig& config) {
#ifdef PZDEBUG
    {
        std::ofstream out("gmeshSolve.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, out);
        
    }
#endif
    
    
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
        ofstream myfile;
        /*Error on MixedPoisson
           [0] L2 for pressure
           [1] L2 for flux
           [2] L2 for div(flux)
           [3] Grad pressure (Semi H1)
           [4] Hdiv norm
           */

          // Erro
          myfile.open("ErrorMixed.txt", ios::app);
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


TPZGeoMesh* ReadGeometricMesh(struct ProblemConfig& config, bool IsgmeshReader) {
    
    
    TPZGeoMesh* gmesh = nullptr;
    int dim = config.dimension;
    
    
    if (IsgmeshReader) {
        
        
        std::string meshfilename = "../LCircle.msh";
        
        if (dim == 3) {
            meshfilename = "../Cube.msh";
        }
        TPZGmshReader gmsh;
        //  gmsh.GetDimNamePhysical().resize(4);
        //  gmsh.GetDimPhysicalTagName().resize(4);
        if (dim == 2) {
            gmsh.GetDimNamePhysical()[1]["dirichlet"] = 2;
            gmsh.GetDimNamePhysical()[2]["domain"] = 1;
        } else {
            gmsh.GetDimNamePhysical()[2]["dirichlet"] = 2;
            gmsh.GetDimNamePhysical()[3]["domain"] = 1;
        }
        config.materialids.insert(1);
        config.bcmaterialids.insert(2);
        
        
        gmsh.SetFormatVersion("4.1");
        gmesh = gmsh.GeometricGmshMesh(meshfilename);
        gmsh.PrintPartitionSummary(std::cout);
        gmesh->SetDimension(dim);
        config.gmesh = gmesh;
        
    } else {
        
        TPZManVector<int, 4> bcids(4, -1);
        gmesh = CreateGeoMesh(2, bcids);
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        config.gmesh = gmesh;
        gmesh->SetDimension(dim);
        
        
    }
    
    return gmesh;
    
    
}

TPZMultiphysicsCompMesh* HybridSolveProblem(TPZMultiphysicsCompMesh* cmesh_HDiv, struct ProblemConfig& config) {
    
    TPZManVector<TPZCompMesh*, 2> hybridmeshvec;
    hybridmeshvec = cmesh_HDiv->MeshVector();
    
    //cria malha hibrida
    std::cout << "Initializing the hybridization procedure" << std::endl;
    
    TPZHybridizeHDiv hybrid;
    auto HybridMesh = hybrid.Hybridize(cmesh_HDiv);
    HybridMesh->CleanUpUnconnectedNodes();//enumerar adequadamente os connects
    HybridMesh->AdjustBoundaryElements();
    delete cmesh_HDiv;
    delete hybridmeshvec[0];
    delete hybridmeshvec[1];
    
    
    std::cout << "---Original PerifericalMaterialId --- " << std::endl;
    std::cout << " LagrangeInterface = " << hybrid.fLagrangeInterface << std::endl;
    std::cout << " HDivWrapMatid = " << hybrid.fHDivWrapMatid << std::endl;
    std::cout << " InterfaceMatid = " << hybrid.fInterfaceMatid << std::endl;


#ifdef PZDEBUG
    {
        
        std::ofstream out2("OriginalFluxMesh.txt");
        HybridMesh->MeshVector()[0]->Print(out2);
        
        std::ofstream out3("OriginalPotentialMesh.txt");
        HybridMesh->MeshVector()[1]->Print(out3);
        
    }
#endif
    
    
    SolveHybridProblem(HybridMesh, hybrid.fInterfaceMatid, config, false);

#ifdef PZDEBUG
    {
        std::ofstream out("OriginalHybridMesh.txt");
        (HybridMesh)->Print(out);
    }
#endif
    
   // PlotLagrangeMultiplier(HybridMesh->MeshVector()[1], config);
    
    cmesh_HDiv = HybridMesh;
    
    //return HybridMesh;
    return cmesh_HDiv;
}

/// Divide lower dimensional elements
void DivideLowerDimensionalElements(TPZGeoMesh* gmesh) {
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


TPZCompMesh* CMeshH1(ProblemConfig problem) {
    
    TPZCompMesh* cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial* mat = 0;
    
    
    for (auto matid : problem.materialids) {
        TPZMatPoisson3d *mix = new TPZMatPoisson3d(matid, cmesh->Dimension());
        mix->SetForcingFunctionExact(problem.exact.operator*().Exact());
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

void hAdaptivity(TPZCompMesh* postProcessMesh, TPZGeoMesh* gmeshToRefine, ProblemConfig& config) {
    
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


TPZGeoMesh* CreateLCircleGeoMesh() {
    
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


TPZGeoMesh* CreateTrapezoidalMesh(int nelx, int nely, REAL Lx, REAL Ly, TPZVec<int>& bcids) {
    
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    
    TPZManVector<REAL, 3> x0(3, 0.), x1(3, 0.);
    TPZManVector<int, 3> nx(2);
    nx[0] = nelx;
    nx[1] = nely;
    x1[0] = Lx;
    x1[1] = Ly;
    
    TPZGenGrid2D gengrid(nx, x0, x1, 1, 0);
    gengrid.SetDistortion(0.25);
    
    gengrid.Read(gmesh);
    gengrid.SetBC(gmesh, 4, bcids[0]);
    gengrid.SetBC(gmesh, 5, bcids[1]);
    gengrid.SetBC(gmesh, 6, bcids[2]);
    gengrid.SetBC(gmesh, 7, bcids[3]);
    
    return gmesh;
}

TPZGeoMesh* CreateLShapeMesh(TPZVec<int>& bcids) {
    
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

TPZGeoMesh* CreateQuadLShapeMesh(TPZVec<int>& bcids) {
    
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

TPZGeoMesh* CreateSingleTriangleMesh(TPZVec<int>& bcids) {
    
    TPZGeoMesh* gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    int matID = 1;
    
    // Creates matrix with node coordinates
    const int NodeNumber = 3;
    REAL coordinates[NodeNumber][3] = {
            {0.,  0., 0.},
            {1.,  0., 0.},
            {1.,  1., 0.}
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
    TPZManVector<int64_t> nodeIDs(3);
    nodeIDs[0] = 0;
    nodeIDs[1] = 1;
    nodeIDs[2] = 2;
    new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodeIDs, matID, *gmesh);
    
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


TPZGeoMesh* CreateSingleQuadMesh(TPZVec<int>& bcids) {
    
    TPZGeoMesh* gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    int matID = 1;
    
    // Creates matrix with node coordinates
    const int NodeNumber = 4;
    REAL coordinates[NodeNumber][4] = {
            {0.,  0., 0.},
            {1.,  0., 0.},
            {1.,  1., 0.},
            {0.,  1., 0.}
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
    nodeIDs[0] = 0;
    nodeIDs[1] = 1;
    nodeIDs[2] = 2;
    nodeIDs[3] = 3;
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeIDs, matID, *gmesh);
    
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

TPZGeoMesh* CreateQuadMeshRefTriang(TPZVec<int>& bcids) {
    
    TPZGeoMesh* gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    int matID = 1;
    
    // Creates matrix with node coordinates
    const int NodeNumber = 4;
    REAL coordinates[NodeNumber][3] = {
        {0.,  0., 0.},
        {1.,  0., 0.},
        {1.,  1., 0.},
        {0.,  1., 0.}
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
    TPZManVector<int64_t> nodeIDs(3);
    
    nodeIDs[0] = 0;
    nodeIDs[1] = 1;
    nodeIDs[2] = 3;
    new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodeIDs, matID, *gmesh);
    nodeIDs[0] = 2;
    nodeIDs[1] = 3;
    nodeIDs[2] = 1;
    
    new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodeIDs, matID, *gmesh);

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

void VectorEnergyNorm(TPZCompMesh *hdivmesh, std::ostream &out,  const ProblemConfig& problem)
{
    
    
    //
     TPZGeoMesh *gmesh = hdivmesh->Reference();
     gmesh->ResetReference();
     int dim = gmesh->Dimension();

     int64_t nel = hdivmesh->NElements();
     
     
     // loop over the elements
     for (int64_t el = 0; el < nel; el++) {
         TPZCompEl *cel = hdivmesh->Element(el);
         if (!cel) continue;
         TPZGeoEl *gel = cel->Reference();
         if (!gel) continue;
         int matId = gel->MaterialId();
         std::cout<<"matId "<<matId<<"\n";
         
         
         int nsides = gel->NSides();
         for (int side = 0; side < nsides; side++) {
             TPZGeoElSide gelside(gel, nsides - 1);
             TPZStack<TPZCompElSide> equal;
             int onlyinterpolated = 1;
             int removeduplicated = 0;
             gelside.EqualLevelCompElementList(equal, onlyinterpolated, removeduplicated);
             int nequal = equal.size();
             if(nequal==0) continue;

             for (int ieq = 0; ieq < nequal; ieq++) {
                 TPZCompEl *celneigh = equal[ieq].Element();
                 TPZInterpolatedElement *intelneigh = dynamic_cast<TPZInterpolatedElement *>(celneigh);
                 TPZVec<REAL> errors(5,0.);

                intelneigh->EvaluateError(problem.exact->ExactSolution(),errors,false);
   

             }
         }
     }
    
    
    //
    
    
    
//        nkaux[iel] = residuo2;
//        out << "\nErrors associated with flux on element Ek\n";
//        out << "L2 Norm flux = "    << nkaux[iel] << endl;
        
      

}
