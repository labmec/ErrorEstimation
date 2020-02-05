//
//  Tools.cpp
//  ErrorEstimation
//
//  Created by Denise De Siqueira on 28/03/19.
//

#include "Tools.h"
#include "pzgengrid.h"
#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "TPZGeoLinear.h"

#include <tuple>
#include <memory>



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
    
    if(problem.prefine){
        Prefinamento(cmesh, problem.ndivisions, problem.porder);
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
        int bctype = 0;
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(problem.exact.Exact());
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
                intel->SetPreferredOrder(problem.porder+problem.hdivmais);
            }
        }
 
    }
    
    if(problem.prefine){
        Prefinamento(cmesh, problem.ndivisions, problem.porder);
    }
    
    
    cmesh->InitializeBlock();
    return cmesh;
    
}

TPZMultiphysicsCompMesh *CreateHDivMesh(const ProblemConfig &problem) {
    
//    std::ofstream out("gmeshMulti.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(problem.gmesh, out);
   
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(problem.gmesh);
    TPZMaterial *mat = NULL;
    TPZFMatrix<REAL> K(3,3,0),invK(3,3,0);
    K.Identity();
    invK.Identity();
    
    if(problem.TensorNonConst && problem.gmesh->Dimension()==3){
        for(int i=0;i<3;i++){
            for(int j=0; j< 3;j++){
                if(i==j){
                    K(i,j) = 2.;
                    invK(i,j) = 3./4.;
                }
                else{
                    K(i,j) = 1.;
                    invK(i,j) = (-1.)/4.;
                }
            }
        }
        
    }

//    K.Print(std::cout);
//    invK.Print(std::cout);
    
    for (auto matid : problem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        mix->SetForcingFunction(problem.exact.ForcingFunction());
        mix->SetForcingFunctionExact(problem.exact.Exact());
        mix->SetPermeabilityTensor(K, invK);
        
        if (!mat) mat = mix;
        
        cmesh->InsertMaterialObject(mix);
    }
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype;
        if(matid == -1){
            bctype = 0;
        }
        else{
            bctype =1;
        }
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(problem.exact.Exact());
        cmesh->InsertMaterialObject(bc);
    }
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
   // std::set<int> matid;
//    matid.insert(1);
//    matid.insert(-1);
    TPZManVector<int> active(2,1);
    TPZManVector<TPZCompMesh *> meshvector(2,0);
    
    meshvector[0] = CreateFluxHDivMesh(problem);
    meshvector[1] = CreatePressureMesh(problem);
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    bool keepmatrix = false;
    bool keeponelagrangian = true;
    TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian, keepmatrix);
    
    return cmesh;
}

void CloneMeshVec(TPZVec<TPZCompMesh *> &meshvec, TPZVec<TPZCompMesh *> &meshvec_clone) {
    for (int i = 0; i < meshvec.size(); i++) {
        meshvec_clone[i] = meshvec[i]->Clone();
    }
}

/// Increase the approximation orders of the sides of the flux elements


void UniformRefinement(int nDiv, TPZGeoMesh *gmesh) {
    
    TPZManVector<TPZGeoEl*> children;
    for(int division = 0; division < nDiv; division++) {
        
        int64_t nels = gmesh->NElements();
        
        for(int64_t elem = 0; elem < nels; elem++) {
            
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            
            if(!gel || gel->HasSubElement()) continue;
            if(gel->Dimension() == 0) continue;
            gel->Divide(children);
        }
    }
}


TPZGeoMesh* CreateGeoMesh(int nel, TPZVec<int>& bcids) {
    
    TPZManVector<int> nx(2, nel);
    TPZManVector<REAL> x0(3, 0.), x1(3, 1.);
    x1[2] = 0.;
    TPZGenGrid gen(nx, x0, x1);
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


void MultiPhysicsCompel(const ProblemConfig &config){
    
    TPZManVector<TPZCompMesh *,2> MeshesHDiv(2);
    TPZMultiphysicsCompMesh * mixed_cmesh = CreateHDivMesh(config);
    MeshesHDiv = mixed_cmesh->MeshVector();
    
    TPZMultiphysicsCompMesh *mphysicCompMesh = new TPZMultiphysicsCompMesh(config.gmesh);
    std::ofstream outgeo("geometria.txt");
    mphysicCompMesh->Reference()->Print(outgeo);
    
    
    //Have to include the materials. Here we just did a copy of previous materials
    TPZCompMesh * cmesh =  dynamic_cast<TPZCompMesh *>(mphysicCompMesh);
    mixed_cmesh->CopyMaterials(*cmesh);
    
    TPZManVector<TPZCompMesh *,3> mp_meshes_vec(3);
    mp_meshes_vec[0] = mixed_cmesh;
    mp_meshes_vec[1] = MeshesHDiv[0];
    mp_meshes_vec[2] = MeshesHDiv[1];
    
    mphysicCompMesh->SetDimModel(2);
    TPZManVector<int,5>  active_approx_spaces(3,1);//teste usando todos os espaços
    mphysicCompMesh->BuildMultiphysicsSpace( active_approx_spaces, mp_meshes_vec);
    
    {
        std::ofstream out("mixed.txt");
        mphysicCompMesh->MeshVector()[0]->Print(out);
        
        std::ofstream out2("hdiv.txt");
        mphysicCompMesh->MeshVector()[1]->Print(out2);
        
        std::ofstream out3("L2.txt");
        mphysicCompMesh->MeshVector()[2]->Print(out3);
        
    }
    
    
    
}

void MultiPhysicsHybrid(const ProblemConfig &config){
    
    
    TPZManVector<TPZCompMesh*, 2> MeshesHDiv(2, 0);
    TPZMultiphysicsCompMesh *mixed_cmesh = CreateHDivMesh(config);//Hdiv x L2
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

    
    
    TPZMultiphysicsCompMesh *mphysicCompMesh = new TPZMultiphysicsCompMesh(config.gmesh);
    std::ofstream outgeo("geometria.txt");
    mphysicCompMesh->Reference()->Print(outgeo);
    
    
    //Have to include the materials. Here we just did a copy of previous materials
    TPZCompMesh * cmesh =  dynamic_cast<TPZCompMesh *>(mphysicCompMesh);
    mixed_cmesh->CopyMaterials(*cmesh);
    
    TPZManVector<TPZCompMesh *,3> mp_meshes_vec(3);
    mp_meshes_vec[0] = mixed_cmesh;
    mp_meshes_vec[1] = MeshesHDiv[0];
    mp_meshes_vec[2] = MeshesHDiv[1];
    
    mphysicCompMesh->SetDimModel(2);
    TPZManVector<int,5>  active_approx_spaces(3,1);//teste usando todos os espaços
    mphysicCompMesh->BuildMultiphysicsSpace( active_approx_spaces, mp_meshes_vec);
    
    {
        std::ofstream out("mixed.txt");
        mphysicCompMesh->MeshVector()[0]->Print(out);
        
        std::ofstream out2("hdiv.txt");
        mphysicCompMesh->MeshVector()[1]->Print(out2);
        
        std::ofstream out3("L2.txt");
        mphysicCompMesh->MeshVector()[2]->Print(out3);
        
    }
    
    
    
    
    
}

void RandomRefine(ProblemConfig &config,int numelrefine){
    
    int64_t nel = config.gmesh->NElements();
    if (numelrefine > nel/2) {
        numelrefine = 1;
    }
    int count = 0;
    while(count < numelrefine) {
        int64_t elindex = rand()%nel;
        TPZGeoEl *gel = config.gmesh->Element(elindex);
        if(gel && gel->Dimension() == config.gmesh->Dimension() && !gel->Father())
        {
            TPZStack<TPZGeoEl *> subels;
            gel->Divide(subels);
            count++;
        }
    }
    nel = config.gmesh->NElements();
    bool changed = true;
    while(changed)
    {
        changed = false;
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = config.gmesh->Element(el);
            if(gel && gel->Dimension() < config.gmesh->Dimension())
            {
                TPZGeoElSide gelside(gel,gel->NSides()-1);
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

void PrintSolAndDerivate(const ProblemConfig config){
    
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


void FunctionTest(){
        TLaplaceExample1 Denise;
    Denise.fExact = TLaplaceExample1::ESinSinDirNonHom;//TLaplaceExample1::ESinMark;//
        TPZVec<FADFADREAL> x(3);
        FADFADREAL x0 = (FADFADREAL) 0.001;
        FADFADREAL x1 = (FADFADREAL) 0.5;
        FADFADREAL x2 = (FADFADREAL) 0;
        x[0]= x0;
        x[1]= x1;
        x[2]= x2;
        TPZVec<FADFADREAL> disp(1);
        Denise.uxy(x, disp);
    std::cout<< "Pto x[0] "<<x[0]<<std::endl;
    std::cout<< "Pto x[1] "<<x[1]<<std::endl;
    std::cout<< "Pto x[2] "<<x[2]<<std::endl;
    
    std::cout<<"valor de ur0 "<<disp[0]<<std::endl;
    
        TPZVec<REAL> x_r(3);
        x_r[0] = x[0].val().val();
        x_r[1] = x[1].val().val();
        x_r[2] = x[2].val().val();
        TPZManVector<REAL,3> grad(3);
        Denise.graduxy(x_r, grad);
    
        REAL force;
        Denise.DivSigma(x_r, force);
    
}


void Prefinamento(TPZCompMesh * cmesh, int ndiv, int porder){
    if(ndiv<1) return;
    int nel = cmesh->NElements();
    for(int iel = 0; iel < nel; iel++){
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        if(!cel) continue;
        
        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
        if(!sp) continue;
        int level = sp->Reference()->Level();
        TPZGeoEl * gel = sp->Reference();
        if((gel->Dimension()==2) && (iel % 2==0)){
            int ordem= 0;
            ordem=porder + (ndiv-1 ) + (level);
            std::cout<<"level "<< level<<" ordem "<<ordem<<std::endl;
            sp->PRefine(ordem);
        }
    }
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    
        std::stringstream sout;
        sout<<"malha computacional apos pRefinamento\n";
        cmesh->Print(sout);

}

void SolveHybridProblem(TPZCompMesh *Hybridmesh,int InterfaceMatId,const ProblemConfig &problem,bool PostProcessingFEM ){
    

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

    if(PostProcessingFEM){

    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("ExactPressure");
    vecnames.Push("Flux");

    std::stringstream sout;
    sout << problem.dir_name << "/" <<  "OriginalHybrid_Order_"<<problem.porder<<"Nref_"<<problem.ndivisions<<".vtk";
    an.DefineGraphMesh(2, scalnames, vecnames, sout.str());
    int resolution = 2;
    an.PostProcess(resolution,Hybridmesh->Dimension());
    }

    
}
void PlotLagrangeMultiplier(TPZCompMesh *cmesh,const ProblemConfig &problem){
    
    TPZAnalysis an(cmesh,false);
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("State");
    
    int dim = cmesh->Reference()->Dimension()-1;
    std::string plotname;
    {
        std::stringstream out;
        out << problem.dir_name << "/" << "OriginalLagrangeMultiplier" << ".vtk";
        plotname = out.str();
    }
    an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
    an.PostProcess(2,dim);
    
}

void SolveMixedProblem(TPZCompMesh *cmesh_HDiv,const ProblemConfig &config)
{
    #ifdef PZDEBUG
            {
                std::ofstream out("gmeshSolve.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, out);
                
            }
    #endif
    
    
    TPZAnalysis an(cmesh_HDiv,false);
    
    
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
    scalnames.Push("ExactPressure");
    vecnames.Push("Flux");
    vecnames.Push("ExactFlux");
    
    int dim = config.gmesh->Dimension();
    
    std::stringstream sout;
    
    sout << config.dir_name << "/"  "OriginalMixed_Order_"<<config.problemname<<"Order"<< config.porder<<"Nref_"<<config.ndivisions<<".vtk";
    
    an.DefineGraphMesh(dim, scalnames, vecnames, sout.str());
    int resolution=2;
    an.PostProcess(resolution,dim);
    
    
    if(config.exact.Exact())
    {
        TPZManVector<REAL> errors(4,0.);
        an.SetThreadsForError(0);
        an.SetExact(config.exact.ExactSolution());
        an.PostProcessError(errors,false);
        
        //Erro
        
        ofstream myfile;
        myfile.open("MixedError.txt", ios::app);
        myfile << "\n\n Error for Mixed formulation " ;
        myfile << "\n-------------------------------------------------- \n";
        myfile << "Ndiv = " << config.ndivisions << " Order k = " << config.porder <<"\n";
        myfile << "Energy norm = " << errors[0] << "\n";//norma energia
        myfile << "error norm L2 = " << errors[1] << "\n";//norma L2
        myfile << "Semi norm H1 = " << errors[2] << "\n";//norma L2
        myfile.close();
        
    }
}


TPZGeoMesh *ReadGeometricMesh(struct ProblemConfig &config, bool IsgmeshReader){
    
   
     TPZGeoMesh *gmesh = nullptr;
     int dim= config.dimension;
    
    
    
    if(IsgmeshReader){


        std::string meshfilename = "../LCircle.msh";

        if(dim==3)
        {
            meshfilename = "../Cube.msh";
        }
        TPZGmshReader gmsh;
        //  gmsh.GetDimNamePhysical().resize(4);
        //  gmsh.GetDimPhysicalTagName().resize(4);
        if(dim==2)
        {
            gmsh.GetDimNamePhysical()[1]["dirichlet"] =2;
            gmsh.GetDimNamePhysical()[2]["domain"] = 1;
        }
        else
        {
            gmsh.GetDimNamePhysical()[2]["dirichlet"] =2;
            gmsh.GetDimNamePhysical()[3]["domain"] = 1;
        }
        config.materialids.insert(1);
        config.bcmaterialids.insert(2);


        gmsh.SetFormatVersion("4.1");
        gmesh = gmsh.GeometricGmshMesh(meshfilename);
        gmsh.PrintPartitionSummary(std::cout);
        gmesh->SetDimension(dim);
        config.gmesh = gmesh;

    }

    else{

        TPZManVector<int,4> bcids(4,-1);
        gmesh = CreateGeoMesh(2, bcids);
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        config.gmesh = gmesh;
        gmesh->SetDimension(dim);



    }
    
    return gmesh;
    
    
}

 TPZMultiphysicsCompMesh * HybridSolveProblem(TPZMultiphysicsCompMesh *cmesh_HDiv, struct ProblemConfig &config){
        
     TPZManVector<TPZCompMesh *,2> hybridmeshvec;
    hybridmeshvec = cmesh_HDiv->MeshVector();
    
    //cria malha hibrida
    std::cout<<"Initializing the hybridization procedure"<<std::endl;
    
    TPZHybridizeHDiv hybrid;
    auto HybridMesh = hybrid.Hybridize(cmesh_HDiv);
    HybridMesh->CleanUpUnconnectedNodes();//enumerar adequadamente os connects
    HybridMesh->AdjustBoundaryElements();
    delete cmesh_HDiv;
    delete hybridmeshvec[0];
    delete hybridmeshvec[1];
    
    
    std::cout<<"---Original PerifericalMaterialId --- "<<std::endl;
    std::cout <<" LagrangeInterface = "<<hybrid.fLagrangeInterface<<std::endl;
    std::cout <<" HDivWrapMatid = "<<hybrid.fHDivWrapMatid<<std::endl;
    std::cout <<" InterfaceMatid = "<<hybrid.fInterfaceMatid<<std::endl;
    
    
    
       #ifdef PZDEBUG
        {
    
            std::ofstream out2("OriginalFluxMesh.txt");
            HybridMesh->MeshVector()[0]->Print(out2);
    
            std::ofstream out3("OriginalPotentialMesh.txt");
            HybridMesh->MeshVector()[1]->Print(out3);
    
        }
    #endif
    
    
        SolveHybridProblem(HybridMesh,hybrid.fInterfaceMatid,config,false);
    
    #ifdef PZDEBUG
        {
            std::ofstream out("OriginalHybridMesh.txt");
            (HybridMesh)->Print(out);
        }
    #endif
    
       PlotLagrangeMultiplier(HybridMesh->MeshVector()[1],config);
     
     cmesh_HDiv = HybridMesh;
    
     //return HybridMesh;
     return cmesh_HDiv;
}

/// Divide lower dimensional elements
void DivideLowerDimensionalElements(TPZGeoMesh *gmesh)
{
    bool haschanged = true;
    int dim = gmesh->Dimension();
    while(haschanged)
    {
        haschanged = false;
        int64_t nel = gmesh->NElements();
        TPZStack<TPZGeoEl *> geldivide;
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel || gel->Dimension() == dim || gel->Dimension() == 0)
            {
                continue;
            }
            if(gel->HasSubElement()) continue;
            int nsides = gel->NSides();
            int ncorner = gel->NCornerNodes();
            for (int side = ncorner; side < nsides; side++) {
                TPZGeoElSide gelside(gel,side);
                TPZGeoElSide neighbour(gelside.Neighbour());
                bool found = false;
                while (neighbour != gelside) {
                    if(neighbour.HasSubElement())
                    {
                        geldivide.Push(gel);
                        found = true;
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                }
                if(found) break;
            }
        }
        if(geldivide.size())
        {
            haschanged = true;
            for (int64_t i = 0; i<geldivide.size(); i++) {
                TPZManVector<TPZGeoEl *> sub;
                geldivide[i]->Divide(sub);
            }
        }
    }
}


TPZCompMesh *CMeshH1( ProblemConfig problem){
    
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial *mat = 0;
    
    
    for (auto matid : problem.materialids) {
        TPZMatPoisson3d *mix = new TPZMatPoisson3d(matid, cmesh->Dimension());
        mix->SetForcingFunctionExact(problem.exact.Exact());
        mix->SetForcingFunction(problem.exact.ForcingFunction());
        
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
        
    }
    
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        val2.Zero();
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(problem.exact.Exact());
        
        cmesh->InsertMaterialObject(bc);
    }
    
    cmesh->SetDefaultOrder(problem.porder);//ordem
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    
    cmesh->AutoBuild();
    
    
    return cmesh;
}
void hAdaptivity(TPZCompMesh *postProcessMesh, TPZGeoMesh *gmeshToRefine) {
    
    // Column of the flux error estimate on the element solution matrix
    const int fluxErrorEstimateCol = 3;
    
    int64_t nelem = postProcessMesh->ElementSolution().Rows();
    
   // postProcessMesh->ElementSolution().Print("ElSolutionForAdaptivity",std::cout);
    
    // Iterates through element errors to get the maximum value
    REAL maxError = 0.;
    for (int64_t iel = 0; iel < nelem; iel++) {
        TPZCompEl *cel = postProcessMesh->ElementVec()[iel];
        if (!cel) continue;
        if (cel->Dimension() != postProcessMesh->Dimension()) continue;
        REAL elementError = postProcessMesh->ElementSolution()(iel, fluxErrorEstimateCol);
    
        
        if (elementError > maxError) {
            maxError = elementError;
        }
    }
    
    // Refines elements which error are bigger than 30% of the maximum error
    REAL threshold = 0.3 * maxError;
    
    for (int64_t iel = 0; iel < nelem; iel++) {
        TPZCompEl *cel = postProcessMesh->ElementVec()[iel];
        if (!cel) continue;
        if (cel->Dimension() != postProcessMesh->Dimension()) continue;
        REAL elementError = postProcessMesh->ElementSolution()(iel, fluxErrorEstimateCol);
        if (elementError > threshold) {
            TPZGeoEl *gel = cel->Reference();
            int iel = gel->Id();
            
            TPZVec<TPZGeoEl *> sons;
            TPZGeoEl *gelToRefine = gmeshToRefine->ElementVec()[iel];
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


TPZGeoMesh *CreateTrapezoidalMesh(int nelx, int nely, REAL Lx, REAL Ly, TPZVec<int> &bcids){

    TPZGeoMesh * gmesh = new TPZGeoMesh;
    
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.);
    TPZManVector<int,3> nx(2);
    nx[0] = nelx;
    nx[1] = nely;
    x1[0] = Lx;
    x1[1] = Ly;
    TPZGenGrid gengrid(nx,x0,x1);
    
    gengrid.SetDistortion(0.25);
    //        gengrid.SetZigZagPattern();

    gengrid.Read(gmesh);
    gengrid.SetBC(gmesh, 4, bcids[0]);
    gengrid.SetBC(gmesh, 5, bcids[1]);
    gengrid.SetBC(gmesh, 6, bcids[2]);
    gengrid.SetBC(gmesh, 7, bcids[3]);
//    x1[0] = Lx;
//    x1[1] = 0.;
//    gengrid.SetBC(gmesh, x0, x1, BC0);
//    x0 = x1;
//    x1[0] = Lx;
//    x1[1] = Ly;
//    gengrid.SetBC(gmesh, x0, x1, BC1);
//    x0 = x1;
//    x1[0] = 0.;
//    x1[1] = Ly;
//    gengrid.SetBC(gmesh, x0, x1, BC2);
//    x0 = x1;
//    x1[0] = 0.;
//    x1[1] = 0.;
//    gengrid.SetBC(gmesh, x0, x1, BC3);
    
    return gmesh;
}
