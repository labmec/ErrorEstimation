//
//  Tools.cpp
//  ErrorEstimation
//
//  Created by Denise De Siqueira on 28/03/19.
//

#include "Tools.h"
#include "pzgengrid.h"



TPZCompMesh *CreatePressureMesh(const ProblemConfig &problem) {
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial *mat = 0;
    for (auto matid : problem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    cmesh->SetDefaultOrder(problem.porder+problem.hdivmais);//ordem + hdivmais
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
        TPZBndCond *bc = mat->CreateBC(mat, matid, 0, val1, val2);
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
        
        if(problem.prefine){
            Prefinamento(cmesh, problem.ndivisions, problem.porder);
        }
        
    }
    cmesh->InitializeBlock();
    return cmesh;
    
}

TPZMultiphysicsCompMesh *CreateHDivMesh(const ProblemConfig &problem) {
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(problem.gmesh);
    TPZMaterial *mat = NULL;
    for (auto matid : problem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        mix->SetForcingFunction(problem.exact.ForcingFunction());
        
        mix->SetForcingFunctionExact(problem.exact.Exact());
        
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
    std::set<int> matid;
    matid.insert(1);
    matid.insert(-1);
    TPZManVector<int> active(2,1);
    TPZManVector<TPZCompMesh *> meshvector(2,0);
    
    meshvector[0] = CreateFluxHDivMesh(problem);
    meshvector[1] = CreatePressureMesh(problem);
    cmesh->BuildMultiphysicsSpace(active, meshvector);
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
            //nao poderia denomiar de left uma vez que esta pegando os elementos a diretira e esqueda do lado e
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
            celstack[0].Element()->Solution(pt1, 0, sol1);//rigth pressure
            celstack[1].Element()->Solution(pt2, 0, sol2);//lefth pressure
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


TPZGeoMesh *CreateGeoMesh(int nel) {
    
    TPZManVector<int> nx(2,nel);
    TPZManVector<REAL> x0(3,0.),x1(3,1.);
    x1[2] = 0.;
    TPZGenGrid gen(nx,x0,x1);
    int matID=1;

    TPZGeoMesh* gmesh = new TPZGeoMesh;
    gen.Read(gmesh);
    gen.SetBC(gmesh, 4, -1);
    gen.SetBC(gmesh, 5, -1);
    gen.SetBC(gmesh, 6, -1);
    gen.SetBC(gmesh, 7, -1);

    
    

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

void RandonRefine(ProblemConfig &config,int numelrefine){
    
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
    Denise.fExact = TLaplaceExample1::ESinMark;//ESinSinDirNonHom;
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


void ExataOmega1(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
    flux.Resize(3, 1);
    
    REAL x,y,z;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    solp[0] = x*x-y*y;
    flux(0,0) = 2.*x;
    flux(1,0) = -2*y;
    flux(2,0) = 0.;
std::cout<<"funcao em omega1 "<<solp[0]<<std::endl;
}

void ExataOmega2(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
    flux.Resize(3, 1);
    
    REAL x,y,z,alpha;
    
    alpha=1.;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    solp[0] = (x*x-y*y)/(alpha*alpha) ;//+ (4./3.);
    
    flux(0,0) = (2.)*x/(alpha*alpha);
    flux(1,0) = (-2.)*y/(alpha*alpha);
    flux(2,0) = 0.;
    
    std::cout<<"funcao em omega2 "<<solp[0]<<std::endl;
    
}

void ExataOmega3(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    ExataOmega1(pt, solp, flux);
    std::cout<<"funcao em omega3 "<<solp[0]<<std::endl;

    
}

void ExataOmega4(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
    flux.Resize(3, 1);
    
    REAL x,y,z,alpha;
    
    alpha=1.;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    REAL alpha4 = pow(alpha, 4);
    
    solp[0] = (x*x-y*y)/(alpha4);
    
    flux(0,0) = 2.*x/(alpha4);
    flux(1,0) = (-2.)*y/(alpha4);
    flux(2,0) = 0.;
    std::cout<<"funcao em omega4 "<<solp[0]<<std::endl;
    
}

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &ff)
{
    ff[0]=0.;
}

void Neumann1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff)
{
   
    TPZFMatrix<STATE> flux;
    ExataOmega1(pt, ff, flux);
     ff[0]=0.;
    
    TPZFMatrix<STATE> normal(3,1);
    normal.Zero();
    
    normal(0,0)=1.;
    
    ff[0]= flux(0,0)*normal(0,0)+flux(1,0)*normal(1,0)+flux(2,0)*normal(2,0);
    std::cout<<"derivada normal em omega1 "<<ff[0]<<std::endl;
    
}

void Neumann2(const TPZVec<REAL> &pt, TPZVec<STATE> &ff)
{
    
    TPZFMatrix<STATE> flux;
    ExataOmega2(pt, ff, flux);
    ff[0]=0.;
    
    TPZFMatrix<STATE> normal(3,1);
    normal.Zero();
    
    normal(1,0) = 1.;
    
    ff[0] = flux(0,0)*normal(0,0)+flux(1,0)*normal(1,0)+flux(2,0)*normal(2,0);
    std::cout<<"derivada normal em omega2 "<<ff[0]<<std::endl;
    
}

void Neumann3(const TPZVec<REAL> &pt, TPZVec<STATE> &ff)
{
    
    TPZFMatrix<STATE> flux;
    
    ExataOmega3(pt, ff, flux);
    ff[0]=0.;
    
    TPZFMatrix<STATE> normal(3,1);
    normal.Zero();
    
    normal(0,0) = (-1.);
    
    ff[0] = flux(0,0)*normal(0,0)+flux(1,0)*normal(1,0)+flux(2,0)*normal(2,0);
    std::cout<<"derivada normal em omega3 "<<ff[0]<<std::endl;
    
}

void Neumann4(const TPZVec<REAL> &pt, TPZVec<STATE> &ff)
{
   
    TPZFMatrix<STATE> flux;
    
    ExataOmega4(pt, ff, flux);
     ff[0]=0.;
    
    TPZFMatrix<STATE> normal(3,1);
    normal.Zero();
    
    normal(1,0) = (-1.);
    
    ff[0] = flux(0,0)*normal(0,0)+flux(1,0)*normal(1,0)+flux(2,0)*normal(2,0);
    std::cout<<"derivada normal em omega4 "<<ff[0]<<std::endl;
    
}




TPZMultiphysicsCompMesh *CreateNeumannHDivMesh(const ProblemConfig &problem) {
    
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(problem.gmesh);
    TPZMaterial *mat = NULL;
    for (auto matid : problem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        
        TPZAutoPointer<TPZFunction<STATE> > solexata;
   //o termo do lado direito ff=0
        mix->SetInternalFlux(0);
    
        if(matid==1){
            solexata = new TPZDummyFunction<STATE>(ExataOmega1,10);
         //   mix->SetForcingFunctionExact(solexata);
             mix->SetPermeability(1.);
        
        }
        
        if(matid==2){
            solexata = new TPZDummyFunction<STATE>(ExataOmega2,10);
      //      mix->SetForcingFunctionExact(solexata);
            STATE alpha2=(problem.alpha)*(problem.alpha);
            
            mix->SetPermeability(alpha2);
            
            
            
        }
        
        if(matid==3){
            solexata = new TPZDummyFunction<STATE>(ExataOmega3,10);
   //         mix->SetForcingFunctionExact(solexata);
            mix->SetPermeability(1.);
            
        }
        if(matid==4){
            solexata = new TPZDummyFunction<STATE>(ExataOmega4,10);
    //        mix->SetForcingFunctionExact(solexata);
            STATE alpha4=pow(problem.alpha, 4);
            mix->SetPermeability(alpha4);
            
        }
 

        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    
 
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 1;
        TPZAutoPointer<TPZFunction<STATE> > bcfunction;
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        
        if(matid==5){
            bcfunction=new TPZDummyFunction<STATE>(ExataOmega1,5);
            int bctype=0;
//            val2(0,0) = 0.;
//            TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
            bc->SetType(bctype);
            TPZAutoPointer<TPZFunction<STATE>> func(bcfunction);
            bc->TPZMaterial::SetForcingFunction(func);
            cmesh->InsertMaterialObject(bc);
        }
        
        if(matid==6){
            bcfunction=new TPZDummyFunction<STATE>(Neumann2,5);
            int bctype = 1;
//            val2(0,0) = 0.;
//            TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
            bc->SetType(bctype);
         
            TPZAutoPointer<TPZFunction<STATE>> func(bcfunction);
            bc->TPZMaterial::SetForcingFunction(func);
            cmesh->InsertMaterialObject(bc);
        }
        if(matid==7){
            bcfunction=new TPZDummyFunction<STATE>(ExataOmega3,5);
            int  bctype=0;
//            val2(0,0) = 0.;
//            TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
             bc->SetType(bctype);
            TPZAutoPointer<TPZFunction<STATE>> func(bcfunction);
            bc->TPZMaterial::SetForcingFunction(func);
            cmesh->InsertMaterialObject(bc);
        }
        
        if(matid==8){
            bcfunction=new TPZDummyFunction<STATE>(Neumann4,5);
            int bctype = 1;
//            val2(0,0)=0;
//            TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
            bc->SetType(bctype);
            TPZAutoPointer<TPZFunction<STATE>> func(bcfunction);
            bc->TPZMaterial::SetForcingFunction(func);
            cmesh->InsertMaterialObject(bc);
        }
        
        
    }
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    
    TPZManVector<int> active(2,1);
    TPZManVector<TPZCompMesh *> meshvector(2,0);
    
    meshvector[0] = CreateNeumannFluxHDivMesh(problem);
    meshvector[1] = CreatePressureMesh(problem);
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    bool keepmatrix = false;
    TPZCompMeshTools::CreatedCondensedElements(cmesh, true, keepmatrix);
    
    return cmesh;
}

TPZCompMesh *CreateNeumannFluxHDivMesh(const ProblemConfig &problem) {
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
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 1;
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
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
        
        if(problem.prefine){
            Prefinamento(cmesh, problem.ndivisions, problem.porder);
        }
        
    }
    cmesh->InitializeBlock();
    return cmesh;
    
}

void SolveHybridProblem(TPZCompMesh *Hybridmesh,int InterfaceMatId,const ProblemConfig &problem){
    
    TPZAnalysis an(Hybridmesh);
    
#ifdef USING_MKL2
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
    //    matIds.insert(1);
    //    matIds.insert(-1);
    //    matIds.insert(-2);
    //
    //   // matIds.insert(4);
    
    
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
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    //   scalnames.Push("ExactPressure");
    vecnames.Push("Flux");
    //  vecnames.Push("ExactFlux");
    // scalnames.Push("Divergence");
    an.DefineGraphMesh(2, scalnames, vecnames, "OriginalHybrid_Problem.vtk");
    //        meshvec_Hybrid[1]->Solution().Print("Press");
    // Post processing
    an.PostProcess(2,2);
    
}
void PlotLagrangreMultiplier(TPZCompMesh *cmesh){
    
    TPZAnalysis an(cmesh,false);
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("State");
    
    int dim = cmesh->Reference()->Dimension()-1;
    std::string plotname;
    {
        std::stringstream out;
        out << "OriginalLagrangeMultiplier" << ".vtk";
        plotname = out.str();
    }
    an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
    an.PostProcess(2,dim);
    
}


