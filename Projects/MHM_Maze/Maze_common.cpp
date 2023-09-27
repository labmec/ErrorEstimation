
#include "Maze_common.h"
#include "TPZDarcyMHMHDivErrorEstimator.h"

TPZCompMesh* MixedTest(ConfigCasesMaze &Conf, int nx, int ny){
  
    TPZGeoMesh *gmesh = GenerateGeoMesh(Conf.GetImageName(),nx,ny);
    int flux_order = Conf.GetFluxOrder();
    int p_order = Conf.GetPressureOrder();
    
    {
#ifdef ERRORESTIMATION_DEBUG
        std::ofstream file("maze.txt");
        gmesh->Print(file);
        
        std::ofstream out("maze.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
#endif
    }
    
    TPZCompMesh *cmesh_flux =CMeshFlux(gmesh,flux_order);
    TPZCompMesh *cmesh_presure =CMeshPressure(gmesh,p_order,Conf);
    
    
    TPZVec<TPZCompMesh *> fmeshvec(2);
    fmeshvec[0]=cmesh_flux;
    fmeshvec[1]=cmesh_presure;
    gmesh->ResetReference();
    
    TPZCompMesh *MixedMesh = CMeshMultphysics(gmesh,fmeshvec,Conf);
    
    std::ofstream file("MixedCMesh.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(MixedMesh, file);
    
    std::ofstream out("MixedCMesh.txt");
    MixedMesh->Print(out);
    std::cout << "number of equations = " << MixedMesh->NEquations() << std::endl;

    //Solving the system:
    bool optimizeBandwidth = true;
    MixedMesh->InitializeBlock();
    
//    TPZCompMesh * cmesh_m_Hybrid;
//    TPZManVector<TPZCompMesh*, 3> meshvector_Hybrid(3);
//    TPZHybridizeHDiv hybridizer;
//    tie(cmesh_m_Hybrid, meshvector_Hybrid) = hybridizer.Hybridize(MixedMesh, fmeshvec, true, -1.);
//    cmesh_m_Hybrid->InitializeBlock();
    
    bool must_opt_band_width_Q = true;
    int number_threads = 4;
    TPZLinearAnalysis *an = new TPZLinearAnalysis(MixedMesh,must_opt_band_width_Q);
    
    //
    TPZSkylineStructMatrix<> skyl_mat(MixedMesh);
    TPZSSpStructMatrix<> sparse_matrix(MixedMesh);
    TPZStepSolver<STATE> step;
    sparse_matrix.SetNumThreads(number_threads);
    step.SetDirect(ELDLt);
//    an->SetStructuralMatrix(sparse_matrix);
    an->SetStructuralMatrix(skyl_mat);
    an->SetSolver(step);
    an->Assemble();
    std::cout << "Norm of rhs " << Norm(an->Rhs()) << std::endl;
    an->Solve();
    std::cout << "Norm of solution " << Norm(an->Solution()) << std::endl;

    std::cout << "Norm of solution in flux mesh " << Norm(cmesh_flux->Solution()) << std::endl;
    std::cout << "Norm of solution in pressure mesh " << Norm(cmesh_presure->Solution()) << std::endl;
    //POS
//    std::cout << "Norm of solution in hybrid flux mesh " << Norm(meshvector_Hybrid[0]->Solution()) << std::endl;
//    std::cout << "Norm of solution in hybrid pressure mesh " << Norm(meshvector_Hybrid[1]->Solution()) << std::endl;
    TPZManVector<std::string,10> scalnames(2), vecnames(1);
    vecnames[0]  = "Flux";
    scalnames[0] = "Pressure";
    scalnames[1] = "Permeability";

    const int dim = an->Mesh()->Dimension();
    int div = 0;
    
    an->DefineGraphMesh(dim,scalnames,vecnames,Conf.GetVTKName());
    an->PostProcess(div,dim);
    std::cout << "Standard post-processing finished." << std::endl;
//    {
//        std::ofstream out("bybrid_mesh.txt");
//        cmesh_m_Hybrid->Print(out);
//        std::ofstream out2("flux_hybrid.txt");
//        meshvector_Hybrid[0]->Print(out2);
//    }
    return cmesh_flux;
}

TPZCompMesh *CMeshPressure(TPZGeoMesh * gmesh, int pOrder,ConfigCasesMaze Conf){
    
    int impervious_mat = 1;
    int permeable_mat = 2;
    int dim = gmesh->Dimension();
    
    REAL perm_0 = Conf.GetImperviousMatPermeability();
    REAL perm_1 = Conf.GetPermeableMatPermeability();
    
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    typedef TPZDarcyFlow TPZMatPoisson3d;
    TPZMatPoisson3d *mat_0 = new TPZMatPoisson3d(impervious_mat,dim);
    TPZMatPoisson3d *mat_1 = new TPZMatPoisson3d(permeable_mat,dim);
    mat_0->SetConstantPermeability(perm_0);
    mat_1->SetConstantPermeability(perm_1);
    
    //  inserting volumetric materials objects
    cmesh->InsertMaterialObject(mat_0);
    cmesh->InsertMaterialObject(mat_1);
    
    
    cmesh->SetName("Pressure");
    cmesh->AutoBuild();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
#ifdef ERRORESTIMATION_DEBUG
    std::ofstream out("cmeshPress.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
    
}
TPZGeoMesh *GeoMeshFromPng(string name, double &l, double &h){
    const int bcDL = -1;
    const int bcB = -2;
    const int bcDR = -3;
    const int bcDT = -4;


#ifdef MACOSX
    Mat image = imread(name,IMREAD_GRAYSCALE);
#else
    Mat image = imread(name,IMREAD_GRAYSCALE);
#endif

    int k=0;
    int px=image.size[0];
    int py=image.size[1];
    l=px;
    h=py;
    int p =px*py;
    if(p==0){
        DebugStop();
    }
    vector<int> vec(p,0);
    
    for (int i = 0; i<px; ++i) {
        for (int j = py  ; j>0; --j) {
            int val =(int)image.at<uchar>(Point(j, i));
            if (val>200){
                val=255;
            }
            int pix =val/255;
            vec[p-k]=pix;
            
            
            
            k++;
        }
    }
    
    
    
    // Creating the Geo mesh
    TPZManVector<REAL,3> x0(3,0.),x1(3,px);
    x1[2] = 0.;
    TPZManVector<int,2> nelx(2,py);
    nelx[0] = px;
    TPZGenGrid2D gengrid(nelx,x0,x1);
    gengrid.SetElementType(MMeshType::EQuadrilateral);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gengrid.Read(gmesh);
    //gengrid.Read(gmesh,2);
    
    //MatsID
    int nels = gmesh->NElements();
    TPZGeoEl *gel_in;
    TPZGeoEl *gel_out;
    TPZGeoEl *gel_in1D;
    TPZGeoEl *gel_out1D;
    
    
    for (int i=0; i<nels; i++) {
        TPZGeoEl *gel =gmesh->Element(i);
        gel->SetMaterialId(vec[i]+ 1 );
        
        if (i<= px) {
            if((vec[i]+1)==2){
                gel_in =gel;
            }
        }
        
        if (i >= (px)*(py-1)) {
            if((vec[i]+1)==2){
                gel_out=gel;
            }
        }
        
    }
    
    //gengrid.SetBC(TPZGeoMesh *gr, int side, int bc)
    gengrid.SetBC(gmesh, 4, bcDL);
    gengrid.SetBC(gmesh, 5, bcB);
    gengrid.SetBC(gmesh, 6, bcDR);
    gengrid.SetBC(gmesh, 7, bcDT);
    
    
    int gel_in_index = gel_in->Index();
    gel_in1D = gmesh->Element(gel_in_index)->Neighbour(4).Element();
    gel_in1D->SetMaterialId(-5);
    
    int gel_out_index = gel_out->Index();
    gel_out1D = gmesh->Element(gel_out_index)->Neighbour(6).Element();
    gel_out1D->SetMaterialId(-6);
    
    
    gmesh->BuildConnectivity();
    return gmesh;
}
TPZCompMesh *CMeshMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, ConfigCasesMaze &Conf ){
    
    //Creating computational mesh for multiphysic elements
    TPZMultiphysicsCompMesh *mphysics = new TPZMultiphysicsCompMesh(gmesh);
    
    int impervious_mat = 1;
    int permeable_mat = 2;
    int dim = gmesh->Dimension();
    
    REAL perm_0 = Conf.GetImperviousMatPermeability();
    REAL perm_1 = Conf.GetPermeableMatPermeability();
    
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    std::cout<<mphysics->NMaterials();
    
    typedef TPZMixedDarcyFlow TPZMixedPoisson;
    
    TPZMixedPoisson *mat_0 = new TPZMixedPoisson(impervious_mat,dim);
    mat_0->SetConstantPermeability(perm_0);
    
    
    TPZMixedPoisson *mat_1 = new TPZMixedPoisson(permeable_mat,dim);
    mat_1->SetConstantPermeability(perm_1);
    
    
    mphysics->InsertMaterialObject(mat_0);
    mphysics->InsertMaterialObject(mat_1);
    
    //Inserir condicoes de contorno
    int type_D = 0;
    int type_N = 1;
    TPZFMatrix<STATE> val1(1, 1, 0.);
    TPZManVector<STATE> val2(1, 0.);
    
    // Insert boundary conditions
    //Neumann boundary conditions (flux = 0)
    int right_bc_id = -2;
    val2[0] = 0.0;
    TPZBndCondT<STATE> * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_N, val1, val2);
    mphysics->InsertMaterialObject(right_bc);
    
    
    int left_bc_id = -4;
    val2[0] = 0.0;
    TPZBndCondT<STATE> * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_N, val1, val2);
    mphysics->InsertMaterialObject(left_bc);
    
    int bottom_bc_1id = -1;
    val2[0] = 0;
    TPZBndCondT<STATE> * bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_N, val1, val2);
    mphysics->InsertMaterialObject(bottom_bc_1);
    
    int top_bc_1id = -3;
    val2[0] = 0.0;
    TPZBndCondT<STATE> * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_N, val1, val2);
    mphysics->InsertMaterialObject(top_bc_1);
    
    //Dirichlet Conditions (p=1 in, p=0 out)
    int bottom_bc_id = -5;
    val2[0] = Conf.GetCCPressureIn();
    TPZBndCondT<STATE> * bottom_bc = mat_0->CreateBC(mat_0, bottom_bc_id, type_D, val1, val2);
    mphysics->InsertMaterialObject(bottom_bc);
    
    int top_bc_id = -6;
    val2[0] = Conf.GetCCPressureOut();
    TPZBndCondT<STATE> * top_bc = mat_0->CreateBC(mat_0, top_bc_id, type_D, val1, val2);
    mphysics->InsertMaterialObject(top_bc);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->SetDimModel(gmesh->Dimension());
//    mphysics->AutoBuild();
    
    mphysics->BuildMultiphysicsSpace(meshvec);
    
//    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
//    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
//    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
#ifdef ERRORESTIMATION_DEBUG
    std::ofstream file("cmesh_mphysics.txt");
    mphysics->Print(file);
#endif
    
    return mphysics;
}

void EstimateError(TPZDarcyMHMHDivErrorEstimator &errorEstimator, ProblemConfig &config) {
    cout << "Error Estimation processing for MHM-Hdiv problem " << endl;

    DebugStop();
    
//    errorEstimator.SetProblemConfig(config);
    errorEstimator.PotentialReconstruction();

    {
        string command = "mkdir -p " + config.dir_name;
        system(command.c_str());

        TPZManVector<REAL, 6> errors;
        TPZManVector<REAL> element_errors;
        std::string vtk_path = "mhm_maze_results.vtk";
        errorEstimator.ComputeErrors(errors, element_errors, vtk_path);
    }
}

void InsertMaterialObjects(TPZMHMixedMeshControl &control)
{
    TPZCompMesh &cmesh = control.CMesh();

    TPZGeoMesh &gmesh = control.GMesh();
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2Flux(1,0.), val2Pressure(1,10.);


    int dim = gmesh.Dimension();
    cmesh.SetDimModel(dim);

    TPZCompMesh *MixedFluxPressureCmesh = &cmesh;

    typedef TPZMixedDarcyFlow TPZMixedPoisson;
    // Material medio poroso
    TPZMixedPoisson * mat = new TPZMixedPoisson(1,dim);
    mat->SetConstantPermeability(1.);
    //    mat->SetForcingFunction(One);
    MixedFluxPressureCmesh->InsertMaterialObject(mat);

    TPZMixedPoisson * mat_2 = new TPZMixedPoisson(2,dim);
    //mat_2->SetPermeability(1000000.0);
    mat_2->SetConstantPermeability(500000.0);
    //    mat->SetForcingFunction(One);
    MixedFluxPressureCmesh->InsertMaterialObject(mat_2);

   // Bc N
    TPZBndCondT<STATE> * bcN = mat->CreateBC(mat, -1, typeFlux, val1, val2Flux);
    //    bcN->SetForcingFunction(0, force);

    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    bcN = mat->CreateBC(mat, -3, typeFlux, val1, val2Flux);
    //    bcN->SetForcingFunction(0, force);

    MixedFluxPressureCmesh->InsertMaterialObject(bcN);

    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, -2, typeFlux, val1, val2Flux);

    MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    bcS = mat->CreateBC(mat, -4, typeFlux, val1, val2Flux);
    MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    val2Pressure[0] = 100.;
    TPZBndCond * bcIn = mat->CreateBC(mat, -5, typePressure, val1, val2Pressure);

    MixedFluxPressureCmesh->InsertMaterialObject(bcIn);
    val2Pressure[0] = -100.;
    TPZBndCond * bcOut = mat->CreateBC(mat, -6, typePressure, val1, val2Pressure);

    MixedFluxPressureCmesh->InsertMaterialObject(bcOut);

}

TPZGeoMesh *GenerateGeoMesh(string name, int nx, int ny){

    double l;
    double h;
    TPZGeoMesh *FineMesh = GeoMeshFromPng(name,l,h);

    std::ofstream file_base_vtk("base.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(FineMesh, file_base_vtk);

    std::ofstream file_base_txt("base.txt");
    FineMesh->Print(file_base_txt);

    // Creating the Geo mesh
    TPZManVector<REAL,3> x0(3,0.),x1(3,l);
    x1[1] = h;
    x1[2] = 0.;
    TPZManVector<int,2> nelx(2,ny);
    nelx[0] = nx;
    TPZGenGrid2D gengrid(nelx,x0,x1);
    gengrid.SetElementType(MMeshType::EQuadrilateral);
    TPZGeoMesh *gmeshcoarse = new TPZGeoMesh;
    gmeshcoarse->SetDimension(2);
    gengrid.SetRefpatternElements(true);
    gengrid.Read(gmeshcoarse);
    //gengrid.Read(gmesh,2);


    //gengrid.SetBC(TPZGeoMesh *gr, int side, int bc)
    gengrid.SetBC(gmeshcoarse, 4, -1);
    gengrid.SetBC(gmeshcoarse, 5, -2);
    gengrid.SetBC(gmeshcoarse, 6, -3);
    gengrid.SetBC(gmeshcoarse, 7, -4);

    gmeshcoarse->BuildConnectivity();

    std::ofstream file_base_c_vtk("base_c.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmeshcoarse, file_base_c_vtk);

    std::ofstream file_base_c_txt("base_c.txt");
    gmeshcoarse->Print(file_base_c_txt);

    //Refine
    TPZVec<REAL> qsi(3,0);
    TPZVec<REAL> result(3,0);
    TPZStack<TPZVec<int64_t>> vecs;
    TPZVec<TPZGeoEl*> indexf;
//    int nref=6; para 128x128 com coarse 2x2
    int nref = log2(l) - log2(nx);
    for(int i=0; i<nref; i++){
        int nel = gmeshcoarse->NElements();
        for(int i=0; i<nel; i++){
            TPZGeoEl * gel = gmeshcoarse->Element(i);
            if (!gel || gel->HasSubElement()) {
                continue;
            }
            gel->Divide(indexf);
        }
    }

    //
    int nel = gmeshcoarse->NElements();

    for(int i=0; i<nel; i++){
        TPZGeoEl *gel = gmeshcoarse->Element(i);
        if (!gel || gel->HasSubElement()) {
            continue;
        }

        if(gel->Dimension()==2){
            TPZFMatrix<REAL> cooridnates1(3,4);
            TPZVec<REAL> qsi(3,0);
            TPZVec<REAL> result(3,0);
            gel->X(qsi,result);
            int flor =floor(result[0]);
            int y =floor(result[1])*l;
            int pos = flor + y;
            TPZGeoEl *gel2 = FineMesh->Element(pos);
            if(!gel2){
                DebugStop();
            }
            int matid= gel2->MaterialId();
            gel->SetMaterialId(matid);

            if(y==0 && matid==2){
                TPZGeoEl *el1D = gel->Neighbour(4).Element();
                el1D->SetMaterialId(-5);
            }
            int niv =y/l;
            if(niv==(l-1) && matid==2){
                TPZGeoEl *el1D = gel->Neighbour(6).Element();
                el1D->SetMaterialId(-6);
            }

        }
    }

    std::ofstream out("mazefine.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(FineMesh, out, true);

    std::ofstream out2("mazehcoarse.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmeshcoarse, out2, true);

    return gmeshcoarse;
}

void LocateElementsToAdapt(TPZDarcyMHMHDivErrorEstimator &errorEstimator, ProblemConfig &config) {

    TPZMultiphysicsCompMesh *postProcMesh = errorEstimator.PostProcMesh();
    postProcMesh->LoadReferences();
    TPZGeoMesh * gmesh = postProcMesh->Reference();

    {
        std::ofstream out("gmeshToAdapt2.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    }

    // This variable stores the difference ||grad u_rec - sigma_fem|| over an skeleton element and the index
    // of this element
    std::set<std::pair<REAL, int64_t>> fluxDiffPerSkeleton;

    int64_t nel = gmesh->NElements();
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZGeoEl *skel_gel = gmesh->Element(iel);
        if (!skel_gel) DebugStop();
        // Filters skeleton gels
        if (skel_gel->MaterialId() != errorEstimator.PressureSkeletonMatId()) continue;

        TPZGeoElSide skel_side(skel_gel);
        if (skel_side.NNeighbours() != 2) DebugStop(); // The skeleton element should have only right/left neighbours

        // Store right and left neighbours
        TPZGeoElSide right_side = skel_side.Neighbour();
        TPZGeoElSide left_side = right_side.Neighbour();

        TPZGeoEl *right_gel = right_side.Element();
        TPZGeoEl *left_gel = left_side.Element();

        if (!right_gel || !left_gel) DebugStop();
        if (right_gel->Dimension() != 2 || left_gel->Dimension() != 2) DebugStop();

        TPZCompEl * right_cel = right_gel->Reference();
        TPZCompEl * left_cel = left_gel->Reference();
        if (!right_cel || !left_cel) DebugStop();

        TPZMultiphysicsElement *m_cel_right = dynamic_cast<TPZMultiphysicsElement *>(right_cel);
        TPZMultiphysicsElement *m_cel_left = dynamic_cast<TPZMultiphysicsElement *>(left_cel);
        if (!m_cel_right || !m_cel_left) DebugStop();

        TPZVec<TPZMaterialDataT<STATE>> right_mat_data(4);
        TPZVec<TPZMaterialDataT<STATE>> left_mat_data(4);
        TPZManVector<int64_t, 4> indexes(4, 0);

        TPZManVector<int> active_spaces = postProcMesh->GetActiveApproximationSpaces();

        indexes[1] = 1;
        indexes[2] = 1;
        m_cel_right->InitMaterialData(right_mat_data, &indexes);
        m_cel_left->InitMaterialData(left_mat_data, &indexes);

        // Transformation of the skeleton to the linear side of the right element
        TPZTransform<REAL> skel_to_right_trans = skel_side.NeighbourSideTransform(right_side);
        // Transformation of the linear side of the right element to its face
        TPZTransform<REAL> tmp = right_gel->SideToSideTransform(right_side.Side(), right_gel->NSides() - 1);
        // Combine the transformations
        skel_to_right_trans = tmp.Multiply(skel_to_right_trans);

        // Analogously for the left neighbour
        TPZTransform<REAL> skel_to_left_trans = skel_side.NeighbourSideTransform(left_side);
        tmp = left_gel->SideToSideTransform(left_side.Side(), left_gel->NSides() - 1);
        skel_to_left_trans = tmp.Multiply(skel_to_left_trans);

        int order = 2; // TODO Gustavo: think if 2 is a good idea
        std::unique_ptr<TPZIntPoints> intRule(skel_gel->CreateSideIntegrationRule(skel_side.Side(), order));
        const int npts = intRule->NPoints();
        TPZManVector<REAL, 3> xi_skel(1, 0);
        TPZManVector<REAL, 3> xi_right(2, 0);
        TPZManVector<REAL, 3> xi_left(2, 0);
        REAL w;

        REAL diff = 0;
        for (auto ipt = 0; ipt < npts; ipt++) {
            intRule->Point(ipt, xi_skel, w);

            skel_to_right_trans.Apply(xi_skel, xi_right);
            skel_to_left_trans.Apply(xi_skel, xi_left);

            if (false)
            {
                TPZManVector<REAL> right_x(3);
                TPZManVector<REAL> left_x(3);
                right_gel->X(xi_right, right_x);
                left_gel->X(xi_left, left_x);

                std::cout << "R: (" << right_x[0] << ", " << right_x[1] << ", " << right_x[2] << ")\n";
                std::cout << "L: (" << left_x[0] << ", " << left_x[1] << ", " << left_x[2] << ")\n";
            }
            DebugStop();
//            TPZVec<TPZTransform<>> tr_vec(0);
//            m_cel_right->ComputeRequiredData(xi_right, tr_vec, right_mat_data, indexes);
            //m_cel_left->ComputeRequiredData(left_mat_data, xi_left);
        }
    }
}
