/**
 * @file
 * @brief Contains the TPZMatLaplacian class.
 */

#include "ConfigCasesMaze.h"
#include "TPZGenGrid2D.h"
#include "TPZMultiphysicsCompMesh.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZNullMaterial.h"
#include "TPZVTKGeoMesh.h"

#include <cmath>
//#include <opencv2/opencv.hpp>
#include <set>
#include <string>
//using namespace cv;
#include "lodepng.h"

TPZGeoMesh *ConfigCasesMaze::GeoMeshFromPng(const std::string &name, int &l, int &h){
//    {
//        std::ofstream test("total4.txt");
//        test<<"oiii"<<"\n";
//
//    }
    const int bcDL = -1;
    const int bcB = -2;
    const int bcDR = -3;
    const int bcDT = -4;



    std::vector<int> vec;
    unsigned width, height;

    {
        std::vector<unsigned char> image; //the raw pixels

//        std::string filename = name;
//        filename = "maze8x8.png";
//        const char *cstr = name.c_str();
//        for(int i=0; i<strlen(cstr); i++) std::cout << cstr[i] << ' '<< (int)cstr[i] << std::endl;
        //decode
        unsigned error = lodepng::decode(image, width, height, name);
        l = width;
        h = height;
        vec.resize(width*height);
        
        //if there's an error, display it
        if(error) {
            std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
            DebugStop();
        }
        for(int i=height-1; i>=0; i--)
        {
            for (int j=0; j<width; j++) {
                int pos = i*width+j;
//                std::cout << i << " " << j << " " << (int)image[pos*4] << std::endl;
                (int)image[pos*4] == 0 ? vec[pos]=0 : vec[pos] = 1;
            }
        }
    }

//    if(0)
//    {
//        Mat image = imread(name,IMREAD_GRAYSCALE);
//
//        int k=0;
//        int px=image.size[0];
//        int py=image.size[1];
//        l=px;
//        h=py;
//        int p =px*py;
//        if(p==0){
//            DebugStop();
//        }
//
//        for (int i = 0; i<px; ++i) {
//            for (int j = py  ; j>0; --j) {
//                int val =(int)image.at<uchar>(Point(j, i));
//                if (val>200){
//                    val=255;
//                }
//                int pix =val/255;
//                vec[p-k]=pix;
//
//
//
//                k++;
//            }
//        }
//    
//    }
    
    // Creating the Geo mesh
    TPZManVector<REAL,3> x0(3,0.),x1(3,width);
    x1[2] = 0.;
    TPZManVector<int,2> nelx(2,height);
    nelx[0] = width;
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
        
        if (i<= width) {
            if((vec[i]+1)==2){
                gel_in =gel;
            }
        }
        
        if (i >= (width)*(height-1)) {
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

TPZMultiphysicsCompMesh *ConfigCasesMaze::CMeshMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> &meshvec){
    
    ConfigCasesMaze &Conf = *this;
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

TPZCompMesh *ConfigCasesMaze::CMeshPressure(TPZGeoMesh * gmesh, int pOrder){
    
    ConfigCasesMaze &Conf = *this;
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
    cmesh->SetAllCreateFunctionsContinuous();
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
    cmesh->ExpandSolution();
#ifdef ERRORESTIMATION_DEBUG
    std::ofstream out("cmeshPress.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
    
}

TPZCompMesh *ConfigCasesMaze::CMeshFlux(TPZGeoMesh * gmesh,int pOrder){
    
    int impervious_mat = 1;
    int permeable_mat = 2;
    int dim = gmesh->Dimension();
    
    TPZVec<REAL> convdir(dim,0.0);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    
    //Definition of the approximation space:
    
    TPZNullMaterial<> *mat_0 = new TPZNullMaterial<>(impervious_mat,dim);
    TPZNullMaterial<> *mat_1 = new TPZNullMaterial<>(permeable_mat,dim);
    
    //  inserting volumetric materials objects
    cmesh->InsertMaterialObject(mat_0);
    cmesh->InsertMaterialObject(mat_1);
    
    cmesh->SetAllCreateFunctionsHDiv(); //Creating H(div) functions
    
    
    int type_D = 0;
    int type_N = 1;
    TPZFMatrix<STATE> val1(1, 1, 0.);
    TPZManVector<REAL,1> val2(1, 0.);
    
    // Insert boundary conditions
    //Neumann boundary conditions (flux = 0)
    int right_bc_id = -2;
    TPZBndCondT<STATE> * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_N, val1, val2);
    cmesh->InsertMaterialObject(right_bc);
    
    int left_bc_id = -4;
    TPZBndCondT<STATE> * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_N, val1, val2);
    cmesh->InsertMaterialObject(left_bc);
    
    int bottom_bc_1id = -1;
    TPZBndCondT<STATE>* bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_N, val1, val2);
    cmesh->InsertMaterialObject(bottom_bc_1);
    
    int top_bc_1id = -3;
    TPZBndCondT<STATE> * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_N, val1, val2);
    cmesh->InsertMaterialObject(top_bc_1);
    
    
    //Dirichlet Conditions (p=1 in, p=0 out)
    int bottom_bc_id = -5;
    TPZBndCondT<STATE> * bottom_bc = mat_0->CreateBC(mat_0, bottom_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(bottom_bc);
    
    int top_bc_id = -6;
    TPZBndCondT<STATE> * top_bc = mat_0->CreateBC(mat_0, top_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(top_bc);
    
    cmesh->SetName("LaberintoFlux");
    cmesh->AutoBuild();
    cmesh->ExpandSolution();
#ifdef ERRORESTIMATION_DEBUG
    std::ofstream file("cmesh_flux.txt");
    cmesh->Print(file);
#endif
    
    return cmesh;
    
}

#include "pzinterpolationspace.h"
/// create a constant value mesh. The singleconnect parameter indicates the mesh is associated
/// with a subdomain and should have a single connect
TPZCompMesh *ConfigCasesMaze::CMeshAverage(TPZGeoMesh *gmesh, bool singleconnect)
{
    int pOrder = 0;
    int impervious_mat = 1;
    int permeable_mat = 2;
    int dim = gmesh->Dimension();
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    
    //Definition of the approximation space:
    
    TPZNullMaterial<> *mat_0 = new TPZNullMaterial<>(impervious_mat,dim);
    TPZNullMaterial<> *mat_1 = new TPZNullMaterial<>(permeable_mat,dim);
    
    //  inserting volumetric materials objects
    cmesh->InsertMaterialObject(mat_0);
    cmesh->InsertMaterialObject(mat_1);
    
    cmesh->SetAllCreateFunctionsDiscontinuous(); //Creating H(div) functions
    cmesh->AutoBuild();
    if(singleconnect)
    {
        int64_t nelem = cmesh->NElements();
        for (int64_t el = 0; el<nelem; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            if(!intel) continue;
            intel->SetConnectIndex(0, 0);
        }
    }
    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    
    
#ifdef ERRORESTIMATION_DEBUG
    std::ofstream file("cmesh_constant.txt");
    cmesh->Print(file);
#endif
    
    return cmesh;
    
}

TPZGeoMesh *ConfigCasesMaze::GenerateGeoMesh(const std::string &name, int nx, int ny){

    int l;
    int h;
    TPZGeoMesh *FineMesh = this->GeoMeshFromPng(name,l,h);

    {
        std::ofstream file_base_vtk("base.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(FineMesh, file_base_vtk);

        std::ofstream file_base_txt("base.txt");
        FineMesh->Print(file_base_txt);
    }
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

TPZMultiphysicsCompMesh *ConfigCasesMaze::BuildMesh(int nx, int ny)
{
    TPZGeoMesh *gmesh = GenerateGeoMesh(GetImageName(), nx, ny);
    //GenerateGeoMesh(Conf.GetImageName(),nx,ny);
    int flux_order = GetFluxOrder();
    int p_order = GetPressureOrder();
    
    {
#ifdef ERRORESTIMATION_DEBUG
        std::ofstream file("maze.txt");
        gmesh->Print(file);
        
        std::ofstream out("maze.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
#endif
    }
    
    TPZCompMesh *cmesh_flux = CMeshFlux(gmesh,flux_order);
    TPZCompMesh *cmesh_presure = CMeshPressure(gmesh,p_order);
    
    
    TPZVec<TPZCompMesh *> fmeshvec(6);
    fmeshvec[0]=cmesh_flux;
    fmeshvec[1]=cmesh_presure;
    fmeshvec[2]=CMeshAverage(gmesh, false);
    fmeshvec[3]=CMeshAverage(gmesh, false);
    fmeshvec[4]=CMeshAverage(gmesh, true);
    fmeshvec[5]=CMeshAverage(gmesh, true);
    gmesh->ResetReference();
    
    TPZMultiphysicsCompMesh *MixedMesh = CMeshMultphysics(gmesh,fmeshvec);

    ConfigureLagrangeLevels(MixedMesh);
    return MixedMesh;
}

/// set the lagrange multipliers so that the global stiffness can be inverted without pivoting
void ConfigCasesMaze::ConfigureLagrangeLevels(TPZMultiphysicsCompMesh *cmesh)
{
    // separate the connects into sets
    std::set<int64_t> internalflux, lambdaflux, onepressure, elpressure,
     gelement, avpressureel, boundflux;
    int64_t gdomain, avpressuredomain, oneelaverage; //, onelambdaflux;
    
    auto meshvec = cmesh->MeshVector();
    auto fluxmesh = meshvec[0];
    auto pressuremesh = meshvec[1];
    
    // identify the first connect index in the multiphysics mesh
    TPZManVector<int64_t,7> firstconnect(7,0);
    for (int mesh=0; mesh < meshvec.size(); mesh++) {
        firstconnect[mesh+1] = firstconnect[mesh]+meshvec[mesh]->NConnects();
    }
    // for each flux element
    // if the element is boundary include the connect in boundflux
    // if the element is internal include internal flux and lambdaflux
    // remove boundflux from lambdaflux
    // identify a onelambdaflux and remove from lambdaflux
    
    // for each pressure element
    // include in onepressure
    // include other connects in elpressure
    
    // for each distr fluxel element
    // include in gelement
    
    // for each av pressure element
    // include in avpressureel
    
    // identify gdomain index - mesh 4
    // identify av pressure domain index mesh 5

    {
        // for each flux element
        // if the element is boundary include the connect in boundflux
        // if the element is internal include internal flux and lambdaflux
        // remove boundflux from lambdaflux
        // identify a onelambdaflux and remove from lambdaflux
        int64_t nelem = fluxmesh->NElements();
        std::set<int64_t> alllambda;
        for (int64_t el = 0; el<nelem; el++) {
            TPZCompEl *cel = fluxmesh->Element(el);
            int nc = cel->NConnects();
            if(nc == 1)
            {
                // we have a boundary element
                boundflux.insert(cel->ConnectIndex(0));
            }
            else
            {
                for (int ic=0; ic<nc-1; ic++) {
                    alllambda.insert(cel->ConnectIndex(ic));
                }
                internalflux.insert(cel->ConnectIndex(nc-1));
            }
        }
        // remove the boundflux elements from the lambdaflux elements
        std::set_difference(alllambda.begin(), alllambda.end(), boundflux.begin(), boundflux.end(), inserter(lambdaflux, lambdaflux.end()));
//        onelambdaflux = *lambdaflux.begin();
//        lambdaflux.erase(onelambdaflux);
#ifdef PZDEBUG
        if(boundflux.size()+lambdaflux.size()+internalflux.size() != fluxmesh->NConnects())
        {
            DebugStop();
        }
#endif
    }
    {
        // for each pressure element
        // include in onepressure
        // include other connects in elpressure
        int64_t nelem = pressuremesh->NElements();
        for (int64_t el = 0; el<nelem ; el++) {
            TPZCompEl *cel = pressuremesh->Element(el);
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace*>(cel);
            if(!intel) DebugStop();
            int nc = intel->NConnects();
            if(nc < 2) DebugStop();
            onepressure.insert(intel->ConnectIndex(0));
            for (int ic=1; ic<nc; ic++) {
                elpressure.insert(intel->ConnectIndex(ic));
            }
        }
    }
    {
        // for each distr fluxel element
        // include in gelement
        TPZCompMesh *mesh = meshvec[2];
        int64_t nelem = mesh->NElements();
        for (int el = 0; el<nelem; el++) {
            TPZCompEl *cel = mesh->Element(el);
#ifdef PZDEBUG
            int nc = cel->NConnects();
            if(nc != 1) DebugStop();
#endif
            gelement.insert(cel->ConnectIndex(0));
        }
    }
    {
        // for each av pressure element
        // include in avpressureel
        TPZCompMesh *mesh = meshvec[3];
        int64_t nelem = mesh->NElements();
        for (int el = 0; el<nelem; el++) {
            TPZCompEl *cel = mesh->Element(el);
#ifdef PZDEBUG
            int nc = cel->NConnects();
            if(nc != 1) DebugStop();
#endif
            avpressureel.insert(cel->ConnectIndex(0));
        }
        oneelaverage = *avpressureel.begin();
        avpressureel.erase(oneelaverage);
    }
    {
    // identify gdomain index - mesh 4
        TPZCompMesh *mesh = meshvec[4];
#ifdef PZDEBUG
            int neq = mesh->NEquations();
            if(neq != 1) DebugStop();
#endif
        gdomain = 0;
    }
    {
        // identify av pressure domain index mesh 5
        TPZCompMesh *mesh = meshvec[5];
#ifdef PZDEBUG
            int neq = mesh->NEquations();
            if(neq != 1) DebugStop();
#endif
        avpressuredomain = 0;
    }
    
    // ADJUST THE LAGRANGE MULTIPLIER LEVELS
    
    // adjust lagrange multiplier level of fluxes
    for(auto it:internalflux)
    {
        TPZConnect &catomic = fluxmesh->ConnectVec()[it];
        catomic.SetLagrangeMultiplier(0);
        TPZConnect &cmulti = cmesh->ConnectVec()[it];
        cmulti.SetLagrangeMultiplier(0);
    }
    for(auto it:lambdaflux)
    {
        TPZConnect &catomic = fluxmesh->ConnectVec()[it];
        catomic.SetLagrangeMultiplier(4);
        TPZConnect &cmulti = cmesh->ConnectVec()[it];
        cmulti.SetLagrangeMultiplier(4);
    }
    for(auto it:boundflux)
    {
        TPZConnect &catomic = fluxmesh->ConnectVec()[it];
        catomic.SetLagrangeMultiplier(8);
        TPZConnect &cmulti = cmesh->ConnectVec()[it];
        cmulti.SetLagrangeMultiplier(8);
    }
//    {
//        TPZConnect &catomic = fluxmesh->ConnectVec()[onelambdaflux];
//        catomic.SetLagrangeMultiplier(7);
//        TPZConnect &cmulti = cmesh->ConnectVec()[onelambdaflux];
//        cmulti.SetLagrangeMultiplier(7);
//    }
    
    // set lagrange multiplier for the pressure mesh
    for(auto it:onepressure)
    {
        TPZConnect &catomic = pressuremesh->ConnectVec()[it];
        catomic.SetLagrangeMultiplier(1);
        TPZConnect &cmulti = cmesh->ConnectVec()[it+firstconnect[1]];
        cmulti.SetLagrangeMultiplier(1);
    }
    for(auto it:elpressure)
    {
        TPZConnect &catomic = pressuremesh->ConnectVec()[it];
        catomic.SetLagrangeMultiplier(3);
        TPZConnect &cmulti = cmesh->ConnectVec()[it+firstconnect[1]];
        cmulti.SetLagrangeMultiplier(3);
    }

    // set the lagrange level of g element
    for(auto it:gelement)
    {
        TPZCompMesh *mesh = meshvec[2];
        TPZConnect &catomic = mesh->ConnectVec()[it];
        catomic.SetLagrangeMultiplier(2);
        TPZConnect &cmulti = cmesh->ConnectVec()[it+firstconnect[2]];
        cmulti.SetLagrangeMultiplier(2);
    }

    // set the lagrange level of average pressure element
    for(auto it:avpressureel)
    {
        TPZCompMesh *mesh = meshvec[3];
        TPZConnect &catomic = mesh->ConnectVec()[it];
        catomic.SetLagrangeMultiplier(5);
        TPZConnect &cmulti = cmesh->ConnectVec()[it+firstconnect[3]];
        cmulti.SetLagrangeMultiplier(5);
    }
    // setting the lagrange level of one element average
    {
        TPZCompMesh *mesh = meshvec[3];
        TPZConnect &catomic = mesh->ConnectVec()[oneelaverage];
        catomic.SetLagrangeMultiplier(7);
        TPZConnect &cmulti = cmesh->ConnectVec()[oneelaverage+firstconnect[3]];
        cmulti.SetLagrangeMultiplier(7);

    }
    // set the lagrange level of the distributed flux of the domain
    {
        TPZCompMesh *mesh = meshvec[4];
        TPZConnect &catomic = mesh->ConnectVec()[gdomain];
        catomic.SetLagrangeMultiplier(6);
        TPZConnect &cmulti = cmesh->ConnectVec()[gdomain+firstconnect[4]];
        cmulti.SetLagrangeMultiplier(6);
    }
    // set the lagrange level of the average pressure of the domain
    {
        TPZCompMesh *mesh = meshvec[5];
        TPZConnect &catomic = mesh->ConnectVec()[avpressuredomain];
        catomic.SetLagrangeMultiplier(9);
        TPZConnect &cmulti = cmesh->ConnectVec()[avpressuredomain+firstconnect[5]];
        cmulti.SetLagrangeMultiplier(9);
    }

}
