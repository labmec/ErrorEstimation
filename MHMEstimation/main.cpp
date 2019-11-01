#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <fstream>
#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"
#include "pzreal.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZParSkylineStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "TPZGeoLinear.h"

#include "pzpoisson3d.h"
//#include "pzhybridpoisson.h"
#include "pzpoisson3dreferred.h"
#include "mixedpoisson.h"
#include "pzelasmat.h"
#include "pzelasthybrid.h"
#include "pzmat1dlin.h"
#include "TPZVecL2.h"
#include "TPZMatLaplacianHybrid.h"
#include "TPZLagrangeMultiplier.h"


#include "pzbuildmultiphysicsmesh.h"
#include "pzelementgroup.h"
#include "TPZCompMeshTools.h"
#include "pzcondensedcompel.h"
#include "pzfunction.h"
#include "pzgraphmesh.h"
#include "pzfmatrix.h"

#include "pzlog.h"

#include "TPZVTKGeoMesh.h"
#include "pzvisualmatrix.h"
#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "pzcheckgeom.h"

#include "TPZMHMeshControl.h"
#include "TPZMHMixedMeshControl.h"
#include "TPZMHMixedMeshChannelControl.h"

#include "TPZMHMixedHybridMeshControl.h"
#include "TPZHybridizeHDiv.h"
#include "meshgen.h"
//#include "ConfigCasesMaze.h"
#include <iostream>
#include <string>
//#include <opencv2/opencv.hpp>
#include <math.h>
#include <set>
#include "pzsolve.h"

#include "TPZPersistenceManager.h"

#include "TPZMHMHDivErrorEstimator.h"
#include "Tools.h"

#define new_identifier


using namespace std;
//using namespace cv;

// Creating the computational H1 mesh
TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int p_order, ProblemConfig &Conf);

TPZCompMesh *CMesh_H1(struct ProblemConfig &problem);

// Creating the computational flux mesh
TPZCompMesh *CMeshFlux(TPZGeoMesh * gmesh,int pOrder);

// Creating the computational pressure mesh
TPZCompMesh *CMeshPressure(TPZGeoMesh * gmesh, int pOrder,ProblemConfig &Conf);

// Creating the computational multphysics mesh
TPZCompMesh *CMeshMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec,ProblemConfig &Conf);

// Read a mesh from a png file. The size of the domain will be npix_x by npix_y (read from the image) . Return l=nx; h=ny.
// on entry image_size_x and image_size_y are equal to the number of MHM elements in x and y
TPZGeoMesh *GeoMeshFromPng(string name, int &image_size_x, int &image_size_y);

// Create a geometric mesh with the given parameters, nx and ny are the coarse elements number. The total number of elements are defined by the image read.
TPZGeoMesh *GenerateGeoMesh(string name, int nx, int ny);

// Create a geoElSide map to be open because it has flux.
// cmesh is the flux mesh obtained from a Hdiv previus simulation.
std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>> IdentifyChanel (TPZCompMesh *cmesh);


// Compute the geometric mesh coarse indices
void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices);

// Insert the necessary objects material in the computational mesh
void InsertMaterialObjects(TPZMHMixedMeshControl &control,const ProblemConfig &config);

// Solve the H1 problem with "Conf" configuration
// Conf contains the maze information and the problem boundary conditions
int H1Test(ProblemConfig &Conf);

// Solve the mixed problem with "Conf" configuration
// Conf contains the maze information and the problem boundary conditions
TPZCompMesh* MixedTest(ProblemConfig &Conf);

// Solve the maze using MHM. By default (2x2 coarse elements)
// Conf contains the maze information and the problem boundary conditions
int MHMTest(ProblemConfig &Conf);

TPZGeoMesh *CreateCircleGeoMesh();
TPZGeoMesh *NonConvexMesh(int ndiv);

TPZGeoMesh *CreateLMHMMesh(int nDiv, TPZVec<int64_t>& coarseIndexes);


bool IsgmeshReader = false;


int main() {
    InitializePZLOG();

    
//    int ndiv = 0;
//    TPZGeoMesh *nonconvexMesh = NonConvexMesh(ndiv);
//
//
//#ifdef PZDEBUG
//    {
//        std::ofstream out("NonConvexMesh.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(nonconvexMesh, out);
//        std::ofstream out2("NonConvexMesh.txt");
//        nonconvexMesh->Print(out2);
//
//    }
//#endif
//
//
//
//    return EXIT_SUCCESS;
    
    
    for(int ndiv=1 ; ndiv<2 ; ndiv++) {
    ProblemConfig config;
    
    config.porder = 1;
    config.hdivmais = 1;
    config.ndivisions = ndiv;
    config.dimension = 2;
    config.prefine=false;
   
    TLaplaceExample1 example;

        config.exact.fExact = example.ESinSin;//example.EBubble;
     
    config.problemname = "3DProblem";
    
    config.dir_name= "MHM_Mixed3D";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());
        
    TPZGeoMesh *gmesh = nullptr;

    if(IsgmeshReader){

       // std::string meshfilename = "../Circular.msh";
       // std::string meshfilename = "../Quad.msh";
        std::string meshfilename = "../Cube.msh";

        TPZGmshReader gmsh;


     //   gmsh.GetDimNamePhysical()[1]["boundary"] = 2;
        gmsh.GetDimNamePhysical()[2]["dirichlet"] = 2;
//        gmsh.GetDimNamePhysical()[1]["boundary"] = 2;
//        gmsh.GetDimNamePhysical()[1]["neumann"] = 10;
        gmsh.GetDimNamePhysical()[3]["domain"] = 1;

 
        config.materialids.insert(1);
        config.bcmaterialids.insert(2);
        //config.bcmaterialids.insert(3);
//        config.bcmaterialids.insert(10);


        gmsh.SetFormatVersion("4.1");
        gmesh = gmsh.GeometricGmshMesh(meshfilename);
        gmsh.PrintPartitionSummary(std::cout);
        gmesh->SetDimension(config.dimension);
        config.gmesh = gmesh;

    }
        
    else{

        TPZManVector<int,4> bcids(4,-1);


        TPZVec<int64_t> coarseindices;

        TPZManVector<REAL,3> x0(3,0.), x1(3,1.);
        x1[2] = 0.;


        int nx = 2;//pow(2, ndiv);
        gmesh = CreateLMHMMesh(ndiv,coarseindices);//CreateGeoMesh(nx, bcids);//MalhaGeomFredQuadrada(nx, nx,x0, x1, coarseindices, 1);// CreateCircleGeoMesh();//CreateGeoMesh(nx, bcids);
        config.gmesh = gmesh;
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
//        config.bcmaterialids.insert(-2);
//        config.bcmaterialids.insert(-3);
//        config.bcmaterialids.insert(-4);
//        config.bcmaterialids.insert(2);
//        config.bcmaterialids.insert(3);
        gmesh->SetDimension(config.dimension);

    }
        
#ifdef PZDEBUG
        {
            std::ofstream out("GmeshMHM.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
            std::ofstream out2("gmeshInitial.txt");
            gmesh->Print(out2);
            
        }
#endif
        
        UniformRefinement(ndiv, gmesh);
        DivideLowerDimensionalElements(gmesh);
        
        TPZGeoMesh *gmesh2 = new TPZGeoMesh();
        *gmesh2 = *gmesh;

    MHMTest(config);
        
  
        
    }
    
   // return 0;
    
    
}


TPZCompMesh* MixedTest(ProblemConfig &Conf){
    
    TPZGeoMesh *gmesh = Conf.gmesh;

    
    TPZCompMesh *cmesh_flux = CreateFluxHDivMesh(Conf);//CMeshFlux(gmesh,flux_order);
    TPZCompMesh *cmesh_presure = CreatePressureMesh(Conf);//CMeshPressure(gmesh,p_order,Conf);
    
    
    TPZVec<TPZCompMesh *> fmeshvec(2);
    fmeshvec[0]=cmesh_flux;
    fmeshvec[1]=cmesh_presure;
    gmesh->ResetReference();
    
    TPZCompMesh *MixedMesh = CMeshMultphysics(gmesh,fmeshvec,Conf);
    
    {
    std::ofstream file("MixedCMesh.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(MixedMesh, file);
    
    std::ofstream out("CMeshFlux.txt");
    cmesh_flux->Print(out);
    
    }
    
    std::cout << "number of equations = " << MixedMesh->NEquations() << std::endl;
    
    
    
    //Solving the system:
  
    MixedMesh->InitializeBlock();
    
    bool must_opt_band_width_Q = true;
    int number_threads = 4;

    
        TPZCompMesh * cmesh_m_Hybrid;
        TPZManVector<TPZCompMesh*, 3> meshvector_Hybrid(3);
        TPZHybridizeHDiv hybridizer;
        tie(cmesh_m_Hybrid, meshvector_Hybrid) = hybridizer.Hybridize(MixedMesh, fmeshvec, true, -1.);
        cmesh_m_Hybrid->InitializeBlock();


        TPZAnalysis *an = new TPZAnalysis(cmesh_m_Hybrid,must_opt_band_width_Q);
    
        TPZSymetricSpStructMatrix sparse_matrix(cmesh_m_Hybrid);
            
            
        TPZStepSolver<STATE> step;
        sparse_matrix.SetNumThreads(number_threads);
        step.SetDirect(ELDLt);
        an->SetStructuralMatrix(sparse_matrix);
        an->SetSolver(step);
        an->Assemble();
        an->Solve();

    
    

    
    
    
    //POS
    TPZManVector<std::string,10> scalnames(2), vecnames(1);
    vecnames[0]  = "Flux";
    scalnames[0] = "Pressure";
    scalnames[1] = "Permeability";
    
    
    const int dim = an->Mesh()->Dimension();
    int div = 0;
    
    
    std::stringstream sout;
    sout << Conf.dir_name << "/" <<  "HybridMixed_"<<Conf.porder<<"Nref_"<<Conf.ndivisions<<".vtk";

    an->DefineGraphMesh(dim,scalnames,vecnames,sout.str());
    an->PostProcess(div,dim);
    std::cout << "Standard post-processing finished." << std::endl;
    std::cout << "Post processing errors " << std::endl;
    
    
    if(Conf.exact.Exact())
    {
        TPZManVector<REAL> errors(4,0.);
        an->SetThreadsForError(2);
        an->SetExact(Conf.exact.ExactSolution());
        an->PostProcessError(errors,false);
        
        //Erro
        
        int nequationHybrid = cmesh_m_Hybrid->NEquations();
        int nequationMixes = cmesh_presure->NEquations() +cmesh_flux->NEquations();
        
        ofstream myfile;
        myfile.open("ArquivosErrosMixedHybrid.txt", ios::app);
        myfile << "\n\n Error for Hybrid Mixed formulation " ;
        myfile << "\n-------------------------------------------------- \n";
        myfile << "Ndiv = " << Conf.ndivisions << " Order k = " << Conf.porder << " Order n = " << Conf.hdivmais << "\n";
        myfile << "DOF Total = " << nequationMixes << "\n";
        myfile << "Energy norm (flux) = " << errors[0] << "\n";//norma energia
        myfile << "error norm L2 (pressure)= " << errors[1] << "\n";//norma L2
        myfile << "Semi norm H1 = " << errors[2] << "\n";//norma L2
        myfile.close();
        
    }
    
    
   // delete cmesh_flux;
    delete  cmesh_m_Hybrid;
    
    return 0;
}


int H1Test(ProblemConfig &Conf)
{

    TPZGeoMesh *gmesh = Conf.gmesh;
    {
#ifdef PZDEBUG
        std::ofstream file("mazeh1.txt");
        gmesh->Print(file);
        
        std::ofstream out("mazeh1.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
        
#endif
    }
    
    //Creando a malla computacional
    int p_order = Conf.porder ;
    int number_threads = 4;
    bool must_opt_band_width_Q = true;
    TPZCompMesh *cmesh = CMesh_H1(Conf);
    TPZAnalysis *an = new TPZAnalysis(cmesh,must_opt_band_width_Q);
    
#ifdef USING_MKL
    TPZSymetricSpStructMatrix struct_mat(cmesh);
    struct_mat.SetNumThreads(number_threads);
    an->SetStructuralMatrix(struct_mat);
#else
    TPZSkylineStructMatrix struct_mat(cmesh);
    struct_mat.SetNumThreads(number_threads);
    an->SetStructuralMatrix(struct_mat);
#endif
    // TPZStepSolver<STATE> step;
    //step.SetDirect(ECholesky);
    
    
    
    int neq = cmesh->NEquations();
    TPZFMatrix<STATE>  cheia(neq,neq, 1.);
    TPZStepSolver<STATE> step(&cheia);
    TPZStepSolver<STATE> precond(step);
    TPZStepSolver<STATE> iterSolver;
    iterSolver.SetCG(100, precond, 0.0000000000001, 0);
    
    an->SetSolver(iterSolver);
    an->Run();
    
    an->Assemble();
    // Solving the LS
    an->Solve();
  //  an->Solve();
    // post-processing step
    {
        const int dim = an->Mesh()->Dimension();
        int div = 0;
        TPZStack<std::string> scalnames, vecnames;
        
        scalnames.push_back("Solution");
        scalnames.push_back("ExactSolution");
        vecnames.push_back("MinusKGradU");
        
        std::stringstream sout;
        sout << Conf.dir_name << "/" <<  "H1Formulation_"<<Conf.porder<<"Nref_"<<Conf.ndivisions<<".vtk";
        
        an->DefineGraphMesh(dim,scalnames,vecnames,sout.str());
        an->PostProcess(div,dim);
        std::cout << "Standard post-processing finished." << std::endl;
    }
    
    return 0;
}



TPZCompMesh *CMeshFlux(TPZGeoMesh * gmesh,int pOrder){
    
    int impervious_mat = 1;
    int permeable_mat = 2;
    int dim = gmesh->Dimension();
 
    TPZVec<REAL> convdir(dim,0.0);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    
    //Definition of the approximation space:
    
    TPZVecL2 *mat_0 = new TPZVecL2(impervious_mat);
    TPZVecL2 *mat_1 = new TPZVecL2(permeable_mat);
    
    //  inserting volumetric materials objects
    cmesh->InsertMaterialObject(mat_0);
    cmesh->InsertMaterialObject(mat_1);
    
    cmesh->SetAllCreateFunctionsHDiv(); //Creating H(div) functions
    
    
    int type_D = 0;
    int type_N = 1;
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    
    // Insert boundary conditions
    //Neumann boundary conditions (flux = 0)
    int right_bc_id = -2;
    TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_N, val1, val2);
    cmesh->InsertMaterialObject(right_bc);
    
    int left_bc_id = -4;
    TPZMaterial * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_N, val1, val2);
    cmesh->InsertMaterialObject(left_bc);
    
    int bottom_bc_1id = -1;
    TPZMaterial * bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_N, val1, val2);
    cmesh->InsertMaterialObject(bottom_bc_1);
    
    int top_bc_1id = -3;
    TPZMaterial * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_N, val1, val2);
    cmesh->InsertMaterialObject(top_bc_1);
    
    
    //Dirichlet Conditions (p=1 in, p=0 out)
    int bottom_bc_id = -5;
    TPZMaterial * bottom_bc = mat_0->CreateBC(mat_0, bottom_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(bottom_bc);
    
    int top_bc_id = -6;
    TPZMaterial * top_bc = mat_0->CreateBC(mat_0, top_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(top_bc);
    
    cmesh->SetName("LaberintoTest");
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream file("cmesh_flux.txt");
    cmesh->Print(file);
#endif
    
    return cmesh;
    
}
TPZCompMesh *CMeshPressure(TPZGeoMesh * gmesh, int pOrder,ProblemConfig &Conf){
    
    int impervious_mat = 1;
    int permeable_mat = 2;
    int dim = gmesh->Dimension();
    
    REAL perm_0 = 1;//Conf.GetImperviousMatPermeability();
    REAL perm_1 = 1;//Conf.GetPermeableMatPermeability();
    
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    
    TPZMatPoisson3d *mat_0 = new TPZMatPoisson3d(impervious_mat,dim);
    TPZMatPoisson3d *mat_1 = new TPZMatPoisson3d(permeable_mat,dim);
    mat_0->SetParameters(perm_0, conv, convdir);
    mat_1->SetParameters(perm_1, conv, convdir);
    
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
    
#ifdef PZDEBUG
    std::ofstream out("cmeshPress.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
    
}

TPZCompMesh *CMeshMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, ProblemConfig &Conf ){
    
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(Conf.gmesh);
    TPZMaterial *mat = NULL;
    TPZFMatrix<REAL> K(3,3,0),invK(3,3,0);
    K.Identity();
    invK.Identity();
    
    
    for (auto matid : Conf.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        mix->SetForcingFunction(Conf.exact.ForcingFunction());
        mix->SetForcingFunctionExact(Conf.exact.Exact());
        mix->SetPermeabilityTensor(K, invK);
        
        if (!mat) mat = mix;
        
        cmesh->InsertMaterialObject(mix);
    }
    for (auto matid : Conf.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(Conf.exact.Exact());
        cmesh->InsertMaterialObject(bc);
    }
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

    TPZManVector<int> active(2,1);
    TPZManVector<TPZCompMesh *> meshvector(2,0);
    
    meshvector[0] = CreateFluxHDivMesh(Conf);
    meshvector[1] = CreatePressureMesh(Conf);
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    bool keepmatrix = false;
    bool keeponelagrangian = true;
    TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian, keepmatrix);
    
    return cmesh;
  
    /*TPZManVector<TPZCompMesh *,2> MeshesHDiv(2);
    TPZMultiphysicsCompMesh * mixed_cmesh = CreateHDivMesh(Conf);//Cria uma malha mista HdivxL2
    MeshesHDiv = mixed_cmesh->MeshVector();
    
    TPZMultiphysicsCompMesh *mphysicCompMesh = new TPZMultiphysicsCompMesh(Conf.gmesh);
    std::ofstream outgeo("geometria.txt");
    mphysicCompMesh->Reference()->Print(outgeo);
    
    
    //Have to include the materials. Here we just did a copy of previous materials
    TPZCompMesh * cmesh =  dynamic_cast<TPZCompMesh *>(mphysicCompMesh);
    mixed_cmesh->CopyMaterials(*cmesh);
    
    TPZManVector<TPZCompMesh *,3> mp_meshes_vec(3);
    mp_meshes_vec[0] = mixed_cmesh;//HdivxL2
    mp_meshes_vec[1] = MeshesHDiv[0];//Hdiv
    mp_meshes_vec[2] = MeshesHDiv[1];//L2
    
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
    
#ifdef PZDEBUG
    std::ofstream file("cmesh_mphysics.txt");
    mphysicCompMesh->Print(file);
#endif
    
    return mphysicCompMesh;*/
    
}

int MHMTest(ProblemConfig &Conf){
    
    TRunConfig Configuration;
    
    Configuration.pOrderInternal = Conf.porder;
    Configuration.pOrderSkeleton = Conf.porder;
    Configuration.numHDivisions = Conf.ndivisions;
    Configuration.hdivmaismais = Conf.hdivmais;
    
    // number of coarse elements in the x and y direction
    //TPZGeoMesh *gmeshcoarse = Conf.gmesh;
    TPZVec<int64_t> coarseindices;
    TPZGeoMesh *gmeshcoarse = CreateLMHMMesh(2, coarseindices);
    std::ofstream file("FineMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmeshcoarse, file);

    
    TPZAutoPointer<TPZMHMixedMeshControl> MHMixed;
    
    {
        TPZAutoPointer<TPZGeoMesh> gmeshauto = new TPZGeoMesh(*gmeshcoarse);
        {
            std::ofstream out("gmeshauto.txt");
            gmeshauto->Print(out);
        }
        TPZMHMixedMeshControl *mhm = new TPZMHMixedMeshControl(gmeshauto);
        // compute for each element the coarse index to which it will belong
        //TPZVec<int64_t> coarseindices;
        //ComputeCoarseIndices(gmeshauto.operator->(), coarseindices);
        
        for(int i =0; i < coarseindices.size(); i++){
                    std::cout << "coarse index i = " << coarseindices[i] << std::endl;
            }
        
        
        // criam-se apenas elementos geometricos
        //mhm->DefinePartitionbyCoarseIndices(coarseindices);
        mhm->DefinePartition(coarseindices);
        //        MHMMixedPref << "MHMixed";
        MHMixed = mhm;
        
        // indicate the boundary material indices to  the MHM control structure
        TPZMHMixedMeshControl &meshcontrol = *mhm;
        std::set<int> matids;
       for (auto matid : Conf.materialids) {
           matids.insert(matid);
           
       }
        mhm->fMaterialIds = matids;
        matids.clear();
        
        for (auto matid : Conf.bcmaterialids) {
            matids.insert(matid);
       
        }
        mhm->fMaterialBCIds = matids;
        
        // insert the material objects in the multiphysics mesh
        InsertMaterialObjects(*mhm,Conf);
        
#ifdef PZDEBUG
        if(1)
        {
            std::ofstream out("MixedMeshControlHDiv.txt");
            meshcontrol.Print(out);
        }
#endif
        meshcontrol.SetInternalPOrder(Configuration.pOrderInternal);
        meshcontrol.SetSkeletonPOrder(Configuration.pOrderSkeleton);
       // meshcontrol.SetHdivmaismais(Configuration.hdivmaismais);
        
        meshcontrol.DivideSkeletonElements(Configuration.numDivSkeleton);
        meshcontrol.DivideBoundarySkeletonElements();

        bool substructure = true;
        meshcontrol.BuildComputationalMesh(substructure);
        
        
        
        
#ifdef PZDEBUG
        {
            std::ofstream file("GMeshControlHDiv.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(meshcontrol.GMesh().operator->(), file);
            std::ofstream out("MixedMeshControlHDiv.txt");
            meshcontrol.Print(out);
        }
#endif
        
        
        std::cout << "MHM Hdiv Computational meshes created\n";
#ifdef PZDEBUG

        {
            std::ofstream gfile("geometryMHMHdiv.txt");
            gmeshauto->Print(gfile);
            std::ofstream out_mhm("MHM_hdiv.txt");
            meshcontrol.CMesh()->Print(out_mhm);
            
        }
#endif
        
        std::cout << "Number of equations MHMixed " << MHMixed->CMesh()->NEquations() << std::endl;
        
        
    }
    
    
    std::string prefix;
    std::cout<<"Solving MHM problem"<<std::endl;
    
    
    SolveProblem(MHMixed->CMesh(), MHMixed->GetMeshes(),Conf.exact, prefix, Configuration);
    
 
    
    std::cout<<"Error Estimation processing for MHM-Hdiv problem "<<std::endl;
    
    // Error estimation
    TPZMultiphysicsCompMesh *InputMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(MHMixed->CMesh().operator->());
    if(!InputMesh) DebugStop();
    
    TPZMHMHDivErrorEstimator ErrorEstimator(*InputMesh, MHMixed.operator->());
    ErrorEstimator.fOriginalIsHybridized = false;
    ErrorEstimator.SetAnalyticSolution(Conf.exact);
    
    ErrorEstimator.CreatePressureSkeleton();
    ErrorEstimator.PotentialReconstruction();
    
    {
        ErrorEstimator.fProblemConfig.problemname = Conf.problemname;
        ErrorEstimator.fProblemConfig.dir_name = Conf.dir_name;
        std::string command = "mkdir " + Conf.dir_name;
        system(command.c_str());
        
        ErrorEstimator.fProblemConfig.porder =Conf.porder;
        ErrorEstimator.fProblemConfig.ndivisions = Conf.ndivisions;
        ErrorEstimator.fProblemConfig.hdivmais = Conf.hdivmais;
        ErrorEstimator.fProblemConfig.adaptivityStep = Conf.adaptivityStep;
        
    TPZManVector<REAL> errors;

    ErrorEstimator.ComputeErrors(errors);
//    hAdaptivity(&ErrorEstimator.fPostProcMesh,gmeshcoarse);
        
    }
    
    return 0;
}

void InsertMaterialObjects(TPZMHMixedMeshControl &control,const ProblemConfig &config)
{
    TPZCompMesh &cmesh = control.CMesh();
    
    TPZGeoMesh &gmesh = control.GMesh();
   
    TPZFMatrix<STATE> val1(1,1,0.), val2Flux(1,1,0.), val2Pressure(1,1,1.);
    
    
    int dim = gmesh.Dimension();
    cmesh.SetDimModel(dim);
    
    TPZCompMesh *MixedFluxPressureCmesh = &cmesh;
    
    // Material meio poroso
    TPZMixedPoisson * mat = new TPZMixedPoisson(1,dim);
//    mat->SetSymmetric();
//    mat->SetPermeability(1.);
    
    TPZFMatrix<REAL> K(3,3,0),invK(3,3,0);
//    K.Identity();
//    invK.Identity();
    for(int i=0;i<3;i++){
        for(int j=0; j<3;j++){
            if (i != j){
                K(i,j)=1.;
                invK(i,j) = -1./4.;
            }
            else{
                K(i,j) = 2.;
                invK(i,i) = 3./4.;
                
            }
        }
    }
    
  //  K.Print(std::cout);
  //  invK.Print(std::cout);
    
    mat->SetForcingFunctionExact(config.exact.Exact());
    mat->SetForcingFunction(config.exact.ForcingFunction());
    mat->SetPermeabilityTensor(K, invK);
    
    
    MixedFluxPressureCmesh->InsertMaterialObject(mat);
    
    for (auto matid : config.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        val2.Zero();
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
       // bc->TPZMaterial::SetForcingFunction(config.exact.Exact());
        MixedFluxPressureCmesh->InsertMaterialObject(bc);
    }
    
}

// compute the coarse indices of the geometric mesh
void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices)
{
    //    {
    //        std::ofstream out("gmeshref.txt");
    //        gmesh->Print(out);
    //    }
     coarseindices.Resize(gmesh->NElements());
    int count = 0;
    for (int64_t el=0; el<gmesh->NElements(); el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel || gel->Dimension() != gmesh->Dimension()) continue;
        if(gel->Father()) continue;
        coarseindices[count] = el;
        count++;
    }
    coarseindices.Resize(count);
}


std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>> IdentifyChanel (TPZCompMesh *cmesh){
    
    int nelements = cmesh->NElements();
    TPZGeoElSide first_gelside;
    for(int iel=0; iel<=nelements; iel++){
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel){continue;}
        TPZGeoEl *gel = cel->Reference();
        TPZGeoElSide gelside(gel,4);
        TPZGeoEl *neig = gelside.Neighbour().Element();
        if(neig->MaterialId()==-5){
            first_gelside.SetElement(gel);
            first_gelside.SetSide(4);
            break;
        }
    }
    
    int count=0;
    int count_chain=0;
    double Flux_Max = 0.0;
    std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>> chain;
    bool exit=false;
    
    while(exit==false){
        
        TPZGeoElSide exit_test = first_gelside;
        exit_test.SetSide(6);
        TPZGeoElSide candidate_exist =exit_test.Neighbour();
        if(candidate_exist.Element()->MaterialId()==-6){
            std::cout<<"Cadena encontrada con exito";
            exit = true;
            break;
        }
        
        
        //    while(candidate.Element()->Dimension()!=2){
        //
        //            candidate=candidate.Neighbour();
        //
        //    }
        int side_in = first_gelside.Side();
        for(int ican=4; ican<8; ican++){
            first_gelside.SetSide(ican);
            TPZGeoElSide candidate = first_gelside.Neighbour();
            while(candidate.Element()->Dimension()!=2){
                
                candidate=candidate.Neighbour();
            }
            if(side_in == ican ){continue;}
            
            //calcula los 3 geoelement side y el fluxo
            TPZVec<REAL> qsi(2);
            qsi[0]=0;
            qsi[1]=0;
            int var=31;
            TPZVec<STATE> sol;
            TPZCompEl *cel = candidate.Element()->Reference();
            cel->Solution(qsi, var,sol);
            double Flux_can_Mag = sqrt(sol[0]*sol[0] + sol[1]*sol[1]);
            
            if(Flux_can_Mag > Flux_Max){
                Flux_Max =Flux_can_Mag;
                first_gelside.SetSide(ican);
                chain[count].first= first_gelside;
                chain[count].second =candidate;
            }
        }
        Flux_Max=0.0;
        first_gelside =chain[count].second;
        count++;
    }
    
    //here
    
    cmesh->LoadReferences();
    TPZGeoMesh *gmesh =cmesh->Reference();
    for(auto it:gmesh->ElementVec()){
        int father_index = it->FatherIndex();
        int element = it->Index();
        std::cout<<"Element: "<<element<<" father: "<<father_index<<std::endl;
        
    }
    
    
    int n_el_chain = chain.size();
    std::map<int,pair<TPZGeoElSide, TPZGeoElSide>> skelchanel;
    int count_skel_chanel=0;
    int matId_skel_chanel = 10;
    count =0;
    for(auto it : chain){
        TPZGeoEl *first_element = it.second.first.Element();
        TPZGeoEl *second_element = it.second.second.Element();
        
        int first_father_index = first_element->LowestFather()->Index();
        int second_father_index = second_element->LowestFather()->Index();
        
        if(first_father_index!=second_father_index){
            TPZGeoElBC(it.second.first, matId_skel_chanel);
            skelchanel[count_skel_chanel] = it.second;
            count_skel_chanel++;
        }
        
    }
    std::ofstream out("TestMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return skelchanel;
    
}


TPZCompMesh *CMesh_H1(struct ProblemConfig &problem) {
    
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
        if(matid==10){
            int bctype = 1;
            TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
             bc->TPZMaterial::SetForcingFunction(problem.exact.Exact());
             cmesh->InsertMaterialObject(bc);
        }
        else{
            TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
            bc->TPZMaterial::SetForcingFunction(problem.exact.Exact());
             cmesh->InsertMaterialObject(bc);
        }
       
    }
    
    cmesh->SetDefaultOrder(problem.porder);//ordem
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    
    cmesh->AutoBuild();
    
    
    return cmesh;
}


TPZGeoMesh *CreateCircleGeoMesh() {
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    
    TPZVec<REAL> coord(3, 0.);
    
    // Inserts node at origin
    gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[0].Initialize(coord, *gmesh);
    
    // Inserts circumference nodes
    for (int64_t i = 0; i < 16; i++) {
        const REAL step = M_PI / 8;
        coord[0] = cos(i * step);
        coord[1] = sin(i * step);
        std::cout<<"{"<<coord[0]<<"," <<coord[1]<<"},"<<std::endl;
        
    
        
        const int64_t newID = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[newID].Initialize(coord, *gmesh);
    }
    
    int matIdTriangle = 1, matIdArc = 2;
    
    // Inserts triangle elements
    TPZManVector<int64_t, 3> nodesIdVec(3);
    for (int64_t i = 0; i < 7; i++) {
        nodesIdVec[0] = 0;
        nodesIdVec[1] = 1 + 2 * i;
        nodesIdVec[2] = 3 + 2 * i;
        std::cout<<"{"<<nodesIdVec[0] <<"," <<nodesIdVec[1] <<", "<<nodesIdVec[2] <<"},"<<std::endl;
        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle>>(nodesIdVec, matIdTriangle, *gmesh);
    }
    nodesIdVec[0] = 0;
    nodesIdVec[1] = 15;
    nodesIdVec[2] = 1;

    new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle>>(nodesIdVec, matIdTriangle, *gmesh);
    // Inserts arc elements
    for (int64_t i = 0; i < 7; i++) {
        nodesIdVec[0] = 1 + 2 * i;
        nodesIdVec[1] = 3 + 2 * i;
        nodesIdVec[2] = 2 + 2 * i;
        new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, matIdArc, *gmesh);
    }
    
    nodesIdVec[0] = 15;
    nodesIdVec[1] = 1;
    nodesIdVec[2] = 16;
    //para o broblema do douglas matIdArc=3
    new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, 3, *gmesh);
//    // Finally, inserts line elements to complete boundary
//    nodesIdVec.Resize(2);
//    nodesIdVec[0] = 0;
//    nodesIdVec[1] = 1;
//    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, matIdArc, *gmesh);
//
//    nodesIdVec[0] = 0;
//    nodesIdVec[1] = 14;
//    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, matIdArc, *gmesh);

    gmesh->BuildConnectivity();
    
    return gmesh;
}

TPZGeoMesh *NonConvexMesh(int ndiv) {
    
    int nx = 6*(ndiv+1);
    int ny = 4*(ndiv+1);
    REAL hx = 1./nx;
    REAL hy = 1./ny;
   
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    
    TPZVec<REAL> coord(3, 0.);
    //1o nivel de nos
    std::cout<<"Ptos={";
    for (int64_t ix = 0; ix <= nx; ix++) {
        
        for (int64_t iy = 0; iy <= ny; iy++) {
            
            coord[0] = hx*ix;
            coord[1] = hy*iy;
            const int64_t newID = gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[newID].Initialize(coord, *gmesh);
            
            std::cout<<"{"<<coord[0]<<","<<coord[1]<<"},"<<std::endl;

        }
    }
    
    std::cout<<"};ListPlot[Ptos, PlotStyle -> {Red}]"<<std::endl;
    
    //2o nivel de nos
    
    std::cout<<"Ptos2={";
    
    for (int64_t ix = 1; ix < nx; ix++) {
        
        for (int64_t iy = 1; iy < ny; iy++) {
            
            coord[0] = (hx/2.)*ix ;
            coord[1] = (hx/2.)*iy ;
            std::cout<<"{"<<coord[0]<<","<<coord[1]<<"},"<<std::endl;
            const int64_t newID = gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[newID].Initialize(coord, *gmesh);
            
        }
    }
    
    std::cout<<"};ListPlot[Ptos, PlotStyle -> {Blue}]"<<std::endl;
    gmesh->BuildConnectivity();
    return gmesh;
}
    

TPZGeoMesh *CreateLMHMMesh(int nDiv, TPZVec<int64_t>& coarseIndexes) {

    int factor = (int)(pow(2, nDiv) + 0.5);

    int xElements = factor * 6;
    int yElements = factor * 4;

    TPZManVector<int> nx(2);
    nx[0] = xElements;
    nx[1] = yElements;

    TPZManVector<REAL> x0(3,0.), x1(3,1.);
    x1[2] = 0.;

    TPZGenGrid gen(nx,x0,x1);

    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gen.Read(gmesh);
    gen.SetBC(gmesh, 4, -1);
    gen.SetBC(gmesh, 5, -1);
    gen.SetBC(gmesh, 6, -1);
    gen.SetBC(gmesh, 7, -1);

    gmesh->SetDimension(2);
    gmesh->BuildConnectivity();

    // Assigns matIDs to create L elements
    {
        int coarseIndex = 0;
        int64_t nelem = gmesh->NElements();
        coarseIndexes.Resize(nelem, -1);
        for (int64_t elem = 0; elem < nelem; elem++) {
            TPZGeoEl *gel = gmesh->ElementVec()[elem];
            if (gel->Dimension() != 2) continue;

            int lineInPattern = elem / nx[0] % 4;
            int colInPattern = elem % nx[0] % 6;

            // IDs of elements in the neighbourhood to which a coarse index has been already assigned
            int leftEl = elem - 1;
            int bottomEl = elem - nx[0];

            if (lineInPattern == 0) {
                if (colInPattern % 2 == 0) {
                    coarseIndexes[elem] = coarseIndex;
                    coarseIndex++;
                } else {
                    coarseIndexes[elem] = coarseIndexes[leftEl];
                }
            } else if (lineInPattern == 1) {
                if (colInPattern == 0 || colInPattern == 2 || colInPattern == 5) {
                    coarseIndexes[elem] = coarseIndexes[bottomEl];
                } else if (colInPattern == 4) {
                    coarseIndexes[elem] = coarseIndexes[leftEl];
                } else {
                    coarseIndexes[elem] = coarseIndex;
                    coarseIndex++;
                }
            } else if (lineInPattern == 2) {
                if (colInPattern == 1 || colInPattern == 4) {
                    coarseIndexes[elem] = coarseIndexes[bottomEl];
                } else if (colInPattern == 2) {
                    coarseIndexes[elem] = coarseIndexes[leftEl];
                } else {
                    coarseIndexes[elem] = coarseIndex;
                    coarseIndex++;
                }
            } else if (lineInPattern == 3) {
                if (colInPattern == 0) {
                    coarseIndexes[elem] = coarseIndexes[bottomEl];
                } else if (colInPattern % 2 == 1) {
                    coarseIndexes[elem] = coarseIndexes[leftEl];
                } else {
                    coarseIndexes[elem] = coarseIndexes[leftEl] + 1;
                }
            }
        }

        for(int64_t elem = 0; elem < nelem; elem++) {
            gmesh->ElementVec()[elem]->SetMaterialId(coarseIndexes[elem]);
        }
    }
    return gmesh;
}
