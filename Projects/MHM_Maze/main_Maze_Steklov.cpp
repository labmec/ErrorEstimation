#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzlog.h"
#include "pzreal.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoelside.h"
#include "tpzgeoelrefpattern.h"
#include "tpzautopointer.h"
#include "TPZLinearAnalysis.h"

#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZStructMatrixT.h"
#include "TPZElementMatrixT.h"

#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZVTKGeoMesh.h"
#include "TPZGenGrid2D.h"
#include "TPZMHMixedMeshControl.h"
#include "TPZMHMixedMeshChannelControl.h"
#include "TPZHybridizeHDiv.h"
#include "ConfigCasesMaze.h"
#include <ToolsMHM.h>
#include <cmath>
#include <opencv2/opencv.hpp>
#include <set>
#include <string>

#include "TPZPersistenceManager.h"

using namespace std;
using namespace cv;

// Creating the computational flux mesh
TPZCompMesh *CMeshFlux(TPZGeoMesh * gmesh,int pOrder);

// Creating the computational pressure mesh
TPZCompMesh *CMeshPressure(TPZGeoMesh * gmesh, int pOrder,ConfigCasesMaze Conf);

// Creating the computational multphysics mesh
TPZCompMesh *CMeshMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec,ConfigCasesMaze &Conf);

// Read a mesh from a png file. The size of the domain will be npix_x by npix_y (read from the image) . Return l=nx; h=ny.
TPZGeoMesh *GeoMeshFromPng(string name, double &l, double &h);

// Create a geometric mesh with the given parameters, nx and ny are the coarse elements number. The total number of elements are defined by the image read.
TPZGeoMesh *GenerateGeoMesh(string name, int nx, int ny);

// Compute the geometric mesh coarse indices
void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices);

// Insert the necessary objects material in the computational mesh
void InsertMaterialObjects(TPZMHMixedMeshControl &control);

// Solve the mixed problem with "Conf" configuration
// Conf contains the maze information and the problem boundary conditions
TPZCompMesh* MixedTest(ConfigCasesMaze &Conf, int nx, int ny);

// Solve the maze using MHM. By default (2x2 coarse elements)
// Conf contains the maze information and the problem boundary conditions
int MHMTest(ConfigCasesMaze &Conf);

// compute the eigenvalues/eigenvectors of the Steklov problems associated with the subdomains
int SteklovTest(ConfigCasesMaze &Conf);

std::map<int,int> matextend;
int matid1BC = 8;
int matid2BC = 9;


void EstimateError(TPZMHMHDivErrorEstimator &errorEstimator, ProblemConfig &config);

void LocateElementsToAdapt(TPZMHMHDivErrorEstimator &errorEstimator, ProblemConfig &config);

int main(){
    TPZLogger::InitializePZLOG();
    
    ConfigCasesMaze ConfCasesMeze;
//    ConfCasesMeze.SetImageName("⁨../Mazes/maze128x128.png");
    ConfCasesMeze.SetImageName("../Mazes/maze8x8.png");
    ConfCasesMeze.SetImperviousMatPermeability(1);//pouco permeavel
    ConfCasesMeze.SetPermeableMatPermeability(100000);//dentro do labirinto
    ConfCasesMeze.SetFluxOrder(1);
    ConfCasesMeze.SetPressureOrder(1);
    ConfCasesMeze.SetCCPressureIn(100);//pressao na entrada
    ConfCasesMeze.SetCCPressureOut(1);//pressao na saida
    ConfCasesMeze.SetMHMOpenChannel(false);
    ConfCasesMeze.SetVTKName("maze8x8.vtk");

    SteklovTest(ConfCasesMeze);
    return 0;
}


int MHMTest(ConfigCasesMaze &Conf){

    TRunConfig Configuration;

    TPZGeoMesh *gmeshcoarse = GenerateGeoMesh(Conf.GetImageName(), 32, 32);
    {
        std::ofstream file(Conf.GetVTKName());
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshcoarse, file);
    }

    int interface_mat_id = 600;
    Conf.SetMHMOpenChannel(true);
    bool OpenChannel = Conf.GetMHMOpenChannel();

    TPZAutoPointer<TPZMHMixedMeshChannelControl> MHMixed;

    {
        TPZAutoPointer<TPZGeoMesh> gmeshauto = new TPZGeoMesh(*gmeshcoarse);
        {
            std::ofstream out("gmeshauto.txt");
            gmeshauto->Print(out);
        }
        TPZMHMixedMeshChannelControl *mhm = new TPZMHMixedMeshChannelControl(gmeshauto);
        TPZVec<int64_t> coarseindices;
        ComputeCoarseIndices(gmeshauto.operator->(), coarseindices);
        gmeshauto->AddInterfaceMaterial(1, 2, interface_mat_id);
        gmeshauto->AddInterfaceMaterial(2, 1, interface_mat_id);


        // criam-se apenas elementos geometricos
        mhm->DefinePartitionbyCoarseIndices(coarseindices);
        MHMixed = mhm;

        TPZMHMixedMeshChannelControl &meshcontrol = *mhm;
        {
            std::set<int> matids;
            matids.insert(1);
            matids.insert(2);
            mhm->fMaterialIds = matids;
            matids.clear();
            matids.insert(-1);
            matids.insert(-2);
            matids.insert(-3);
            matids.insert(-4);
            matids.insert(-5);
            matids.insert(-6);
            mhm->fMaterialBCIds = matids;
        }

        InsertMaterialObjects(*mhm);

        meshcontrol.SetInternalPOrder(1);
        meshcontrol.SetSkeletonPOrder(1);

//        meshcontrol.DivideSkeletonElements(2);
        meshcontrol.DivideBoundarySkeletonElements();

        bool substructure = true;
        std::map<int, std::pair<TPZGeoElSide, TPZGeoElSide>> test;
        if (OpenChannel) {
            TPZCompMesh *flux_temp = MixedTest(Conf,32,32);
            std::cout << "flux_temp norm of solution " << Norm(flux_temp->Solution()) << std::endl;
            test = IdentifyChanel(flux_temp);
            flux_temp->Reference()->ResetReference();
//            delete flux_temp;
        }

        meshcontrol.BuildComputationalMesh(substructure, OpenChannel, test);

#ifdef ERRORESTIMATION_DEBUG
        if (1) {
            std::ofstream file("GMeshControlHDiv.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(meshcontrol.GMesh().operator->(), file);
        }
#endif

        std::cout << "MHM Hdiv Computational meshes created\n";

        std::cout << "Number of equations MHMixed " << MHMixed->CMesh()->NEquations() << std::endl;

    }
    
    TPZCompMesh *MixedMesh = MHMixed->CMesh().operator->();
    SolveProblem(MHMixed->CMesh(), MHMixed->GetMeshes(), Conf.GetExactSolution(),  Conf.GetVTKName(), Configuration);
    ProblemConfig config;
    config.dimension = 2;
    config.exact = nullptr;
    config.problemname = "MazeHdiv128x128";
    config.dir_name = "Results128x128";
    config.porder = 3;
    config.hdivmais = 3;
    config.materialids = {1, 2};
    config.bcmaterialids = {-1, -2, -3, -4, -5, -6};
    config.makepressurecontinuous = true;
    config.ndivisions = 0;
    config.gmesh = MixedMesh->Reference();

    TPZMultiphysicsCompMesh *originalMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(MHMixed->CMesh().operator->());
    bool postProcWithHdiv = false;
//    TPZMHMHDivErrorEstimator ErrorEstimator(*originalMesh, MHMixed.operator->(), postProcWithHdiv);
//    EstimateError(ErrorEstimator, config);
    //LocateElementsToAdapt(ErrorEstimator, config);

    return 0;
}

void AnalyseSteklov(TPZSubCompMesh *sub, int count);

// add geometric elements between macro domains with material ids determined by wrapids
void AddDomainWrapElements(TPZMHMeshControl &mhm, std::map<int,int> &wrapids);

int SteklovTest(ConfigCasesMaze &Conf){

    TRunConfig Configuration;

    TPZGeoMesh *gmeshcoarse = GenerateGeoMesh(Conf.GetImageName(), 2, 2);
    {
        std::ofstream file(Conf.GetVTKName());
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshcoarse, file);
    }

    int interface_mat_id = 600;
    matextend[1] = matid1BC;
    matextend[2] = matid2BC;
    Conf.SetMHMOpenChannel(true);
    bool OpenChannel = Conf.GetMHMOpenChannel();

    TPZAutoPointer<TPZMHMixedMeshChannelControl> MHMixed;

    // configure the MHM computational mesh
    {
        TPZAutoPointer<TPZGeoMesh> gmeshauto = new TPZGeoMesh(*gmeshcoarse);
        {
            std::ofstream out("gmeshauto.txt");
            gmeshauto->Print(out);
        }
        TPZMHMixedMeshChannelControl *mhm = new TPZMHMixedMeshChannelControl(gmeshauto);
        TPZVec<int64_t> coarseindices;
        ComputeCoarseIndices(gmeshauto.operator->(), coarseindices);
        gmeshauto->AddInterfaceMaterial(1, 2, interface_mat_id);
        gmeshauto->AddInterfaceMaterial(2, 1, interface_mat_id);


        // criam-se apenas elementos geometricos
        mhm->DefinePartitionbyCoarseIndices(coarseindices);
        MHMixed = mhm;

        AddDomainWrapElements(*mhm, matextend);
        
        TPZMHMixedMeshChannelControl &meshcontrol = *mhm;
        {
            std::set<int> matids;
            matids.insert(1);
            matids.insert(2);
            mhm->fMaterialIds = matids;
            matids.clear();
            matids.insert(-1);
            matids.insert(-2);
            matids.insert(-3);
            matids.insert(-4);
            matids.insert(-5);
            matids.insert(-6);
            matids.insert(matid1BC);
            matids.insert(matid2BC);
            mhm->fMaterialBCIds = matids;
        }

        InsertMaterialObjects(*mhm);
        
        {
            TPZCompMesh &cmesh = mhm->CMesh();
            auto *mat = dynamic_cast<TPZMixedDarcyFlow *>(cmesh.FindMaterial(1));
            TPZFNMatrix<2,REAL> val1(1,1,0.);
            TPZManVector<REAL> val2(1,0.);
            auto *bnd1 = mat->CreateBC(mat, matid1BC, 0, val1, val2);
            cmesh.InsertMaterialObject(bnd1);
            auto *bnd2 = mat->CreateBC(mat, matid2BC, 0, val1, val2);
            cmesh.InsertMaterialObject(bnd2);
        }

        meshcontrol.SetInternalPOrder(1);
        meshcontrol.SetSkeletonPOrder(1);

//        meshcontrol.DivideSkeletonElements(2);
        meshcontrol.DivideBoundarySkeletonElements();

        bool substructure = true;
        std::map<int, std::pair<TPZGeoElSide, TPZGeoElSide>> test;
        if (OpenChannel) {
            TPZCompMesh *flux_temp = MixedTest(Conf,2,2);
            std::cout << "flux_temp norm of solution " << Norm(flux_temp->Solution()) << std::endl;
            test = IdentifyChanel(flux_temp);
            flux_temp->Reference()->ResetReference();
//            delete flux_temp;
        }

        meshcontrol.BuildComputationalMesh(substructure, OpenChannel, test);

#ifdef ERRORESTIMATION_DEBUG
        if (1) {
            std::ofstream file("GMeshControlHDiv.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(meshcontrol.GMesh().operator->(), file);
        }
#endif

        std::cout << "MHM Hdiv Computational meshes created\n";

        std::cout << "Number of equations MHMixed " << MHMixed->CMesh()->NEquations() << std::endl;

    }
    
    TPZCompMesh *MixedMesh = MHMixed->CMesh().operator->();
    {
        MixedMesh->ComputeNodElCon();
        std::ofstream out("EntryMesh.txt");
        MixedMesh->Print(out);
    }
    
    int64_t nelem = MixedMesh->NElements();
    int64_t count = 0;
    for (int64_t el = 0; el < nelem; el++) {
        TPZCompEl *cel = MixedMesh->Element(el);
        auto *sub = dynamic_cast<TPZSubCompMesh *>(cel);
        if(sub) {
            AnalyseSteklov(sub,count);
            count++;
        }
    }
    



    SolveProblem(MHMixed->CMesh(), MHMixed->GetMeshes(), Conf.GetExactSolution(),  Conf.GetVTKName(), Configuration);
    ProblemConfig config;
    config.dimension = 2;
    config.exact = nullptr;
    config.problemname = "MazeHdiv128x128";
    config.dir_name = "Results128x128";
    config.porder = 3;
    config.hdivmais = 3;
    config.materialids = {1, 2};
    config.bcmaterialids = {-1, -2, -3, -4, -5, -6};
    config.makepressurecontinuous = true;
    config.ndivisions = 0;
    config.gmesh = MixedMesh->Reference();

    TPZMultiphysicsCompMesh *originalMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(MHMixed->CMesh().operator->());
    bool postProcWithHdiv = false;
//    TPZMHMHDivErrorEstimator ErrorEstimator(*originalMesh, MHMixed.operator->(), postProcWithHdiv);
//    EstimateError(ErrorEstimator, config);
    //LocateElementsToAdapt(ErrorEstimator, config);

    return 0;
}

void AnalyseSteklov(TPZSubCompMesh *sub, int count){
    // Identify element/sides that belong to a different mesh
    int64_t nel = sub->NElements();
    TPZCompMesh *father = sub->Mesh();
    father->LoadReferences();
    TPZGeoMesh *gmesh = father->Reference();
    // plot the geometric elements
    {
        std::set<int64_t> elindices;
        int64_t nel = sub->NElements();
        for (int64_t el = 0; el<nel; el++) {
            auto cel = sub->Element(el);
            if(!cel) continue;
            TPZStack<TPZCompEl *> celstack;
            cel->GetCompElList(celstack);
            for(auto el : celstack) {
                auto gel = el->Reference();
                if(!gel) DebugStop();
                elindices.insert(gel->Index());
            }
        }
        std::stringstream sout;
        sout << "SubMesh_geo_" << count << ".vtk";
        std::ofstream out(sout.str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,elindices,out);
    }
    // compute the stiffness matrix
    TPZElementMatrixT<STATE> ek,ef;
    sub->CalcStiff(ek, ef);
    {
        std::stringstream sout;
        sout << "SubMesh_matrix_" << count << ".txt";
        std::ofstream out(sout.str());
        ek.Print(out);
    }
    // compute the mass matrix (how?)
    // solve the eigenvalue problem
    // plot the eigenvectors
    // print the eigenvalues
    
}

// add geometric elements between macro domains with material ids determined by wrapids
void AddDomainWrapElements(TPZMHMeshControl &mhm, std::map<int,int> &wrapids)
{
    auto gmesh = mhm.GMesh();
    auto &geoToMHM = mhm.GetGeoToMHMDomain();
    int64_t nel = gmesh->NElements();
    int meshdim = gmesh->Dimension();
    if(gmesh->NElements() != geoToMHM.size()) DebugStop();
    // Identify element/sides that belong to a different mesh
    TPZStack<std::pair<TPZGeoElSide,TPZGeoElSide>> bounds;
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        int64_t geldomain = geoToMHM[el];
        if(!gel || gel->HasSubElement()) continue;
        int geldim = gel->Dimension();
        if(geldim != meshdim) continue;
        int firstside = gel->FirstSide(geldim-1);
        int lastside = gel->NSides()-1;
        for (int side = firstside; side < lastside; side++) {
            TPZGeoElSide gelside(gel,side);
            for (auto neigh = gelside.Neighbour(); neigh != gelside; neigh++) {
                auto neighgel = neigh.Element();
                // we look for a neighbour of same dimension
                if(neighgel->Dimension() != meshdim) {
                    continue;
                }
                // if the mesh of the neighbour is different then we have a boundary element
                auto neighdomain = geoToMHM[neighgel->Index()];
                if(neighdomain != geldomain) {
                    bounds.Push({gelside,neigh});
                }
            }
        }
    }
    int64_t nbounds = bounds.size();
    geoToMHM.Resize(nel+nbounds, -1);
    for (int el = 0; el<nbounds; el++) {
        auto gelside = bounds[el].first;
        auto neighside = bounds[el].second;
        int64_t geldomain = geoToMHM[gelside.Element()->Index()];
        if(geldomain == -1) DebugStop();
        int neighmatid = neighside.Element()->MaterialId();
        if(wrapids.find(neighmatid) == wrapids.end()) DebugStop();
        int bcid = wrapids[neighmatid];
        TPZGeoElBC gbc(gelside,bcid);
        int64_t gbcindex = gbc.CreatedElement()->Index();
        geoToMHM[gbcindex] = geldomain;
    }
    nel = gmesh->NElements();
    geoToMHM.Resize(nel, -1);
}

