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
#include <libInterpolate/Interpolate.hpp>

typedef _2D::BicubicInterpolator<REAL> Interpolator;

// Global variables
Interpolator interpolator;


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
#include "TPZVTKGenerator.h"
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
TPZCompMesh *CMeshPressure(TPZGeoMesh * gmesh, int pOrder, ConfigCasesMaze Conf);

// Creating the computational multphysics mesh
TPZCompMesh *CMeshMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec,ConfigCasesMaze &Conf);

// Read a mesh from a png file. The size of the domain will be npix_x by npix_y (read from the image) . Return l=nx; h=ny.
TPZGeoMesh *GeoMeshFromPng(string name, double &l, double &h);

// Create a geometric mesh with the given parameters, nx and ny are the coarse elements number. The total number of elements are defined by the image read.
TPZGeoMesh *CreateSPE10GeoMesh();

// Compute the geometric mesh coarse indices
void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices);

// Insert the necessary objects material in the computational mesh
void InsertMaterials(TPZCompMesh *cmesh);

// Solve the mixed problem with "Conf" configuration
// Conf contains the maze information and the problem boundary conditions
TPZCompMesh* MixedTest(ConfigCasesMaze &Conf, int nx, int ny);

// Solve the maze using MHM. By default (2x2 coarse elements)
// Conf contains the maze information and the problem boundary conditions
int MHMTest(ConfigCasesMaze &Conf, std::map<int,std::pair<int64_t,int64_t>> &intersectGeoElIndex, std::map<int64_t,int> &indexToSide);

// compute the eigenvalues/eigenvectors of the Steklov problems associated with the subdomains
int SteklovTest(ConfigCasesMaze &Conf, std::map<int,std::pair<int64_t,int64_t>> &intersectGeoElIndex, std::map<int64_t,int> &indexToSide);

std::map<int,int> matextend;
int matid1BC = 8;
int matid2BC = 9;

constexpr int layer = 36;
constexpr int nx = 220;
constexpr int ny = 60;
constexpr int n_cells = nx * ny;

constexpr int n_subdomainsx = 11;
constexpr int n_subdomainsy = 3;

STATE PermeabilityFunction(const TPZVec<REAL> &x);

void ReadSPE10CellPermeabilities(TPZVec<REAL>*perm_vec, int layer);


void AssociateGeoElSides(TPZVec<std::set<TPZGeoElSide>> &eigGeoElSides, std::map<int,std::pair<int64_t,int64_t>> &intersectGeoElIndex, std::map<int64_t,int> &indexToSide);

std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>> IdentifyIntersections(TPZCompMesh *cmesh, std::map<int,std::pair<int64_t,int64_t>> &intersectGeoElIndex, std::map<int64_t,int> &indexToSide);

int main(){
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    ConfigCasesMaze ConfCasesMaze;
    // ConfCasesMaze.SetImageName("Mazes/maze128x128.png");
    ConfCasesMaze.SetImageName("Mazes/maze8x8.png");
    ConfCasesMaze.SetImperviousMatPermeability(1);//pouco permeavel
    ConfCasesMaze.SetPermeableMatPermeability(100000);//dentro do labirinto
    ConfCasesMaze.SetFluxOrder(1);
    ConfCasesMaze.SetPressureOrder(0);
    ConfCasesMaze.SetCCPressureIn(100);//pressao na entrada
    ConfCasesMaze.SetCCPressureOut(1);//pressao na saida
    ConfCasesMaze.SetMHMOpenChannel(false);
    ConfCasesMaze.SetVTKName("maze128x128.vtk");
    ConfCasesMaze.SetNumberOfSubdomains(2);
    ConfCasesMaze.SetSkeletonDivision(2);

    std::map<int,std::pair<int64_t,int64_t>> intersectGeoElIndex;
    std::map<int64_t,int> indexToSide;
    SteklovTest(ConfCasesMaze, intersectGeoElIndex, indexToSide);

    std::cout << "intersectGeoElIndex = ";
    for (const auto &it:intersectGeoElIndex)
    {
        std::cout << it.second.first << " ";
    }

    // ConfCasesMaze.SetMHMOpenChannel(true);
    // MHMTest(ConfCasesMaze, intersectGeoElIndex, indexToSide);

    return 0;
}


int MHMTest(ConfigCasesMaze &Conf, std::map<int,std::pair<int64_t,int64_t>> &intersectGeoElIndex, std::map<int64_t,int> &indexToSide){

    TRunConfig Configuration;

    auto perm_vec = TPZManVector<REAL, n_cells>(n_cells, 1);
    ReadSPE10CellPermeabilities(&perm_vec, layer);
    
    TPZGeoMesh *gmesh = CreateSPE10GeoMesh();
    {
        std::ofstream file(Conf.GetVTKName());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    }

    std::vector<REAL> x, y, perm;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            const int cell_id = ny * i + j;
            const double cell_perm = perm_vec[cell_id];
            x.push_back(0.5 + i);
            y.push_back(0.5 + j);
            perm.push_back(cell_perm);
        }
    }

    interpolator.setData(x.size(), x.data(), y.data(), perm.data());


    int interface_mat_id = 600; // it will not be used here...
    Conf.SetMHMOpenChannel(true);
    bool OpenChannel = false;// Conf.GetMHMOpenChannel();

    TPZAutoPointer<TPZMHMixedMeshChannelControl> mhm = new TPZMHMixedMeshChannelControl(gmesh);
    {
        std::ofstream out("gmeshauto.txt");
        gmesh->Print(out);
    }

    TPZVec<int64_t> coarseindices;
    ComputeCoarseIndices(gmesh, coarseindices);
    // gmesh->AddInterfaceMaterial(1, 2, interface_mat_id);
    // gmesh->AddInterfaceMaterial(2, 1, interface_mat_id);

    int nInternalRef = 3;
    Tools::UniformRefinement(nInternalRef, 2, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    // criam-se apenas elementos geometricos
    mhm->DefinePartitionbyCoarseIndices(coarseindices);
    mhm->SetHDivFamily(HDivFamily::EHDivConstant);

    mhm->fMaterialIds = {1};
    mhm->fMaterialBCIds = {-1,-2,-3,-4};
    
    InsertMaterials(mhm->CMesh().operator->());

    mhm->SetInternalPOrder(1);
    mhm->SetSkeletonPOrder(1);

    mhm->DivideSkeletonElements(0);
    mhm->DivideBoundarySkeletonElements();

    bool substructure = true;

    std::map<int, std::pair<TPZGeoElSide, TPZGeoElSide>> intersectGeoElSide;
    mhm->BuildComputationalMesh(substructure,OpenChannel, intersectGeoElSide);

    //Jeferson: The following code is needed to enrich the MHM mesh with the Steklov eigenvectors. The code is commented because it may need refactor for the SPE10 problem
    // std::map<int, std::pair<TPZGeoElSide, TPZGeoElSide>> intersectGeoElSide;
    // intersectGeoElSide = IdentifyIntersections(mhm->FluxMesh().operator->(),intersectGeoElIndex,indexToSide);

    // std::cout << "intersectGeoElIndex = ";
    // for (const auto &it:intersectGeoElIndex)
    // {
    //     std::cout << it.second.first << " ";
    // }
    // std::cout << std::endl;
    // std::cout << "intersectGeoElSide = ";
    // for (const auto &it:intersectGeoElSide)
    // {
    //     std::cout << it.second.first.Element()->Index() << " ";
    // }
    // std::cout << std::endl;
    

    // // meshcontrol.BuildComputationalMesh(substructure, OpenChannel, test);
    // meshcontrol.BuildComputationalMesh(substructure, OpenChannel, intersectGeoElSide);

#ifdef ERRORESTIMATION_DEBUG
    if (1) {
        std::ofstream file("GMeshControlHDiv.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(mhm->GMesh().operator->(), file);
    }
#endif

    std::cout << "MHM Hdiv Computational meshes created\n";

    std::cout << "Number of equations MHMixed " << mhm->CMesh()->NEquations() << std::endl;
    
    TPZCompMesh *MixedMesh = mhm->CMesh().operator->();
    SolveProblem(mhm->CMesh(), mhm->GetMeshes(), Conf.GetExactSolution(),  Conf.GetVTKName(), Configuration);
    ProblemConfig config;
    config.dimension = 2;
    config.exact = nullptr;
    config.problemname = "MazeHdiv128x128";
    config.dir_name = "Results128x128";
    config.porder = 1;
    config.hdivmais = 1;
    config.materialids = {1};
    config.bcmaterialids = {-1, -2, -3, -4};
    config.makepressurecontinuous = true;
    config.ndivisions = 0;
    config.gmesh = MixedMesh->Reference();


    //open solution file
    std::ofstream solutionFile(Conf.solutionFileName, std::ios::app);
    auto solutionMHM = mhm->CMesh()->Solution();
    solutionMHM.Print("solution MHM", solutionFile,EMathematicaInput);
 

    TPZMultiphysicsCompMesh *originalMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm->CMesh().operator->());
    bool postProcWithHdiv = false;
//    TPZMHMHDivErrorEstimator ErrorEstimator(*originalMesh, mhm.operator->(), postProcWithHdiv);

    return 0;
}

void AnalyseSteklov(TPZSubCompMesh *sub, int count, int skelmat, TPZVec<std::set<TPZGeoElSide>> &eigGeoElSides);

// add geometric elements between macro domains with material ids determined by wrapids
void AddDomainWrapElements(TPZMHMeshControl &mhm, std::map<int,int> &wrapids);

int SteklovTest(ConfigCasesMaze &Conf, std::map<int,std::pair<int64_t,int64_t>> &intersectGeoElIndex, std::map<int64_t,int> &indexToSide){

    TRunConfig Configuration;

    //Creates geomesh and reads the permeability values from the SPE10 data
    auto perm_vec = TPZManVector<REAL, n_cells>(n_cells, 1);
    ReadSPE10CellPermeabilities(&perm_vec, layer);
    
    TPZGeoMesh *gmesh = CreateSPE10GeoMesh();
    {
        std::ofstream file(Conf.GetVTKName());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    }

    std::vector<REAL> x, y, perm;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            const int cell_id = ny * i + j;
            const double cell_perm = perm_vec[cell_id];
            x.push_back(0.5 + i);
            y.push_back(0.5 + j);
            perm.push_back(cell_perm);
        }
    }

    interpolator.setData(x.size(), x.data(), y.data(), perm.data());
    int interface_mat_id = 600; // it will not be used here...
    matextend[1] = matid1BC;
    matextend[2] = matid2BC;
    Conf.SetMHMOpenChannel(true);
    bool OpenChannel = Conf.GetMHMOpenChannel();

    TPZAutoPointer<TPZMHMixedMeshChannelControl> mhm = new TPZMHMixedMeshChannelControl(gmesh);
    {
        std::ofstream out("gmeshauto.txt");
        gmesh->Print(out);
    }

    TPZVec<int64_t> coarseindices;
    ComputeCoarseIndices(gmesh, coarseindices);
    // gmesh->AddInterfaceMaterial(1, 2, interface_mat_id);
    // gmesh->AddInterfaceMaterial(2, 1, interface_mat_id);

    int nInternalRef = 3;
    Tools::UniformRefinement(nInternalRef, 2, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    // criam-se apenas elementos geometricos
    mhm->DefinePartitionbyCoarseIndices(coarseindices);
    mhm->SetHDivFamily(HDivFamily::EHDivConstant);
    // mhm->SetHDivFamily(HDivFamily::EHDivStandard);

    AddDomainWrapElements(*mhm, matextend);
    
    mhm->fMaterialIds = {1};
    mhm->fMaterialBCIds = {-1,-2,-3,-4, matid1BC, matid2BC};

    InsertMaterials(mhm->CMesh().operator->());
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

    mhm->SetInternalPOrder(Conf.GetFluxOrder());
    mhm->SetSkeletonPOrder(Conf.GetFluxOrder());

    mhm->DivideSkeletonElements(Conf.GetSkeletonDivision());
    OpenChannel = false;
    mhm->DivideBoundarySkeletonElements();

    bool substructure = true;
    std::map<int, std::pair<TPZGeoElSide, TPZGeoElSide>> test;
    if (OpenChannel) {
        TPZCompMesh *flux_temp = MixedTest(Conf,2,2);
        std::cout << "flux_temp norm of solution " << Norm(flux_temp->Solution()) << std::endl;
        test = IdentifyChanel(flux_temp);
        flux_temp->Reference()->ResetReference();
//            delete flux_temp;
    }

    mhm->BuildComputationalMesh(substructure, OpenChannel, test);

    if (1) {
        std::ofstream file("GMeshControlHDiv.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(mhm->GMesh().operator->(), file);
    }

    std::cout << "MHM Hdiv Computational meshes created\n";

    std::cout << "Number of equations MHMixed " << mhm->CMesh()->NEquations() << std::endl;

    
    TPZCompMesh *MixedMesh = mhm->CMesh().operator->();
    {
        MixedMesh->ComputeNodElCon();
        std::ofstream out("EntryMesh.txt");
        MixedMesh->Print(out);
    }
    
    TPZCompMesh &cmesh = mhm->CMesh();
    // we have to zero the Neumann condition in order to identify the eigenvectors
    TPZBndCondT<STATE> *bc5 = dynamic_cast<TPZBndCondT<STATE> *>(cmesh.FindMaterial(-2));
    TPZVec<REAL> val2(1,0.),val25;
    val25 = bc5->Val2();
    bc5->SetVal2(val2);
    TPZBndCondT<STATE> *bc6 = dynamic_cast<TPZBndCondT<STATE> *>(cmesh.FindMaterial(-4));
    TPZVec<REAL> val26;
    val26 = bc5->Val2();
    bc6->SetVal2(val2);
    
    std::ofstream myfile("MixedMesh.txt");
    MixedMesh->Print(myfile);
    
    int64_t nelem = MixedMesh->NElements();
    int64_t count = 0;

    // TPZCompMesh *MHMCopy = MixedMesh;

    auto nsub = Conf.GetNumberOfSubdomains();
    TPZVec<std::set<TPZGeoElSide>> eigGeoElSides(nsub*nsub);

    for (int64_t el = 0; el < nelem; el++) {
        TPZCompEl *cel = MixedMesh->Element(el);
        auto *sub = dynamic_cast<TPZSubCompMesh *>(cel);
        if(sub) {
            AnalyseSteklov(sub,count,mhm->fSkeletonMatId,eigGeoElSides);
            count++;
        }
    }

    // for (int i = 0; i < eigGeoElSides.size(); i++)
    // {
    //     std::cout << "eigGeoElSides[" << i<< "]= " ;
    //     for (auto it:eigGeoElSides[i])
    //     {
    //         std::cout << it.Element()->Index() << " ";
    //     }
    //     std::cout << std::endl;
    // }


    // for (int i = 0; i < eigGeoElSides.size(); i++)
    // {
    //     std::cout << "eigGeoElSides[" << i << "] = " ;
    //     for (auto it:eigGeoElSides[i])
    //     {
    //         std::cout << it << " ";
    //     }
    //     std::cout << std::endl;
    // }

    count = 0;
    for (int64_t el = 0; el < nelem; el++) {
        TPZCompEl *cel = MixedMesh->Element(el);
        auto *sub = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!sub) continue;
        TPZMultiphysicsCompMesh *father = dynamic_cast<TPZMultiphysicsCompMesh *>(sub->Mesh());
        father->LoadReferences();
        TPZGeoMesh *gmesh = father->Reference();
        int64_t nel = gmesh->NElements();
        for (int i = 0; i < nel; i++){
            TPZGeoEl *gel = gmesh->ElementVec()[i];
            if (!gel) continue;
            int nSides = gel->NSides();
            int nCorder = gel->NCornerNodes();
            for (int iside = nCorder; iside < nSides; iside++){
                TPZGeoElSide gelside(gel,iside);
                if (eigGeoElSides[count].find(gelside) != eigGeoElSides[count].end()){
                    TPZGeoElBC gbc(gelside,100*(count+2));
                }
            }
        }
        count++;
        
    }
    std::ofstream file("GMeshAux.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(MixedMesh->Reference(), file); 

    AssociateGeoElSides(eigGeoElSides,intersectGeoElIndex,indexToSide);
    
    std::cout << "intersectGeoElIndex = ";
    for (const auto &it:intersectGeoElIndex)
    {
        std::cout << it.second.first << " ";
    }
    std::cout << std::endl;
    
    
    count = 0;
    for (int64_t el = 0; el < nelem; el++) {
        TPZCompEl *cel = MixedMesh->Element(el);
        auto *sub = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!sub) continue;
        TPZMultiphysicsCompMesh *father = dynamic_cast<TPZMultiphysicsCompMesh *>(sub->Mesh());
        father->LoadReferences();
        TPZGeoMesh *gmesh = father->Reference();
        int64_t nel = gmesh->NElements();
        for (int i = 0; i < nel; i++){
            TPZGeoEl *gel = gmesh->ElementVec()[i];
            if (!gel) continue;
            
            for (const auto &gelind : intersectGeoElIndex){
                if (gel->Index() != gelind.second.first) continue;

                TPZGeoElSide gelside(gel,indexToSide[gel->Index()]);
                TPZGeoElBC gbc(gelside,100*(count+6));
            }
        }
        count++;
        
    }
    std::ofstream file2("GMeshAuxNew.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(MixedMesh->Reference(), file2); 


    bc5->SetVal2(val25);
    bc6->SetVal2(val26);

    std::cout << "Number of equations MHMixed enriched " << mhm->CMesh()->NEquations() << std::endl;

    SolveProblem(mhm->CMesh(), mhm->GetMeshes(), Conf.GetExactSolution(),  Conf.GetVTKName(), Configuration);

    //open solution file
    std::ofstream solutionFile(Conf.solutionFileName, std::ios::app);
    auto solutionenriched = mhm->CMesh()->Solution();
    solutionenriched.Print("solution enriched", solutionFile,EMathematicaInput);

    return 0;
}

void AnalyseSteklov(TPZSubCompMesh *sub, int count, int skelmat, TPZVec<std::set<TPZGeoElSide>> &eigGeoElSides){
    

    // Identify element/sides that belong to a different mesh
    int64_t nel = sub->NElements();
    TPZMultiphysicsCompMesh *father = dynamic_cast<TPZMultiphysicsCompMesh *>(sub->Mesh());
    father->LoadReferences();
    auto &matvec = father->MaterialVec();
    auto &matvecsub = sub->MaterialVec();
        // pointers to the original materials
    auto matvecorig = matvec;
    auto *mat8 = matvecsub[8];
    // In the case of SPE 10, we do not sepparate into permeable or impermeable materials
    auto *mat9 = matvecsub[9];// should be a null pointer
    std::set<int> bndmat = {skelmat};
    std::map<int64_t,TPZGeoElSide> permeableconnects;
    TPZMixedDarcyFlow *darcy = dynamic_cast<TPZMixedDarcyFlow *> (matvec[1]);
    PermeabilityFunctionType a = PermeabilityFunction;
    darcy->SetPermeabilityFunction(a);
    TPZFNMatrix<2,REAL> val1(1,1,1.);
    TPZManVector<REAL> val2(1,0.);
    TPZGeoMesh *gmesh = father->Reference();
    int dim = gmesh->Dimension();
    // create a connection between the connects of the submesh and the skeleton elements
    std::map<int64_t,TPZCompEl *> connectToSkel;
    TPZCompEl *subcel = sub;
    int64_t ncon = subcel->NConnects();
    std::set<int64_t> activecon;
    // keep track of the connects of the submesh
    int64_t neq = 0;
    TPZManVector<int64_t,50> connectindexes(ncon);
    for (int64_t ic = 0; ic<ncon; ic++) {
        int64_t cindex = subcel->ConnectIndex(ic);
        TPZConnect &c = subcel->Connect(ic);
        neq += c.NShape()*c.NState();
        activecon.insert(cindex);
        connectindexes[ic] = cindex;
    }
    // plot the geometric elements
    {
        std::set<int64_t> elindices;
        int64_t nel = sub->NElements();
        for (int64_t el = 0; el<nel; el++) {
            auto cel = sub->Element(el);
            if(!cel) continue;
            TPZStack<TPZCompEl *> celstack;
            cel->GetCompElList(celstack);
            // identify the domain of the volumetric element
            // matid 2 = permeable
            // matid 1 = impermeable
            // matid 9 = permeable
            // matid 8 = impermeable
            int dommainmat = -1;
            for(auto el : celstack) {
                auto gel = el->Reference();
                if(!gel) DebugStop();
                if(gel->Dimension() == dim) {
                    dommainmat = gel->MaterialId();
                    break;
                }
            }
            if(dommainmat == -1) DebugStop();
            for(auto el : celstack) {
                auto gel = el->Reference();
                if(!gel) DebugStop();
                if(gel->Dimension() == dim-1)
                {
                    TPZGeoElSide gelside(gel);
                    elindices.insert(gel->Index());
                    auto neigh = gelside.HasNeighbour(bndmat);
                    if(!neigh) continue;
                    TPZCompEl *celskel = neigh.Element()->Reference();
                    if(!celskel) DebugStop();
                    if(celskel->Mesh() != father) continue;
                    int64_t celskelcindex = celskel->ConnectIndex(0);
                    if(activecon.find(celskelcindex) == activecon.end()) DebugStop();
                    connectToSkel[celskelcindex] = celskel;
                    
                    if(dommainmat == 2) {
                        if(gelside.HasNeighbour(9)) {
                            TPZGeoElSide neigh = gelside.HasNeighbour(dommainmat);
                            permeableconnects[celskelcindex] = neigh;
                        } else if(gelside.HasNeighbour(8)) {
                            
                        }
                        else
                        {
                            DebugStop();
                        }
                    }
                }
            }
        }
        int64_t nelfather = father->NElements();
        for(int64_t el = 0; el<nelfather; el++)
        {
            TPZCompEl *cel = father->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if(!gel) continue;
            if(gel->Dimension() != dim-1) DebugStop();
            if(gel->MaterialId() == skelmat) continue;
            int cindex = cel->ConnectIndex(0);
            if(activecon.find(cindex) != activecon.end()) {
                if(connectToSkel.find(cindex) != connectToSkel.end()) DebugStop();
                connectToSkel[cindex] = cel;
            }
        }
        std::stringstream sout;
        sout << "SubMesh_geo_" << count << ".vtk";
        std::ofstream out(sout.str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,elindices,out);
    }
    {
        TPZFMatrix GK,GM;
        std::stringstream sout;
        sout << "SubMesh_matrix_" << count << ".txt";
        std::ofstream out(sout.str());
        {
            int addlayer = 1;
            
            // switch the material objects
            // 9 permeable, 8 impermeable
            val1(0,0) = 1.;
            auto *bnd8 = darcy->CreateBC(darcy, 8, 2, val1, val2);
            val1(0,0) = 250000.;
            auto *bnd9 = darcy->CreateBC(darcy, 9, 2, val1, val2);
            if(addlayer) {
                matvecsub.erase(8);
                matvecsub.erase(9);
                sub->InsertMaterialObject(bnd8);
                sub->InsertMaterialObject(bnd9);
            }
            // compute the stiffness matrix
            TPZElementMatrixT<STATE> ek,ef;
            sub->CalcStiff(ek, ef);
            GK = ek.fMat;
            ek.Print(out);
            // switch the material objects back
            if(addlayer) {
                matvecsub.erase(8);
                matvecsub.erase(9);
                sub->InsertMaterialObject(mat8);
                sub->InsertMaterialObject(mat9);
                delete bnd8;
                delete bnd9;
            }
        }
        {
            // change the material of the skeleton to represent an L2 projection
            matvec.erase(skelmat);
            for(auto it : matvecorig) {
                if(it.first < 0) matvec.erase(it.first);
            }
            val1(0,0) = 1.;
            auto *bnd = darcy->CreateBC(darcy, skelmat, 2, val1, val2);
            father->InsertMaterialObject(bnd);
            for(auto it : matvecorig) {
                if(it.first < 0) {
                    auto *bnd = darcy->CreateBC(darcy, it.first, 2, val1, val2);
                    father->InsertMaterialObject(bnd);
                }
            }
            TPZElementGroup *celgr = new TPZElementGroup(*father);
            for(auto it : connectToSkel) celgr->AddElement(it.second);
            celgr->ReorderConnects(connectindexes);
            // compute the stiffness matrix
            TPZElementMatrixT<STATE> ek,ef;
            celgr->CalcStiff(ek, ef);
            GM = ek.fMat;
            ek.Print(out);
            celgr->Unwrap();
            // switch the material object back
            matvec.erase(skelmat);
            father->InsertMaterialObject(matvecorig[skelmat]);
            for(auto it : matvecorig) {
                if(it.first < 0) {
                    matvec.erase(it.first);
                    father->InsertMaterialObject(it.second);
                }
            }
            delete bnd;
        }
        int neq = GK.Rows();
        TPZVec<std::complex<double>> Lambda;
        TPZFMatrix<std::complex<double>> EigenVector;
        TPZFMatrix<STATE> GKcopy(GK),GMcopy(GM);
        GKcopy.SolveGeneralisedEigenProblem(GMcopy, Lambda, EigenVector);
        // std::cout << Lambda << std::endl;
        // EigenVector.Print(std::cout);
        {
            TPZFMatrix<STATE> Q(neq,neq);
            TPZFMatrix<STATE> eigvals(neq,1);
            for(int i=0; i<neq; i++) for(int j=0; j<neq; j++) Q(i,j) = EigenVector(i,j).real();
            for(int i=0; i<neq; i++) eigvals(i,0) = Lambda[i].real();
            Q.Print("Q = ",out,EMathematicaInput);
            eigvals.Print("Lambda = ",out,EMathematicaInput);
        }
        int neig = Lambda.size();
        std::stringstream filename;
        filename << "EigSub." << count;
        TPZVTKGenerator gen(sub, {"Flux","Pressure","DivFlux"}, filename.str(), 0);

        REAL tol = 1.e-2;
        for (int i = 0; i<neig; i++) {
            if (Lambda[i].real() > tol || Lambda[i].real()==0) continue;
            // std::cout << "Plot sequence " << i << " eigenvalue " << Lambda[i].real() << std::endl;
            TPZFMatrix<STATE> sol(neq,1);
            for(int ieq = 0; ieq<neq; ieq++) sol(ieq,0) = EigenVector(ieq,i).real();
            auto &block = father->Block();
            TPZFMatrix<STATE> &solmesh = father->Solution();
            solmesh.Zero();
            int loccount = 0;
            for(int ic = 0; ic < connectindexes.size(); ic++)
            {
                int64_t cindex = connectindexes[ic];
                bool isPermeable = false;
                // if(permeableconnects.find(cindex) != permeableconnects.end())
                if(permeableconnects[cindex])
                {
                    // std::cout << "permeable connect " << cindex << " sol ";
                    isPermeable = true;
                }
                else {
                    // std::cout << "impermeable connect " << cindex << " sol ";
                }
                TPZConnect &c = father->ConnectVec()[cindex];
                int64_t seqnum = c.SequenceNumber();
                int blsize = c.NShape()*c.NState();
                for (int ibl = 0; ibl<blsize; ibl++) {
                    int64_t pos = block.Index(seqnum, ibl);
                    solmesh(pos,0) = sol(loccount,0);
                    // cout << sol(loccount) << " ";
                    if(isPermeable && fabs(sol(loccount)) >= 1.e-5) {
                        eigGeoElSides[count].insert(permeableconnects[cindex]);
                    }
                    loccount++;
                }
                // std::cout << endl;
            }
            
            TPZFMatrix<STATE> residual;
            auto a = GK*sol;
            auto b = GM*sol;
//            a.Print("GK * sol",std::cout);
//            b.Print("GM * sol",std::cout);
            residual = a-Lambda[i].real()*b;
            // std::cout << "Residual of eigenvector " << Norm(residual) << std::endl;
//            solmesh.Print(std::cout);
            father->TPZCompMesh::LoadSolution(solmesh);
            father->TransferMultiphysicsSolution();
            
            gen.Do();
        }
        
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


void AssociateGeoElSides(TPZVec<std::set<TPZGeoElSide>> &eigGeoElSides, std::map<int,std::pair<int64_t,int64_t>> &intersectGeoElIndex, std::map<int64_t,int> &indexToSide){
    
    int nsub = eigGeoElSides.size();
    int count = 0;
    for (int isub = 0; isub < nsub; isub++)
    {  
        for (TPZGeoElSide gelside:eigGeoElSides[isub])
        {
            int thisMatId = gelside.Element()->MaterialId();
            for (int jsub = 0; jsub < nsub; jsub++)
            {  
                if (isub == jsub) continue;
                for (auto neighsides : eigGeoElSides[jsub])
                {
                    auto isNeigh = gelside.IsNeighbour(neighsides);
                    if (isNeigh){
                        //Check if the gelside already exists in the map
                        bool intersectExists = false;
                        for (const auto &gside:intersectGeoElIndex)
                        {
                            if(gside.second.first == neighsides.Element()->Index()){
                                intersectExists = true;
                                break;
                            }
                        }
                        if (!intersectExists){
                            intersectGeoElIndex[count]=std::make_pair(gelside.Element()->Index(),neighsides.Element()->Index());
                            indexToSide[gelside.Element()->Index()] = gelside.Side();
                            indexToSide[neighsides.Element()->Index()] = neighsides.Side();
                            count++;
                            break;
                        }
                    }    
                }
            }
        }
    }    
}


std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>> IdentifyIntersections(TPZCompMesh *cmesh, std::map<int,std::pair<int64_t,int64_t>> &intersectGeoElIndex, std::map<int64_t,int> &indexToSide){

    std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>> intersectGeoElSide;
    int count = 0;
    int64_t nel = cmesh->Reference()->NElements();
    for (int iel = 0; iel < nel; iel++)
    {
        TPZGeoEl *gel = cmesh->Reference()->ElementVec()[iel];
        if (!gel) continue;

        for (const auto &intIndex:intersectGeoElIndex)
        {
            if (gel->Index() == intIndex.second.first){
                TPZGeoElSide gelside(gel,indexToSide[intIndex.second.first]);
                TPZStack<TPZGeoElSide> allneigh;
                gelside.AllNeighbours(allneigh);
                for (const auto &ineigh : allneigh)
                {
                    if (ineigh.Element()->Index() == intIndex.second.second){
                        intersectGeoElSide[count] = std::make_pair(gelside,ineigh);
                        count++;
                        break;
                    }
                }
            }
            // if (gel->Index() == intIndex.second.second){
            //     TPZGeoElSide gelside(gel,indexToSide[intIndex.second.second]);
            //     TPZStack<TPZGeoElSide> allneigh;
            //     gelside.AllNeighbours(allneigh);
            //     for (const auto &ineigh : allneigh)
            //     {
            //         if (ineigh.Element()->Index() == intIndex.second.first){
            //             intersectGeoElSide[count] = std::make_pair(gelside,ineigh);
            //             count++;
            //             break;
            //         }
            //     }
            // }
        }
        
    }
    

    return intersectGeoElSide;
}


TPZGeoMesh *CreateSPE10GeoMesh() {
    std::cout << "Creating SPE10 initial grid...\n";

    const TPZManVector<REAL, 3> x0 = {0, 0, 0};
    const TPZManVector<REAL, 3> x1 = {nx, ny, 0.};
    const TPZManVector<int, 3> ndiv = {n_subdomainsx, n_subdomainsy, 0};

    TPZGenGrid2D gen(ndiv, x0, x1);

    gen.SetRefpatternElements(true);
    auto gmesh = new TPZGeoMesh;
    gen.Read(gmesh);

    gen.SetBC(gmesh, 4, -1);
    gen.SetBC(gmesh, 5, -2);
    gen.SetBC(gmesh, 6, -3);
    gen.SetBC(gmesh, 7, -4);

    std::cout << "SPE10 initial grid created. NElem: " << gmesh->NElements() << "\n";

    return gmesh;
}

void ReadSPE10CellPermeabilities(TPZVec<REAL> *perm_vec, const int layer) {

    std::cout << "Reading permeability data...\n";

    std::ifstream perm_file("../../../Projects/SPE10/InputData/spe_perm.dat", std::ios::in);
    if (!perm_file) {
        std::cerr << "Unable to open input file\n";
        DebugStop();
    }

    int cell_id = 0;
    const int n_cells = perm_vec->size();
    const int start_line = 1 + n_cells * (layer - 1) / 6;

    int line_num = 0;
    int line_num2 = 0;
    while (perm_file) {
        line_num++;
        line_num2++;
        std::string line;
        std::getline(perm_file, line, '\n');

        if (line_num < start_line) continue;

        std::stringstream stream(line);
        for (int i = 0; i < 6; i++) {
            stream >> perm_vec->operator[](cell_id);
            cell_id++;
        }
        if (cell_id == n_cells) break;
    }
    std::cout << "Finished reading permeability data from input file!\n";
}

void InsertMaterials(TPZCompMesh *cmesh) {

    typedef TPZMixedDarcyFlow TPZMixedPoisson;
    auto *mix = new TPZMixedPoisson(1, cmesh->Dimension());
    PermeabilityFunctionType a = PermeabilityFunction;
    mix->SetPermeabilityFunction(a);

    TPZFNMatrix<1, REAL> val1(1, 1, 0.);
    TPZManVector<REAL> val2(1, 0.);
    constexpr int dirichlet_bc = 0;

    // Pressure at reservoir boundary
    val2[0] = 1;
    TPZBndCond *pressure_left = mix->CreateBC(mix, -4, dirichlet_bc, val1, val2);

    val2[0] = 0;
    TPZBndCond *pressure_right = mix->CreateBC(mix, -2, dirichlet_bc, val1, val2);

    //Zero flux
    val2.Fill(0);
    TPZBndCond *pressure_bottom = mix->CreateBC(mix, -1, 1, val1, val2);

    TPZBndCond *pressure_top = mix->CreateBC(mix, -3, 1, val1, val2);

    cmesh->InsertMaterialObject(mix);
    cmesh->InsertMaterialObject(pressure_left);
    cmesh->InsertMaterialObject(pressure_right);
    cmesh->InsertMaterialObject(pressure_bottom);
    cmesh->InsertMaterialObject(pressure_top);
}

STATE PermeabilityFunction(const TPZVec<REAL> &x) {

    STATE perm;
    for (int i = 0; i < 2; i++) {
        perm = interpolator(x[0], x[1]);
        if (perm <= 1) {
            perm = 1;
        } else {
            perm += 1;
        }
    }
    // std::cout << "[" << x[0] << ", " << x[1] << "]\n";
    // std::cout << " perm = " << perm << std::endl;
    return perm;
    //std::cout << "[" << x[0] << ", " << x[1] << "]\n";
    //std::cout << "[" << res_mat(0, 0) << " " << res_mat(0, 1) << "\n";
    //std::cout        << res_mat(1, 0) << " " << res_mat(1, 1) << "\n";
    //std::cout        << res_mat(2, 0) << " " << res_mat(2, 1) << "\n";
    //std::cout        << res_mat(3, 0) << " " << res_mat(3, 1) << "]\n\n";
}