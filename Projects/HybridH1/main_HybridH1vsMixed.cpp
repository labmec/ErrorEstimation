// Comparison between HybridH1 and Mixed Poisson
// Modified from main_HybridH1.cpp
//
// Created by victor on 28/04/2020.
//

#include "pzlog.h"
#include "TPZRefPatternDataBase.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"

#include "tpzquadraticquad.h"
#include "tpzquadratictrig.h"
#include "tpzgeoelrefpattern.h"
#include "tpzarc3d.h"
#include "TPZGenGrid2D.h"

#include "ProblemConfig.h"

#include "mixedpoisson.h"
#include "TPZMatLaplacianHybrid.h"
#include "TPZNullMaterial.h"
#include "TPZVecL2.h"
#include "pzbndcond.h"

#include "pzintel.h"

#include "TPZMultiphysicsCompMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZHybridizeHDiv.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZMultiphysicsCompMesh.h"

#include "TPZHDivErrorEstimatorH1.h"

#include "Tools.h"

#include "TPZBFileStream.h"

#include "TPZCreateMultiphysicsSpace.h"
#include <tuple>
#include <memory>
#include <stack>
#include "TPZHybridH1ErrorEstimator.h"

bool neumann = true;

struct ErrorData
{
    std::ofstream ErroH1,ErroHybridH1,ErroMixed,Erro, timer;
    TPZVec<REAL> *LogH1,*LogHybridH1,*LogMixed, *rate, *Log;
    int maxdiv = 3;

    bool isMultiK = 0;
    REAL perm_Q1 = 5;
    REAL perm_Q2 = 1;

    REAL hLog = -1, h = -1000;
    int numErrors = 4;

    std::string plotfile;
    int mode = 3;           // 1 = "ALL"; 2 = "H1"; 3 = "HybridH1"; 4 = "Mixed"; 5 = "HybridSquared;
    int argc = 1;

    bool last = false, post_proc = true;
    int exp = 2; // Initial exponent of mesh refinement (numElem = 2*2^exp)
};
////Insert materials
void InsertMaterialObjectsH1Hybrid(TPZMultiphysicsCompMesh *cmesh, ProblemConfig &config,ErrorData &eData);
void BuildFluxMesh(TPZCompMesh *cmesh_flux, ProblemConfig &config);
void BuildFluxMesh_MultiK(TPZCompMesh *cmesh_flux, ProblemConfig &config);
void BuildPotentialMesh(TPZCompMesh *cmesh_p, ProblemConfig &config);
void BuildPotentialMesh_MultiK(TPZCompMesh *cmesh_p, ProblemConfig &config);
void InsertMaterialMixed(TPZMultiphysicsCompMesh *cmesh_mixed, ProblemConfig config,ErrorData &eData);
void SetMultiPermeMaterials(TPZGeoMesh* gmesh);
void CreateMaterialMultiK_Hybrid(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, REAL permQ1, REAL permQ2,ProblemConfig &config);
void CreateMaterialMultiK_Mixed(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, REAL permQ1, REAL permQ2,ProblemConfig &config);
TPZGeoMesh* CreateGeoMesh_OriginCentered(int nel, TPZVec<int>& bcids);
////Computational mesh and FEM solvers
void CreateHybridH1ComputationalMesh(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, int &InterfaceMatId,ErrorData &eData, ProblemConfig &config,int hybridLevel);
void CreateMixedComputationalMesh(TPZMultiphysicsCompMesh *cmesh_H1Mixed,ErrorData &eData, ProblemConfig &config);
void SolveH1Problem(TPZCompMesh *cmeshH1,struct ProblemConfig &config, struct ErrorData &eData);
void SolveHybridH1Problem(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, int InterfaceMatId, struct ProblemConfig config, struct ErrorData &eData,int hybridLevel);
void SolveMixedProblem(TPZMultiphysicsCompMesh *cmesh_Mixed,struct ProblemConfig config,struct ErrorData &eData);
void StockErrorsH1(TPZAnalysis &an,TPZCompMesh *cmesh,ofstream &Erro, TPZVec<REAL> *Log, ErrorData &eData);
void StockErrors(TPZAnalysis &an,TPZMultiphysicsCompMesh *cmesh,ofstream &Erro, TPZVec<REAL> *Log, ErrorData &eData);
void nElementsPerMaterial(TPZMultiphysicsCompMesh *cmesh);
////Output setup
void FlushTable(ErrorData &eData,char *argv[]);
void CleanErrors(string file);
void InvertError(string file);
void FillErrors(ofstream &table,string f1,string f2,string f3);
void FillErrors(ofstream &table,string f,int mode);
void FillLegend(ofstream &table,int hash_count, int it_count);
void FlushTime(ErrorData &eData, clock_t start);
////Command Line Entry
////The command line entry is useful for consecutive runs
////One can build a .sh documment to run all of them.
void EvaluateEntry(int argc, char *argv[],ErrorData &eData);
void InitializeOutstream(ErrorData &eData,char *argv[]);
void IsInteger(char *argv);

void Configure(ProblemConfig &config,int ndiv,ErrorData &eData,char *argv[]){
    config.porder = 3;         // Potential and internal flux order
    config.hdivmais = 2;       // p_order - hdivmais = External flux order
    config.H1Hybridminus = 2;  // p_order - H1HybridMinus = Flux order
    config.ndivisions = ndiv;
    config.dimension = 2;
    config.prefine = false;

    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin/*EArcTan*//*ESinSin*/;
    config.exact.operator*().fSignConvention = 1;

    config.problemname = "ESinSin"/*"EArcTan"*//*"ESinSin"*/;

    config.dir_name = "Pressure_Reconstruction_VTK";
    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    // geometric mesh
    TPZManVector<int, 4> bcids(4, -1);
    TPZGeoMesh *gmesh;
    if(eData.isMultiK == false) gmesh = Tools::CreateGeoMesh(1, bcids); //rectangular mesh [0,1]x[0,1], matID = 1;
    else {
        gmesh = CreateGeoMesh_OriginCentered(1, bcids); //rectangular mesh [-1,1]x[-1,1], matID_Q1-Q3 = alpha, matID_Q2-Q4 = beta
    }
    Tools::UniformRefinement(config.ndivisions, gmesh);

    config.gmesh = gmesh;
    config.materialids.insert(1);
    if(eData.isMultiK == false) config.bcmaterialids.insert(-1);
    if(eData.isMultiK == true) {
        config.materialids.insert(2);config.materialids.insert(3);
    }
    if(eData.isMultiK == true) {
        config.bcmaterialids.insert(-5);config.bcmaterialids.insert(-6);config.bcmaterialids.insert(-8);config.bcmaterialids.insert(-9);
    }

    int refinement_depth = 2;
    if(config.ndivisions > 0) {
        /*Tools::RandomRefinement(config.gmesh, 1, refinement_depth);*/
        config.gmesh->SetDimension(2);
        Tools::RefineElements(config.gmesh, {5, 6});
    }

    if(eData.argc != 1) {
        config.porder = atoi(argv[2]);
        config.hdivmais = atoi(argv[4]);
        config.H1Hybridminus = atoi(argv[3]);
    }
}

/**
 * @brief This "main" solves a given mesh using the H1, HybridH1 and Mixed approximations;
 * @brief One can run this "main" using Command Line or manually within your programming IDE;
 * @brief Running manually: Set p-order, H1Hybridmais and hdivmais within "Configure",
 * @brief maximum number of refinements (maxdiv) and "mode" within ErrorData;
 * @brief mode stands for 1 = "ALL", 2 = "H1", 3 = "HybridH1" and 4 = "Mixed";
 * @brief Running via Command Line: call "HybridH1vsMixed mode porder H1Hybridminus hdivmais"
 * @brief at the executable directory. One may define a bash file running multiple simulations
 * @brief using the Command Line. Check BashFile.sh within HybridH1 directory.
 */
int main(int argc, char *argv[]) {

    //Testing random stuff (begin)
    TPZManVector<int, 4> bcids(4, -1);
    TPZGeoMesh *geomesh = Tools::CreateGeoMesh(1, bcids);
    Tools::UniformRefinement(1, geomesh);
    TPZGeoEl *gel = geomesh->Element(0);
    TPZGeoElSide gelside(gel, gel->NSides() - 1);
    gelside.LowerLevelCompElementList2(1);


    //Testing random stuff (end)
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    ErrorData eData;
    EvaluateEntry(argc,argv,eData);

    InitializeOutstream(eData,argv);

    const clock_t begin_iter = clock();

    for (int ndiv = 1; ndiv < /*eData.maxdiv+2*/2; ndiv++) { //ndiv = 1 corresponds to a 2x2 mesh.
        if (ndiv == eData.maxdiv+1) eData.last = true;
        eData.h = 1./eData.exp;

        ProblemConfig config;
        Configure(config,ndiv,eData,argv);
        std::cout <<"HyybridH1Minus = " << config.H1Hybridminus << "\n\n\n\n";
        if(eData.isMultiK){
            SetMultiPermeMaterials(config.gmesh);
        }


        if(eData.last && eData.post_proc) {
            std::ofstream out(eData.plotfile + "/gmesh.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, out);
            std::ofstream out2(eData.plotfile + "/gmeshInitial.txt");
            config.gmesh->Print(out2);
        }

        const clock_t start = clock();

        //H1
        TPZCompMesh *cmeshH1 = Tools::CMeshH1(config);
        TPZCompMeshTools::CreatedCondensedElements(cmeshH1, false, false);
        if(eData.mode == 1 || eData.mode == 2) {
            SolveH1Problem(cmeshH1, config, eData);
        }

        //HybridH1
        TPZMultiphysicsCompMesh *cmesh_H1Hybrid = new TPZMultiphysicsCompMesh(config.gmesh);
        if(eData.mode == 1 || eData.mode == 3){
            clock_t start = clock();
            int interfaceMatID = -10;
            int hybridLevel = 1;
            CreateHybridH1ComputationalMesh(cmesh_H1Hybrid, interfaceMatID,eData, config,hybridLevel);
            SolveHybridH1Problem(cmesh_H1Hybrid, interfaceMatID, config, eData,hybridLevel);
            TPZHybridH1ErrorEstimator HybridH1Estimate(*cmesh_H1Hybrid);
            HybridH1Estimate.fProblemConfig = config;
            HybridH1Estimate.SetAnalyticSolution(config.exact);
            HybridH1Estimate.PotentialReconstruction();

            TPZManVector<REAL> elementerrors;
            TPZVec<REAL> errorVec;
            HybridH1Estimate.ComputeErrors(errorVec, elementerrors, true);

            FlushTime(eData,start);
        }

        //Mixed
        TPZMultiphysicsCompMesh *cmesh_mixed = new TPZMultiphysicsCompMesh(config.gmesh);
        if(eData.mode == 1 || eData.mode == 4) {
            clock_t start = clock();
            CreateMixedComputationalMesh(cmesh_mixed, eData, config);
            SolveMixedProblem(cmesh_mixed, config, eData);
            TPZHybridHDivErrorEstimator test(*cmesh_mixed);
            //std::set<int> bcmatID;
            //bcmatID.insert(-5);bcmatID.insert(-6);bcmatID.insert(-8);bcmatID.insert(-9);
            //config.bcmaterialids = bcmatID;
            test.fProblemConfig = config;
            test.fOriginalIsHybridized = true;
            test.PotentialReconstruction();
            FlushTime(eData,start);
        }

        //HybridSquared
        TPZMultiphysicsCompMesh *cmesh_HybridSquared = new TPZMultiphysicsCompMesh(config.gmesh);
        if(eData.mode == 1 || eData.mode == 5){
            clock_t start = clock();
            int interfaceMatID = -10;
            int hybridLevel = 2;
            CreateHybridH1ComputationalMesh(cmesh_HybridSquared, interfaceMatID,eData, config,hybridLevel);
            SolveHybridH1Problem(cmesh_HybridSquared, interfaceMatID, config, eData, hybridLevel);
            TPZHybridH1ErrorEstimator test(*cmesh_HybridSquared);
            test.fProblemConfig = config;
            test.PotentialReconstruction();
            FlushTime(eData,start);
        }
        //finish clock

        if(eData.last && eData.post_proc) {
            const clock_t begin_time = clock();
            if(eData.mode == 1 || eData.mode == 3){
                std::ofstream out2(eData.plotfile + "/H1HybridMesh.txt");
                cmesh_H1Hybrid->Print(out2);
            }
            if(eData.mode == 1 || eData.mode == 2) {
                std::ofstream out3(eData.plotfile + "/H1Mesh.txt");
                cmeshH1->Print(out3);
            }
            if(eData.mode == 1 || eData.mode == 4) {
                std::ofstream out4(eData.plotfile + "/MixedMesh.txt");
                cmesh_mixed->Print(out4);
            }
            std::cout << "printing comp mesh time: " << float( clock () - begin_time )/CLOCKS_PER_SEC <<endl;
        }

        eData.hLog = eData.h;
        eData.exp *=2;
    }
    std::string command = "cp Erro.txt " + eData.plotfile + "/Erro.txt";
    system(command.c_str());

    std::cout << "Total Time: " <<float( clock () - begin_iter )/CLOCKS_PER_SEC <<"\n";
    FlushTable(eData,argv);

    return 0.;
}

void nElementsPerMaterial(TPZMultiphysicsCompMesh *cmesh){
    cmesh->LoadReferences();
    std::map<int,int> matelem;
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
//                TPZManVector<REAL,3> center(3);
//                TPZGeoElSide gelside(gel);
//                gelside.CenterX(center);
//                std::cout << "Matid " << gel->MaterialId() << " center " << center << std::endl;
        matelem[gel->MaterialId()]++;
    }
    std::cout << __PRETTY_FUNCTION__ << " number of computational elements by material \n";
    for (auto it : matelem) {
        std::cout << "Material id " << it.first << " number of elements " << it.second << std::endl;
    }
}

void CreateMaterialMultiK_Mixed(TPZMultiphysicsCompMesh *cmesh_mixed, REAL permQ1, REAL permQ2,ProblemConfig &config){
    int matID_Q1 = 2;
    int matID_Q2 = 3;
    int dim = config.gmesh->Dimension();
    int dirichlet = 0;

    cmesh_mixed -> SetDefaultOrder(config.porder);
    cmesh_mixed -> SetDimModel(dim);
    cmesh_mixed->SetAllCreateFunctionsMultiphysicElem();

    TLaplaceExample1 *mat1 = new TLaplaceExample1,*mat2 = new TLaplaceExample1;
    mat1->fExact = TLaplaceExample1::ESinSin/*EArcTan*/;
    mat2->fExact = TLaplaceExample1::ESinSin/*EArcTan*/;
    mat1->fSignConvention = 1;
    mat2->fSignConvention = 1;
    TPZFNMatrix<9,REAL> K, invK;
    K.Resize(3,3); invK.Resize(3,3);
    K.Identity();invK.Identity();
    K(0,0) = K(1,1) = permQ1;
    invK(0,0) = invK(1,1) = 1./permQ1;
    mat1->setPermeabilyTensor(K,invK);
    K(0,0) = K(1,1) = permQ2;
    invK(0,0) = invK(1,1) =  1./permQ2;
    mat2->setPermeabilyTensor(K,invK);

    TPZMixedPoisson *material_Q1 = new TPZMixedPoisson(matID_Q1,dim); //Using standard PermealityTensor = Identity.
    TPZMixedPoisson *material_Q2 = new TPZMixedPoisson(matID_Q2,dim);
    material_Q1->SetForcingFunction(config.exact.operator*().ForcingFunction());
    material_Q1->SetForcingFunctionExact(config.exact.operator*().Exact());
    material_Q2->SetForcingFunction(config.exact.operator*().ForcingFunction());
    material_Q2->SetForcingFunctionExact(config.exact.operator*().Exact());

    material_Q1->SetPermeability(permQ1);
    material_Q2->SetPermeability(permQ2);

    cmesh_mixed->InsertMaterialObject(material_Q1);
    cmesh_mixed->InsertMaterialObject(material_Q2);

    //Boundary Conditions
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);

    TPZMaterial *BCond0_Q1 = material_Q1->CreateBC(material_Q1, -5, dirichlet, val1, val2);
    BCond0_Q1->SetForcingFunction(config.exact.operator*().Exact());

    TPZMaterial *BCond0_Q2 = material_Q2->CreateBC(material_Q2, -6, dirichlet, val1, val2);
    BCond0_Q2->SetForcingFunction(config.exact.operator*().Exact());

    TPZMaterial *BCond1_Q1 = material_Q1->CreateBC(material_Q1, -8, neumann, val1, val2);
    TPZMaterial *BCond1_Q2 = material_Q1->CreateBC(material_Q2, -9, neumann, val1, val2);

    cmesh_mixed->InsertMaterialObject(BCond0_Q1);
    cmesh_mixed->InsertMaterialObject(BCond0_Q2);
    cmesh_mixed->InsertMaterialObject(BCond1_Q1);
    cmesh_mixed->InsertMaterialObject(BCond1_Q2);
}

void CreateMaterialMultiK_Hybrid(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, REAL permQ1, REAL permQ2,ProblemConfig &config){

    TPZGeoMesh* gmesh = cmesh_H1Hybrid->Reference();
    int matID_Q1 = 2;
    int matID_Q2 = 3;
    int dim = gmesh->Dimension();
    int dirichlet = 0;

    TLaplaceExample1 *mat1 = new TLaplaceExample1,*mat2 = new TLaplaceExample1;
    mat1->fExact = TLaplaceExample1::ESinSin/*EArcTan*/;
    mat2->fExact = TLaplaceExample1::ESinSin/*EArcTan*/;
    mat1->fSignConvention = 1;
    mat2->fSignConvention = 1;
    TPZFNMatrix<9,REAL> K, invK;
    K.Resize(3,3); invK.Resize(3,3);
    K.Identity();invK.Identity();
    K(0,0) = K(1,1) = permQ1;
    invK(0,0) = invK(1,1) = 1./permQ1;
    mat1->setPermeabilyTensor(K,invK);
    K(0,0) = K(1,1) = permQ2;
    invK(0,0) = invK(1,1) =  1./permQ2;
    mat2->setPermeabilyTensor(K,invK);

    TPZMatLaplacianHybrid *material_Q1 = new TPZMatLaplacianHybrid(matID_Q1, dim);
    TPZMatLaplacianHybrid *material_Q2 = new TPZMatLaplacianHybrid(matID_Q2, dim);

    material_Q1->SetPermeability(permQ1);
    material_Q2->SetPermeability(permQ2);

    cmesh_H1Hybrid->InsertMaterialObject(material_Q1);
    cmesh_H1Hybrid->InsertMaterialObject(material_Q2);

    if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {

        //material_Q1->SetForcingFunction(config.exact.operator*().ForcingFunction());
        //material_Q1->SetForcingFunctionExact(config.exact.operator*().Exact());
        material_Q1->SetForcingFunction(mat1->ForcingFunction());
        material_Q1->SetForcingFunctionExact(mat1->Exact());

        //material_Q2->SetForcingFunction(config.exact.operator*().ForcingFunction());
        //material_Q2->SetForcingFunctionExact(config.exact.operator*().Exact());
        material_Q2->SetForcingFunction(mat2->ForcingFunction());
        material_Q2->SetForcingFunctionExact(mat2->Exact());
    }

    // Inserts boundary conditions
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 1.);
    TPZMaterial *BCond0_Q1 = material_Q1->CreateBC(material_Q1, -5, dirichlet, val1, val2);
    TPZMaterial *BCond0_Q2 = material_Q2->CreateBC(material_Q2, -6, dirichlet, val1, val2);
    if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
        //BCond0_Q1->SetForcingFunction(config.exact.operator*().Exact());
        //BCond0_Q2->SetForcingFunction(config.exact.operator*().Exact());
        BCond0_Q1->SetForcingFunction(mat1->Exact());
        BCond0_Q2->SetForcingFunction(mat2->Exact());
    }
    val2.Zero();
    TPZMaterial *BCond1_Q1 = material_Q1->CreateBC(material_Q1, -8, neumann, val1, val2);
    TPZMaterial *BCond1_Q2 = material_Q1->CreateBC(material_Q1, -9, neumann, val1, val2);

    cmesh_H1Hybrid->InsertMaterialObject(BCond0_Q1);
    cmesh_H1Hybrid->InsertMaterialObject(BCond0_Q2);
    cmesh_H1Hybrid->InsertMaterialObject(BCond1_Q1);
    cmesh_H1Hybrid->InsertMaterialObject(BCond1_Q2);

    /*config.materialids.insert(matID_Q1);
    config.materialids.insert(matID_Q2);
    config.bcmaterialids.insert(-5);
    config.bcmaterialids.insert(-6);
    config.bcmaterialids.insert(-8);
    config.bcmaterialids.insert(-9);*/
}

void SetMultiPermeMaterials(TPZGeoMesh* gmesh){

    //TPZGeoMesh* gmesh = cmesh_H1Hybrid->Reference();
    TPZAdmChunkVector<TPZGeoEl *> &elvec = gmesh->ElementVec();
    int numEl = elvec.NElements();
    TPZVec<int64_t> nodeInd;
    int numberNodes;
    TPZManVector<REAL,2> elCenterCoord = {0.,0.};

    //Getting center coordinates of the element
    for(int ind = 0; ind < numEl; ind++) {
        TPZGeoEl *gel = elvec[ind];
        gel->GetNodeIndices(nodeInd);
        numberNodes = gel->NNodes();
        elCenterCoord[0] = elCenterCoord[1] = 0.;
        for (int nodeIter = 0; nodeIter < numberNodes; nodeIter++) {
            TPZGeoNode nodo = gmesh->NodeVec()[nodeInd[nodeIter]];
            elCenterCoord[0] += nodo.Coord(0);
            elCenterCoord[1] += nodo.Coord(1);
        }
        elCenterCoord[0] /= numberNodes;
        elCenterCoord[1] /= numberNodes;

        // Changing the permeability
        //std::cout <<"el: " << ind  << "x_ave: " << elCenterCoord[0] << "y_ave: " << elCenterCoord[1] << std::endl;
        if((elCenterCoord[0] >=0. && elCenterCoord[1] >= 0.) || (elCenterCoord[0] < 0. && elCenterCoord[1] <= 0.)) { // first/third quadrants
            if(gel->Dimension() == gmesh->Dimension()) { /*if(gel->MaterialId() == 1)*/ gel->SetMaterialId(2);}
            else {/*if (gel->MaterialId() < 0 && gel->MaterialId() > -5)*/ gel->SetMaterialId(-5);} // Assumes dirichlet boundary condition
        }
        else{
            if(gel->Dimension() == gmesh->Dimension())  { /*if(gel->MaterialId() == 1)*/ gel->SetMaterialId(3);}
            else {/*if (gel->MaterialId() < 0 && gel->MaterialId() > -5)*/ gel->SetMaterialId(-6);}// Assumes dirichlet boundary condition
        }
    }
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
    //gmesh->Print(std::cout);
}

TPZGeoMesh* CreateGeoMesh_OriginCentered(int nel, TPZVec<int>& bcids) {

    TPZManVector<int> nx(2, nel);
    TPZManVector<REAL> x0(3, -1.), x1(3, 1.);
    x1[2] = x0[2] = 0.;
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

void FlushTime(ErrorData &eData, clock_t start){
    float timer = float( clock () - start )/CLOCKS_PER_SEC;
    eData.timer << "Simulation time (" << eData.h << "x" << eData.h <<"): " << timer << "\n";
    eData.timer.flush();
}

void CreateMixedComputationalMesh(TPZMultiphysicsCompMesh *cmesh_Mixed, ErrorData &eData, ProblemConfig &config){

    int matID = 1;
    int dim = config.gmesh->Dimension();

    //Flux mesh creation
    TPZCompMesh *cmesh_flux = new TPZCompMesh(config.gmesh);
    BuildFluxMesh(cmesh_flux, config);

    //Potential mesh creation
    TPZCompMesh *cmesh_p = new TPZCompMesh(config.gmesh);
    BuildPotentialMesh(cmesh_p, config);

    //Multiphysics mesh build
    InsertMaterialMixed(cmesh_Mixed, config,eData);
    TPZManVector<int> active(2, 1);
    TPZManVector<TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmesh_flux;
    meshvector[1] = cmesh_p;
    TPZCompMeshTools::AdjustFluxPolynomialOrders(meshvector[0], config.hdivmais); //Increases internal flux order by "hdivmais"
    TPZCompMeshTools::SetPressureOrders(meshvector[0], meshvector[1]);//Set the pressure order the same as the internal flux

    cmesh_Mixed->BuildMultiphysicsSpace(active,meshvector);
    bool keeponelagrangian = true, keepmatrix = false;
    TPZCompMeshTools::CreatedCondensedElements(cmesh_Mixed, keeponelagrangian, keepmatrix);
    cmesh_Mixed->LoadReferences();
    cmesh_Mixed->InitializeBlock();
}

void CreateHybridH1ComputationalMesh(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,int &interFaceMatID , ErrorData &eData, ProblemConfig &config,int hybridLevel){
    auto spaceType = TPZCreateMultiphysicsSpace::EH1Hybrid;
    if(hybridLevel == 2) {
        spaceType = TPZCreateMultiphysicsSpace::EH1HybridSquared;
    }
    else if(hybridLevel != 1) {
        DebugStop();
    }

    TPZCreateMultiphysicsSpace createspace(config.gmesh, spaceType);
    //TPZCreateMultiphysicsSpace createspace(config.gmesh);
    std::cout << cmesh_H1Hybrid->NEquations();

    createspace.SetMaterialIds({1,2,3}, {-6,-5,-2,-1});
    createspace.fH1Hybrid.fHybridizeBCLevel = 1;//opcao de hibridizar o contorno
    createspace.ComputePeriferalMaterialIds();

    TPZManVector<TPZCompMesh *> meshvec;

    int lagrangeOrder = config.porder - config.H1Hybridminus;
    createspace.CreateAtomicMeshes(meshvec,config.porder,lagrangeOrder);

    InsertMaterialObjectsH1Hybrid(cmesh_H1Hybrid, config,eData);
    createspace.InsertPeriferalMaterialObjects(cmesh_H1Hybrid);
    cmesh_H1Hybrid->BuildMultiphysicsSpace(meshvec);
    createspace.InsertLagranceMaterialObjects(cmesh_H1Hybrid);

    createspace.AddInterfaceElements(cmesh_H1Hybrid);
    createspace.GroupandCondenseElements(cmesh_H1Hybrid);

    cmesh_H1Hybrid->InitializeBlock();
    cmesh_H1Hybrid->ComputeNodElCon();

    interFaceMatID = createspace.fH1Hybrid.fLagrangeMatid.first;

}

void InitializeOutstream(ErrorData &eData, char *argv[]){
    //Error buffer
    if(eData.mode == 1) {

        if( remove( "ErroH1.txt" ) != 0 || remove( "ErroHybridH1.txt" ) != 0 || remove( "ErroMixed.txt" ) != 0 )
            perror( "Error deleting file" );
        else
            puts( "Error log successfully deleted" );

        eData.ErroH1.open("ErroH1.txt", std::ofstream::app);
        eData.ErroHybridH1.open("ErroHybridH1.txt", std::ofstream::app);
        eData.ErroMixed.open("ErroMixed.txt", std::ofstream::app);

        eData.ErroH1 << "----------COMPUTED ERRORS (H1)----------\n";
        eData.ErroH1.flush();
        eData.ErroHybridH1 << "-------COMPUTED ERRORS (HYBRID_H1)------\n";
        eData.ErroHybridH1.flush();
        eData.ErroMixed << "--------COMPUTED ERRORS (MIXED)---------\n";
        eData.ErroMixed.flush();

        // Log Initialization
        eData.LogH1 = new TPZVec<REAL>(eData.numErrors, -1);
        eData.LogHybridH1 = new TPZVec<REAL>(eData.numErrors, -1);
        eData.LogMixed = new TPZVec<REAL>(eData.numErrors, -1);
        eData.rate = new TPZVec<REAL>(eData.numErrors, -1);
    }
    else{
        if( remove( "Erro.txt" ) != 0) perror( "Error deleting file" );
        else puts( "Error log successfully deleted" );

        eData.Erro.open("Erro.txt",std::ofstream::app);
        eData.Erro << "----------COMPUTED ERRORS----------\n";

        eData.Log = new TPZVec<REAL>(eData.numErrors, -1);
        eData.rate = new TPZVec<REAL>(eData.numErrors, -1);
    }



    ProblemConfig config;
    Configure(config,0,eData,argv);

    std::stringstream out;
    if(eData.mode == 1) {
        out << "Quad_" << config.problemname << "___porder_"
            << config.porder << "___hdivmais_" << config.hdivmais
            << "___H1Hybridminus_" << config.H1Hybridminus;
        eData.plotfile = out.str();
    }
    else{
        switch(eData.mode) {
            case 2:
                out << "Quad_" << config.problemname << "___porder_"
                    << config.porder;
                eData.plotfile = out.str();
                break;
            case 3:
                out << "Quad_" << config.problemname << "___porder_"
                    << config.porder << "___H1Hybridminus_" << config.H1Hybridminus;
                eData.plotfile = out.str();
                break;
            case 4:
                out << "Quad_" << config.problemname << "___porder_"
                    << config.porder << "___hdivmais_" << config.hdivmais;
                eData.plotfile = out.str();
                break;
            case 5:
                out << "Quad_" << config.problemname << "___porder_"
                    << config.porder << "___HybridSquaredminus_" << config.H1Hybridminus;
                eData.plotfile = out.str();
                break;
            default:
                std::cout << "Invalid mode number";
                break;
        }
    }

    std::string command = "mkdir " + eData.plotfile;
    system(command.c_str());

    std::string timer_name = eData.plotfile + "/timer.txt";

    if( remove(timer_name.c_str()) != 0)
        perror( "Error deleting file" );
    else
        puts( "Error log successfully deleted" );

    eData.timer.open(timer_name, std::ofstream::app);
}

void BuildFluxMesh(TPZCompMesh *cmesh_flux, ProblemConfig &config){

    int dim = config.gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    int flux_order = config.porder - config.hdivmais;
    cmesh_flux->SetDefaultOrder(flux_order);
    cmesh_flux->SetDimModel(config.gmesh->Dimension());

    cmesh_flux->SetAllCreateFunctionsHDiv();

    TPZNullMaterial *material = new TPZNullMaterial(matID);
    material->SetDimension(dim);
    cmesh_flux->InsertMaterialObject(material);

    //(EE)Create Boundary conditions
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);

    //(PROPRIETARY)
    TPZMaterial *BCond0 = material->CreateBC(material,-1,dirichlet,val1,val2);
    cmesh_flux->InsertMaterialObject(BCond0);

    TPZMaterial *BCond1 = material->CreateBC(material,-2,neumann,val1,val2);
    cmesh_flux->InsertMaterialObject(BCond1);

    //This part is for multi-K problemas
    BuildFluxMesh_MultiK(cmesh_flux,config);

    cmesh_flux->AutoBuild();
    cmesh_flux->InitializeBlock();
}

void BuildFluxMesh_MultiK(TPZCompMesh *cmesh_flux, ProblemConfig &config) {

    int dim = config.gmesh->Dimension();
    int matID_Q1 = 2;
    int matID_Q2 = 3;
    int dirichlet = 0;
    int neumann = 1;

    TPZNullMaterial *material_Q1 = new TPZNullMaterial(matID_Q1);
    TPZNullMaterial *material_Q2 = new TPZNullMaterial(matID_Q2);
    material_Q1->SetDimension(dim);
    material_Q2->SetDimension(dim);
    cmesh_flux->InsertMaterialObject(material_Q1);
    cmesh_flux->InsertMaterialObject(material_Q2);

    //(EE)Create Boundary conditions
    TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);

    //(PROPRIETARY)
    TPZMaterial *BCond0_Q1 = material_Q1->CreateBC(material_Q1, -5, dirichlet, val1, val2);
    TPZMaterial *BCond0_Q2 = material_Q1->CreateBC(material_Q2, -6, dirichlet, val1, val2);
    cmesh_flux->InsertMaterialObject(BCond0_Q1);
    cmesh_flux->InsertMaterialObject(BCond0_Q2);

    TPZMaterial *BCond1_Q1 = material_Q1->CreateBC(material_Q1, -8, neumann, val1, val2);
    TPZMaterial *BCond1_Q2 = material_Q2->CreateBC(material_Q2, -9, neumann, val1, val2);
    cmesh_flux->InsertMaterialObject(BCond1_Q1);
    cmesh_flux->InsertMaterialObject(BCond1_Q2);
}


void BuildPotentialMesh(TPZCompMesh *cmesh_p, ProblemConfig &config){
    int matID = 1;
    int dim = config.gmesh->Dimension();

    int potential_order = config.porder - config.hdivmais;
    cmesh_p->SetDefaultOrder(potential_order);
    cmesh_p->SetDimModel(dim);

    cmesh_p->SetAllCreateFunctionsContinuous(); //H1 functions
    cmesh_p->ApproxSpace().CreateDisconnectedElements(true);

    TPZNullMaterial *material = new TPZNullMaterial(matID); material->SetDimension(dim);
    cmesh_p->InsertMaterialObject(material);

    BuildPotentialMesh_MultiK(cmesh_p, config);

    cmesh_p->AutoBuild();
    cmesh_p->ExpandSolution();

    TPZAdmChunkVector<TPZConnect> &nodeIt = cmesh_p->ConnectVec();
    for(auto &nodo : nodeIt){
        nodo.SetLagrangeMultiplier(1);
    }
}

void BuildPotentialMesh_MultiK(TPZCompMesh *cmesh_p, ProblemConfig &config){
    int matID_Q1 = 2;
    int matID_Q2 = 3;
    int dim = config.gmesh->Dimension();

    TPZNullMaterial *material_Q1 = new TPZNullMaterial(matID_Q1); material_Q1->SetDimension(dim);
    cmesh_p->InsertMaterialObject(material_Q1);

    TPZNullMaterial *material_Q2 = new TPZNullMaterial(matID_Q2); material_Q2->SetDimension(dim);
    cmesh_p->InsertMaterialObject(material_Q2);
}

void InsertMaterialMixed(TPZMultiphysicsCompMesh *cmesh_mixed, ProblemConfig config, ErrorData &eData){

    int dim = config.gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    if(!eData.isMultiK) {
        cmesh_mixed->SetDefaultOrder(config.porder);
        cmesh_mixed->SetDimModel(dim);
        cmesh_mixed->SetAllCreateFunctionsMultiphysicElem();

        TPZMixedPoisson *material = new TPZMixedPoisson(matID, dim); //Using standard PermealityTensor = Identity.
        material->SetForcingFunction(config.exact.operator*().ForcingFunction());
        material->SetForcingFunctionExact(config.exact.operator*().Exact());
        cmesh_mixed->InsertMaterialObject(material);

        //Boundary Conditions
        TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);

        TPZMaterial *BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2);
        BCond0->SetForcingFunction(config.exact.operator*().Exact());

        TPZMaterial *BCond1 = material->CreateBC(material, -2, neumann, val1, val2);

        cmesh_mixed->InsertMaterialObject(BCond0);
        cmesh_mixed->InsertMaterialObject(BCond1);
    }

    else {
        CreateMaterialMultiK_Mixed(cmesh_mixed, eData.perm_Q1, eData.perm_Q2, config);
    }
}

void InsertMaterialObjectsH1Hybrid(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, ProblemConfig &config, ErrorData &eData)
{
    TPZGeoMesh *gmesh = cmesh_H1Hybrid->Reference();
    int dim = gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    // Creates Poisson material
    if(!eData.isMultiK) {
        TPZMatLaplacianHybrid *material = new TPZMatLaplacianHybrid(matID, dim);

        material->SetPermeability(1.);

        cmesh_H1Hybrid->InsertMaterialObject(material);
        if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
            material->SetForcingFunction(
                    config.exact.operator*().ForcingFunction());
            material->SetForcingFunctionExact(config.exact.operator*().Exact());
        }
        //    TPZMaterial * mat(material);
        //    cmesh->InsertMaterialObject(mat);

        // Inserts boundary conditions
        TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 1.);
        TPZMaterial *BCond0 =
                material->CreateBC(material, -1, dirichlet, val1, val2);
        if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
            BCond0->SetForcingFunction(config.exact.operator*().Exact());
        }
        val2.Zero();
        TPZMaterial *BCond1 = material->CreateBC(material, -2, neumann, val1, val2);

        cmesh_H1Hybrid->InsertMaterialObject(BCond0);
        cmesh_H1Hybrid->InsertMaterialObject(BCond1);
    }
    else {
        CreateMaterialMultiK_Hybrid(cmesh_H1Hybrid, eData.perm_Q1, eData.perm_Q2, config);
    }
}

void SolveH1Problem(TPZCompMesh *cmeshH1,struct ProblemConfig &config, struct ErrorData &eData){

    config.exact.operator*().fSignConvention = -1;

    std::cout << "Solving H1 " << std::endl;

    TPZAnalysis an(cmeshH1);

#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmeshH1);
    strmat.SetNumThreads(8);
    //        strmat.SetDecomposeType(ELDLt);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmeshH1);
    strmat.SetNumThreads(0);
    //        TPZSkylineStructMatrix strmat3(cmesh_HDiv);
    //        strmat3.SetNumThreads(8);
#endif

    std::set<int> matids;
    matids.insert(1);

    for(auto mat:config.bcmaterialids){
        matids.insert(mat);
    }

    strmat.SetMaterialIds(matids);
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();//resolve o problema misto ate aqui

    int64_t nelem = cmeshH1->NElements();
    cmeshH1->LoadSolution(cmeshH1->Solution());
    cmeshH1->ExpandSolution();
    cmeshH1->ElementSolution().Redim(nelem, 10);

    ////Calculo do erro
    std::cout << "Computing Error H1 " << std::endl;

    an.SetExact(config.exact.operator*().ExactSolution());

    if(eData.mode == 1) StockErrorsH1(an,cmeshH1,eData.ErroH1,eData.LogH1,eData);
    else StockErrorsH1(an,cmeshH1,eData.Erro,eData.Log,eData);

    ////PostProcess
    if(eData.last && eData.post_proc) {
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Solution");
        vecnames.Push("Derivative");
        scalnames.Push("ExactSolution");

        int dim = cmeshH1->Reference()->Dimension();

        std::string plotname;
        {
            std::stringstream out;
            out << eData.plotfile /* << config.dir_name*/ << "/" << "H1_Problem" << config.porder << "_" << dim
                << "D_" << config.problemname << "Ndiv_ " << config.ndivisions << ".vtk";
            plotname = out.str();
        }
        int resolution=0;
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(resolution,dim);
    }

    std::cout << "FINISHED!" << std::endl;
}
void SolveHybridH1Problem(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,int InterfaceMatId, struct ProblemConfig config,struct ErrorData &eData,int hybridLevel){

    config.exact.operator*().fSignConvention = 1;

    std::cout << "Solving HYBRID_H1 " << std::endl;

    TPZAnalysis an(cmesh_H1Hybrid);

#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh_H1Hybrid);
    strmat.SetNumThreads(0);
    //        strmat.SetDecomposeType(ELDLt);
#else
    //    TPZFrontStructMatrix<TPZFrontSym<STATE> > strmat(Hybridmesh);
        //    strmat.SetNumThreads(2);
        //    strmat.SetDecomposeType(ELDLt);
        TPZSkylineStructMatrix strmat(cmesh_H1Hybrid);
        strmat.SetNumThreads(0);
#endif
    std::set<int> matIds;
    for (auto matid : config.materialids) matIds.insert(matid);
    for (auto matidbc : config.bcmaterialids) matIds.insert(matidbc);

    matIds.insert(InterfaceMatId);
    strmat.SetMaterialIds(matIds);
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE>* direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();

    //TPZSymetricSpStructMatrix sparse(cmesh_H1Hybrid);
    /*TPZSkylineStructMatrix skylstr(cmesh_H1Hybrid);
    an.SetStructuralMatrix(skylstr);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Run();*/

    int64_t nelem = cmesh_H1Hybrid->NElements();
    cmesh_H1Hybrid->LoadSolution(cmesh_H1Hybrid->Solution());
    cmesh_H1Hybrid->ExpandSolution();
    cmesh_H1Hybrid->ElementSolution().Redim(nelem, 5);

    ////Calculo do erro
    std::cout << "Computing Error HYBRID_H1 " << std::endl;

    an.SetExact(config.exact.operator*().ExactSolution());

    std::cout << "DOF = " << cmesh_H1Hybrid->NEquations() << std::endl;

    if(eData.mode == 1) StockErrors(an,cmesh_H1Hybrid,eData.ErroHybridH1,eData.LogHybridH1,eData);
    else StockErrors(an,cmesh_H1Hybrid,eData.Erro,eData.Log,eData);

    ////PostProcess

    if(/*eData.last &&*/ eData.post_proc) {
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Pressure");
        scalnames.Push("PressureExact");
        vecnames.Push("Flux");
        //vecnames.Push(/*"ExactFlux"*/"ExactFluxShiftedOrigin");

        int dim = 2;
        std::string plotname;
        std::string hybridType;
        if(hybridLevel == 1) hybridType = "_HybridH1_k-";
        else hybridType = "_HybridSquared_k-";
        {
            std::stringstream out;
            out << eData.plotfile /* << config.dir_name*/ << "/"
                << config.problemname << hybridType << config.porder-config.H1Hybridminus
                << "_Mais-" << config.H1Hybridminus << "_ref_" << 1/eData.h << " x " << 1/eData.h <<".vtk";
            plotname = out.str();
        }
        int resolution = 0;
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(resolution, dim);
    }
}

void SolveMixedProblem(TPZMultiphysicsCompMesh *cmesh_Mixed,struct ProblemConfig config,struct ErrorData &eData) {

    config.exact.operator*().fSignConvention = 1;
    bool optBW = true;

    std::cout << "Solving Mixed " << std::endl;
    TPZAnalysis an(cmesh_Mixed, optBW); //Cria objeto de análise que gerenciará a analise do problema

    //MKL solver
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh_Mixed);
    strmat.SetNumThreads(8);
    //        strmat.SetDecomposeType(ELDLt);
#else
    //    TPZFrontStructMatrix<TPZFrontSym<STATE> > strmat(Hybridmesh);
        //    strmat.SetNumThreads(2);
        //    strmat.SetDecomposeType(ELDLt);
        TPZSkylineStructMatrix strmat(cmesh_H1Hybrid);
        strmat.SetNumThreads(0);
#endif
    //std::set<int> matIds;
    //for (auto matid : config.materialids) matIds.insert(matid);
    //for (auto matidbc : config.bcmaterialids) matIds.insert(matidbc);

    //matIds.insert(InterfaceMatId);
    //strmat.SetMaterialIds(matIds);
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE>* direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();

    //Previous solver
    /*TPZSkylineStructMatrix matskl(cmesh_Mixed); //caso simetrico ***
    int numthreads = 0;
    matskl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Assemble();//Assembla a matriz de rigidez (e o vetor de carga) global
    an.Solve();*/

    ////Calculo do erro
    std::cout << "Computing Error MIXED " << std::endl;

    an.SetExact(config.exact.operator*().ExactSolution());

    std::cout << "DOF = " << cmesh_Mixed->NEquations() << std::endl;

    if(eData.mode == 1) StockErrors(an,cmesh_Mixed,eData.ErroMixed,eData.LogMixed,eData);
    else StockErrors(an,cmesh_Mixed,eData.Erro,eData.Log,eData);

    ////PostProcess
    if(/*eData.last &&*/ eData.post_proc) {

        int dim = config.gmesh->Dimension();
        std::string plotname;
        {
            std::stringstream out;
            out << eData.plotfile /* << config.dir_name*/ << "/"
                << config.problemname << "_Mixed_k-" << config.porder - config.hdivmais
                << "_Mais-" << config.hdivmais << "_ref-" << 1/eData.h <<" x " << 1/eData.h << ".vtk";
            plotname = out.str();
        }

        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Pressure");
        scalnames.Push("ExactPressure");
        vecnames.Push("Flux");
        vecnames.Push(/*"ExactFlux"*/"ExactFluxShiftedOrigin");

        int resolution = 0;
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(resolution, dim);
    }
}

void FlushTable(ErrorData &eData, char *argv[]){

    std::string plotname;
    plotname = eData.plotfile + "/" + eData.plotfile + ".csv";

    remove(plotname.c_str());
    ofstream table(plotname.c_str(), ios::app);

    ProblemConfig config;
    Configure(config, 0, eData, argv);

    if(eData.mode == 1) {
        eData.ErroH1.close();
        eData.ErroMixed.close();
        eData.ErroHybridH1.close();
        string file1 = "ErroH1.txt";
        CleanErrors(file1);
        string file2 = "ErroHybridH1.txt";
        CleanErrors(file2);
        string file3 = "ErroMixed.txt";
        CleanErrors(file3);

        InvertError(file2);
        InvertError(file3);

        table << "Geometry" << "," << "Quadrilateral" << "\n";
        table << "Refinement" << "," << "Uniform" << "\n";
        table << "domain" << "," << "[0 1]x[0 1]" << "\n";
        table << "Case" << "," << config.problemname << "\n";
        table << "Approximation" << "," << "H1" << "," << "HybridH1" << "," << "Mixed" << "\n";
        table << "Internal order" << "," << config.porder << "," << config.porder << "," << config.porder<< "\n";
        table << "External flux" << "," << "---" << "," << config.porder - config.H1Hybridminus << "," << config.porder - config.hdivmais << "\n\n";
        table << "Norm" << "," << "H1" << "," << "HybridH1" << "," << "Mixed" << "\n";

        FillErrors(table, file1, file2, file3);
    }

    else{
        eData.Erro.close();
        string file = "Erro.txt";
        CleanErrors(file);

        if(eData.mode == 3 || eData.mode == 4 || eData.mode == 5) InvertError(file);

        table << "Geometry" << "," << "Quadrilateral" << "\n";
        table << "Refinement" << "," << "Uniform" << "\n";
        table << "domain" << "," << "[0 1]x[0 1]" << "\n";
        table << "Case" << "," << config.problemname << "\n";
        switch(eData.mode) {
            case 2:
                table << "Approximation" << "," << "H1" << "\n";
                table << "Internal order"  << "," << config.porder << "\n";
                table << "---" << "," << "---" <<  "\n\n";
                table << "Norm" << "," << "H1" << "\n";
                break;
            case 3:
                table << "Approximation" << "," << "HybridH1" << "\n";
                table << "Internal order"  << "," << config.porder << "\n";
                table << "External flux" << "," << config.porder - config.H1Hybridminus <<  "\n\n";
                table << "Norm" << "," << "HybridH1" << "\n";
                break;
            case 4:
                table << "Approximation" << "," << "Mixed" << "\n";
                table << "Internal order" << "," << config.porder << "\n";
                table << "External flux" << "," << config.porder - config.hdivmais <<  "\n\n";
                table << "Norm" << "," << "Mixed" << "\n";
                break;
            case 5:
                table << "Approximation" << "," << "HybridSquared" << "\n";
                table << "Internal order"  << "," << config.porder << "\n";
                table << "External flux" << "," << config.porder - config.H1Hybridminus <<  "\n\n";
                table << "Norm" << "," << "HybridSquared" << "\n";
                break;
        }

        FillErrors(table, file, eData.mode);

    }

    table.close();
}

void StockErrorsH1(TPZAnalysis &an,TPZCompMesh *cmesh, ofstream &Erro, TPZVec<REAL> *Log,ErrorData &eData){

    TPZManVector<REAL,6> Errors;
    Errors.resize(eData.numErrors);
    bool store_errors = false;

    an.PostProcessError(Errors, store_errors, Erro);

    if ((*Log)[0] != -1) {
        for (int j = 0; j < 3; j++) {
            (*eData.rate)[j] =
                    (log10(Errors[j]) - log10((*Log)[j])) /
                    (log10(eData.h) - log10(eData.hLog));
            Erro << "rate " << j << ": " << (*eData.rate)[j] << std::endl;
        }
    }

    Erro << "h = " << eData.h << std::endl;
    Erro << "DOF = " << cmesh->NEquations() << std::endl;
    for (int i = 0; i < eData.numErrors; i++)
        (*Log)[i] = Errors[i];
    Errors.clear();
}

void StockErrors(TPZAnalysis &an,TPZMultiphysicsCompMesh *cmesh, ofstream &Erro, TPZVec<REAL> *Log,ErrorData &eData){

    TPZManVector<REAL,6> Errors;
    Errors.resize(eData.numErrors);
    bool store_errors = false;

    an.PostProcessError(Errors, store_errors, Erro);

    if ((*Log)[0] != -1) {
        for (int j = 0; j < 3; j++) {
            (*eData.rate)[j] =
                    (log10(Errors[j]) - log10((*Log)[j])) /
                    (log10(eData.h) - log10(eData.hLog));
            Erro << "rate " << j << ": " << (*eData.rate)[j] << std::endl;
        }
    }

    Erro << "h = " << eData.h << std::endl;
    Erro << "DOF = " << cmesh->NEquations() << std::endl;
    for (int i = 0; i < eData.numErrors; i++)
        (*Log)[i] = Errors[i];
    Errors.clear();
}

void InvertError(string file){
    std::vector<std::string> erro,rate;
    std::string sErro,sRate;
    int size = 0;

    std::ifstream iErro(file);
    std::ofstream temp("temp.txt");

    int it_count = -1, hash_count = 0;
    string Line;
    while(getline(iErro,Line)){
        it_count++;
        if (it_count == 0) {
            temp << Line << endl;
            continue;
        }

        hash_count++;
        if (Line.find("#") != string::npos) hash_count = 0;

        switch (hash_count) {
            case (1):
                sErro = Line;
                temp << Line << endl;
                continue;
            case (2):
                temp << sErro << endl;
                erro.push_back(Line);
                break;
            case (4):
                if (it_count == 5)
                    temp << Line << endl;
                else{
                    sRate = Line;
                    temp << Line << endl;
                }
                break;
            case (5):
                if (it_count == 6)
                    temp << Line << endl;
                else {
                    temp << sRate << endl;
                    rate.push_back(Line);
                }
                break;
            default:
                temp << Line << endl;
                break;
        }
    }
    temp.close();

    remove(file.c_str());
    std::ifstream itemp("temp.txt");
    std::ofstream Erro(file.c_str());

    it_count = -1; hash_count =0;
    int erro_counter =0, rate_counter =0;
    while(getline(itemp,Line)){
        it_count++;
        if (it_count == 0) {
            Erro << Line << endl;
            continue;
        }

        hash_count++;
        if (Line.find("#") != string::npos) hash_count = 0;

        switch (hash_count) {
            case (1):
                Erro << erro[erro_counter] << endl;
                erro_counter++;
                break;
            case (4):
                if (it_count == 5)
                    Erro << Line << endl;
                else{
                    Erro << rate[rate_counter] << endl;
                    rate_counter++;
                }
                break;
            default:
                Erro << Line << endl;
                continue;
        }
    }
    itemp.close();
    remove("temp.txt");
}

void FillErrors(ofstream &table,string f,int mode){
    std::ifstream iErro(f.c_str());

    int it_count = -1, hash_count = 0;
    string Line;
    while(getline(iErro,Line)) {
        it_count++;

        if (it_count == 0) {
            continue;
        }
        hash_count++;
        if (Line.find("#") != string::npos) hash_count = 0;

        FillLegend(table,hash_count,it_count);

        if(hash_count == 0) table <<"\n";
        else table << Line << "\n";
    }
}

void FillErrors(ofstream &table,string f1,string f2,string f3){

    std::ifstream iErroH1(f1.c_str());
    std::ifstream iErroHybridH1(f2.c_str());
    std::ifstream iErroMixed(f3.c_str());

    int it_count = -1, hash_count = 0;
    string Line0,Line1,Line2;
    while(getline(iErroH1,Line0)) {
        it_count++;
        getline(iErroHybridH1, Line1);
        getline(iErroMixed, Line2);

        if (it_count == 0) {
            continue;
        }
        hash_count++;
        if (Line0.find("#") != string::npos) hash_count = 0;

        FillLegend(table,hash_count,it_count);

        if(hash_count == 0) table <<"\n";
        else table << Line0 << "," << Line1 << "," << Line2 << "\n";
    }
}

void FillLegend(ofstream &table,int hash_count,int it_count){
    switch (hash_count) {
        case (0):
            table << "\n";
            break;
        case (1):
            table << "semiH1-error" << ",";
            break;
        case (2):
            table << "L2-error" << ",";
            break;
        case (3):
            table << "3rd-error" << ",";
            break;
        case (4):
            if (it_count == 5)
                table << "h" << ",";
            else
                table << "semiH1-rate" << ",";
            break;
        case (5):
            if (it_count == 6)
                table << "DOF" << ",";
            else
                table << "L2-rate" << ",";
            break;
        case (6):
            table << "3rd-rate" << ",";
            break;
        case (7):
            table << "h" << ",";
            break;
        case (8):
            table << "DOF" << ",";
            break;
        default:
            break;
    }
}

void CleanErrors(string file){
    std::ifstream iErro(file);
    std::ofstream temp("temp.txt");
    size_t last_index;
    string Line, residue, other = "other";


    int counter = -1;
    while (getline(iErro, Line)) {
        size_t found = Line.find(other);
        if (found != string::npos) continue;

        last_index = Line.find("=");

        if (last_index == string::npos) {
            last_index = Line.find(":");
            if(last_index == string::npos)
                temp << Line << endl;
            else{
                residue = Line.substr(last_index+2);
                temp << residue <<endl;
            }
        }
        else {
            residue = Line.substr(last_index+2);
            temp << residue <<endl;
        }
    }
    iErro.close(); temp.close();
    remove(file.c_str());
    rename("temp.txt", file.c_str());
}

void EvaluateEntry(int argc, char *argv[],ErrorData &eData){
    if(argc != 1 && argc != 5){
        std::cout << "Invalid entry";
        DebugStop();
    }
    if(argc == 5){
        eData.argc = argc;
        for(int i = 2; i < 5 ; i++)
            IsInteger(argv[i]);
        if(std::strcmp(argv[1], "All") == 0)
            eData.mode = 1;
        if(std::strcmp(argv[1], "H1") == 0)
            eData.mode = 2;
        if(std::strcmp(argv[1], "HybridH1") == 0)
            eData.mode = 3;
        if(std::strcmp(argv[1], "Mixed") == 0)
            eData.mode = 4;
    }
    if(argc == 1){
        std::cout << "The polynomial order used here is defined within the code.\n"
                     "One can also define the polynomial order at the command line.\n"
                     "For that, insert 3 integer values after the executable command:\n"
                     "mode --- Polynomial Order --- H1Hybridminus --- HdivMais\n"
                     "mode = 'All', 'H1', 'HybridH1' or 'Mixed'\n"
                     "The command line entry is useful for consecutive runs as one can\n"
                     "build a .sh document to run all of them.";
    }
}

void IsInteger(char *argv){
    std::istringstream ss(argv);
    int x;
    if (!(ss >> x)) {
        std::cerr << "Invalid number: " << argv << '\n';
    } else if (!ss.eof()) {
        std::cerr << "Trailing characters after number: " << argv << '\n';
    }
}