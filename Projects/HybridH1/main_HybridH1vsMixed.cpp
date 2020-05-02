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

bool neumann = true;

struct ErrorData
{
    std::ofstream ErroH1,ErroHybridH1,ErroMixed;
    TPZVec<REAL> *LogH1,*LogHybridH1,*LogMixed, *rate;
    int maxdiv = 5;
    int orderlagrange = 1;
    REAL hLog = -1, h = -1000;
    int numErrors = 4;

    bool last = false, post_proc = true;
    int exp = 2; // Initial exponent of mesh refinement (numElem = 2*2^exp)
};

void InsertMaterialObjectsH1Hybrid(TPZMultiphysicsCompMesh *cmesh, ProblemConfig &config);
void BuildANDinsertBCflux(TPZCompMesh *cmesh_flux, ProblemConfig &config);
void BuildANDinsertBCpotential(TPZCompMesh *cmesh_p, ProblemConfig &config);
void InsertMaterialMixed(TPZMultiphysicsCompMesh *cmesh_mixed, ProblemConfig config);

int CreateHybridH1ComputationalMesh(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,ErrorData &eData, ProblemConfig &config);
void CreateMixedComputationalMesh(TPZMultiphysicsCompMesh *cmesh_H1Mixed,ErrorData &eData, ProblemConfig &config);
void SolveH1Problem(TPZCompMesh *cmeshH1,struct ProblemConfig &config, struct ErrorData &eData);
void SolveHybridH1Problem(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,int InterfaceMatId,struct ProblemConfig config, struct ErrorData &eData);
void SolveMixedProblem(TPZMultiphysicsCompMesh *cmesh_Mixed,struct ProblemConfig config,struct ErrorData &eData);

void InitializeOutstream(ErrorData &eData);
void FlushTable(ErrorData &eData);
void CleanErrors(string file);
void InvertError(string file);
void FillErrors(ofstream &table,string f1,string f2,string f3);

void Configure(ProblemConfig &config,int ndiv){
    config.porder = 2;
    config.hdivmais = 1;
    config.ndivisions = ndiv;
    config.dimension = 2;
    config.prefine = false;

    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
    config.exact.operator*().fSignConvention = 1;

    config.problemname = "ESinSin";

    config.dir_name = "HybridH1_ESinSin";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    // geometric mesh
    TPZManVector<int, 4> bcids(4, -1);
    TPZGeoMesh *gmesh = CreateGeoMesh(1, bcids);
    UniformRefinement(config.ndivisions, gmesh);

    config.gmesh = gmesh;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
}

int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    ErrorData eData;
    InitializeOutstream(eData);

    for (int ndiv = 1; ndiv < eData.maxdiv+2; ndiv++) {
        if (ndiv == eData.maxdiv+1) eData.last = true;
        eData.h = 1./eData.exp;

        ProblemConfig config;
        Configure(config,ndiv);

#ifdef PZDEBUG
        if(eData.last && eData.post_proc) {
            std::ofstream out("gmesh.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, out);
            std::ofstream out2("gmeshInitial.txt");
            config.gmesh->Print(out2);
        }
#endif
        //H1
        TPZCompMesh *cmeshH1 = CMeshH1(config);
        if(eData.last && eData.post_proc) {ofstream arg1("CompMeshH1.txt"); cmeshH1->Print(arg1);}
        SolveH1Problem(cmeshH1, config, eData);

        //HybridH1
        TPZMultiphysicsCompMesh *cmesh_H1Hybrid = new TPZMultiphysicsCompMesh(config.gmesh);
        int lagrangeId = CreateHybridH1ComputationalMesh(cmesh_H1Hybrid, eData, config);
        SolveHybridH1Problem(cmesh_H1Hybrid,lagrangeId,config,eData);

        //Mixed
        TPZMultiphysicsCompMesh *cmesh_mixed = new TPZMultiphysicsCompMesh(config.gmesh);
        CreateMixedComputationalMesh(cmesh_mixed,eData,config);
        SolveMixedProblem(cmesh_mixed,config,eData);

#ifdef PZDEBUG
        {
            std::ofstream out2("H1HybridMesh.txt");
            cmesh_H1Hybrid->Print(out2);
            std::ofstream out3("gmeshHybridH1.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, out3);
        }
#endif

    eData.hLog = eData.h;
    eData.exp *=2;
    }

    //Not ready yet
    FlushTable(eData);

    return 0.;
}

void CreateMixedComputationalMesh(TPZMultiphysicsCompMesh *cmesh_Mixed, ErrorData &eData, ProblemConfig &config){

    int matID = 1;
    int dim = config.gmesh->Dimension();

    //Flux mesh creation
    TPZCompMesh *cmesh_flux = new TPZCompMesh(config.gmesh);
    BuildANDinsertBCflux(cmesh_flux, config);

    //Potential mesh creation
    TPZCompMesh *cmesh_p = new TPZCompMesh(config.gmesh);
    BuildANDinsertBCpotential(cmesh_p, config);

    //Multiphysics mesh build
    InsertMaterialMixed(cmesh_Mixed, config);
    TPZManVector<int> active(2, 1);
    TPZManVector<TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmesh_flux;
    meshvector[1] = cmesh_p;
    TPZCompMeshTools::AdjustFluxPolynomialOrders(meshvector[0], config.hdivmais);

    cmesh_Mixed->BuildMultiphysicsSpace(active,meshvector);
    cmesh_Mixed->LoadReferences();
    cmesh_Mixed->InitializeBlock();
}

int CreateHybridH1ComputationalMesh(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, ErrorData &eData, ProblemConfig &config){
    TPZCreateMultiphysicsSpace createspace(config.gmesh);

    createspace.SetMaterialIds({1}, {-2,-1});
    createspace.fH1Hybrid.fHybridizeBC = false;//opcao de hibridizar o contorno
    createspace.ComputePeriferalMaterialIds();

    TPZManVector<TPZCompMesh *> meshvec;

    createspace.CreateAtomicMeshes(meshvec,config.porder,eData.orderlagrange);

    InsertMaterialObjectsH1Hybrid(cmesh_H1Hybrid, config);
    createspace.InsertPeriferalMaterialObjects(cmesh_H1Hybrid);
    cmesh_H1Hybrid->BuildMultiphysicsSpace(meshvec);

    createspace.AddInterfaceElements(cmesh_H1Hybrid);
    createspace.GroupandCondenseElements(cmesh_H1Hybrid);

    cmesh_H1Hybrid->InitializeBlock();
    cmesh_H1Hybrid->ComputeNodElCon();

    return createspace.fH1Hybrid.fLagrangeMatid.first;
}

void InitializeOutstream(ErrorData &eData){
    //Error buffer
    if( remove( "ErroH1.txt" ) != 0 || remove( "ErroHybridH1.txt" ) != 0 || remove( "ErroMixed.txt" ) != 0 )
        perror( "Error deleting file" );
    else
        puts( "Error log successfully deleted" );

    //Outstream file
    eData.ErroH1.open("ErroH1.txt",std::ofstream::app);
    eData.ErroHybridH1.open("ErroHybridH1.txt",std::ofstream::app);
    eData.ErroMixed.open("ErroMixed.txt",std::ofstream::app);

    eData.ErroH1 <<       "----------COMPUTED ERRORS (H1)----------\n"; eData.ErroH1.flush();
    eData.ErroHybridH1 << "-------COMPUTED ERRORS (HYBRID_H1)------\n"; eData.ErroHybridH1.flush();
    eData.ErroMixed <<    "--------COMPUTED ERRORS (MIXED)---------\n"; eData.ErroMixed.flush();

    //Log Initialization
    eData.LogH1 = new TPZVec<REAL>(eData.numErrors, -1);
    eData.LogHybridH1 = new TPZVec<REAL>(eData.numErrors, -1);
    eData.LogMixed = new TPZVec<REAL>(eData.numErrors, -1);
    eData.rate = new TPZVec<REAL>(eData.numErrors, -1);
}

void BuildANDinsertBCflux(TPZCompMesh *cmesh_flux, ProblemConfig &config){

    int dim = config.gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    cmesh_flux->SetDefaultOrder(config.porder);
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

    cmesh_flux->AutoBuild();
    cmesh_flux->InitializeBlock();
}

void BuildANDinsertBCpotential(TPZCompMesh *cmesh_p, ProblemConfig &config){
    int matID = 1;
    int dim = config.gmesh->Dimension();

    cmesh_p->SetDefaultOrder(config.porder);
    cmesh_p->SetDimModel(dim);

    cmesh_p->SetAllCreateFunctionsContinuous(); //H1 functions
    cmesh_p->ApproxSpace().CreateDisconnectedElements(true);

    TPZNullMaterial *material = new TPZNullMaterial(matID); material->SetDimension(dim);
    cmesh_p->InsertMaterialObject(material);

    //cmesh_p->SetDefaultOrder(config.porder + config.hdivmais);

    cmesh_p->AutoBuild();
    cmesh_p->ExpandSolution();

    TPZAdmChunkVector<TPZConnect> &nodeIt = cmesh_p->ConnectVec();
    for(auto &nodo : nodeIt){
        nodo.SetLagrangeMultiplier(1);
    }
}

void InsertMaterialMixed(TPZMultiphysicsCompMesh *cmesh_mixed, ProblemConfig config){

    int dim = config.gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    cmesh_mixed -> SetDefaultOrder(config.porder);
    cmesh_mixed -> SetDimModel(dim);
    cmesh_mixed->SetAllCreateFunctionsMultiphysicElem();

    TPZMixedPoisson *material = new TPZMixedPoisson(matID,dim); //Using standard PermealityTensor = Identity.
    material->SetForcingFunction(config.exact.operator*().ForcingFunction());
    material->SetForcingFunctionExact(config.exact.operator*().Exact());
    cmesh_mixed->InsertMaterialObject(material);

    //Boundary Conditions
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);

    TPZMaterial *BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2);
    BCond0->SetForcingFunction(config.exact.operator*().Exact());

    TPZMaterial *BCond1 = material->CreateBC(material, -2, neumann, val1, val2);

    cmesh_mixed->InsertMaterialObject(BCond0);
    cmesh_mixed->InsertMaterialObject(BCond1);
}

void InsertMaterialObjectsH1Hybrid(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, ProblemConfig &config)
{
    TPZGeoMesh *gmesh = cmesh_H1Hybrid->Reference();
    int dim = gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    // Creates Poisson material
    TPZMatLaplacianHybrid *material = new TPZMatLaplacianHybrid(matID, dim);

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

void SolveH1Problem(TPZCompMesh *cmeshH1,struct ProblemConfig &config, struct ErrorData &eData){

    config.exact.operator*().fSignConvention = -1;

    TPZAnalysis an(cmeshH1);

#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmeshH1);
    strmat.SetNumThreads(0);
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
    TPZManVector<REAL,4> Errors;
    Errors.resize(eData.numErrors);

    an.SetExact(config.exact.operator*().ExactSolution());

    bool store_errors = false;
    an.PostProcessError(Errors, store_errors, eData.ErroH1);

    eData.LogH1->Print(std::cout);
    if((*eData.LogH1)[0] !=-1){
        for(int j =0; j < 3 ; j++) {
            (*eData.rate)[j] = (log10(Errors[j]) - log10((*eData.LogH1)[j]))/(log10(eData.h)-log10(eData.hLog));
            eData.ErroH1 << "rate " << j << ": " << (*eData.rate)[j] << std::endl;
        }
    }

    eData.ErroH1 << "h = " << eData.h << std::endl;
    eData.ErroH1 << "DOF = " << cmeshH1->NEquations() << std::endl;
    for(int i =0; i < eData.numErrors;i++) (*eData.LogH1)[i] = Errors[i];
    Errors.clear();

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
            out << config.dir_name << "/" << "H1_Problem" << config.porder << "_" << dim
                << "D_" << config.problemname << "Ndiv_ " << config.ndivisions << ".vtk";
            plotname = out.str();
        }
        int resolution=0;
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(resolution,dim);
    }

    std::cout << "FINISHED!" << std::endl;
}
void SolveHybridH1Problem(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,int InterfaceMatId, struct ProblemConfig config,struct ErrorData &eData){

    config.exact.operator*().fSignConvention = 1;

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
    TPZManVector<REAL,4> Errors;
    Errors.resize(eData.numErrors);

    an.SetExact(config.exact.operator*().ExactSolution());

    bool store_errors = true;
    an.PostProcessError(Errors, store_errors, eData.ErroHybridH1); std::cout << Errors[0] << "   " << Errors[0] << "   " << Errors[0] << "   ";

    if((*eData.LogHybridH1)[0] !=-1){
        for(int j =0; j < 3 ; j++) {
            (*eData.rate)[j] = (log10(Errors[j]) - log10((*eData.LogHybridH1)[j]))/(log10(eData.h)-log10(eData.hLog));
            eData.ErroHybridH1 << "rate " << j << ": " << (*eData.rate)[j] << std::endl;
        }
    }

    eData.ErroHybridH1 << "h = " << eData.h << std::endl;
    eData.ErroHybridH1 << "DOF = " << cmesh_H1Hybrid->NEquations() << std::endl;
    for(int i =0; i < eData.numErrors;i++) (*eData.LogHybridH1)[i] = Errors[i];
    Errors.clear();

    ////PostProcess

    if(eData.last && eData.post_proc) {
        TPZStack<std::string> scalnames, vecnames;
        //scalnames.Push("Pressure");
        scalnames.Push("PressureExact");

        int dim = 2;
        std::string plotname;
        {
            std::stringstream out;
            out << config.dir_name << "/"
                << "HybridH1" << config.porder << "_" << dim << "D_"
                << config.problemname << "Ndiv_ " << config.ndivisions
                << "HdivMais" << config.hdivmais << ".vtk";
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

    TPZAnalysis an(cmesh_Mixed, optBW); //Cria objeto de análise que gerenciará a analise do problema
    TPZSkylineStructMatrix matskl(cmesh_Mixed); //caso simetrico ***
    int numthreads = 0;
    matskl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Assemble();//Assembla a matriz de rigidez (e o vetor de carga) global
    an.Solve();

    ////Calculo do erro
    std::cout << "Computing Error MIXED " << std::endl;
    TPZManVector<REAL,6> Errors;
    Errors.resize(eData.numErrors);

    an.SetExact(config.exact.operator*().ExactSolution());

    bool store_errors = false;
    an.PostProcessError(Errors, store_errors, eData.ErroMixed);

    if((*eData.LogMixed)[0] !=-1){
        for(int j =0; j < 3 ; j++) {
            (*eData.rate)[j] = (log10(Errors[j]) - log10((*eData.LogMixed)[j]))/(log10(eData.h)-log10(eData.hLog));
            eData.ErroMixed << "rate " << j << ": " << (*eData.rate)[j] << std::endl;
        }
    }

    eData.ErroMixed << "h = " << eData.h << std::endl;
    eData.ErroMixed << "DOF = " << cmesh_Mixed->NEquations() << std::endl;
    for(int i =0; i < eData.numErrors;i++) (*eData.LogMixed)[i] = Errors[i];
    Errors.clear();

    ////PostProcess
    if(eData.last && eData.post_proc) {

        int dim = config.gmesh->Dimension();
        std::string plotname;
        {
            std::stringstream out;
            out << config.dir_name << "/"
                << "Mixed" << config.porder << "_" << dim << "D_"
                << config.problemname << "Ndiv_ " << config.ndivisions
                << "HdivMais" << config.hdivmais << ".vtk";
            plotname = out.str();
        }

        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Pressure");
        scalnames.Push("ExactPressure");
        vecnames.Push("Flux");
        vecnames.Push("ExactFlux");

        int postProcessResolution = 4; //  keep low as possible

        an.DefineGraphMesh(dim,scalnames,vecnames,plotname);
        an.PostProcess(postProcessResolution,dim);

        int resolution = 0;
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(resolution, dim);
    }
}

//Incomplete
void FlushTable(ErrorData &eData){
    eData.ErroH1.close(); eData.ErroMixed.close(); eData.ErroHybridH1.close();
    string file1 = "ErroH1.txt"; CleanErrors(file1);
    string file2 = "ErroHybridH1.txt"; CleanErrors(file2);
    string file3 = "ErroMixed.txt"; CleanErrors(file3);

    InvertError(file2);
    InvertError(file3);

    ProblemConfig config;
    Configure(config,0);
    std::string plotname;
    {
        std::stringstream out;
        std::string command = "mkdir ComparisonTable";
        system(command.c_str());

        out << "ComparisonTable" << "/" << "Quad_"<<config.problemname <<"___porder_" << config.porder
            << "___hdivmais_" << config.hdivmais << "___lagrangeOrder_" <<
            eData.orderlagrange <<".csv";
        plotname = out.str();
    }
    std::cout << plotname;
    remove(plotname.c_str());
    ofstream table(plotname.c_str(),ios::app);

    table << "Geometry" << ","<< "Quadrilateral" <<"\n";
    table << "Refinement" << "," << "Uniform" <<"\n";
    table << "domain" << "," << "[0 1]x[0 1]" <<"\n";
    table << "Case" << "," << config.problemname << "\n";
    table <<"Approximation" << "," << "H1" << "," << "HybridH1" << "," << "Mixed" <<"\n";
    table << "p-order" << "," << config.porder << "," << config.porder << "," << config.porder << "\n";
    table << "lagrange/hdivmais" << "," << "---" << "," << eData.orderlagrange << "," << config.hdivmais << "\n\n";
    table << "Norm" << "," << "H1" << "," << "HybridH1" << "," << "Mixed" << "\n";

    FillErrors(table,file1,file2,file3);

    table.close();
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

        switch (hash_count) {
        case (0):
            table << "\n";
            continue;
        case (1):
            table << "H1-error" << ",";
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
                table << "H1-rate" << ",";
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
        table << Line0 << "," << Line1 << "," << Line2 << "\n";
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