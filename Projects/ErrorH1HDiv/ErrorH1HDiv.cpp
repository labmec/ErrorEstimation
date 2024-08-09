/**
 * @file This file implements an error estimator in space H1.
 */

#include "pzlog.h"
#include "TPZLinearAnalysis.h"
#include "pzstepsolver.h"
#include "TPZBndCondT.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include <TPZNullMaterial.h>
#include <pzbuildmultiphysicsmesh.h>
#include <TPZMultiphysicsCompMesh.h>
#include "pzgeoelbc.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSYSMPMatrix.h"
#include "TPZRefPatternDataBase.h"
#include "TPZVTKGeoMesh.h"
#include "TPZGmshReader.h"
#include "TPZAnalyticSolution.h"
#include "TPZFrontSym.h"
#include "pzcheckgeom.h"
#include <TPZSimpleTimer.h>
#include "TPZVTKGenerator.h"
#include "Projection/TPZL2ProjectionCS.h"            
#include "TPZHDivApproxCreator.h"       

#include <fstream>
#include <ctime>
#include <cstdio>
#include <cmath>

//#include "pzgengrid.h"
#include "TPZGenGrid2D.h"


#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "TPZGeoLinear.h"

#include "tpzblendnaca.h"
#include "tpznacaprofile.h"

#include <iostream>

using namespace pzgeom;
/// @brief verify is the derivative of the NACA coordinate is correct
/// @param naca profile object
/// @param point parametric coordinate around which the derivative will be verified
void VerifyDerivative(TPZBlendNACA &naca, REAL point);

/// @brief Read a gmsh file and return a geometric mesh
/// @param filename 
/// @return a geometric mesh
TPZGeoMesh *ReadGmshSimple(const std::string &filename);

/// @brief Create a computational mesh with H1 elements
TPZCompMesh *CreateH1CompMeshQuadMesh(TPZGeoMesh *gmesh, int64_t porder);

/// @brief Create a computational mesh with L2 elements
TPZCompMesh *CreateL2CompMesh(TPZGeoMesh *gmesh);

/// @brief Create the computational mesh with HDiv elements
TPZCompMesh *CreateHDivCompMesh(TPZGeoMesh *gmesh);

TPZMultiphysicsCompMesh* CreateHDivMultiphysicsCompQuadMesh(TPZGeoMesh* gmesh, const int porder);


/// @brief print the results of the analysis
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh);


void SolveSyst(TPZLinearAnalysis &an, TPZCompMesh *cmesh, DecomposeType dtype);


enum EMatId {ENone = 0, EDOMAIN = 1, EBCL = 2, EBCR = 3, EBCT = 4, EBCD = 5};

REAL ComputeDifference(TPZMultiphysicsCompMesh *cmeshHDiv, TPZCompMesh *cmeshH1);

TLaplaceExample1 laplaceex1;

struct ErrorData {
    int grau;
    REAL erro;
};


int main() {

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    TPZGeoMesh *gmesh = ReadGmshSimple("quadmesh.msh");
    // {
    //     std::ofstream out("gmesh.txt");
    //     gmesh->Print(out);

    //     std::ofstream out2("gmesh.vtk"); 
    //     TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
    // }
    std::map<int, TLaplaceExample1::EExactSol> type_exact_sol;
    type_exact_sol[0] = TLaplaceExample1::EBubble2D;
    type_exact_sol[1] = TLaplaceExample1::ESinSin;
    type_exact_sol[2] = TLaplaceExample1::EArcTan;

    std::map<int, std::string> type_exact_sol_name;
    type_exact_sol_name[0] = "EBubble2D";
    type_exact_sol_name[1] = "ESinSin";
    type_exact_sol_name[2] = "EArcTan";
    
    
    // int refinements[] = {1}; 
    TPZVec<int> refinements = {1,2,3};
    TPZVec<int> porders = {1}; 
    // std::cout << "Size: " <<  porders.size() << std::endl;
    std::ofstream outfile("errestimates.txt");

    for (int i = 0; i < 1 ; i++) {
        laplaceex1.fExact = type_exact_sol[i];
        int porder = porders[0];

        outfile << type_exact_sol_name[i]  << std::endl;
            for (int j=0; j < refinements.size();j++) {
                TPZCheckGeom check(gmesh);
                check.UniformRefine(refinements[j]);

                // {
                //     std::ofstream out4("gmeshfine.txt");
                //     gmesh->Print(out4);

                //     std::ofstream out5("gmeshfine.vtk"); 
                //     TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out5);
                // }

            
                TPZMultiphysicsCompMesh *cmeshHDiv = nullptr;
                {
                    DecomposeType dtype = ELDLt;
                    cmeshHDiv = CreateHDivMultiphysicsCompQuadMesh(gmesh,porder);
                    TPZLinearAnalysis an(cmeshHDiv,RenumType::EDefault);
                    SolveSyst(an,cmeshHDiv,dtype);
                    PrintResults(an,cmeshHDiv);
                }
                TPZCompMesh *cmeshH1 = nullptr;
                {
                    DecomposeType dtype = ECholesky;
                    cmeshH1 = CreateH1CompMeshQuadMesh(gmesh,porder);
                    TPZLinearAnalysis an(cmeshH1,RenumType::EDefault);
                    SolveSyst(an,cmeshH1,dtype);
                    PrintResults(an,cmeshH1);
                }
                REAL errestimate = ComputeDifference(cmeshHDiv,cmeshH1);
                std::cout << "Error estimate " << errestimate << std::endl;
            // errestimates.push_back({refinements[i], errestimate});

            
                outfile << "Refinement Level: " << refinements[j] << ", Error: " << errestimate << std::endl;
                delete cmeshHDiv;
                delete cmeshH1;
                }
    }

    outfile.close();
    delete gmesh;
    return 0;
}


TPZGeoMesh *ReadGmshSimple(const std::string &meshfilename)
{
    TPZGmshReader gmsh;
    gmsh.SetVerbose(1);
    gmsh.GetDimNamePhysical()[2]["dom"] = EDOMAIN;
    gmsh.GetDimNamePhysical()[1]["bcL"] = EBCL;
    gmsh.GetDimNamePhysical()[1]["bcR"] = EBCR;
    gmsh.GetDimNamePhysical()[1]["bcT"] = EBCT;  
    gmsh.GetDimNamePhysical()[1]["bcD"] = EBCD;
    auto gmesh = gmsh.GeometricGmshMesh(meshfilename);

    return gmesh;

}



void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
//Define a função do tipo void chamada PrintResults, que recebe como parâmetros TPZLinearAnalysis &an e  TPZCompMesh *cmesh
{
 
    std::cout << "--------- Post Process ---------" << std::endl;
    //printa para o usuário "--------- Post Process ---------" indicando o início da fase de pós-processamento. 
    TPZSimpleTimer postProc("Post processing time");
    //declara uma variável chamada postProc, do tipo TPZSimpleTimer, chamando um construtor com uma string como argumento, igual a "Post processing time".
    //inicializa um temporizador chamado postProc que será usado para medir o tempo gasto no pós-processamento.
    TPZMultiphysicsCompMesh *mphysics = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);
    std::stringstream sout;
    sout << "Darcy_" << laplaceex1.Name() << "_p" << cmesh->GetDefaultOrder();
    if(mphysics) {
        sout << "_Hdiv";
    } else {
        sout << "_H1";
    }
    const std::string plotfile = sout.str();
    //define o nome base do arquivo de saída para o pós-processamento. O nome base é "postprocess".
    constexpr int vtkRes{2};
    //define a variável do tipo inteiro denominada vtkRes, do tipo constexpr, que significa que é uma expressão constante, ou seja,  vtkRes é um valor constante e não pode ser alterado. Ainda, {0} indica o valor associado a essa constante, e portanto não será alterado, com valor determinado na hora de compilação.
    //define a resolução para o formato de arquivo VTK. Neste caso, a resolução é definida como 0, o que geralmente significa que a resolução será automática.
    TPZVec<std::string> fields = {
        "Pressure",
        "Flux",
        "EstimatedError",
        "TrueError",
        "EffectivityIndex"
    };
    //nesse conjunto de linhas de código, temos que TPZVec é uma estrutura do tipo vetor que contém como argumento uma variável chamda "fields" que é uma lista de strings, que, pelo que se chamam, são relacionadas à pressão e ao fluxo.
    //cria um vetor de strings chamado fields que contém os nomes dos campos que serão pós-processados. Neste caso, os campos incluem "Pressure" (pressão) e "Flux" (fluxo). Esses campos representam propriedades do problema que desejamos visualizar após a simulação.
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    //essa linha de código declara uma variável chamada vtk do tipo auto, o que significa que o compilador irá deduzir o tipo que ela terá a depender do que ela é igual. No caso, ela é igual a função TPZVTKGenerator, de parâmetros cmesh, fields, plotfile, vtkRes.
    //cria um objeto vtk da classe TPZVTKGenerator, que é usado para gerar arquivos VTK a partir dos dados da malha computacional cmesh. Os argumentos passados para o construtor incluem a malha computacional, os campos a serem pós-processados, o nome base do arquivo de saída (plotfile) e a resolução VTK (vtkRes).
    vtk.SetNThreads(0);
    //define o número de threads a serem usadas durante o pós-processamento. A variável global_nthread provavelmente contém o número desejado de threads.
    vtk.Do();
    //inicia o processo de geração dos arquivos VTK. Esta função gera arquivos de saída contendo informações sobre os campos especificados na malha computacional.
    std::cout << "Total time = " << postProc.ReturnTimeDouble()/1000. << " s" << std::endl;
    //imprime o tempo gasto no pós-processamento, convertido para segundos.
    
    return;
    //a função é concluída e retorna.
}



void SolveSyst(TPZLinearAnalysis &an, TPZCompMesh *cmesh, DecomposeType dtype)
{
    TPZSkylineStructMatrix<STATE> strmat(cmesh);
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;

    step.SetDirect(dtype);
    an.SetSolver(step);
    an.Run();
    TPZMultiphysicsCompMesh *mphysics = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);
    int nerrors = 5;
    if(!mphysics) nerrors = 6;
    TPZManVector<REAL,6> errorsum(nerrors,0.);
    int64_t nelem = cmesh->NElements();
    cmesh->ElementSolution().Redim(nelem,nerrors);
    cmesh->EvaluateError(true, errorsum);
    std::cout << "Errors associated with H1 and HDiv space\n" << errorsum << std::endl;
    int errindex = 2;
    if(mphysics) errindex = 1;
    TPZFMatrix<STATE> &elsol = cmesh->ElementSolution();
    for(int64_t el = 0; el < nelem; el++)
    {
        elsol(el,1) = elsol(el,errindex);
    }
}

TPZMultiphysicsCompMesh* CreateHDivMultiphysicsCompQuadMesh(TPZGeoMesh* gmesh, const int porder) {
    TPZHDivApproxCreator approxCreator(gmesh);
    approxCreator.HdivFamily() = HDivFamily::EHDivStandard;
    approxCreator.ProbType() = ProblemType::EDarcy;
    approxCreator.SetDefaultOrder(porder);

    const int dim = approxCreator.GeoMesh()->Dimension();

    // Creating domain material
    TPZMixedDarcyFlow* matdarcy = nullptr;
    matdarcy = new TPZMixedDarcyFlow(EDOMAIN,dim);
    matdarcy->SetConstantPermeability(1.);
    matdarcy->SetForcingFunction(laplaceex1.ForceFunc(),3);
    matdarcy->SetExactSol(laplaceex1.ExactSolution(),3);
    approxCreator.InsertMaterialObject(matdarcy);

    // ========> Boundary Conditions
    // -----------------------------
    TPZBndCondT<STATE> *BCond1 = nullptr, *BCond2 = nullptr, *BCond3 = nullptr, *BCond4 = nullptr;
    const int dirType = 0, neuType = 1;

    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,1.);
    BCond1 = matdarcy->CreateBC(matdarcy, EBCL, dirType, val1, val2);
    BCond1->SetForcingFunctionBC(laplaceex1.ExactSolution(),3);
    val2[0] = 0.;
    BCond2 = matdarcy->CreateBC(matdarcy, EBCR, dirType, val1, val2);
    BCond2->SetForcingFunctionBC(laplaceex1.ExactSolution(),3);
    BCond3 = matdarcy->CreateBC(matdarcy, EBCD, neuType, val1, val2);
    BCond3->SetForcingFunctionBC(laplaceex1.ExactSolution(),3);
    BCond4 = matdarcy->CreateBC(matdarcy, EBCT, neuType, val1, val2);
    BCond4->SetForcingFunctionBC(laplaceex1.ExactSolution(),3);
    
    
    if(BCond1) approxCreator.InsertMaterialObject(BCond1);
    if(BCond2) approxCreator.InsertMaterialObject(BCond2);
    if(BCond3) approxCreator.InsertMaterialObject(BCond3);
    if(BCond4) approxCreator.InsertMaterialObject(BCond4);
    
    TPZMultiphysicsCompMesh *cmesh = approxCreator.CreateApproximationSpace();
    return cmesh;
}

TPZCompMesh *CreateH1CompMeshQuadMesh(TPZGeoMesh *gmesh, int64_t porder) {

    // Mesh creation
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    int dim = 2;
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(porder);

    // Inserting darcy flow material in the mesh
    TPZDarcyFlow *material = new TPZDarcyFlow(EDOMAIN,dim);
    material->SetExactSol(laplaceex1.ExactSolution(),5);
    material->SetForcingFunction(laplaceex1.ForceFunc(),3);
    cmesh->InsertMaterialObject(material);

    // Data structure for bcs
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);

    // Inserting boundary conditions
    // Left of domain
    const int dirichlet = 0, neumann = 1;
    val2[0] = 1;
    auto bndL = material->CreateBC(material, EBCL, neumann, val1, val2);
    bndL->SetForcingFunctionBC(laplaceex1.ExactSolution(),1);
    cmesh->InsertMaterialObject(bndL);

    // Right of domain
    val2[0] = 0;
    auto bndR = material->CreateBC(material, EBCR, dirichlet, val1, val2);
    bndR->SetForcingFunctionBC(laplaceex1.ExactSolution(),1);
    cmesh->InsertMaterialObject(bndR);

    // Down (bottom) of domain
    auto bndD = material->CreateBC(material, EBCD, dirichlet, val1, val2);
    bndD->SetForcingFunctionBC(laplaceex1.ExactSolution(),1);
    cmesh->InsertMaterialObject(bndD);

    // Top of domain
    auto bndT = material->CreateBC(material, EBCT, dirichlet, val1, val2);
    bndT->SetForcingFunctionBC(laplaceex1.ExactSolution(),1);
    cmesh->InsertMaterialObject(bndT);  

    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();

    return cmesh;
}

#include "pzintel.h"
#include "pzvec_extras.h"
#include "pzcondensedcompel.h"

REAL ComputeDifference(TPZMultiphysicsCompMesh *cmeshHDiv, TPZCompMesh *cmeshH1)
{
    REAL error = 0.;
    REAL errorH1 = 0.;
    REAL errorHdiv = 0.;
    TPZGeoMesh *gmesh = cmeshHDiv->Reference();
    int dim = gmesh->Dimension();
    gmesh->ResetReference();
    cmeshH1->LoadReferences();
    int64_t nel = cmeshHDiv->NElements();

    for(int64_t el = 0; el < nel; el++)
    {
        TPZCompEl *cel = cmeshHDiv->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() != gmesh->Dimension()) continue;
        TPZCompEl *celH1 = gel->Reference();
        if(!cel || !celH1) continue;
        TPZCondensedCompEl *condense = dynamic_cast<TPZCondensedCompEl *>(cel);
        if(!condense) DebugStop();
        TPZMultiphysicsElement *intelHDiv = dynamic_cast<TPZMultiphysicsElement *>(condense->ReferenceCompEl());
        TPZInterpolatedElement *intelH1 = dynamic_cast<TPZInterpolatedElement *>(celH1);
        if( !intelH1) DebugStop();
        TPZManVector<STATE,3> dataHDiv(3,0.), dataH1(3,0.),diff(3,0.),x(3,0.),sol(1,0.);
        TPZFNMatrix<9,REAL> gradx(3,dim),gradsol(3,1);
        TPZFNMatrix<9,REAL> jac(dim,dim), jacinv(dim,dim), axes(3,dim);
        REAL detjac;
        TPZIntPoints &intpoints = intelHDiv->GetIntegrationRule();
        int np = intpoints.NPoints();
        REAL w;
        TPZManVector<REAL,3> pos(gmesh->Dimension(),0.);
        REAL errorEstimate = 0., ElerrorH1 = 0., ElerrorHDiv = 0.;

        for(int ip = 0; ip < np; ip++) {
            intpoints.Point(ip,pos,w);
            gel->X(pos,x);
            gel->GradX(pos,gradx);
            gel->Jacobian(gradx,jac,axes,detjac,jacinv);
            laplaceex1.ExactSolution()(x,sol,gradsol);
            intelHDiv->Solution(pos,1,dataHDiv);
            intelH1->Solution(pos,7,dataH1);
            std::cout << "HDiv " << dataHDiv << " H1 " << dataH1 << " gradsol ";
            for(int i=0; i<3; i++) {
                std::cout << gradsol(i,0) << " ";
            }
            std::cout << std::endl;
            REAL diffnorm = 0., diffH1 = 0., diffHDiv = 0.;
            for (int i=0; i<3; i++) {
                diff[i] = dataHDiv[i] - dataH1[i];
                diffnorm += diff[i]*diff[i];
                diff[i] = dataHDiv[i] + gradsol(i,0);
                diffHDiv += diff[i]*diff[i];
                diff[i] = dataH1[i] + gradsol(i,0);
                diffH1 += diff[i]*diff[i];
            }
            errorEstimate += w*detjac*diffnorm;
            ElerrorH1 += w*detjac*diffH1;
            ElerrorHDiv += w*detjac*diffHDiv;
        }
        error += errorEstimate;
        errorH1 += ElerrorH1;
        errorHdiv += ElerrorHDiv;
        TPZFMatrix<REAL> &hdivelsol = cmeshHDiv->ElementSolution();
        hdivelsol(el,0) = sqrt(errorEstimate);
        TPZFMatrix<REAL> &h1elsol = cmeshH1->ElementSolution();
        h1elsol(celH1->Index(),0) = sqrt(errorEstimate);
        REAL effHDiv = hdivelsol(el,0)/hdivelsol(el,1);
        if(effHDiv > 10) effHDiv = 10;
        REAL effH1 = h1elsol(celH1->Index(),0)/h1elsol(celH1->Index(),1);
        if(effH1 > 10) effH1 = 10;
        hdivelsol(el,2) = effHDiv;
        h1elsol(celH1->Index(),2) = effH1;
    }
    std::cout << "Error estimate " << sqrt(error) << " H1 error " << sqrt(errorH1) << " HDiv error " << sqrt(errorHdiv) << std::endl;
    return sqrt(error);
}
