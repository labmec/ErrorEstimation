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


TLaplaceExample1 laplaceex1;

int main() {

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    const int porder = 1;
    laplaceex1.fExact = TLaplaceExample1::EBubble2D;
    laplaceex1.fExact = TLaplaceExample1::ESinSin;

    TPZGeoMesh *gmesh = ReadGmshSimple("quadmesh.msh");

    {
        std::ofstream out("gmesh.txt");
        gmesh->Print(out);

        std::ofstream out2("gmesh.vtk"); 
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
    }

    // TPZCheckGeom check(gmesh);
    // check.UniformRefine(4);

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


    // {
    //     std::ofstream out("cmesh.txt");
    //     cmesh->Print(out);
    // }
    delete cmeshHDiv;
    delete cmeshH1;

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
        "TrueError"
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

