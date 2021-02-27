//
//  main.cpp
//  Steklov problem
//
//  Created by Denise De Siqueira on 25/11/20.
//

#include "TPZGmshReader.h"
#include "TPZRefPatternDataBase.h"
#include "pzlog.h"
#include "tpzgeoelrefpattern.h"
#include "ProblemConfig.h"
#include "pzbndcond.h"
#include "TPZHybridizeHDiv.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZHybridHDivErrorEstimator.h"
#include "Tools.h"
#include "TPZBFileStream.h"
#include "TPZSteklovMaterial.h"


TPZMultiphysicsCompMesh * CreateMixedMesh(const ProblemConfig& problem);

void RunSteklovProblem();


int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeAllUniformRefPatterns();

  
    RunSteklovProblem();
  

    return 0;
}

void RunSteklovProblem() {
    ProblemConfig config;
    config.dimension = 2;
    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinMark;
    config.problemname = "SteklovOnSquare";
    config.dir_name = "SteklovResults";
    config.porder = 1;
    config.hdivmais = 1;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-2);
    config.makepressurecontinuous = true;

    int nRef = 0;
    int nElem = 1 ;

    config.ndivisions = nRef;
    TPZManVector<int, 4> bcIDs(4, -2);
    bcIDs[0] = -1;
    config.gmesh = Tools::CreateNewGeoMesh(nElem, bcIDs);
    
   
    string command = "mkdir " + config.dir_name;
    system(command.c_str());
    {
         std::string fileName = config.dir_name + "/" + config.problemname + "GeoSquareMesh.vtk";
         std::ofstream file(fileName);
         TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    
     }
    
    TPZMultiphysicsCompMesh *MixedCMesh = CreateMixedMesh(config);
    
    MixedCMesh->InitializeBlock();
    
    {
         std::string fileName = config.dir_name + "/" + config.problemname + "CMesh.txt";
         std::ofstream file(fileName);
        
        
        MixedCMesh->Print(file);
    
     }
    
       
    Tools::SolveMixedProblem(MixedCMesh,config);


}


TPZMultiphysicsCompMesh * CreateMixedMesh(const ProblemConfig& problem) {

    TPZMultiphysicsCompMesh* cmesh = new TPZMultiphysicsCompMesh(problem.gmesh);
    TPZMaterial* mat = NULL;
    TPZFMatrix<REAL> K(3, 3, 0), invK(3, 3, 0);
    K.Identity();
    invK.Identity();
    
    STATE Km = problem.Km;

    
 

//    K.Print(std::cout);
//    invK.Print(std::cout);
    
    for (auto matid : problem.materialids) {
        TPZSteklovMaterial *mix = new TPZSteklovMaterial(matid, cmesh->Dimension());
        mix->SetPermeabilityTensor(K, invK);

        if (!mat) mat = mix;

        cmesh->InsertMaterialObject(mix);
    }
        
    
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype;
    
        switch (matid) {
            case -1 :{
            bctype = 0;
                break;
            }
                
                
            case -2:{
            bctype = 1;
    
            break;
            }
            case -3:{
            bctype = 4;// different from mixed (bctype 2) already implemented on TPZMixedPoisson3d
            val1(0,0) = Km ;

                
            break;
            }
        }
        TPZBndCond* bc = mat->CreateBC(mat, matid, bctype, val1, val2);
       // bc->TPZMaterial::SetForcingFunction(problem.exact.operator*().Exact());
        cmesh->InsertMaterialObject(bc);
    }
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    
    TPZManVector<int> active(2, 1);
    TPZManVector<TPZCompMesh*> meshvector(2, 0);
    
    meshvector[0] = Tools::CreateFluxHDivMesh(problem);
    meshvector[1] = Tools::CreatePressureMesh(problem);
    
    TPZCompMeshTools::AdjustFluxPolynomialOrders(meshvector[0], problem.hdivmais);
    TPZCompMeshTools::SetPressureOrders(meshvector[0], meshvector[1]);
    
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    bool keepmatrix = false;
    bool keeponelagrangian = true;
    TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian, keepmatrix);
    cmesh->InitializeBlock();
    return cmesh;
}
