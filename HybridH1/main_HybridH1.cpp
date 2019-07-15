//
//  New_main.cpp
//  ErrorEstimateHDiv
//
//  Created by Denise De Siqueira on 01/04/19.
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
#include "TPZVecL2.h"
#include "pzbndcond.h"

#include "pzintel.h"

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




bool IsgmeshReader = true;
bool neumann = true;

void InsertMaterialObjectsH1Hybrid(TPZMultiphysicsCompMesh *cmesh, ProblemConfig &config);

int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    // Initializing uniform refinements for reference elements
//    gRefDBase.InitializeUniformRefPattern(EOned);
//    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
//    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    ProblemConfig config;
    
    config.porder = 1;
    config.hdivmais = 0;
    config.ndivisions = 0;
    config.dimension = 2;
    config.prefine=false;
    config.makepressurecontinuous = true;
    
    
    config.exact.fExact = TLaplaceExample1::EX;//ESinSinDirNonHom;//ESinSin;//ESinMark;//EX;//EConst;//EArcTanSingular;//EArcTan;//
    config.problemname = "Constant k=1 e n=0 Up=2";//"EConst";//"ESinSinDirNonHom";//"ESinSin";//" //"EArcTanSingular_PRef";//""ArcTang";//
    
    config.dir_name= "HybridH1_Const";
    //config.dir_name= "ESinSin";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());
    
    
    //geometric mesh
    TPZManVector<int,4> bcids(4,-1);
    TPZGeoMesh *gmesh = CreateGeoMesh(2, bcids);
    if(1)
    {
        TPZManVector<TPZGeoEl *> sub;
        gmesh->Element(0)->Divide(sub);
        DivideLowerDimensionalElements(gmesh);
    }
    
    UniformRefinement(config.ndivisions, gmesh);
    // RandomRefine(config, config.ndivisions);
    
#ifdef PZDEBUG
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmeshInitial.txt");
        gmesh->Print(out2);
        
    }
#endif
    
    TPZCreateMultiphysicsSpace createspace(gmesh);
    
    createspace.SetMaterialIds({1}, {-2,-1});
    createspace.fH1Hybrid.fHybridizeBC = false;
    createspace.ComputePeriferalMaterialIds();
    
    TPZManVector<TPZCompMesh *> meshvec;
    createspace.CreateAtomicMeshes(meshvec);
    
    TPZMultiphysicsCompMesh *cmesh_H1Hybrid = new TPZMultiphysicsCompMesh(gmesh);
    InsertMaterialObjectsH1Hybrid(cmesh_H1Hybrid, config);
    createspace.InsertPeriferalMaterialObjects(cmesh_H1Hybrid);
    cmesh_H1Hybrid->BuildMultiphysicsSpace(meshvec);
    
    createspace.AddInterfaceElements(cmesh_H1Hybrid);
    createspace.GroupandCondenseElements(cmesh_H1Hybrid);
    
    cmesh_H1Hybrid->InitializeBlock();
    cmesh_H1Hybrid->ComputeNodElCon();
    
    TPZAnalysis an(cmesh_H1Hybrid);
    TPZSymetricSpStructMatrix sparse(cmesh_H1Hybrid);
    TPZSkylineStructMatrix skylstr(cmesh_H1Hybrid);
    an.SetStructuralMatrix(skylstr);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    
    an.Run();
    
#ifdef PZDEBUG
    {
        
        std::ofstream out2("H1HybridMesh.txt");
        cmesh_H1Hybrid->Print(out2);
        std::ofstream out3("gmeshHybridH1.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out3);
        
    }
#endif
    
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
    if(config.exact.fExact != TLaplaceExample1::ENone)
    {
        material->SetForcingFunction(config.exact.ForcingFunction());
    }
    //    TPZMaterial * mat(material);
    //    cmesh->InsertMaterialObject(mat);
    
    // Inserts boundary conditions
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 1.);
    TPZMaterial * BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2);
    if(config.exact.fExact != TLaplaceExample1::ENone)
    {
        BCond0->SetForcingFunction(config.exact.Exact());
    }
    val2.Zero();
    TPZMaterial * BCond1 = material->CreateBC(material, -2, neumann, val1, val2);
    
    cmesh_H1Hybrid->InsertMaterialObject(BCond0);
    cmesh_H1Hybrid->InsertMaterialObject(BCond1);
    

}
