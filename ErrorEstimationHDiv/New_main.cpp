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

#include "TPZHybridHDivErrorEstimator.h"
#include "Tools.h"

#include "TPZBFileStream.h"
#include <tuple>
#include <memory>


bool IsgmeshReader=true;
bool neumann=true;


int main5(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    
    ProblemConfig config;
    config.porder = 2;
    config.hdivmais = 0;
    
    config.ndivisions=1;
    config.prefine=false;
    config.makepressurecontinuous = true;
    
    if(!neumann){
    config.exact.fExact = TLaplaceExample1::ESinSinDirNonHom;//ESinSin;//EArcTanSingular;//ESinMark;//EArcTan;//
    config.problemname = "ESinMark";//"ESinSinDirNonHom";//"ESinSin";////"EArcTanSingular_PRef";//""ArcTang";//
    }
    
 //   FunctionTest();
    
    int dim=2;
    
    //malha geometrica
    TPZGeoMesh *gmesh = nullptr;
    
    if(IsgmeshReader){
        
//        if(neumann){
//
//            config.alpha=1;
//            std::string meshfilename = "../MeshHetero.msh";
//            TPZGmshReader gmsh;
//            gmsh.GetDimNamePhysical().resize(8);
//            gmsh.GetDimPhysicalTagName().resize(8);
//
//
//            gmsh.GetDimNamePhysical()[2]["Omega1"] = 1;
//            gmsh.GetDimNamePhysical()[2]["Omega2"] = 2;
//            gmsh.GetDimNamePhysical()[2]["Omega3"] = 3;
//            gmsh.GetDimNamePhysical()[2]["Omega4"] = 4;
//
//            gmsh.GetDimNamePhysical()[1]["neumann1"] =5;
//            gmsh.GetDimNamePhysical()[1]["neumann2"] =6;
//            gmsh.GetDimNamePhysical()[1]["neumann3"] =7;
//            gmsh.GetDimNamePhysical()[1]["neumann4"] =8;
//
//            config.materialids.insert(1);
//            config.materialids.insert(2);
//            config.materialids.insert(3);
//            config.materialids.insert(4);
//
//            config.bcmaterialids.insert(5);
//            config.bcmaterialids.insert(6);
//            config.bcmaterialids.insert(7);
//            config.bcmaterialids.insert(8);
//
//
//            gmsh.SetFormatVersion("4.1");
//            gmesh = gmsh.GeometricGmshMesh(meshfilename);
//            gmsh.PrintPartitionSummary(std::cout);
//            gmesh->SetDimension(dim);
//            config.gmesh = gmesh;
//
//        }
//
//        else{
        
            std::string meshfilename = "../LCircle.msh";
            TPZGmshReader gmsh;
            gmsh.GetDimNamePhysical()[1]["dirichlet"] =2;
            gmsh.GetDimNamePhysical()[2]["domain"] = 1;
            
            config.materialids.insert(1);
            config.bcmaterialids.insert(2);
        
            
            gmsh.SetFormatVersion("4.1");
            gmesh = gmsh.GeometricGmshMesh(meshfilename);
            gmsh.PrintPartitionSummary(std::cout);
            gmesh->SetDimension(dim);
            config.gmesh = gmesh;
  //      }
    }
    
    else{
    
    
//    int nDiv=6;
//    int nelem = 1<<nDiv;
    gmesh = CreateGeoMesh(1);
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.bcmaterialids.insert(-2);
    config.gmesh = gmesh;
    }
    
  //  UniformRefinement(config.ndivisions, gmesh);
    
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmeshInitial.txt");
        gmesh->Print(out2);
        
    }
    
    TPZManVector<TPZCompMesh*, 2> meshvec_HDiv(2, 0);
    
    TPZMultiphysicsCompMesh *cmesh_HDiv=nullptr;
    
//    if(neumann){
//
//        cmesh_HDiv = CreateNeumannHDivMesh(config);//Hdiv x L2
//
//        TPZAnalysis an(cmesh_HDiv);
//
//#ifdef USING_MKL2
//        TPZSymetricSpStructMatrix strmat(cmesh_HDiv);
//        strmat.SetNumThreads(0);
//        //        strmat.SetDecomposeType(ELDLt);
//        an.SetStructuralMatrix(strmat);
//#else
//        TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmesh_HDiv);
//        strmat.SetNumThreads(0);
//        //        TPZSkylineStructMatrix strmat3(cmesh_HDiv);
//        //        strmat3.SetNumThreads(8);
//#endif
//
//        TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
//        direct->SetDirect(ELDLt);
//        an.SetSolver(*direct);
//        delete direct;
//        direct = 0;
//        an.Assemble();
//        an.Solve();//resolve o problema misto ate aqui
//        TPZStack<std::string> scalnames, vecnames;
//        scalnames.Push("Pressure");
//        vecnames.Push("Flux");
//        an.DefineGraphMesh(2, scalnames, vecnames, "Original_Misto.vtk");
//        //        meshvec_Hybrid[1]->Solution().Print("Press");
//        // Post processing
//        an.PostProcess(0,2);
//
//
//        {
//            std::ofstream out("OriginalMixedMesh.txt");
//            cmesh_HDiv->Print(out);
//        }
//
//
//
//
//
//    }
//    else {
        cmesh_HDiv = CreateHDivMesh(config);//Hdiv x L2
 //   }
    
    cmesh_HDiv->InitializeBlock();
    meshvec_HDiv = cmesh_HDiv->MeshVector();
    
    //cria malha hibrida
    
    TPZHybridizeHDiv hybrid;
    auto HybridMesh = hybrid.Hybridize(cmesh_HDiv);
    HybridMesh->CleanUpUnconnectedNodes();//enumerar adequadamente os connects
    HybridMesh->AdjustBoundaryElements();
    delete cmesh_HDiv;
    delete meshvec_HDiv[0];
    delete meshvec_HDiv[1];
    
    int n=hybrid.fLagrangeInterface;
    int n1=hybrid.fHDivWrapMatid;
    int n2=hybrid.fInterfaceMatid;
    
    std::cout<<"---Original PerifericalMaterialId --- "<<std::endl;
    std::cout <<" LagrangeInterface = "<<n<<std::endl;
    std::cout <<" HDivWrapMatid = "<<n1<<std::endl;
    std::cout <<" InterfaceMatid = "<<n2<<std::endl;
    
    
    cmesh_HDiv=(HybridMesh);//malha hribrida
    meshvec_HDiv[0] = (HybridMesh)->MeshVector()[0];//malha Hdiv
    meshvec_HDiv[1] = (HybridMesh)->MeshVector()[1];//malha L2
    
    {
        
        //                std::ofstream outgeo("HrybridGeometria.txt");
        //                std::get<0>(HybridMesh)->Reference()->Print(outgeo);
                        std::ofstream out("OriginalHybridMesh.txt");
                        (HybridMesh)->Print(out);
        //
        std::ofstream out2("OriginalFluxMesh.txt");
        meshvec_HDiv[0]->Print(out2);
        
        std::ofstream out3("OriginalPotentialMesh.txt");
        meshvec_HDiv[1]->Print(out3);
        
        
    }
    
    
    SolveHybridProblem(cmesh_HDiv,n2,config);
    

    
    {
        std::ofstream out("OriginalHybridMesh.txt");
        (HybridMesh)->Print(out);
    }

        return 0;
    PlotLagrangreMultiplier(meshvec_HDiv[1]);
    
    
    //reconstroi potencial e calcula o erro
    {
        TPZHybridHDivErrorEstimator HDivEstimate(*cmesh_HDiv);
        
        HDivEstimate.fProblemConfig = config;
        HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
        HDivEstimate.SetAnalyticSolution(config.exact);
        
        HDivEstimate.PotentialReconstruction();
        
        TPZManVector<REAL> elementerrors;
        HDivEstimate.ComputeErrors(elementerrors);
        
    }
    
    delete cmesh_HDiv;
    delete meshvec_HDiv[0];
    delete meshvec_HDiv[1];
    return 0;
    
    
}

