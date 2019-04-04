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

#include <tuple>
#include <memory>

void PlotLagrangreMultiplier(TPZCompMesh *cmesh);
void SolveHybridProblem(TPZCompMesh *Hybridmesh);


int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    
    ProblemConfig config;
    config.porder = 1;
    config.hdivmais = 1;
    config.makepressurecontinuous = true;
    config.exact.fExact = TLaplaceExample1::ESinSinDirNonHom;//ESinSin;//EArcTanSingular;//EArcTan;//ESinSinDirNonHom;//
    config.problemname = "ESinSin";//"ArcTang";//"SinSin";//"SinSinNonHom";//
    
    TPZGeoMesh *gmesh = CreateGeoMesh();
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        config.bcmaterialids.insert(-2);
        config.gmesh = gmesh;
    
    
    int nDiv=1;
    
    UniformRefinement( nDiv,gmesh);

    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmesh.txt");
        gmesh->Print(out2);
        
    }
    
    
    TPZManVector<TPZCompMesh*, 2> meshvec_HDiv(2, 0);
    TPZCompMesh *cmesh_HDiv = CreateHDivMesh(config, meshvec_HDiv);//Hdiv x L2
    cmesh_HDiv->InitializeBlock();

    //cria malha hibrida
        TPZHybridizeHDiv hybrid;
        auto HybridMesh = hybrid.Hybridize(cmesh_HDiv, meshvec_HDiv);
        std::get<0>(HybridMesh)->CleanUpUnconnectedNodes();//enumerar adequadamente os connects
        delete cmesh_HDiv;
        delete meshvec_HDiv[0];
        delete meshvec_HDiv[1];
    
        cmesh_HDiv=std::get<0>(HybridMesh);//malha hribrida
        meshvec_HDiv[0] = std::get<1>(HybridMesh)[0];//malha Hdiv
        meshvec_HDiv[1] = std::get<1>(HybridMesh)[1];//malha L2
    
        {
            std::ofstream out("HybridMesh.txt");
            std::get<0>(HybridMesh)->Print(out);
            
            std::ofstream out2("FluxMesh.txt");
           meshvec_HDiv[0]->Print(out2);
            
            std::ofstream out3("PotentialMesh.txt");
            meshvec_HDiv[1]->Print(out3);
            
           
        }

    
    
    SolveHybridProblem(cmesh_HDiv);
    
    PlotLagrangreMultiplier(meshvec_HDiv[1]);
    
 
    //reconstroi potencial e calcula o erro
    {
        TPZManVector<TPZCompMesh *> MeshesHDiv(3);
        MeshesHDiv[0] = cmesh_HDiv;//malha hibrida
        MeshesHDiv[1] = meshvec_HDiv[0];//Hdiv
        MeshesHDiv[2] = meshvec_HDiv[1];//L2
        TPZHybridHDivErrorEstimator HDivEstimate(MeshesHDiv);

        HDivEstimate.fProblemConfig = config;
        HDivEstimate.fUpliftPostProcessMesh = 0;
        HDivEstimate.SetAnalyticSolution(&config.exact);
        
        HDivEstimate.PotentialReconstruction();
        
        TPZManVector<REAL> elementerrors;
        HDivEstimate.ComputeErrors(elementerrors);
    }
    
    return 0;
    
    
}
    

void SolveHybridProblem(TPZCompMesh *Hybridmesh){
    TPZAnalysis an(Hybridmesh);
    
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(Hybridmesh);
    strmat.SetNumThreads(0);
    //        strmat.SetDecomposeType(ELDLt);
    an.SetStructuralMatrix(strmat);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(Hybridmesh);
    strmat.SetNumThreads(0);
    //        TPZSkylineStructMatrix strmat3(cmesh_HDiv);
    //        strmat3.SetNumThreads(8);
#endif
    
    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    vecnames.Push("Flux");
    an.DefineGraphMesh(2, scalnames, vecnames, "Original.vtk");
    //        meshvec_Hybrid[1]->Solution().Print("Press");
    // Post processing
    an.PostProcess(2,2);
    
    


    
    
    
    
}
void PlotLagrangreMultiplier(TPZCompMesh *cmesh){
    
    TPZAnalysis an(cmesh,false);
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("State");

    int dim = cmesh->Reference()->Dimension()-1;
    std::string plotname;
    {
        std::stringstream out;
        out << "LagrangeMultiplier" << ".vtk";
        plotname = out.str();
    }
    an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
    an.PostProcess(2,dim);
    
}