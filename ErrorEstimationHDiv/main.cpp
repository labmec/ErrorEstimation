/**
 * @file Poisson 3D in hexahedra with shock problem
 */
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

bool testeCompelMultPhysics=false;
bool testeSolandDerivate=false;


bool gmeshreader=false;
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
    bool random_refine = false;
    
    if(testeSolandDerivate) PrintSolAndDerivate(config);
    
    
    TPZGeoMesh *gmesh = 0;
    if(gmeshreader){
        std::string meshfilename = "../BasicMesh.msh";
        TPZGmshReader gmsh;
        gmsh.GetDimNamePhysical()[1]["dirichlet"] = -1;
        gmsh.GetDimNamePhysical()[1]["neuman"] = -2;
        gmsh.GetDimNamePhysical()[2]["domain"] = 1;
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        config.bcmaterialids.insert(-2);
        
        
        gmsh.SetFormatVersion("3.0");
        
#ifdef MACOSX
       gmesh = gmsh.GeometricGmshMesh(meshfilename);
       //gmesh = gmsh.GeometricGmshMesh("../BasicMesh.msh");
#else
        gmesh = gmsh.GeometricGmshMesh("BasicMesh.msh");
#endif
        gmesh->SetDimension(2);
        config.gmesh = gmesh;
    }
    
    else{
        int nel = 1;
        TPZManVector<int,4> bcids(4,-1);
        gmesh=CreateGeoMesh(nel, bcids);
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        config.gmesh = gmesh;
    }
    
        
        int nDiv=1;

        UniformRefinement( nDiv,gmesh);
        
        {
            
            std::ofstream out("gmesh.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        }
    
    if(random_refine){
        int numelrefine=4;
        RandomRefine(config,numelrefine);
  
    }
    
    
    if(testeCompelMultPhysics){
        
        MultiPhysicsCompel(config);
        return 0;
        
    }
    
    
    
    
    TPZManVector<TPZCompMesh*, 2> meshvec_HDiv(2, 0);
    TPZMultiphysicsCompMesh *cmesh_HDiv = CreateHDivMesh(config);//Hdiv x L2
    cmesh_HDiv->InitializeBlock();
    meshvec_HDiv = cmesh_HDiv->MeshVector();
    /// clones the atomic meshes in meshvec_HDiv and creates a multiphysics hybrid mesh
//    std::tuple<TPZCompMesh *, TPZVec<TPZCompMesh *> > Hybridize(TPZCompMesh *cmesh_Multiphysics, TPZVec<TPZCompMesh *> &meshvec_HDiv, bool group_elements=true, double Lagrange_term_multiplier = 1.);

    //cria malha hibrida
//    TPZHybridizeHDiv hybrid;
//    auto meshAuto = hybrid.Hybridize(cmesh_HDiv, meshvec_HDiv);
//    std::get<0>(meshAuto)->CleanUpUnconnectedNodes();//enumerar adequadamente os connects
//    delete cmesh_HDiv;
//    delete meshvec_HDiv[0];
//    delete meshvec_HDiv[1];
    
//    cmesh_HDiv=std::get<0>(meshAuto);
//    meshvec_HDiv[0] = std::get<1>(meshAuto)[0];
//    meshvec_HDiv[1] = std::get<1>(meshAuto)[1];
    

//    {
//        std::ofstream out("meshvec_HDiv_flux.txt");
//        meshvec_HDiv[0]->Print(out);
//    }
//    {
//        std::ofstream out("meshvec_HDiv_pres.txt");
//        meshvec_HDiv[1]->Print(out);
//    }

    {
//        {
//            std::ofstream out("MeshHyb.txt");
//            std::get<0>(meshAuto)->Print(out);
//        }
        
        TPZAnalysis an(cmesh_HDiv);
        
#ifdef USING_MKL
        TPZSymetricSpStructMatrix strmat(cmesh_HDiv);
        strmat.SetNumThreads(0);
//        strmat.SetDecomposeType(ELDLt);
        an.SetStructuralMatrix(strmat);
#else
        TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmesh_HDiv);
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
        an.Solve();//resolve o problema misto ate aqui
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Pressure");
        vecnames.Push("Flux");
        an.DefineGraphMesh(2, scalnames, vecnames, "Original.vtk");
        //        meshvec_Hybrid[1]->Solution().Print("Press");
        // Post processing
        an.PostProcess(0,2);

//        {
//            {
//
//
//
//                TPZAnalysis an(meshvec_HDiv[1],false);
//                TPZStack<std::string> scalnames, vecnames;
//                scalnames.Push("State");
//
//                int dim = meshvec_HDiv[1]->Reference()->Dimension()-1;
//                std::string plotname;
//                {
//                    std::stringstream out;
//                    out << "LagrangeMultiplier" << ".vtk";
//                    plotname = out.str();
//                }
//                an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
//                an.PostProcess(2,dim);
//
//            }
//        }
    }
    
    //reconstroi potencial e calcula o erro
    {
        TPZHybridHDivErrorEstimator HDivEstimate(*cmesh_HDiv);


        
        HDivEstimate.fProblemConfig = config;
        HDivEstimate.fUpliftPostProcessMesh = 0;
        HDivEstimate.SetAnalyticSolution(config.exact);
        TPZManVector<REAL> elementerrors;
        HDivEstimate.ComputeErrors(elementerrors);
    }

    return 0;
}

