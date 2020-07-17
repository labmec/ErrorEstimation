//
//  New_main.cpp
//  ErrorEstimateHDiv
//
//  Created by Denise De Siqueira on 01/04/19.
//

#include "TPZGmshReader.h"
#include "TPZRefPatternDataBase.h"
#include "TPZVTKGeoMesh.h"
#include "pzlog.h"

#include "tpzarc3d.h"
#include "tpzgeoelrefpattern.h"
#include "tpzquadraticquad.h"
#include "tpzquadratictrig.h"

#include "ProblemConfig.h"

#include "TPZVecL2.h"
#include "mixedpoisson.h"
#include "pzbndcond.h"

#include "pzintel.h"

#include "TPZCompMeshTools.h"
#include "TPZHybridizeHDiv.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZMultiphysicsCompMesh.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "TPZHDivErrorEstimatorH1.h"
#include "TPZHybridHDivErrorEstimator.h"

#include "Tools.h"

#include "TPZBFileStream.h"
#include <memory>
#include <tuple>

bool IsgmeshReader = false; // para ler a malha
bool neumann = true; // para o problema local de neumann da forlmulacao Mark

bool mixedsolution = true; // se quiser rodar o problema misto

bool PostProcessingFEM = true; // para graficos da solucao FEM

int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(ECube);

    ProblemConfig config;
    

    config.dimension = 2;
    config.prefine = false;
    config.makepressurecontinuous = true;
    config.TensorNonConst = false; // para problem 3d com tensor nao constante

    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;

    bool RunMark = false;
    config.problemname = "ESinSin";

    config.dir_name = "RobinBC";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    int dim = config.dimension;
    
    TPZManVector<int,4> bcids(4,-3);
    //bcids[0] = -1;
    config.bcmaterialids.insert(-3);//Robin
    config.bcmaterialids.insert(-1);//Dirichlet

    config.Km = 100;
    config.hdivmais = 3;
    
for (int p = 1; p <2; p ++){
        
config.porder = p;
    for (int ndiv = 1; ndiv < 2; ndiv++) {

        config.ndivisions = ndiv;

        config.adaptivityStep = ndiv;

        // malha geometrica
        TPZGeoMesh *gmesh = nullptr;

        if (IsgmeshReader) {

            std::string meshfilename = "../LCircle.msh";
            // std::string meshfilename = "../esfera2.msh";

            if (dim == 3) {
                meshfilename = "../Cube.msh";
            }
            TPZGmshReader gmsh;
            //  gmsh.GetDimNamePhysical().resize(4);
            //  gmsh.GetDimPhysicalTagName().resize(4);
            if (dim == 2) {
                gmsh.GetDimNamePhysical()[1]["dirichlet"] = 2;
                // gmsh.GetDimNamePhysical()[1]["boundary"] =2;
                gmsh.GetDimNamePhysical()[2]["domain"] = 1;
                //  gmsh.GetDimNamePhysical()[1]["boundary2"] =3;

            } else {
                gmsh.GetDimNamePhysical()[2]["dirichlet"] = 2;
                gmsh.GetDimNamePhysical()[3]["domain"] = 1;
            }
            config.materialids.insert(1);
            config.bcmaterialids.insert(2);
            // config.bcmaterialids.insert(3);

            gmsh.SetFormatVersion("4.1");
            gmesh = gmsh.GeometricGmshMesh(meshfilename);
            gmsh.PrintPartitionSummary(std::cout);
            gmesh->SetDimension(dim);
            config.gmesh = gmesh;

        }

        else {
            
          //  TPZManVector<int,4> bcids(4,-3);
            //TPZManVector<int,4> bcids(4,-3);
            //TPZManVector<int,4> bcids(3,-3);
          //  bcids[1] = -1;
            //constants for Robin boundary conditions
            // sigma.n=Km(u-u_d)-g
            //Particular cases: 1) Km=0---> Neumann, 2) Km=infinity-->Dirichlet
            //config.coefG = 0.;//nao passar mais isso
            //config.Km = 1.e12;//pow(10,2);
            

            int nel = 1;//pow(2, ndiv);

            gmesh = CreateGeoMesh(nel, bcids); //CreateLShapeMesh(bcids);//CreateQuadMeshRefTriang(bcids); //CreateSingleTriangleMesh(bcids);// CreateTrapezoidalMesh(nelT,
                             // nelT, 1.,1.,bcids);//CreateLCircleGeoMesh();//
            config.materialids.insert(1);
            //config.bcmaterialids.insert(-1); // dirichlet
            //config.bcmaterialids.insert(-2); // neumann
           // config.bcmaterialids.insert(-3); // Robin
            config.gmesh = new TPZGeoMesh;
            *config.gmesh = *gmesh;
            gmesh->SetDimension(dim);
        }

#ifdef PZDEBUG
        {
            std::ofstream out("gmeshInit.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
            std::ofstream out2("gmeshTo.txt");
            gmesh->Print(out2);
        }
#endif
        

        

        TPZGeoMesh *hybridEstimatorMesh = new TPZGeoMesh();
        *hybridEstimatorMesh = *gmesh;

        UniformRefinement(config.ndivisions, gmesh);
        DivideLowerDimensionalElements(gmesh);

        *config.gmesh = *gmesh;
        
        

#ifdef PZDEBUG
        {
            std::ofstream out("gmesh.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
            std::ofstream out2("gmeshInitial.txt");
            gmesh->Print(out2);
        }
#endif

        TPZManVector<TPZCompMesh *, 2> meshvec_HDiv(2, 0);

        TPZMultiphysicsCompMesh *cmesh_HDiv = nullptr;

        cmesh_HDiv = CreateHDivMesh(config); // Hdiv x L2
        cmesh_HDiv->InitializeBlock();

#ifdef PZDEBUG
        {

            std::ofstream out2("MalhaMista.txt");
            cmesh_HDiv->Print(out2);
        }
#endif

       //  SolveMixedProblem(cmesh_HDiv,config);

        meshvec_HDiv = cmesh_HDiv->MeshVector();

        // cria malha hibrida

        TPZHybridizeHDiv hybrid;
        auto HybridMesh = hybrid.Hybridize(cmesh_HDiv);
        HybridMesh
            ->CleanUpUnconnectedNodes(); // enumerar adequadamente os connects
        HybridMesh->AdjustBoundaryElements();
        delete cmesh_HDiv;
        delete meshvec_HDiv[0];
        delete meshvec_HDiv[1];

        std::cout << "---Original PerifericalMaterialId --- " << std::endl;
        std::cout << " LagrangeInterface = " << hybrid.fLagrangeInterface
                  << std::endl;
        std::cout << " HDivWrapMatid = " << hybrid.fHDivWrapMatid << std::endl;
        std::cout << " InterfaceMatid = " << hybrid.fInterfaceMatid
                  << std::endl;

        cmesh_HDiv = (HybridMesh);                       // malha hribrida
        meshvec_HDiv[0] = (HybridMesh)->MeshVector()[0]; // malha Hdiv
        meshvec_HDiv[1] = (HybridMesh)->MeshVector()[1]; // malha L2

#ifdef PZDEBUG2
        {

            std::ofstream out2("OriginalFluxMesh.txt");
            meshvec_HDiv[0]->Print(out2);

            std::ofstream out3("OriginalPotentialMesh.txt");
            meshvec_HDiv[1]->Print(out3);
        }
#endif

        SolveHybridProblem(cmesh_HDiv, hybrid.fInterfaceMatid, config, false);

        

#ifdef PZDEBUG2
        {
            std::ofstream out("OriginalHybridMesh.txt");
            (HybridMesh)->Print(out);
        }
#endif
      //  PlotLagrangeMultiplier(meshvec_HDiv[1],config);

        // reconstroi potencial e calcula o erro
        {

            if (RunMark) {
                TPZHDivErrorEstimatorH1 HDivEstimate(*cmesh_HDiv);
                HDivEstimate.fProblemConfig = config;
                HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
                HDivEstimate.SetAnalyticSolution(config.exact);
                HDivEstimate.fperformUplift = true;
                HDivEstimate.fUpliftOrder = config.hdivmais;

                HDivEstimate.PotentialReconstruction();

                TPZManVector<REAL> elementerrors;
                TPZVec<REAL> errorVec;
                HDivEstimate.ComputeErrors(errorVec,elementerrors,true);

            }

            else {

                TPZHybridHDivErrorEstimator HDivEstimate(*cmesh_HDiv);
                HDivEstimate.SetHybridizer(hybrid);
                HDivEstimate.fProblemConfig = config;
                HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
                HDivEstimate.SetAnalyticSolution(config.exact);
                HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;

                HDivEstimate.fPostProcesswithHDiv = false;

                HDivEstimate.PotentialReconstruction();

                TPZManVector<REAL> elementerrors;
                TPZVec<REAL> errorVec;
                HDivEstimate.ComputeErrors(errorVec,elementerrors,true);
       //         hAdaptivity(&HDivEstimate.fPostProcMesh, hybridEstimatorMesh,config);
            }
        }

        delete cmesh_HDiv;
        delete meshvec_HDiv[0];
        delete meshvec_HDiv[1];
        // return 0;
 
    }
}
}
