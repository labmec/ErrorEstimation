//
// Created by Gustavo Batistela on 12/07/19.
//

#include "TPZRefPatternDataBase.h"
#include "TPZGmshReader.h"

#include "tpzarc3d.h"
#include <math.h>
#include <tuple>

#include "ProblemConfig.h"
#include "Tools.h"
#include "TPZHybridHDivErrorEstimator.h"
#include "TPZHDivErrorEstimatorH1.h"

bool IsgmeshReader = true;
bool mixedsolution = false;

int main(int argc, char *argv[]) {

#ifdef LOG4CXX
    InitializePZLOG();
#endif

    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);


    for (int ndiv = 0; ndiv < 1; ndiv++) {
        ProblemConfig config;

        config.porder = 2;
        config.hdivmais = 1;
        config.ndivisions = ndiv;
        config.dimension = 2;
        config.prefine = false;
        config.makepressurecontinuous = true;

        config.exact.fExact = TLaplaceExample1::ESinMark;
        config.problemname = "AdaptivityMark";

        config.dir_name = "ESinMark";
        std::string command = "mkdir " + config.dir_name;
        system(command.c_str());


        //malha geometrica
        TPZGeoMesh *gmesh = nullptr;

        if (IsgmeshReader) {
            std::string meshfilename = "LCircle.msh";
            //std::string meshfilename = "../esfera2.msh";
#ifdef MACOSX
            meshfilename = "../" + meshfilename;
#endif
            TPZGmshReader gmsh;

            gmsh.GetDimNamePhysical()[1]["dirichlet"] = 2;
            gmsh.GetDimNamePhysical()[2]["domain"] = 1;

            config.materialids.insert(1);
            config.bcmaterialids.insert(2);
            config.bcmaterialids.insert(3);

            gmsh.SetFormatVersion("4.1");

            gmesh = gmsh.GeometricGmshMesh(meshfilename);
            int dim = gmsh.Dimension();
            gmesh->SetDimension(dim);

            gmsh.PrintPartitionSummary(std::cout);

            config.gmesh = gmesh;
        }

        UniformRefinement(config.ndivisions, gmesh);
        // RandomRefine(config, config.ndivisions);

#ifdef PZDEBUG
        {
            std::ofstream out("GeometricMesh.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
            std::ofstream out2("gmeshInitial.txt");
            gmesh->Print(out2);
        }
#endif

        TPZMultiphysicsCompMesh *cmesh_HDiv = nullptr;

        cmesh_HDiv = CreateHDivMesh(config); //Hdiv x L2
        cmesh_HDiv->InitializeBlock();

        TPZManVector<TPZCompMesh *, 2> meshvec_HDiv(2, 0);
        meshvec_HDiv = cmesh_HDiv->MeshVector();

        //cria malha hibrida
        TPZHybridizeHDiv hybrid;
        auto HybridMesh = hybrid.Hybridize(cmesh_HDiv);
        HybridMesh->CleanUpUnconnectedNodes();//enumerar adequadamente os connects
        HybridMesh->AdjustBoundaryElements();
        delete cmesh_HDiv;
        delete meshvec_HDiv[0];
        delete meshvec_HDiv[1];

        std::cout << "---Original PerifericalMaterialId --- " << std::endl;
        std::cout << " LagrangeInterface = " << hybrid.fLagrangeInterface << std::endl;
        std::cout << " HDivWrapMatid = " << hybrid.fHDivWrapMatid << std::endl;
        std::cout << " InterfaceMatid = " << hybrid.fInterfaceMatid << std::endl;

        cmesh_HDiv = (HybridMesh);//malha hribrida
        meshvec_HDiv[0] = (HybridMesh)->MeshVector()[0];//malha Hdiv
        meshvec_HDiv[1] = (HybridMesh)->MeshVector()[1];//malha L2

        SolveHybridProblem(cmesh_HDiv, hybrid.fInterfaceMatid, config);

        // reconstroi potencial e calcula o erro
        {
            TPZHDivErrorEstimatorH1 HDivEstimate(*cmesh_HDiv);
            HDivEstimate.fProblemConfig = config;
            HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
            HDivEstimate.SetAnalyticSolution(config.exact);

            HDivEstimate.fperformUplift = true;
            HDivEstimate.fUpliftOrder = 1;

            HDivEstimate.PotentialReconstruction();

            TPZManVector<REAL> elementerrors;
            HDivEstimate.ComputeErrors(elementerrors);
        }

        delete cmesh_HDiv;
        delete meshvec_HDiv[0];
        delete meshvec_HDiv[1];

    }
    return 0;
}

