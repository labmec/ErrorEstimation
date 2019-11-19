//
// Created by gustavo on 16/11/2019.
//

#include <pybind11/pybind11.h>

#include "TPZHDivErrorEstimatorH1.h"
#include "Tools.h"

namespace py = pybind11;

// ---------------------------------------------------
#include "TPZRefPatternDataBase.h"
#include "pzgengrid.h"
#include "TPZVTKGeoMesh.h"
// ---------------------------------------------------

int add() {

    bool neumann = true; //para o problema local de neumann da forlmulacao Mark
    bool mixedsolution = true;//se quiser rodar o prolbema misto
    bool PostProcessingFEM = true;//para graficos da solucao FEM

    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);

    ProblemConfig config;

    config.porder = 1;
    config.hdivmais = 0;
    config.ndivisions = 2;
    config.dimension = 2;
    config.prefine=false;
    config.makepressurecontinuous = true;
    config.adaptivityStep = 2;
    config.TensorNonConst = false;//para problem 3d com tensor nao constante

    config.exact.fExact = TLaplaceExample1::ESinMark;//ESinSinDirNonHom;//ESinSin;//ESinMark;//EX;//EConst;//EArcTanSingular;//EArcTan;//
    config.problemname = "ESinSin";//"EConst";//"ESinSinDirNonHom";//"ESinSin";//" //"EArcTanSingular_PRef";//""ArcTang";//

    config.dir_name= "ESinMark";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());
    int dim = config.dimension;

    // Geometric mesh
    TPZGeoMesh *gmesh = nullptr;

    TPZManVector<int,4> bcids(4,-1);
    gmesh = CreateGeoMesh(2, bcids);
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.gmesh = gmesh;
    gmesh->SetDimension(dim);

    UniformRefinement(config.ndivisions, gmesh);

    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmeshInitial.txt");
        gmesh->Print(out2);
    }

    TPZManVector<TPZCompMesh*, 2> meshvec_HDiv(2, 0);

    TPZMultiphysicsCompMesh *cmesh_HDiv=nullptr;

    cmesh_HDiv = CreateHDivMesh(config);//Hdiv x L2
    cmesh_HDiv->InitializeBlock();

    if(mixedsolution) SolveMixedProblem(cmesh_HDiv,config);


    meshvec_HDiv = cmesh_HDiv->MeshVector();

    //cria malha hibrida

    TPZHybridizeHDiv hybrid;
    auto HybridMesh = hybrid.Hybridize(cmesh_HDiv);
    HybridMesh->CleanUpUnconnectedNodes();//enumerar adequadamente os connects
    HybridMesh->AdjustBoundaryElements();
    delete cmesh_HDiv;
    delete meshvec_HDiv[0];
    delete meshvec_HDiv[1];


    std::cout<<"---Original PerifericalMaterialId --- "<<std::endl;
    std::cout <<" LagrangeInterface = "<<hybrid.fLagrangeInterface<<std::endl;
    std::cout <<" HDivWrapMatid = "<<hybrid.fHDivWrapMatid<<std::endl;
    std::cout <<" InterfaceMatid = "<<hybrid.fInterfaceMatid<<std::endl;


    cmesh_HDiv=(HybridMesh);//malha hribrida
    meshvec_HDiv[0] = (HybridMesh)->MeshVector()[0];//malha Hdiv
    meshvec_HDiv[1] = (HybridMesh)->MeshVector()[1];//malha L2

    SolveHybridProblem(cmesh_HDiv,hybrid.fInterfaceMatid,config);

    //reconstroi potencial e calcula o erro
    {
        bool RunMark = 0;
        if(RunMark){
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

        else {

            TPZHybridHDivErrorEstimator HDivEstimate(*cmesh_HDiv);
            HDivEstimate.fProblemConfig = config;
            HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
            HDivEstimate.SetAnalyticSolution(config.exact);
            HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;

            HDivEstimate.fPostProcesswithHDiv = true;

            HDivEstimate.PotentialReconstruction();

            TPZManVector<REAL> elementerrors;
            HDivEstimate.ComputeErrors(elementerrors);
            // hAdaptivity(&HDivEstimate.fPostProcMesh, hybridEstimatorMesh);

        }


    }

    return 0;
}

PYBIND11_MODULE(errorestimation, m) {

    m.doc() = R"pbdoc(
        --------------------------------------------------
        Python bindings for NeoPZ Error Estimation classes
        --------------------------------------------------
    )pbdoc";

    py::class_<ProblemConfig>(m, "ProblemConfig")
            .def(py::init<>());

    py::class_<TPZMultiphysicsCompMesh>(m, "TPZMultiphysicsCompMesh")
        .def(py::init<>())
        .def("AutoBuild", &TPZMultiphysicsCompMesh::AutoBuild)
        .def("__repr__",
             [](const TPZMultiphysicsCompMesh &obj) {
                 return "It works, go home.";
             }
        );

    m.def("add", &add, "A function which adds two numbers");

    m.def("Create2DGridMesh",
        [] (const int x_nel, const int y_nel, TPZManVector<REAL> x0, TPZManVector<REAL> x1, const TPZManVector<int> bcIDs) {

            TPZManVector<int> nel(2);
            nel[0] = x_nel;
            nel[1] = y_nel;

            TPZGenGrid gen(nel, x0, x1);
            gen.SetRefpatternElements(true);
            TPZGeoMesh* gmesh = new TPZGeoMesh;
            gen.Read(gmesh);
            gen.SetBC(gmesh, 4, bcIDs[0]);
            gen.SetBC(gmesh, 5, bcIDs[1]);
            gen.SetBC(gmesh, 6, bcIDs[2]);
            gen.SetBC(gmesh, 7, bcIDs[3]);

            return gmesh;
        }, py::return_value_policy::reference, "A function that creates a 2D grid mesh");

    m.def("UniformRefinement", &UniformRefinement, "A function that refines a Geometric Mesh uniformly");

    m.def("PrintGMeshToVTK",
        [] (TPZGeoMesh* gmesh, std::string fileName) {
            std::ofstream out(fileName);
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
          }, "Prints a VTK file with the TPZGeoMesh information");
}