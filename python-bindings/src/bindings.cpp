//
// Created by gustavo on 16/11/2019.
//

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "TPZHDivErrorEstimatorH1.h"
#include "Tools.h"

namespace py = pybind11;

// ---------------------------------------------------
#include "pzlog.h"
#include "TPZRefPatternDataBase.h"
#include "pzgengrid.h"
#include "TPZVTKGeoMesh.h"
// ---------------------------------------------------

void EstimateErrorWithH1Reconstruction(TPZMultiphysicsCompMesh* cmesh_HDiv, ProblemConfig& config) {
    TPZHDivErrorEstimatorH1 HDivEstimate(*cmesh_HDiv);
    HDivEstimate.fProblemConfig = config;
    HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
    HDivEstimate.SetAnalyticSolution(config.exact);
    HDivEstimate.fperformUplift = true;
    HDivEstimate.fUpliftOrder = 1;

    HDivEstimate.PotentialReconstruction();

    TPZManVector<REAL> elementErrors;
    HDivEstimate.ComputeErrors(elementErrors);
}

void EstimateErrorWithHdivReconstruction(TPZMultiphysicsCompMesh* cmeshHDiv, ProblemConfig& config) {
    TPZHybridHDivErrorEstimator HDivEstimate(*cmeshHDiv);
    HDivEstimate.fProblemConfig = config;
    HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
    HDivEstimate.SetAnalyticSolution(config.exact);
    HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;

    HDivEstimate.fPostProcesswithHDiv = true;

    HDivEstimate.PotentialReconstruction();

    TPZManVector<REAL> elementErrors;
    HDivEstimate.ComputeErrors(elementErrors);
}

TPZGeoMesh* Create2DGridMesh(const int x_nel, const int y_nel, TPZManVector<REAL> x0, TPZManVector<REAL> x1,
                             const TPZManVector<int> bcIDs) {

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
    gmesh->SetDimension(2);
    return gmesh;
}

void Initialize2DUniformRefPatterns() {
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
}

TPZMultiphysicsCompMesh* CreateHybridMultiphysicsMesh(ProblemConfig& config) {

    TPZManVector<TPZCompMesh*, 2> meshvecHDiv(2, 0);
    TPZMultiphysicsCompMesh* cmeshHDiv = nullptr;

    cmeshHDiv = CreateHDivMesh(config);//Hdiv x L2
    cmeshHDiv->InitializeBlock();

    TPZHybridizeHDiv hybrid;
    auto HybridMesh = hybrid.Hybridize(cmeshHDiv);
    HybridMesh->CleanUpUnconnectedNodes();
    HybridMesh->AdjustBoundaryElements();

    config.materialids.insert(hybrid.fInterfaceMatid);

    return HybridMesh;
}

int add() {
// ===============================================================
    InitializePZLOG();
// ===============================================================
    Initialize2DUniformRefPatterns();
// ===============================================================

    // Geometric mesh
    TPZGeoMesh* gmesh = nullptr;
    TPZManVector<int, 4> bcids(4, -1);
    TPZManVector<REAL, 3> coord1(3, 0.);
    TPZManVector<REAL, 3> coord2(3, 1.);
    coord2[2] = 0;
    gmesh = Create2DGridMesh(2, 2, coord1, coord2, bcids);
    gmesh->SetDimension(2);


// ===============================================================

    ProblemConfig config;

    config.ndivisions = 2;
    UniformRefinement(config.ndivisions, gmesh);
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmeshAfterAdd.txt");
        gmesh->Print(out2);
    }


    config.gmesh = gmesh;
    config.porder = 1;
    config.hdivmais = 0;
    config.dimension = gmesh->Dimension();
    config.prefine = false;
    config.makepressurecontinuous = true;
    config.adaptivityStep = 2;
    config.TensorNonConst = false;//para problem 3d com tensor nao constante

    config.exact.fExact = TLaplaceExample1::ESinSinDirNonHom;//ESinSinDirNonHom;//ESinSin;//ESinMark;//EX;//EConst;//EArcTanSingular;//EArcTan;//
    config.problemname = "ESinSinDirNomHom";//"EConst";//"ESinSinDirNonHom";//"ESinSin";//" //"EArcTanSingular_PRef";//""ArcTang";//

    config.dir_name = "ESinAdd";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());
    int dim = config.dimension;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);

// ===============================================================
    TPZMultiphysicsCompMesh* cmesh_HDiv = CreateHybridMultiphysicsMesh(config);
// ===============================================================
    SolveHybridProblem(cmesh_HDiv, config);
// ===============================================================

    // Reconstructs the potential and calculate errors
    bool RunMark = 1;
    if (RunMark) {
        EstimateErrorWithH1Reconstruction(cmesh_HDiv, config);
    } else {
        EstimateErrorWithHdivReconstruction(cmesh_HDiv, config);
    }

    return 0;
}

PYBIND11_MODULE(errorestimation, m) {

    m.doc() = R"pbdoc(
        --------------------------------------------------
        Python bindings for NeoPZ Error Estimation classes
        --------------------------------------------------
    )pbdoc";

    m.def("add", &add, "A function which adds two numbers");

    m.def("InitializePZLog", py::overload_cast<>(&InitializePZLOG),
          "Initializes log file for log4cxx with common name log4cxx.cfg");

    m.def("Initialize2DUniformRefPatterns", &Initialize2DUniformRefPatterns,
          "Initializes refinement patterns for topologies up to 2 dimensions");

    m.def("Create2DGridMesh", &Create2DGridMesh, "A function that creates a 2D grid mesh");

    m.def("UniformRefinement", &UniformRefinement, "A function that refines a Geometric Mesh uniformly");

    m.def("PrintGMeshToVTK",
          [](TPZGeoMesh* gmesh, std::string fileName) {
              std::ofstream out(fileName);
              TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
          }, "Prints a VTK file with the TPZGeoMesh information");

    m.def("PrintGMeshToTXT",
          [](TPZGeoMesh* gmesh, std::string fileName) {
              std::ofstream out(fileName);
              gmesh->Print(out);
          }, "Prints a text file with the TPZGeoMesh information");

    m.def("CreateHybridMultiphysicsMesh", &CreateHybridMultiphysicsMesh, py::return_value_policy::reference,
          "Creates a multiphysics computational mesh which is already hybridized");

    m.def("EstimateErrorWithH1Reconstruction", &EstimateErrorWithH1Reconstruction,
          "Estimate approximation error using Ainsworth's proposal");

    m.def("EstimateErrorWithHdivReconstruction", &EstimateErrorWithHdivReconstruction,
          "Estimate approximation error using H(div) reconstruction");

    m.def("SolveHybridProblem", py::overload_cast<TPZMultiphysicsCompMesh*, const ProblemConfig&>(&SolveHybridProblem),
          "Solves a hybrid H(div) finite element problem");

    // ProblemConfig
    py::class_<ProblemConfig>(m, "ProblemConfig")
            .def(py::init())
            .def_property("gmesh", &ProblemConfig::getGmesh, &ProblemConfig::setGmesh)
            .def_property("Porder", &ProblemConfig::getPorder, &ProblemConfig::setPorder)
            .def_property("Hdivmais", &ProblemConfig::getHdivmais, &ProblemConfig::setHdivmais)
            .def_property("Makepressurecontinuous", &ProblemConfig::isMakepressurecontinuous,
                          &ProblemConfig::setMakepressurecontinuous)
            .def_property("NDivisions", &ProblemConfig::getNdivisions, &ProblemConfig::setNdivisions)
            .def_property("AdaptivityStep", &ProblemConfig::getAdaptivityStep, &ProblemConfig::setAdaptivityStep)
            .def_property("Dimension", &ProblemConfig::getDimension, &ProblemConfig::setDimension)
            .def_property("Prefine", &ProblemConfig::isPrefine, &ProblemConfig::setPrefine)
            .def_property("Steklovexample", &ProblemConfig::isSteklovexample, &ProblemConfig::setSteklovexample)
            .def_property("GalvisExample", &ProblemConfig::isGalvisExample, &ProblemConfig::setGalvisExample)
            .def_property("TensorNonConst", &ProblemConfig::isTensorNonConst, &ProblemConfig::setTensorNonConst)
            .def_property("MeshNonConvex", &ProblemConfig::isMeshNonConvex, &ProblemConfig::setMeshNonConvex)
            .def_property("Alpha", &ProblemConfig::getAlpha, &ProblemConfig::setAlpha)
            .def_property("DirName", &ProblemConfig::getDirName, &ProblemConfig::setDirName)
            .def_property("Problemname", &ProblemConfig::getProblemname, &ProblemConfig::setProblemname)
            .def_property("Materialids", &ProblemConfig::getMaterialids, &ProblemConfig::setMaterialids)
            .def_property("Bcmaterialids", &ProblemConfig::getBcmaterialids, &ProblemConfig::setBcmaterialids)
            .def_property("Exact", &ProblemConfig::getExact, &ProblemConfig::setExact);

    py::class_<TLaplaceExample1> laplaceExample (m, "TLaplaceExample");
    laplaceExample.def_readwrite("Exact", &TLaplaceExample1::fExact);

    py::enum_<TLaplaceExample1::EExactSol>(laplaceExample, "ExactSol")
            .value("ENone", TLaplaceExample1::EExactSol::ENone)
            .value("EConst", TLaplaceExample1::EExactSol::EConst)
            .value("EX", TLaplaceExample1::EExactSol::EX)
            .value("ESinSin", TLaplaceExample1::EExactSol::ESinSin)
            .value("ECosCos", TLaplaceExample1::EExactSol::ECosCos)
            .value("EArcTan", TLaplaceExample1::EExactSol::EArcTan)
            .value("EArcTanSingular", TLaplaceExample1::EExactSol::EArcTanSingular)
            .value("ESinDist", TLaplaceExample1::EExactSol::ESinDist)
            .value("E10SinSin", TLaplaceExample1::EExactSol::E10SinSin)
            .value("ESinSinDirNonHom", TLaplaceExample1::EExactSol::ESinSinDirNonHom)
            .value("ESinMark", TLaplaceExample1::EExactSol::ESinMark)
            .value("ESteklovNonConst", TLaplaceExample1::EExactSol::ESteklovNonConst)
            .value("EGalvisNonConst", TLaplaceExample1::EExactSol::EGalvisNonConst)
            .value("EBoundaryLayer", TLaplaceExample1::EExactSol::EBoundaryLayer)
            .value("EBubble", TLaplaceExample1::EExactSol::EBubble)
            .value("ESinCosCircle", TLaplaceExample1::EExactSol::ESinCosCircle)
            .export_values();
}