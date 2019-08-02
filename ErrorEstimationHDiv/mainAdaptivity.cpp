//
// Created by Gustavo Batistela on 12/07/19.
//

#include "TPZRefPatternDataBase.h"
#include "TPZGmshReader.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "TPZGeoLinear.h"

#include <math.h>
#include <tuple>

#include "ProblemConfig.h"
#include "Tools.h"
#include "TPZHybridHDivErrorEstimator.h"
#include "TPZHDivErrorEstimatorH1.h"

bool readGeoMeshFromFile = false;

TPZGeoMesh *CreateGeoMesh();
TPZGeoMesh *CreateLCircleGeoMesh();


int main(int argc, char *argv[]) {

#ifdef LOG4CXX
    InitializePZLOG();
#endif

    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);

    // Creates geometric mesh
    TPZGeoMesh *gmeshOriginal = CreateGeoMesh();
    
#ifdef PZDEBUG
    {
        std::ofstream out("Gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshOriginal, out);
        
    }
#endif
    
    

    int refinementSteps = 5;

    // Copies meshes to be used with each proposal
    TPZGeoMesh *hybridEstimatorMesh = new TPZGeoMesh();
    *hybridEstimatorMesh = *gmeshOriginal;
    TPZGeoMesh *markEstimatorMesh = new TPZGeoMesh();
    *markEstimatorMesh = *gmeshOriginal;
    
    
    
#ifdef PZDEBUG
    {
        std::ofstream out("GmeshPosCopy.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(hybridEstimatorMesh, out);
        
    }
#endif

    // Run tests with hdiv++ proposal
    
    for (int i = 0; i < refinementSteps; i++) {
        ProblemConfig config;
        config.dir_name = "AdaptivityHybridSinSin";
        config.adaptivityStep = i;

//        config.gmesh = new TPZGeoMesh();
        config.gmesh = hybridEstimatorMesh;

        config.materialids.insert(1);
        config.bcmaterialids.insert(2);
        config.bcmaterialids.insert(3);

        config.porder = 1;
        config.hdivmais = 1;
        config.dimension = 2;
        config.prefine = false;
        config.makepressurecontinuous = true;

        config.exact.fExact = TLaplaceExample1::ESinSin;
        config.problemname = "AdaptivityTest";

        std::string command = "mkdir " + config.dir_name;
        system(command.c_str());

        TPZMultiphysicsCompMesh *mixedMesh = nullptr;

        mixedMesh = CreateHDivMesh(config); // H(div) x L2
        mixedMesh->InitializeBlock();

        TPZManVector<TPZCompMesh *, 2> mixedMeshVector(2, 0);
        mixedMeshVector = mixedMesh->MeshVector();

        // Hybridizes mixed mesh
        TPZHybridizeHDiv hybrid;
        auto HybridMesh = hybrid.Hybridize(mixedMesh);
        HybridMesh->CleanUpUnconnectedNodes();
        HybridMesh->AdjustBoundaryElements();
        delete mixedMesh;
        delete mixedMeshVector[0];
        delete mixedMeshVector[1];

        mixedMesh = (HybridMesh); // Substitute mixed by hybrid mesh
        mixedMeshVector[0] = (HybridMesh)->MeshVector()[0]; // Hdiv mesh
        mixedMeshVector[1] = (HybridMesh)->MeshVector()[1]; // L2 mesh

        // Solves problem
        SolveHybridProblem(mixedMesh, hybrid.fInterfaceMatid, config);

        // Reconstructs pressure and calculates error
        TPZHybridHDivErrorEstimator HDivEstimate(*mixedMesh);
        HDivEstimate.fProblemConfig = config;
        HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
        HDivEstimate.SetAnalyticSolution(config.exact);

        HDivEstimate.PotentialReconstruction();

        TPZManVector<REAL> elementErrors;
        HDivEstimate.ComputeErrors(elementErrors);

//        delete mixedMesh;
//        delete mixedMeshVector[0];
//        delete mixedMeshVector[1];

        hAdaptivity(&HDivEstimate.fPostProcMesh, hybridEstimatorMesh);
    }
    return 0;
    // Run tests with Ainsworth's proposal
    for (int i = 0; i < refinementSteps; i++) {

        ProblemConfig config;
        config.dir_name = "AdaptivityMarkSin";
        config.adaptivityStep = i;

//        config.gmesh = new TPZGeoMesh();
        config.gmesh = markEstimatorMesh;

        config.materialids.insert(1);
        config.bcmaterialids.insert(2);
        config.bcmaterialids.insert(3);

        config.porder = 1;
        config.hdivmais = 0;
        config.dimension = 2;
        config.prefine = false;
        config.makepressurecontinuous = true;

        config.exact.fExact = TLaplaceExample1::ESinMark;
        config.problemname = "AdaptivityTest";

        std::string command = "mkdir " + config.dir_name;
        system(command.c_str());
        
//        std::string plotname;
//        {
//            std::stringstream out;
//            out << config.dir_name << "/" << "GMeshT_" << i<<".vtk";
//            plotname = out.str();
//            TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, out);
//
//        }
        TPZMultiphysicsCompMesh *mixedMesh = nullptr;

        mixedMesh = CreateHDivMesh(config); // H(div) x L2
        mixedMesh->InitializeBlock();

        TPZMultiphysicsCompMesh *hybridmesh = HybridSolveProblem(mixedMesh, config);

        // Reconstructs pressure and calculates error
        TPZHDivErrorEstimatorH1 HDivEstimate(*hybridmesh);
        HDivEstimate.fProblemConfig = config;
        HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;

        HDivEstimate.SetAnalyticSolution(config.exact);

        HDivEstimate.fperformUplift = true;
        HDivEstimate.fUpliftOrder = 2;

        HDivEstimate.PotentialReconstruction();

        TPZManVector<REAL> elementErrors;
        HDivEstimate.ComputeErrors(elementErrors);

        hAdaptivity(&HDivEstimate.fPostProcMesh, markEstimatorMesh);
    }

    return 0;
}

TPZGeoMesh *CreateLCircleGeoMesh() {

    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);

    TPZVec<REAL> coord(3, 0.);

    // Inserts node at origin
    gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[0].Initialize(coord, *gmesh);

    // Inserts circumference nodes
    for (int64_t i = 0; i < 13; i++) {
        const REAL step = M_PI / 8;
        coord[0] = cos(i * step);
        coord[1] = sin(i * step);
        const int64_t newID = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[newID].Initialize(coord, *gmesh);
    }

    int matIdTriangle = 1, matIdArc = 2;

    // Inserts triangle elements
    TPZManVector<int64_t, 3> nodesIdVec(3);
    for (int64_t i = 0; i < 6; i++) {
        nodesIdVec[0] = 0;
        nodesIdVec[1] = 1 + 2 * i;
        nodesIdVec[2] = 3 + 2 * i;
        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle>>(nodesIdVec, matIdTriangle, *gmesh);
    }
    // Inserts arc elements
    for (int64_t i = 0; i < 6; i++) {
        nodesIdVec[0] = 1 + 2 * i;
        nodesIdVec[1] = 3 + 2 * i;
        nodesIdVec[2] = 2 + 2 * i;
        new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, matIdArc, *gmesh);
    }
    // Finally, inserts line elements to complete boundary
    nodesIdVec.Resize(2);
    nodesIdVec[0] = 0;
    nodesIdVec[1] = 1;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, matIdArc, *gmesh);

    nodesIdVec[0] = 0;
    nodesIdVec[1] = 13;
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, matIdArc, *gmesh);

    gmesh->BuildConnectivity();

    return gmesh;
}

<<<<<<< HEAD

=======
void hAdaptivity(TPZCompMesh *postProcessMesh, TPZGeoMesh *gmeshToRefine) {

    // Column of the flux error estimate on the element solution matrix
    const int fluxErrorEstimateCol = 3;
//    postProcessMesh->ElementSolution().Print(std::cout);
    int64_t nelem = postProcessMesh->ElementSolution().Rows();

    // Iterates through element errors to get the maximum value
    REAL maxError = 0.;
    for (int64_t iel = 0; iel < nelem; iel++) {

        TPZCompEl *cel = postProcessMesh->Element(iel);
        if (!cel) continue;
        if (cel->Dimension() != 2) continue;
        REAL elementError = postProcessMesh->ElementSolution()(iel, fluxErrorEstimateCol);


        if (elementError > maxError) {
            maxError = elementError;
        }
    }

    // Refines elements which error are bigger than 30% of the maximum error
    REAL threshold = 0.3 * maxError;

    for (int64_t iel = 0; iel < nelem; iel++) {
        TPZCompEl *cel = postProcessMesh->ElementVec()[iel];
        if (!cel) continue;
        if (cel->Dimension() != 2) continue;
        REAL elementError = postProcessMesh->ElementSolution()(iel, fluxErrorEstimateCol);
        if (elementError > threshold) {

            TPZGeoEl *gel = cel->Reference();
//            int eltorefine = gel->Id();

            TPZVec<TPZGeoEl *> sons;
            TPZGeoEl *gelToRefine = gmeshToRefine->ElementVec()[gel->Index()];
            if (gelToRefine && !gelToRefine->HasSubElement()) {
                gelToRefine->Divide(sons);
#ifdef LOG4CXX
                int nsides = gelToRefine->NSides();
                TPZVec<REAL> loccenter(3);
                TPZVec<REAL> center(3);
                gelToRefine->CenterPoint(nsides - 1, loccenter);

                gelToRefine->X(loccenter, center);
                static LoggerPtr logger(Logger::getLogger("HDivErrorEstimator"));
                if (logger->isDebugEnabled()) {
                    std::stringstream sout;
                    sout << "\nCenter coord: = " << center[0] << " " << center[1] << "\n";
                    sout << "Error = " << elementError << "\n\n";
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
            }
        }
    }
    DivideLowerDimensionalElements(gmeshToRefine);
}

TPZGeoMesh *CreateGeoMesh() {
    TPZGeoMesh * gmesh = nullptr;
    if (readGeoMeshFromFile) {
        TPZGmshReader gmsh;
        std::string meshfilename = "LCircle.msh";

        gmsh.GetDimNamePhysical()[1]["dirichlet"] = 2;
        gmsh.GetDimNamePhysical()[2]["domain"] = 1;

        gmsh.SetFormatVersion("4.1");
        gmsh.PrintPartitionSummary(std::cout);

        gmesh = gmsh.GeometricGmshMesh(meshfilename);
        gmesh->SetDimension(2);
    } else {
        gmesh = CreateLCircleGeoMesh();
    }
    int initialRefinement = 1;
    UniformRefinement(initialRefinement, gmesh);

#ifdef PZDEBUG
    {
        std::ofstream out("OriginalGeometricMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
#endif
    return gmesh;
}
