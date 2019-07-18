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

TPZGeoMesh *CreateGeoMesh2D();

void hAdaptivity(TPZCompMesh * postProcessMesh, TPZGeoMesh * gmeshToRefine);

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
    config.hdivmais = 0;
    config.ndivisions = 1;
    config.dimension = 2;
    config.prefine = false;
    config.makepressurecontinuous = true;

    config.exact.fExact = TLaplaceExample1::ESinSin;
    config.problemname = "AdaptivityTest";
    config.dir_name = "ESinSinAdaptivity";

    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    // Creates geometric mesh
    TPZGeoMesh *gmesh = nullptr;

    if (readGeoMeshFromFile) {
        TPZGmshReader gmsh;
        std::string meshfilename = "../LCircle.msh";

        gmsh.GetDimNamePhysical()[1]["dirichlet"] = 2;
        gmsh.GetDimNamePhysical()[2]["domain"] = 1;

        config.materialids.insert(1);
        config.bcmaterialids.insert(2);
        config.bcmaterialids.insert(3);

        gmsh.SetFormatVersion("4.1");
        gmsh.PrintPartitionSummary(std::cout);

        config.gmesh = gmsh.GeometricGmshMesh(meshfilename);
        config.gmesh->SetDimension(2);
    }
    else {
        config.materialids.insert(1);
        config.bcmaterialids.insert(2);
        config.bcmaterialids.insert(3);
        config.gmesh = CreateGeoMesh2D();
    }
    gmesh = config.gmesh;

    UniformRefinement(config.ndivisions, gmesh);
    int nelem = gmesh->NElements();

#ifdef PZDEBUG
    {
        std::ofstream out("GeometricMeshMark.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
#endif

    TPZMultiphysicsCompMesh *mixedMesh = nullptr;

    mixedMesh = CreateHDivMesh(config); //Hdiv x L2
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

    std::cout << "---Original PerifericalMaterialID --- " << std::endl;
    std::cout << " LagrangeInterface = " << hybrid.fLagrangeInterface << std::endl;
    std::cout << " HDivWrapMatID = " << hybrid.fHDivWrapMatid << std::endl;
    std::cout << " InterfaceMatID = " << hybrid.fInterfaceMatid << std::endl;

    mixedMesh = (HybridMesh); // Substitute mixed by hybrid mesh
    mixedMeshVector[0] = (HybridMesh)->MeshVector()[0]; // Hdiv mesh
    mixedMeshVector[1] = (HybridMesh)->MeshVector()[1]; // L2 mesh

    // Solves problem
    SolveHybridProblem(mixedMesh, hybrid.fInterfaceMatid, config);

    // reconstroi potencial e calcula o erro
    {
        TPZHDivErrorEstimatorH1 HDivEstimate(*mixedMesh);
        HDivEstimate.fProblemConfig = config;
        HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
        HDivEstimate.SetAnalyticSolution(config.exact);

        HDivEstimate.fperformUplift = true;
        HDivEstimate.fUpliftOrder = 1;

        HDivEstimate.PotentialReconstruction();

        TPZManVector<REAL> elementerrors;
        HDivEstimate.ComputeErrors(elementerrors);

        hAdaptivity(&HDivEstimate.fPostProcMesh, gmesh);
    }

    delete mixedMesh;
    delete mixedMeshVector[0];
    delete mixedMeshVector[1];

    return 0;
}


TPZGeoMesh *CreateLCircularBlendMesh() {
    TPZGeoMesh *gmesh = nullptr;

    return gmesh;
}

TPZGeoMesh *CreateGeoMesh2D() {

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

void hAdaptivity(TPZCompMesh * postProcessMesh, TPZGeoMesh * gmeshToRefine) {

    // Column of the flux error estimate on the element solution matrix
    const int fluxErrorEstimateCol = 3;

    int64_t nelem = postProcessMesh->ElementSolution().Rows();

    // Iterates through element errors to get the maximum value
    REAL maxError = 0.;
    for (int64_t iel = 0; iel < nelem; iel++) {
        TPZCompEl * cel = postProcessMesh->ElementVec()[iel];
        if (!cel) continue;
        if (cel->Dimension() != 2) continue;
        REAL elementError = postProcessMesh->ElementSolution()(iel, fluxErrorEstimateCol);

        std::cout << "geo element: " << cel->Reference()->Id() << " error: " << elementError << std::endl;

        if (elementError > maxError) {
            maxError = elementError;
        }
    }

    // Refines elements which error are bigger than 30% of the maximum error
    REAL threshold = 0.8 * maxError; // TODO voltar pra 0.3

    for (int64_t iel = 0; iel < nelem; iel++) {
        TPZCompEl * cel = postProcessMesh->ElementVec()[iel];
        if (!cel) continue;
        if (cel->Dimension() != 2) continue;
        REAL elementError = postProcessMesh->ElementSolution()(iel, fluxErrorEstimateCol);
        if (elementError > threshold) {
            TPZGeoEl * gel = cel->Reference();
            int iel  = gel->Id();

            TPZVec<TPZGeoEl *> sons;
            TPZGeoEl *gelToRefine = gmeshToRefine->ElementVec()[iel];
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

    {
        std::string meshFileName = "gmeshHopefullyRefined.vtk";
        std::ofstream outVTK(meshFileName.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshToRefine, outVTK, true);
        outVTK.close();
    }

}
