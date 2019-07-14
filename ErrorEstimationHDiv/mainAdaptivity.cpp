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
    config.ndivisions = 2;
    config.dimension = 2;
    config.prefine = false;
    config.makepressurecontinuous = true;

    config.exact.fExact = TLaplaceExample1::ESinMark;
    config.problemname = "AdaptivityMark";
    config.dir_name = "ESinMarkAdaptivity";

    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    // Creates geometric mesh
    TPZGeoMesh *gmesh = nullptr;

    if (readGeoMeshFromFile) {
        config.meshFileName = "LCircle.msh";
        config.InsertDomainMat("domain", 1, 2);
        config.InsertBCMat("dirichlet", 2, 1);
        config.ReadGeoMeshFromFile();
        gmesh = config.gmesh;
    }
    else {
        config.gmesh = CreateGeoMesh2D();
    }
    gmesh = config.gmesh;

    UniformRefinement(config.ndivisions, gmesh);

#ifdef PZDEBUG
    {
        std::ofstream out("GeometricMeshMark.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
#endif

    return 0;
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
    meshvec_HDiv[0] = (HybridMesh)->MeshVector()[0]; // malha Hdiv
    meshvec_HDiv[1] = (HybridMesh)->MeshVector()[1]; // malha L2

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

    return 0;
}


TPZGeoMesh *CreateLCircularBlendMesh() {
    TPZGeoMesh *gmesh = nullptr;

    return gmesh;
}

TPZGeoMesh *CreateGeoMesh2D() {

    TPZGeoMesh *gmesh = new TPZGeoMesh();

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

    {
        TPZVec<TPZGeoEl *> sons;
        for (int iDiv = 0; iDiv < 0; iDiv++) {
            const int nel = gmesh->NElements();
            for (int iel = 0; iel < nel; iel++) {
                TPZGeoEl *geo = gmesh->ElementVec()[iel];
                if (geo && !geo->HasSubElement()) {
                    geo->Divide(sons);
                }
            }
        }
    }

    {
        std::string meshFileName = "blendmesh2D.vtk";
        std::ofstream outVTK(meshFileName.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        outVTK.close();
    }
    return gmesh;
}
