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

#include "TPZHDivErrorEstimatorH1.h"

#include "Tools.h"

#include "TPZBFileStream.h"
#include <tuple>
#include <memory>




bool IsgmeshReader = false;
bool neumann = true;

bool mixedsolution = false;


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
    config.ndivisions = 1;
    config.dimension = 2;
    config.prefine=false;
    config.makepressurecontinuous = true;
    
    config.exact.fExact = TLaplaceExample1::ESinSinDirNonHom;//EConst;//;//ESinSin;//ESinMark;//EArcTanSingular;//EArcTan;//
    config.problemname = "ESinSinDirNonHom";//"EConst";////"ESinSin";//" ESinMark";////"EArcTanSingular_PRef";//""ArcTang";//
    
    //config.dir_name= "ESinSinDirNonHom";
    config.dir_name= "EConst";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    
    //geometric mesh

     TPZGeoMesh *gmesh = ReadGeometricMesh(config, IsgmeshReader);
    

    UniformRefinement(config.ndivisions, gmesh);
   // RandomRefine(config, config.ndivisions);
    
#ifdef PZDEBUG
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmeshInitial.txt");
        gmesh->Print(out2);
        
    }
#endif


        
    TPZMultiphysicsCompMesh *cmesh_HDiv = CreateHDivMesh(config);//Hdiv x L2
    cmesh_HDiv->InitializeBlock();
    
#ifdef PZDEBUG
    {
        
        std::ofstream out2("MalhaMista.txt");
        cmesh_HDiv->Print(out2);
        
    }
#endif
    
    TPZMultiphysicsCompMesh *hybridmesh= HybridSolveProblem(cmesh_HDiv, config);

    //reconstroi potencial e calcula o erro
    {
        // TPZHybridHDivErrorEstimator HDivEstimate(*hybridmesh);
        TPZHDivErrorEstimatorH1 HDivEstimate(*hybridmesh);
        HDivEstimate.fProblemConfig = config;
        HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
        HDivEstimate.SetAnalyticSolution(config.exact);
        
        HDivEstimate.PotentialReconstruction();
        
        TPZManVector<REAL> elementerrors;
        HDivEstimate.ComputeErrors(elementerrors);
        
    }
    
    TPZManVector<TPZCompMesh *,2> meshvector;
    meshvector = hybridmesh->MeshVector();
    delete hybridmesh;
    delete meshvector[0];
    delete meshvector[1];
    return 0;
}

