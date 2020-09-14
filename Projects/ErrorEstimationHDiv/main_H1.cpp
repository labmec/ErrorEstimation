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




bool IsgmeshReader = true;
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
    
    for(int ndiv=1; ndiv<2; ndiv++) {
        ProblemConfig config;

        config.porder = 1;
        config.hdivmais = 0;
        config.ndivisions = ndiv;
        config.dimension = 2;
        config.prefine = false;
        config.makepressurecontinuous = true;

        config.exact = new TLaplaceExample1;
        config.exact.operator*().fExact = TLaplaceExample1::
            ESinMark; // ESinSinDirNonHom;//ESinSin;//ESinMark;//EX;//EConst;//EArcTanSingular;//EArcTan;//
        config.problemname =
            "ESinSinCircle k=1 e n=0 Up=2"; //"EConst";//"ESinSinDirNonHom";//"ESinSin";//"
                                            ////"EArcTanSingular_PRef";//""ArcTang";//

        config.dir_name = "ESinMark";
        // config.dir_name= "ESinSin";
        std::string command = "mkdir " + config.dir_name;
        system(command.c_str());

        // geometric mesh

        TPZGeoMesh *gmesh = Tools::ReadGeometricMesh(config, IsgmeshReader);


        Tools::UniformRefinement(config.ndivisions, gmesh);
   // RandomRefinement(config, config.ndivisions);
    
#ifdef PZDEBUG
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmeshInitial.txt");
        gmesh->Print(out2);
        
    }
#endif


        
    TPZMultiphysicsCompMesh *cmesh_HDiv = Tools::CreateHDivMesh(config);//Hdiv x L2
    cmesh_HDiv->InitializeBlock();
    
#ifdef PZDEBUG2
    {
        
        std::ofstream out2("MalhaMista.txt");
        cmesh_HDiv->Print(out2);
        
    }
#endif

        TPZMultiphysicsCompMesh *hybridmesh = Tools::HybridSolveProblem(cmesh_HDiv, config);

    //reconstroi potencial e calcula o erro
    {
        // TPZHybridHDivErrorEstimator HDivEstimate(*hybridmesh);
        TPZHDivErrorEstimatorH1 HDivEstimate(*hybridmesh);
        HDivEstimate.fProblemConfig = config;
        HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
        
        HDivEstimate.SetAnalyticSolution(config.exact);
        
        HDivEstimate.fperformUplift = true;
        HDivEstimate.fUpliftOrder = 2;
        
        HDivEstimate.PotentialReconstruction();
        
        TPZManVector<REAL> elementerrors;
        TPZManVector<REAL> errvec;
        bool store = true;
        HDivEstimate.ComputeErrors(errvec,elementerrors,store);
        
    }
    
    TPZManVector<TPZCompMesh *,2> meshvector;
    meshvector = hybridmesh->MeshVector();
    delete hybridmesh;
    delete meshvector[0];
    delete meshvector[1];
    //return 0;
    }
}

