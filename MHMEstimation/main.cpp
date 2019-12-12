#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <fstream>
#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"
#include "pzreal.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZParSkylineStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "TPZGeoLinear.h"

#include "pzpoisson3d.h"
//#include "pzhybridpoisson.h"
#include "pzpoisson3dreferred.h"
#include "mixedpoisson.h"
#include "pzelasmat.h"
#include "pzelasthybrid.h"
#include "pzmat1dlin.h"
#include "TPZVecL2.h"
#include "TPZMatLaplacianHybrid.h"
#include "TPZLagrangeMultiplier.h"

#include "pzbuildmultiphysicsmesh.h"
#include "pzelementgroup.h"
#include "TPZCompMeshTools.h"
#include "pzcondensedcompel.h"
#include "pzfunction.h"
#include "pzgraphmesh.h"
#include "pzfmatrix.h"

#include "pzlog.h"

#include "TPZVTKGeoMesh.h"
#include "pzvisualmatrix.h"
#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "pzcheckgeom.h"

#include "TPZMHMeshControl.h"
#include "TPZMHMixedMeshControl.h"
#include "TPZMHMixedMeshChannelControl.h"

#include "TPZMHMixedHybridMeshControl.h"
#include "TPZHybridizeHDiv.h"
#include "meshgen.h"
//#include "ConfigCasesMaze.h"
#include <iostream>
#include <string>
//#include <opencv2/opencv.hpp>
#include <math.h>
#include <set>
#include "pzsolve.h"

#include "TPZPersistenceManager.h"

#include "TPZMHMHDivErrorEstimator.h"
#include "Tools.h"
#include "ToolsMHM.h"

using namespace std;



bool IsgmeshReader = false;


int main() {
    InitializePZLOG();


    for(int ndiv=1; ndiv<2 ; ndiv++) {
    ProblemConfig config;

    config.porder = 1;
    config.hdivmais = 1;
    config.ndivisions = ndiv;
    config.dimension = 2;
    config.prefine=false;
    config.TensorNonConst = false;
    config.MeshNonConvex = false ;

    TLaplaceExample1 example;

        config.exact.fExact = example.EX;//ESinSin;

    config.problemname = "EXHdiv";

    config.dir_name= "TestesMHM";
    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    TPZGeoMesh *gmesh = nullptr;

    if(IsgmeshReader){
        std::string meshfilename = "../Cube.msh";

        TPZGmshReader gmsh;
        gmsh.GetDimNamePhysical()[2]["dirichlet"] = 2;
        gmsh.GetDimNamePhysical()[3]["domain"] = 1;


        config.materialids.insert(1);
        config.bcmaterialids.insert(2);


        gmsh.SetFormatVersion("4.1");
        gmesh = gmsh.GeometricGmshMesh(meshfilename);
        gmsh.PrintPartitionSummary(std::cout);
        gmesh->SetDimension(config.dimension);
        config.gmesh = gmesh;

    }

    else {

        TPZManVector<int,4> bcids(4,-1);


        TPZVec<int64_t> coarseindices;

        TPZManVector<REAL,3> x0(3,0.), x1(3,1.);
        x1[2] = 0.;


        int nx = pow(2, ndiv);
        gmesh = CreateGeoMesh(nx, bcids);//MalhaGeomFredQuadrada(nx, nx,x0, x1, coarseindices, 1);// CreateCircleGeoMesh();//CreateGeoMesh(nx, bcids);
        config.gmesh = gmesh;
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        gmesh->SetDimension(config.dimension);

    }


#ifdef PZDEBUG
        {
            std::ofstream out("GmeshMHM_Coarse.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
            std::ofstream out2("gmeshInitial.txt");
            gmesh->Print(out2);

        }
#endif

        UniformRefinement(ndiv, gmesh);
        DivideLowerDimensionalElements(gmesh);
        
        TPZManVector<TPZCompMesh*, 2> meshvec_HDiv(2, 0);
        
        TPZMultiphysicsCompMesh *cmesh_HDiv=nullptr;
        
        
        cmesh_HDiv = CreateHDivMesh(config);//Hdiv x L2
        cmesh_HDiv->InitializeBlock();
        
#ifdef PZDEBUG
        {
            
            std::ofstream out2("MalhaMista.txt");
            cmesh_HDiv->Print(out2);
            
        }
#endif
        
        
        
       // SolveMixedProblem(cmesh_HDiv,config);
        
        
        
        MHMTest(config);
        
        
        
    }
    
    // return 0;
    
    
}
