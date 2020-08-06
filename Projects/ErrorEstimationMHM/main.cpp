#include "pzlog.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"

#include <string>
#include <cmath>
#include <set>

#include "Tools.h"
#include "meshgen.h"
#include "ToolsMHM.h"
#include "TPZMHMHDivErrorEstimator.h"

using namespace std;

int main() {
    InitializePZLOG();

    ProblemConfig config;
    config.dimension = 2;
    config.prefine = false;
    bool IsgmeshReader = false;
    config.TensorNonConst = false;
    config.MeshNonConvex = true;

    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;

    config.problemname = "NonConvexMesh";

    config.dir_name = "NonConvexMHM";

    config.problemname = "ESinSin";

    config.dir_name= "MHMTesteCilamce";

    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());
    int dim =2;
    
    
    //    TPZManVector<int,4> bcids(8,-3);
    //    bcids[0] = -1;
        TPZManVector<int,4> bcids(4,-1);
    //config.bcmaterialids.insert(-3);//Robin
    config.bcmaterialids.insert(-1);//Dirichlet
   // config.bcmaterialids.insert(-2);//Newmann
    

    for (int orderp = 1; orderp < 2; orderp++) {
        config.porder = orderp;
//
        for (int hdivmais = 1; hdivmais < 2; hdivmais++) {
            config.hdivmais = hdivmais;
            

            for(int ndiv=2; ndiv<3 ; ndiv++) {

            config.ndivisions = ndiv;
                
            TPZGeoMesh *gmesh = nullptr;

            if(IsgmeshReader){
                std::string meshfilename = "../LMesh.msh";

                TPZGmshReader gmsh;
                gmsh.GetDimNamePhysical()[2]["dirichlet"] = 2;
                gmsh.GetDimNamePhysical()[3]["domain"] = 1;
            }
            else{
                
            int nel =pow(2, ndiv);

             gmesh =CreateGeoMesh(nel, bcids); // CreateLShapeMesh(bcids);//CreateSingleTriangleMesh(bcids);//CreateQuadMeshRefTriang(bcids); // CreateTrapezoidalMesh(nelT,
                              // nelT, 1.,1.,bcids);//CreateLCircleGeoMesh();//
             config.materialids.insert(1);

             config.gmesh = new TPZGeoMesh;
             *config.gmesh = *gmesh;
            gmesh->SetDimension(dim);
                
            #ifdef PZDEBUG
                    {
                        std::ofstream out("CoarseMesh.vtk");
                        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

                    }
            #endif
                
                
            }

//
                UniformRefinement(1, gmesh);
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
       }
    }
}
