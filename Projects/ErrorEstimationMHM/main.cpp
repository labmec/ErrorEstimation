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
    config.TensorNonConst = false;
    config.MeshNonConvex = true;

    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;

<<<<<<< Updated upstream
    config.problemname = "NonConvexMesh";

    config.dir_name = "NonConvexMHM";

    config.problemname = "ESinSin";

    config.dir_name= "MHMTesteCilamce";

    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    for (int orderp = 1; orderp < 2; orderp++) {
        config.porder = orderp;

        for (int hdivmais = 1; hdivmais < 2; hdivmais++) {
            config.hdivmais = hdivmais;

        
        for(int hdivmais=1 ; hdivmais<2; hdivmais++){
        config.hdivmais = hdivmais;
            for(int ndiv=2; ndiv<3 ; ndiv++) {
          


            config.ndivisions = ndiv;




            if(IsgmeshReader){
                std::string meshfilename = "../LMesh.msh";

                TPZGmshReader gmsh;
                gmsh.GetDimNamePhysical()[2]["dirichlet"] = 2;
                gmsh.GetDimNamePhysical()[3]["domain"] = 1;



            for (int ndiv = 0; ndiv < 1; ndiv++) {
                config.ndivisions = ndiv;
                config.materialids.insert(1);
                config.bcmaterialids.insert(-1);

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
