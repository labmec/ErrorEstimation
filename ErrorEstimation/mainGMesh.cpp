 /**
 * @file This file implements an error estimator in space H1.
 */

#include "pzlog.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzbndcond.h"
#include "pzpoisson3d.h"
#include "pzgeoelbc.h"
#include "TPZSSpStructMatrix.h"
#include "TPZRefPatternDataBase.h"
#include "TPZVTKGeoMesh.h"
#include "TPZGmshReader.h"
#include "TPZAnalyticSolution.h"
#include "TPZFrontSym.h"

#include <fstream>
#include <ctime>
#include <cstdio>
#include <cmath>

#include "ProblemConfig.h"
#include "TPZPostProcessError.h"

// Global variables
const int problemDimension = 2;
const bool readGMeshFromFile = true;

const int matID = 1;
const int dirichletMatID = -1;
const int neumannMatID = -2;
const int dirichlet = 0;
const int neumann = 1;

// Functions declarations
TPZGeoMesh *CreateGeoMesh();
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder);
void UniformRefinement(int nDiv, TPZGeoMesh *gmesh);
bool SolvePoissonProblem(struct SimulationCase &sim_case);
bool PostProcessProblem(TPZAnalysis &an, TPZGeoMesh * gmesh, TPZCompMesh * pressuremesh);

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.refine"));
#endif

// Laplace equation on square 1D 2D 3D - Volker John article 2000
int main(int argc, char *argv[]) {

#ifdef LOG4CXX
	InitializePZLOG();
#endif

    // Initializing uniform refinements patterns.
    gRefDBase.InitializeAllUniformRefPatterns();

    TPZGeoMesh *gmesh = nullptr;

	if (readGMeshFromFile) {
        TPZGmshReader gmsh;

        // Assigns IDs of 1D and 2D elements defining boundary conditions.
        gmsh.fPZMaterialId[1]["dirichlet"] = -1;
        gmsh.fPZMaterialId[1]["neumann"] = -2;
        gmsh.fPZMaterialId[2]["dirichlet"] = -1;
        gmsh.fPZMaterialId[2]["neumann"] = -2;

        // Assigns IDs of 2D and 3D elements defining the problem domain.
        gmsh.fPZMaterialId[2]["domain"] = 1;
        gmsh.fPZMaterialId[3]["domain"] = 1;

#ifdef MACOSX
        gmesh = gmsh.GeometricGmshMesh("../Quad.msh");
#else
        gmesh = gmsh.GeometricGmshMesh("Quad.msh");
#endif
	}
	else {
        gmesh = CreateGeoMesh();
	}

    struct SimulationCase Case1;
    Case1.nthreads = 0;
    Case1.numinitialrefine = 0;
    Case1.dir_name = "QuadSimCase1";
    Case1.gmesh = gmesh;

    if(!SolvePoissonProblem(Case1)) {
        return 1;
    }

    return 0;
}

 TPZGeoMesh *CreateGeoMesh() {

     TPZGeoMesh* gmesh = new TPZGeoMesh();

     if (problemDimension == 2) {

         gmesh->SetDimension(2);

         // Creates matrix with quadrilateral node coordinates.
         const int quadNodeNumber = 4;
         REAL coordinates[quadNodeNumber][3] = {
             {0., 0., 0.},
             {1., 0., 0.},
             {1., 1., 0.},
             {0., 1., 0.}
         };

         // Inserts coordinates in the TPZGeoMesh object.
         for(int i = 0; i < quadNodeNumber; i++) {
             int64_t nodeID = gmesh->NodeVec().AllocateNewElement();

             TPZVec<REAL> nodeCoord(3);
             nodeCoord[0] = coordinates[i][0];
             nodeCoord[1] = coordinates[i][1];
             nodeCoord[2] = coordinates[i][2];

             gmesh->NodeVec()[nodeID] = TPZGeoNode(i, nodeCoord, *gmesh);
         }

         // Creates quadrilateral element.
         int64_t index;
         TPZManVector<int64_t> nodeIDs(quadNodeNumber);

         for(int n = 0; n < quadNodeNumber; n++) {
            nodeIDs[n] = n;
         }
         gmesh->CreateGeoElement(EQuadrilateral, nodeIDs, matID, index);

         // Creates line elements where boundary conditions will be inserted.
         nodeIDs.Resize(2);
         for (int i = 0; i < 4; i++) {

             nodeIDs[0] = i % 4;
             nodeIDs[1] = (i + 1) % 4;

             gmesh->CreateGeoElement(EOned, nodeIDs, -1, index);
         }

         gmesh->BuildConnectivity();
         return gmesh;
     }
     else if (problemDimension == 3) {
         // Creates hexahedra element.
         gmesh->SetDimension(3);
         const int hexahedraNodes = 8;
         REAL coordinates[hexahedraNodes][3] = {
             {0., 0., 0.},
             {1., 0., 0.},
             {1., 1., 0.},
             {0., 1., 0.},
             {0., 0., 1.},
             {1., 0., 1.},
             {1., 1., 1.},
             {0., 1., 1.},
         };

         for(int n = 0; n < hexahedraNodes; n++) {
             int64_t nodeID = gmesh->NodeVec().AllocateNewElement();
             TPZVec<REAL> coord(3);
             coord[0] = coordinates[n][0];
             coord[1] = coordinates[n][1];
             coord[2] = coordinates[n][2];
             gmesh->NodeVec()[nodeID] = TPZGeoNode(n, coord, *gmesh);
         }

         TPZVec<int64_t> nodeID(hexahedraNodes);
         for(int n = 0; n < hexahedraNodes; n++) {
             nodeID[n] = n;
         }

         // Inserts Dirichlet BC
// CODE dsahdajshdkaj
         gmesh->BuildConnectivity();
         return gmesh;
     }
     else {
         DebugStop();
     }
 }
    
bool SolvePoissonProblem(struct SimulationCase &sim_case) {

    // Creating the directory
    std::string command = "mkdir " + sim_case.dir_name;
    system(command.c_str());

    // Output files
    std::string file_name = sim_case.dir_name + "/" + "ErrorsHP_Poisson.txt";
    std::ofstream fileerrors(file_name, ios::app);   // To store all errors calculated by TPZAnalysis (PosProcess)

    // Initializing the auto adaptive process
    TPZVec<REAL> ervec, ErrorVec(100, 0.0);
    TPZVec<int64_t> NEquations(100, 0L);
    TPZVec<REAL> ervecbyel;
    TPZVec<REAL> gradervecbyel;

    /** Variable names for post processing */
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("POrder");
    scalnames.Push("Pressure");
    vecnames.Push("Derivative");

    fileerrors.flush();

    TPZGeoMesh *gmesh = sim_case.gmesh;
    {
        // Refines an element
        UniformRefinement(sim_case.numinitialrefine, gmesh);
    }

    // Creates computational mesh (approximation space and materials)
    TPZCompEl::SetgOrder(sim_case.porder);
    gmesh->SetName("Original GeoMesh");

    TPZManVector<TPZCompMesh *> meshvec(0);

    TPZCompMesh *pressuremesh = CMeshPressure(gmesh, sim_case.porder);
    pressuremesh->AdjustBoundaryElements();

    {
        std::ofstream out("CompMesh.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(pressuremesh, out);
        std::ofstream out2("CompMesh.txt");
        pressuremesh->Print(out2);
    }

    TLaplaceExample1 example;
    example.fExact = TLaplaceExample1::ESinSin;
    example.fDimension = gmesh->Dimension();
    example.fSignConvention = -1;

    {
        for (auto it:pressuremesh->MaterialVec()) {
            TPZMaterial *mat = it.second;
            TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
            if (!bc) {
                mat->SetForcingFunction(example.ForcingFunction());
            }
            else {
                bc->SetForcingFunction(0, example.Exact());
            }
        }
    }

    TPZAnalysis an(pressuremesh, true);
    an.SetExact(example.ExactSolution());

    {
        std::stringstream sout;
        sout << sim_case.dir_name << "/" << "Poisson" << gmesh->Dimension() << "numref" << sim_case.numinitialrefine
             << "Porder" << sim_case.porder << ".vtk";
        an.DefineGraphMesh(gmesh->Dimension(), scalnames, vecnames, sout.str());
    }

    pressuremesh->SetName("Adapted CompMesh");

    // Printing geometric and computational mesh
#ifdef PZDEBUG
    {
        std::ofstream out("../PressureGeoMesh.txt");
        pressuremesh->Reference()->Print(out);
    }
#endif

#ifdef USING_MKL
    // Solves using a symmetric matrix then using Cholesky decomposition (direct method)
    TPZSymetricSpStructMatrix strmat(pressuremesh);
    strmat.SetNumThreads(sim_case.nthreads);
    an.SetStructuralMatrix(strmat);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(pressuremesh);
    strmat.SetNumThreads(sim_case.nthreads);
    strmat.SetDecomposeType(ECholesky);
#endif

    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ECholesky);
    an.SetSolver(*direct);
    delete direct;

    an.Run();

    PostProcessProblem(an, gmesh, pressuremesh);

    return true;
}

bool PostProcessProblem(TPZAnalysis &an, TPZGeoMesh * gmesh, TPZCompMesh * pressuremesh) {
    // Post processing
    an.PostProcess(1, gmesh->Dimension());
    TPZPostProcessError error(pressuremesh);
    
    TPZVec<STATE> estimatedelementerror, exactelementerror;
    error.ComputeElementErrors(estimatedelementerror);

    {
        int64_t nels = pressuremesh->ElementVec().NElements();
        pressuremesh->ElementSolution().Redim(nels, 6);
    }

    bool store_errors = true;
    an.PostProcessError(exactelementerror, store_errors);
    std::cout << "Exact error " << exactelementerror << std::endl;
    gmesh->ResetReference();
    error.MultiPhysicsMesh()->LoadReferences();

    {
        TPZFMatrix<STATE> true_elerror(pressuremesh->ElementSolution());
        TPZFMatrix<STATE> estimate_elerror(error.MultiPhysicsMesh()->ElementSolution());
        int64_t nel = true_elerror.Rows();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = pressuremesh->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            TPZCompEl *mphys = gel->Reference();
            int64_t elindex2 = mphys->Index();
            true_elerror(el,0) = estimate_elerror(elindex2,2);
            true_elerror(el,1) = true_elerror(el,2);
            if (true_elerror(el,1) > 1.e-8) {
                true_elerror(el,2) = true_elerror(el,0)/true_elerror(el,1);
            }
        }
        pressuremesh->ElementSolution() = true_elerror;
    }
    
    {
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("State");
        scalnames.Push("Error");
        scalnames.Push("TrueError");
        scalnames.Push("EffectivityIndex");
        an.DefineGraphMesh(pressuremesh->Dimension(), scalnames, vecnames, "Errors.vtk");
        an.PostProcess(1);
    }
    
    if(gmesh) delete gmesh;

    return true;
}

void UniformRefinement(int nDiv, TPZGeoMesh *gmesh) {

	TPZManVector<TPZGeoEl*> children;
    for(int division = 0; division < nDiv; division++) {

        int64_t nels = gmesh->NElements();

        for(int64_t elem = 0; elem < nels; elem++) {

            TPZGeoEl * gel = gmesh->ElementVec()[elem];

			if(!gel || gel->HasSubElement()) continue;
			if(gel->Dimension() == 0) continue;
            gel->Divide(children);
        }
    }
}

TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder) {

    int dim = gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    // Creates Poisson material
    TPZMatPoisson3d *material = new TPZMatPoisson3d(matID, dim);

    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->InsertMaterialObject(material);

//    TPZMaterial * mat(material);
//    cmesh->InsertMaterialObject(mat);
    
    // Inserts boundary conditions
    TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);
    TPZMaterial * BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(material, -2, neumann, val1, val2);

    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);

    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetAllCreateFunctionsContinuous();

    // Adjusts computational data structure
    cmesh->AutoBuild();
    return cmesh;
}
