 /**
 * @file Poisson 3D in hexahedra with shock problem
 */


#include "pzlog.h"
#include "TPZRefPatternDataBase.h"
//#include "pzvec.h"
//#include "pzadmchunk.h"
//#include "pzcmesh.h"
//#include "pzvec_extras.h"
//#include "pzcheckgeom.h"
//#include "pzcheckmesh.h"
//
//#include "pzgeoel.h"
//#include "pzgnode.h"
//#include "pzgeoelside.h"
//#include "pzgeoelbc.h"
//
//#include "pzintel.h"
//#include "pzcompel.h"
//
//#include "pzmatrix.h"

#include "pzanalysis.h"
//#include "pzfstrmatrix.h"
//#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
//#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
//#include "TPZFrontStructMatrix.h"
//#include "TPZParFrontStructMatrix.h"

//#include "TPZParSkylineStructMatrix.h"
//#include "pzsbstrmatrix.h"
//#include "pzfstrmatrix.h"
//
//#include "TPZMaterial.h"
#include "pzbndcond.h"
//#include "pzelasmat.h"
//#include "pzplaca.h"
#include "pzpoisson3d.h"
//#include "pzmathyperelastic.h"
//#include "pzmattest3d.h"
//#include "pzmatplaca2.h"
//
//#include "pzfunction.h"
//
//#include "pzgengrid.h"
//#include "TPZExtendGridDimension.h"
//#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"

//#include "pzshapelinear.h"
//
//#include "TPZRefPatternTools.h"

#include "TPZPostProcessError.h"

#include <time.h>
#include <stdio.h>
#include <math.h>

#include <fstream>
#include <cmath>


//#include "problem.h"
//#include "TPZCreateMesh.h"
//#include "pzbuildmultiphysicsmesh.h"
//#include "pzcondensedcompel.h"

#include "TPZAnalyticSolution.h"

#include "TPZGmshReader.h"


using namespace std;
//using namespace pzshape;
//using namespace pzgeom;


/**  Global variables  */
// Crea malla computacional sem forcingfunction quando hasforcingfunction = 0, ou toma diferentes forcingfuncition para diferentes
// valores de hasforcingfunction
TPZCompMesh *CreateComputationalMesh(TPZGeoMesh *gmesh,int dim,int hasforcingfunction);

TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder);


/** Fucntions to apply refinement. */
void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh);


/** PROBLEMS */
bool SolvePoissonProblem(struct SimulationCase &sim_case);

#ifdef LOG4CXX
static LoggerPtr  logger(Logger::getLogger("pz.refine"));
#endif

// Simulation Case
struct SimulationCase {
    int   nthreads = 0;
    int numinitialrefine = 0;
    int porder = 1;
    std::string  dir_name;
    std::string  file_name;
    
    SimulationCase()
    {
        dir_name = "dump";
    }
    SimulationCase(const SimulationCase &cp) : nthreads(cp.nthreads), numinitialrefine(cp.numinitialrefine), porder(cp.porder),
        dir_name(cp.dir_name), file_name(cp.file_name)
    {
    }
};

// MAIN FUNCTION TO NUMERICAL SOLVE WITH AUTO ADAPTIVE HP REFINEMENTS
/** Laplace equation on square 1D 2D 3D - Volker John article 2000 */
int main(int argc,char *argv[]) {
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Initializing uniform refinements for reference elements
	gRefDBase.InitializeAllUniformRefPatterns();
//    gRefDBase.InitializeRefPatterns();

    
    struct SimulationCase Case_1;
    Case_1.nthreads = 12;
    Case_1.numinitialrefine = 0;
    Case_1.dir_name = "Quad_Case_1";
    Case_1.file_name = "Quad.msh";
        
    struct SimulationCase Case_2;
    Case_2.nthreads = 12;
    Case_2.dir_name = "Cube_Case_1";
    Case_2.file_name = "Cube.msh";
        
    
    if(!SolvePoissonProblem(Case_2)){
        return 1;
    }
        
    

    
    return 0;
}


bool SolvePoissonProblem(struct SimulationCase &sim_case) {
	// Variables

    TPZGeoMesh *gmesh = 0;
    {
        TPZGmshReader gmsh;
        gmsh.fPZMaterialId[1]["dirichlet"] = -1;
        gmsh.fPZMaterialId[1]["neuman"] = -2;
        gmsh.fPZMaterialId[2]["domain"] = 1;
        gmsh.fPZMaterialId[2]["dirichlet"] = -1;
        gmsh.fPZMaterialId[2]["neuman"] = -2;
        gmsh.fPZMaterialId[3]["domain"] = 1;
#ifdef MACOSX
        gmesh = gmsh.GeometricGmshMesh("../"+ sim_case.file_name);
#else
        gmesh = gmsh.GeometricGmshMesh(sim_case.file_name);
#endif
    }

    // Creating the directory
    std::string command = "mkdir " + sim_case.dir_name;
    system(command.c_str());
    
	
	// Output files
    std::string file_name = sim_case.dir_name + "/" + "ErrorsHP_Poisson.txt";
	std::ofstream fileerrors(file_name,ios::app);   // To store all errors calculated by TPZAnalysis (PosProcess)
	// Initial message to print computed errors
	

	// Initializing the auto adaptive process
	TPZVec<REAL> ervec, ErrorVec(100,0.0);
	TPZVec<int64_t> NEquations(100,0L);
	TPZVec<REAL> ervecbyel;
	TPZVec<REAL> gradervecbyel;



	/** Variable names for post processing */
	TPZStack<std::string> scalnames, vecnames;
	scalnames.Push("POrder");
	scalnames.Push("Pressure");
    vecnames.Push("Derivative");

	fileerrors.flush();
    
    {
        UniformRefinement(sim_case.numinitialrefine,gmesh);
    }
    // refinar um elemento

    if(0)
    {
        for (auto gel:gmesh->ElementVec()) {
            if(gel && gel->Dimension() == gmesh->Dimension() && gel->NSubElements()==0)
            {
                TPZStack<TPZGeoEl *> subels;
                gel->Divide(subels);
                break;
            }
        }
    }
	// Creating computational mesh (approximation space and materials)
	TPZCompEl::SetgOrder(sim_case.porder);
	gmesh->SetName("Malha Geometrica original");
    
    
//    TPZManVector<TPZCompMesh *,5> meshvec(5,0);

    
    
    TPZManVector<TPZCompMesh *> meshvec(0);
    TPZCompMesh *pressuremesh = 0;
	// loop solving iteratively
    if(1)
    {
        pressuremesh = CMeshPressure(gmesh, sim_case.porder);
        
        pressuremesh->AdjustBoundaryElements();
    }
    
    {
        std::ofstream out("CompMesh.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(pressuremesh, out);
        std::ofstream out2("CompMesh.txt");
        pressuremesh->Print(out2);

    }
    
    TLaplaceExample1 example;
    example.fExact = TLaplaceExample1::ECosCos;
    example.fDimension = gmesh->Dimension();
    example.fSignConvention = -1;
//    example.fCenter[0] = 0.5;
//    example.fCenter[1] = 0.5;
    {
        for (auto it:pressuremesh->MaterialVec())
        {
            TPZMaterial *mat = it.second;
            TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
            if(!bc)
            {
                mat->SetForcingFunction(example.ForcingFunction());
            }
            else
            {
                bc->SetForcingFunction(0, example.Exact());
            }
        }
    }
    
    
    TPZAnalysis an(pressuremesh,true);
    an.SetExact(example.ExactSolution());
    {
        std::stringstream sout;
        sout << sim_case.dir_name << "/" << "Poisson" << gmesh->Dimension() << "numref" << sim_case.numinitialrefine  << "Porder" << sim_case.porder << ".vtk";
        an.DefineGraphMesh(gmesh->Dimension(),scalnames,vecnames,sout.str());
    }
    
    
    pressuremesh->SetName("Malha computacional adaptada");
    // Printing geometric and computational mesh
#ifdef PZDEBUG
    if(1)
    {
        std::ofstream out("../PressureGeoMesh.txt");
        pressuremesh->Reference()->Print(out);
//        pressuremesh->Print(out);
    }
#endif
    // Solve using symmetric matrix then using Cholesky (direct method)
        
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(pressuremesh);
    strmat.SetNumThreads(sim_case.nthreads);
    an.SetStructuralMatrix(strmat);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(pressuremesh);
    strmat.SetNumThreads(0);
    strmat.SetDecomposeType(ECholesky);
//		TPZSkylineStructMatrix strmat3(cmesh);
//        strmat3.SetNumThreads(8);
#endif
		
    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ECholesky);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
            
    an.Assemble();
    
    

    an.Solve();
    

    // Post processing
    an.PostProcess(1,gmesh->Dimension());
    
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
    
    
//    error.ComputeExactH1SemiNormErrors(*exact, exactelementerror);
    {
        TPZStack<std::string> scalnames,vecnames;
        scalnames.Push("State");
        scalnames.Push("Error");
        scalnames.Push("TrueError");
        scalnames.Push("EffectivityIndex");
        an.DefineGraphMesh(pressuremesh->Dimension(), scalnames, vecnames, "Errors.vtk");
        an.PostProcess(1);
    }
    
				
	if(gmesh)
		delete gmesh;
	gmesh = NULL;
	return true;
}



void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh) {
	TPZManVector<TPZGeoEl*> filhos;
    for(int D=0; D<nDiv; D++)
    {
        int64_t nels = gmesh->NElements();
        for(int64_t elem = 0; elem < nels; elem++)
        {    
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
			if(!gel || gel->HasSubElement())
				continue;
			if(gel->Dimension() == 0) continue;
            gel->Divide(filhos);
        }
    }
}

TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder)
{
    int dim = gmesh->Dimension();
    int matid = 1;
    int dirichlet = 0;
    int neumann = 1;
    /// criar materiais
    //TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    TPZMatPoisson3d *material = new TPZMatPoisson3d(matid,dim);
    material->NStateVariables();
    
    
    TPZCompMesh * cmesh;
    cmesh = new TPZCompMesh(gmesh);
    
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond0 = material->CreateBC(mat, -1,dirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, -2,neumann, val1, val2);

    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->SetDefaultOrder(pOrder);
    //cmesh->SetDimModel(dim);
    
    cmesh->SetAllCreateFunctionsContinuous();
    
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    
    return cmesh;
    
}





