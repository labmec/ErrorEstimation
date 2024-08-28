// Implement the Estimation for Elasticity problem 2D
//
// Created by denise on 06/09/2023.
//

#include "InputTreatment.h"
#include "Solver.h"
#include "Output.h"
#include "tpzgeoelrefpattern.h"
#include "DataStructure.h"
#include "Tools.h"
#include <tuple>
#include <fenv.h>

//%%%%%%

#include "common_files.h"

#include <TPZGeoMeshTools.h>
#include "TPZAnalyticSolution.h"
#include <TPZGmshReader.h>
#include "TPZCompMeshTools.h"
#include "pzlog.h"
#include "TPZLinearAnalysis.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZHDivApproxCreator.h"

#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
#include "TPZTimer.h"
#include "TPZSimpleTimer.h"
#include "TPZVTKGenerator.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "Elasticity/TPZMixedElasticityND.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "pzblockdiag.h"
#include "pzbdstrmatrix.h"
#include "TPZVTKGeoMesh.h"
#include <iostream>
#include <fstream>
#include "TPZElasticityErrorEstimator.h"

std::ofstream printerrors("results_errors.txt",std::ios::app);

enum EMatid  {ENone, EDomain, EBoundary};

/**
   @brief Creates a geometric mesh with elements of a given type on a unit square or cube (depending on the mesh dimension).
   @param[in] nDivs Number of divisions (rows of elements) in x, y and z.
   @param[in] volId Material identifier for the volumetric region.
   @param[in] bcId Material identifier for the boundary.
*/
template<class tshape>
TPZGeoMesh*
CreateGeoMesh(TPZVec<int> &nDivs, int volId, int bcId);

// The test function
template<class tshape>
void SolveFEMProblem(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily,ProblemConfig &config);

void InsertMaterials(const int &dim, TPZHDivApproxCreator& hdivc, TPZAnalyticSolution *fAn);
void EstimateError(ProblemConfig &config, TPZMultiphysicsCompMesh *multimesh);

//%%%%


int main() {
    
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    
//    //----
    ProblemConfig pConfig;
//    pConfig.k = 1;
//  //  pConfig.n = 2;
//  //  pConfig.problem = "ESinSin";                 //// {"ESinSin","EArcTan",ESteklovNonConst", "EBubble2D", "ELaplace","ESing2D, "EProb"}
//    pConfig.integrationorder = 6;
//  //  pConfig.maxIter = 100;                     //// Maximum iterations for computing the exact solution (only for ELaplace)
//    pConfig.approx = "Hybrid";                 //// {"H1","Hybrid", "Mixed"}//
// //   pConfig.topology = "Quadrilateral";        //// Triangular, Quadrilateral, LQuad, Tetrahedral, Hexahedral, Prism
//    pConfig.refLevel = 1;                      //// How many uniform refinements
//    pConfig.numberAdapativitySteps = 1;        //// Maximum number of adapativity refinement steps.
//    pConfig.estimateError = true;              //// Wheater Error Estimation procedure is invoked
//    pConfig.debugger = false;                   //// Print geometric and computational mesh for the simulation (Error estimate not involved).
//    pConfig.vtkResolution = 1;                 //// Vtk resolution. Set 0 to see a paraview mesh equals the  simulation mesh.
//
    //----

    const int xdiv = 3; //Number of elements in each direction
    const int pOrder = 1; // Polynomial degree
    pConfig.porder = pOrder;
    pConfig.ndivisions = xdiv;
    // Family of HDiv approximation spaces.
    // The possible choices are HDivFamily::EHDivStandard, HDivFamily::EHDivConstant and HDivFamily::EHDivKernel
    HDivFamily hdivfam = HDivFamily::EHDivStandard;
    
    //Creates the geometric mesh for the given topology and solve the FEM problem.
    SolveFEMProblem<pzshape::TPZShapeQuad>(xdiv,pOrder,hdivfam, pConfig);
   
    return 0;
}


template<class tshape>
void SolveFEMProblem(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily,ProblemConfig &config)
{

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    const int DIM = tshape::Dimension;
    TPZVec<int> nDivs;
    TPZVec<int> divs = {6};
    // TPZVec<int> divs = {5};
    int pend = 1;
    
    for (int iorder = 1; iorder < pend+1; iorder++){
        for (int idiv = 0; idiv < divs.size(); idiv++){
        int divx = divs[idiv];


    if (DIM == 2) nDivs = {divx,divx};
    if (DIM == 3) nDivs = {divx,divx,divx};
    
    // Creates/import a geometric mesh
    auto gmesh = CreateGeoMesh<tshape>(nDivs, 1, -1);

    // Creates an hdivApproxCreator object. It is an environment developped to
    // help creating H(div)-family possible approximation spaces.
    TPZHDivApproxCreator hdivCreator(gmesh);
    //Set the family of H(div) functions: Standard, Constant or Kernel
    hdivCreator.HdivFamily() = hdivfamily;
    //Set the problem type to be solved: Only EDarcy and EElastic are currently available
    hdivCreator.ProbType() = ProblemType::EElastic;
    //Includes the rigid body spaces (constant flux and pressure) if set as true
    hdivCreator.IsRigidBodySpaces() = false;
    //Set the default polynomial order
    hdivCreator.SetDefaultOrder(iorder);
    //Set the extra polynomial order for the bubble functions. If zero, the polynomial degree
    //of the internal functions are the same as the default order
    hdivCreator.SetExtraInternalOrder(0);
    //Sets if the resulting problem should or not be condensed
    hdivCreator.SetShouldCondense(false);
    // hdivCreator.SetShouldCondense(false);
    
    //Sets the type of hybridizantion desired.
    //The current options are HybridizationType::ENone, HybridizationType::EStandard
    //and HybridizationType::ESemi (the last only works with H(div)-constant spaces)
    hdivCreator.HybridType() = HybridizationType::ENone;
    // hdivCreator.HybridType() = HybridizationType::EStandard;

    // Prints gmesh mesh properties
    std::string vtk_name = "geoMesh.vtk";
    std::ofstream vtkfile(vtk_name.c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
            
    config.gmesh=gmesh;

    //Creates an analytical solution to test
    TPZAnalyticSolution *gAnalytic = 0;
    if (hdivCreator.ProbType() == ProblemType::EDarcy){
        TLaplaceExample1 *lap = new TLaplaceExample1;
        lap->fExact = TLaplaceExample1::EHarmonic;
        gAnalytic = lap;
    } else if (hdivCreator.ProbType() == ProblemType::EElastic){
        if (DIM == 2){
            TElasticity2DAnalytic *elas = new TElasticity2DAnalytic;
            // double lambda = 123.;
            // double mu = 79.3;
            elas->gE = 100.;//mu*(3*lambda+2*mu)/(lambda+mu);
            elas->gPoisson = 0.3;//0.5*lambda/(lambda+mu);
            elas->fProblemType = TElasticity2DAnalytic::EDispx;
            // elas->fProblemType = TElasticity2DAnalytic::EThiago;
            // elas->fPlaneStress = 0;
            gAnalytic = elas;
        } else if (DIM == 3){
            TElasticity3DAnalytic *elas = new TElasticity3DAnalytic;
            elas->fE = 1.;
            elas->fPoisson = 0.0;
            elas->fProblemType = TElasticity3DAnalytic::ELoadedBeam;
            gAnalytic = elas;
        }
    } else {
        DebugStop();
    }
         //   config.exact2 = gAnalytic;
            
    //Insert Materials
    InsertMaterials(DIM,hdivCreator,gAnalytic);

    //Gets the Multiphysics mesh from the HdivApproxCreator
    TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();
    //Here the cmesh is printed
    // std::string txt = "cmesh.txt";
    // std::ofstream myfile(txt);
    // cmesh->Print(myfile);

    //Create the analysis environment
    TPZLinearAnalysis an(cmesh,RenumType::ESloan);
    // TPZLinearAnalysis an(cmesh,RenumType::EMetis);
    an.SetExact(gAnalytic->ExactSolution(),4);
    // if (hdivCreator.ProbType() == ProblemType::EDarcy){
        
    // } else if (hdivCreator.ProbType() == ProblemType::EElastic){
    //     an.SetExact(gAnalytic,solOrder);
    // } else {
    //     DebugStop();
    // }

    // Solve problem
    constexpr int nThreads{20};
            
   // TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> stiffness(cmesh);
            
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<> stiffness(cmesh);
            //stiffness.SetNumThreads(0);
#else
    TPZSkylineStructMatrix<STATE> stiffness(cmesh);
           // stiffness.SetNumThreads(0);
#endif
            
    stiffness.SetNumThreads(nThreads);
    an.SetStructuralMatrix(stiffness);
            
            

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    an.SetSolver(step);
    
    //Assemble and solve the problem.
    an.Run();
  
    //Printing results in a vtk file
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(cmesh->MeshVector(), an.Mesh());
    std::string plotfile = "res_h"+std::to_string(idiv)+"p"+std::to_string(iorder);//sem o .vtk no final
    constexpr int vtkRes{0};//Resolution
    {
        //Fields to be printed
        TPZVec<std::string> fields;
        if (hdivCreator.ProbType() == ProblemType::EDarcy){
            fields = {"Pressure","ExactPressure","Flux","ExactFlux"};
        } else if (hdivCreator.ProbType() == ProblemType::EElastic){
            fields = {"Displacement","SigmaX","SigmaY","TauXY","ExactDisplacement","ExactStress"};
        } else {
            DebugStop();
        }
        
        auto vtk = TPZVTKGenerator( an.Mesh(), fields, plotfile, vtkRes);
        vtk.Do();
    }

    // //Compute error
    std::ofstream anPostProcessFile("postprocess.txt");
    TPZManVector<REAL,7> error;
    an.LoadSolution();
    // cmesh->LoadSolution(cmesh->Solution());
    // cmesh->ExpandSolution();
    int64_t nelem = 0;
    for (int64_t i = 0; i < an.Mesh()->NElements(); i++)
    {
        auto el = an.Mesh()->ElementVec()[i];
        if (el->Dimension() == DIM) nelem++;
    }
        
    an.Mesh()->ElementSolution().Redim(nelem, 7);
    an.SetExact(gAnalytic->ExactSolution(),5);
    an.PostProcessError(error,true,anPostProcessFile);

    std::ofstream postVTK(plotfile+".0.vtk",std::ios::app);
    TPZFMatrix<STATE> SolMatrix = an.Mesh()->ElementSolution();
    postVTK << "CELL_DATA " << nelem << std::endl;
    postVTK << "SCALARS ErrorStress float \nLOOKUP_TABLE default" << std::endl;
    for (int64_t iel = 0; iel < nelem; iel++){
        postVTK << SolMatrix.GetVal(iel,0) << " ";
    }
    postVTK << std::endl;
    postVTK << "SCALARS ErrorEnergy float \nLOOKUP_TABLE default" << std::endl;
    for (int64_t iel = 0; iel < nelem; iel++){
        postVTK << SolMatrix.GetVal(iel,1) << " ";
    }
    postVTK << std::endl;
    postVTK << "SCALARS ErrorDivStress float \nLOOKUP_TABLE default" << std::endl;
    for (int64_t iel = 0; iel < nelem; iel++){
        postVTK << SolMatrix.GetVal(iel,2) << " ";
    }
    postVTK << std::endl;
    postVTK << "SCALARS ErrorDisplacement float \nLOOKUP_TABLE default" << std::endl;
    for (int64_t iel = 0; iel < nelem; iel++){
        postVTK << SolMatrix.GetVal(iel,3) << " ";
    }
    postVTK << std::endl;
    postVTK << "SCALARS ErrorRotation float \nLOOKUP_TABLE default" << std::endl;
    for (int64_t iel = 0; iel < nelem; iel++){
        postVTK << SolMatrix.GetVal(iel,4) << " ";
    }
    postVTK << std::endl;
    postVTK << "SCALARS ErrorSymmetry float \nLOOKUP_TABLE default" << std::endl;
    for (int64_t iel = 0; iel < nelem; iel++){
        postVTK << SolMatrix.GetVal(iel,5) << " ";
    }
    postVTK << std::endl;
    
    
    // //Print Errors
    // std::cout << "POrder = " << iorder << std::endl;
    // std::cout << "h = " << std::fixed <<  1./double(divx) << std::endl;
    // std::cout << "L2 Stress = " << std::scientific << std::setprecision(15) << error[0] << std::endl;
    // std::cout << "En Stress = " << error[1] << std::endl;
    // std::cout << "L2 DivStr = " << error[2] << std::endl;
    // std::cout << "L2 Displa = " << error[3] << std::endl;
    // std::cout << "L2 Rotati = " << error[4] << std::endl;
    // std::cout << "L2 Symmet = " << error[5] << std::endl;
    // std::cout << "En Displa = " << error[6] << std::endl;
    // printerrors << iorder << " " << std::fixed << std::setprecision(5) <<  1./double(divx) << " "
    //             << std::scientific << std::setprecision(15) << error[0] << " "
    //             << error[1] << " " << error[2] << " " << error[3] << " " << error[4]
    //             << " " << error[5] << " " << error[6] << std::endl;
            
            // Estimation Error for elasticity problem
            EstimateError(config,cmesh);
        }
    }
}
//Estimation error

void EstimateError(ProblemConfig &config, TPZMultiphysicsCompMesh *multimesh) {

    std::cout << "\nError Estimation processing for Elasticity problem " << std::endl;

    // Error estimation

    //criar construtor para estimador de elasticidade
    
    auto *mhm = new TPZMHMixedMeshControl(config.gmesh);
    
    TPZVec<int64_t> coarseindices(config.gmesh->NElements());
    int64_t count = 0;
    int coarselevel = 0;
    int64_t nel = config.gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = config.gmesh->Element(el);
        if(gel && gel->Dimension() == config.gmesh->Dimension() && gel->Level()==coarselevel)
        {
            coarseindices[count++] = el;
        }
    }
    coarseindices.resize(count);
    mhm->DefinePartitionbyCoarseIndices(coarseindices);
    // mhm->fMaterialIds = config.materialids;
    // mhm->fMaterialBCIds = config.bcmaterialids;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.exactElast = new TElasticity2DAnalytic;
    config.exactElast.operator*().fProblemType = TElasticity2DAnalytic::EDispx;
    {
        TPZCompMesh &cmesh = mhm->CMesh();
        int dim = mhm->GMesh()->Dimension();
        cmesh.SetDimModel(dim);
        double E= 1000.;
        double nu = 0.3;
        auto matelast = new TPZMixedElasticityND(1,E,nu,0,0,1,dim);
        matelast->SetExactSol(config.exactElast->ExactSolution(),5);
        matelast->SetForcingFunction(config.exactElast->ForceFunc(),5);
        cmesh.InsertMaterialObject(matelast);
        for (auto matid : config.bcmaterialids) {
            TPZFNMatrix<1, REAL> val1(dim, dim, 0.);
            TPZManVector<REAL> val2(dim, 0.);
            int bctype = 0;
            TPZBndCondT<STATE> *bc = matelast->CreateBC(matelast, matid, bctype, val1, val2);
            bc->SetForcingFunctionBC(config.exactElast->ExactSolution(),4);
            cmesh.InsertMaterialObject(bc);
        }
    }
    mhm->fMaterialIds = {1};
    mhm->fMaterialBCIds = {-1};
    bool substructure = true;
    mhm->BuildComputationalMesh(substructure);
    
    bool postProcWithHDiv = false;
    TPZElasticityErrorEstimator ErrorEstimator(*multimesh, mhm, postProcWithHDiv);
    
    ErrorEstimator.SetAnalyticSolution(config.exactElast);
    //create displacement reconstruction
    ErrorEstimator.PrimalReconstruction();
    // ErrorEstimator.DisplacementReconstruction();

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    TPZManVector<REAL, 6> errors;
    TPZManVector<REAL, 6> elementerrors;
    std::stringstream outVTK;
    outVTK << config.dir_name << "/" << config.problemname << "-" << config.ndivisions << "-" << config.ninternalref
           << "-Errors.vtk";
    std::string outVTKstring = outVTK.str();
    ErrorEstimator.ComputeErrors(errors, elementerrors, outVTKstring);

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "-GlobalErrors.txt";
        std::ofstream file(fileName, std::ios::app);
        Tools::PrintElasticityErrors(file, config, errors);
    }
}

//Create
template <class tshape>
TPZGeoMesh*
CreateGeoMesh(TPZVec<int> &nDivs, int volId, int bcId)
{
    
    MMeshType meshType;
    int dim = tshape::Dimension;

    switch (tshape::Type())
    {
    case ETriangle:
        meshType = MMeshType::ETriangular;
        break;
    case EQuadrilateral:
        meshType = MMeshType::EQuadrilateral;
        break;
    case ETetraedro:
        meshType = MMeshType::ETetrahedral;
        break;
    case ECube:
        meshType = MMeshType::EHexahedral;
        break;
        case EPrisma:
        meshType = MMeshType::EPrismatic;
        break;
    default:
        DebugStop();
    }

    TPZManVector<REAL,3> minX = {0,0,0};
    TPZManVector<REAL,3> maxX = {1,1,1};
    int nMats = 2*dim+1;

    //all bcs share the same id
    constexpr bool createBoundEls{true};
    TPZVec<int> matIds(nMats,bcId);
    matIds[0] = volId;
    // matIds[1] = bcId;
    // matIds[2] = -11;
    // matIds[3] = -11;
    // matIds[4] = -11;
    
    TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
                        matIds, nDivs, meshType,createBoundEls);
    // TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshSingleEl(meshType,
    //                     volId,createBoundEls, bcId);
    
    return gmesh;
    
}

void InsertMaterials(const int &dim, TPZHDivApproxCreator& hdivCreator, TPZAnalyticSolution *fAn) {
    if (hdivCreator.ProbType() == ProblemType::EDarcy) {
        TPZMixedDarcyFlow* matdarcy = new TPZMixedDarcyFlow(1, dim);
        // matdarcy->SetConstantPermeability(1.);
        TLaplaceExample1* lapl = dynamic_cast<TLaplaceExample1*> (fAn);
        matdarcy->SetExactSol(lapl->ExactSolution(), 4);
        matdarcy->SetForcingFunction(lapl->ForceFunc(), 4);

        hdivCreator.InsertMaterialObject(matdarcy);

        TPZFMatrix<STATE> val1(1, 1, 0.);
        TPZManVector<STATE> val2(1, 0.);
        TPZBndCondT<STATE> *BCond1 = matdarcy->CreateBC(matdarcy, -1, 0, val1, val2);
        BCond1->SetForcingFunctionBC(lapl->ExactSolution(), 4);
        hdivCreator.InsertMaterialObject(BCond1);
    } else if (hdivCreator.ProbType() == ProblemType::EElastic) {
        if (dim == 2) {
            TElasticity2DAnalytic *elas2D = dynamic_cast<TElasticity2DAnalytic*> (fAn);
            auto matelas = new TPZMixedElasticityND(1, elas2D->gE, elas2D->gPoisson, 0, 0, elas2D->fPlaneStress, dim);
            matelas->SetExactSol(elas2D->ExactSolution(), 4);
            matelas->SetForcingFunction(elas2D->ForceFunc(), 4);
            hdivCreator.InsertMaterialObject(matelas);

            TPZFMatrix<STATE> val1(dim, dim, 0.);
            TPZManVector<STATE> val2(dim, 0.);
            TPZBndCondT<STATE> *BCond1 = matelas->CreateBC(matelas, -1, 0, val1, val2);
            BCond1->SetForcingFunctionBC(elas2D->ExactSolution(), 4);
            hdivCreator.InsertMaterialObject(BCond1);
        }
        if (dim == 3) {
            TElasticity3DAnalytic *elas3D = dynamic_cast<TElasticity3DAnalytic*> (fAn);
            auto matelas = new TPZMixedElasticityND(1, elas3D->fE, elas3D->fPoisson, 0, 0, 0, dim);
            matelas->SetExactSol(elas3D->ExactSolution(), 4);
            matelas->SetForcingFunction(elas3D->ForceFunc(), 4);
            hdivCreator.InsertMaterialObject(matelas);

            TPZFMatrix<STATE> val1(dim, dim, 0.);
            TPZManVector<STATE> val2(dim, 0.);
            TPZBndCondT<STATE> *BCond1 = matelas->CreateBC(matelas, -1, 0, val1, val2);
            BCond1->SetForcingFunctionBC(elas3D->ExactSolution(), 4);
            hdivCreator.InsertMaterialObject(BCond1);
        }
    } else {
        DebugStop(); //Material Not Implemented
    }
}
