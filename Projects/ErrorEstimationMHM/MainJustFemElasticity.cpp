//
//  MainJustFemElasticity.cpp
//  MHM_ElastPhil
//
//  Created by Denise De Siqueira on 12/12/23.
//

#include <stdio.h>

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
#include "TPZH1ApproxCreator.h"
#include "Elasticity/TPZElasticity2D.h"


//

#include "pzgmesh.h"
#include <Pre/TPZGenGrid3D.h>
#include <Pre/TPZMHMixedMeshControl.h>
#include <TPZMFSolutionTransfer.h>
//#include <Tools.h>
#include <ToolsMHM.h>
#include <Util/pzlog.h>
#include "TPZLinearAnalysis.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "Elasticity/TPZMixedElasticityND.h"
#include "TPZBndCondT.h"
#include "Pre/TPZMHMeshControl.h"
#include "TPZElasticityMHMHDivErrorEstimator.h"
#include "TPZElasticityErrorEstimator.h"
//
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);

std::ofstream printerrors("ErrorsH1.txt",std::ios::app);

enum EMatid  {ENone, EDomain, EBoundary};

/**
   @brief Creates a geometric mesh with elements of a given type on a unit square or cube (depending on the mesh dimension).
   @param[in] nDivs Number of divisions (rows of elements) in x, y and z.
   @param[in] volId Material identifier for the volumetric region.
   @param[in] bcId Material identifier for the boundary.
*/
template<class tshape>
TPZGeoMesh*
CreateGeoMesh(TPZVec<int> &nDivs, EMatid volId, EMatid bcId, REAL distortion);

// The test function
template<class tshape>
void SolveFEMProblem(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily,ProblemConfig &config);
template<class tshape>
void SolveH1Problem(const int &xdiv, const int &pOrder, H1Family &h1family,ProblemConfig &config);
template<class tshape>
void SolvingH1Displacement(TPZCompMesh *cH1Mesh,ProblemConfig &config);


void InsertMaterials(int &dim, TPZHDivApproxCreator& hdivc, TPZAnalyticSolution *fAn);
void InsertH1Materials(int &dim, TPZH1ApproxCreator& h1Creator,TPZAnalyticSolution *fAn);
void EstimateError(ProblemConfig &config, TPZMultiphysicsCompMesh *multimesh);
void EstimateErrorElasticity(ProblemConfig &config, TPZMultiphysicsCompMesh *originalMesh, int step);
TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity2DAnalytic *elas,ProblemConfig &config);
template<class tshape>
void RunSmoothProblemSquareMesh(ProblemConfig &pConfig);

template<class tshape>
void RunSmoothProblemTrapMesh(ProblemConfig &pConfig);

template<class tshape>
void RunLShapeProblem(ProblemConfig &pConfig);

template<class tshape>
void RunLambdaTest(ProblemConfig &pConfig);

template<class tshape>
void SolveFEMProblemNew(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily,ProblemConfig &config);


//%%%%
int main() {
    

    ProblemConfig pConfig;
    pConfig.vtkResolution = 0;
   
    pConfig.exactElast = new TElasticity2DAnalytic;
    RunSmoothProblemSquareMesh<pzshape::TPZShapeQuad>(pConfig);
    // RunSmoothProblemTrapMesh<pzshape::TPZShapeQuad>(pConfig);
     //RunLShapeProblem<pzshape::TPZShapeQuad>(pConfig);
   // RunLambdaTest<pzshape::TPZShapeQuad>(pConfig);
   
    return 0;
}


int main2() {
    

    ProblemConfig pConfig;
    pConfig.geometry = ProblemConfig::EGeometry::EQuad;
    pConfig.exactElast = new TElasticity2DAnalytic;
    //pConfig.exactElast.operator*().fProblemType = TElasticity2DAnalytic::EDispy;
    switch (pConfig.geometry){
        case ProblemConfig::EGeometry::ECrack:
            pConfig.problemname = "ECrack";
            pConfig.exactElast->fProblemType = TElasticity2DAnalytic::ECrack;
            pConfig.coefgE = 100.;
            pConfig.coefgPoisson = 0.3;
            break;
        case ProblemConfig::EGeometry::ELShape:
            pConfig.problemname = "ELShape";
            pConfig.exactElast->fProblemType = TElasticity2DAnalytic::ELShape;
            pConfig.mu = 1.;
            pConfig.lambda = 5.0;
            pConfig.dir_name = "LShapeAdapt";
            //pConfig.dir_name = "LShape-Uniform";
            break;
        case ProblemConfig::EGeometry::EQuad:
        case ProblemConfig::EGeometry::ETrap:
        default:
            pConfig.exactElast->fProblemType = TElasticity2DAnalytic::EHarmonic;
            pConfig.lambda= 123.;//1000.;//123.;//100;//10000.;//10000.;
            pConfig.mu= 79.3;//1.;//
            //pConfig.problemname="EHarmonic-ETrap";
            pConfig.dir_name = "ArticleEx01";
            //pConfig.dir_name = "CompareErrors";
            //pConfig.dir_name = "QuadMeshRef-Trap";
            break;
    }
    
    const int xdiv = 10; //Number of elements in each direction
    const int pOrder = 4; // Polynomial degree
    //pConfig.porder = pOrder;
    pConfig.ndivisions = xdiv;
    pConfig.hdivmais = 1;// internal order
    pConfig.isAdaptivity = false;
    pConfig.adaptivityStep = 2;//numero de steps no refinamento
   
   
    // Family of HDiv approximation spaces.
    // The possible choices are HDivFamily::EHDivStandard, HDivFamily::EHDivConstant and HDivFamily::EHDivKernel
    HDivFamily hdivfam = HDivFamily::EHDivStandard;
    
    //Creates the geometric mesh for the given topology and solve the FEM problem.
   // TPZVec<int64_t> values = {1,10,100,1000,10000};
    
   // for (int64_t ilambda=0; ilambda< values.size(); ilambda++) {
        
     //   pConfig.lambda = values[ilambda];
//    for (int64_t iorder=1; iorder< pOrder;iorder++) {
//        pConfig.porder = iorder;
//    
//        
//        SolveFEMProblem<pzshape::TPZShapeQuad>(xdiv,pConfig.porder,hdivfam, pConfig);
//    }

  //  }
    
 //RunSmoothProblemSquareMesh<pzshape::TPZShapeQuad>(pConfig);
   
    return 0;
}
template<class tshape>
void RunSmoothProblemSquareMesh(ProblemConfig &pConfig){
    
    pConfig.geometry = ProblemConfig::EGeometry::EQuad;
    pConfig.exactElast->fProblemType = TElasticity2DAnalytic::EHomoDir;//EThiago;//EDispx;//
    pConfig.lambda= 123.;
    pConfig.mu= 79.3;
    pConfig.problemname="EHomoDir-Problem";
    pConfig.dir_name = "SmoothProb-Quad";
   // pConfig.dir_name = "SymmetricTest";
    
    const int xdiv = 2; //Number of elements in each direction

    const int pOrder = 2;

    pConfig.ndivisions = xdiv;
    pConfig.hdivmais = 1;// internal order
    pConfig.isAdaptivity = false;
    pConfig.adaptivityStep = 1;//numero de steps no refinamento
    HDivFamily hdivfam = HDivFamily::EHDivStandard;
    TPZGeoMesh *gmesh;
    REAL distortion = 0;
    int DIM = tshape::Dimension;
    TPZVec<int> nDivs = {1,1};
   
    
    TPZVec<int> divs = {4,8,16,32};//,64};
    
    for (int64_t iorder=pOrder; iorder< pOrder+1;iorder++) {
        pConfig.porder = iorder;
        
        for (int idiv = 0; idiv < divs.size(); idiv++){
            pConfig.ndivisions =divs[idiv];
            int divx = divs[idiv];
            nDivs = {divx,divx};
            
            gmesh = CreateGeoMesh<tshape>(nDivs, EDomain, EBoundary, distortion);
            
            pConfig.gmesh=gmesh;
            
            SolveFEMProblemNew<pzshape::TPZShapeQuad>(xdiv,pConfig.porder,hdivfam, pConfig);
          //  H1Family h1family=H1Family::EH1Standard;
          // SolveH1Problem<pzshape::TPZShapeQuad>(xdiv, pConfig.porder, h1family,pConfig);
        }
        
        
    }
    
}
template<class tshape>
void RunLambdaTest(ProblemConfig &pConfig){
    
    pConfig.geometry = ProblemConfig::EGeometry::EQuad;
    pConfig.exactElast->fProblemType = TElasticity2DAnalytic::EHarmonic;
    pConfig.mu= 1;
    pConfig.problemname="EHarmonic-EQuad";
    pConfig.dir_name = "SmoothLambdaTest";
    
    const int xdiv = 10; //Number of elements in each direction
    const int pOrder = 4;

    pConfig.ndivisions = xdiv;
    pConfig.hdivmais = 1;// internal order
    pConfig.isAdaptivity = false;
    pConfig.adaptivityStep = 2;//numero de steps no refinamento
    HDivFamily hdivfam = HDivFamily::EHDivStandard;
    TPZGeoMesh *gmesh;
    REAL distortion = 0;
    int DIM = tshape::Dimension;
    TPZVec<int> nDivs = {4,4};
   
    
    TPZVec<int> divs = {4,8,16,32,64,128,256,512};
    
     TPZVec<int64_t> values = {1,10,100,1000,10000};
     
for (int64_t ilambda=0; ilambda< values.size(); ilambda++) {
         
         pConfig.lambda = values[ilambda];
    
    for (int64_t iorder=1; iorder< pOrder;iorder++) {
        pConfig.porder = iorder;
        
        for (int idiv = 0; idiv < divs.size(); idiv++){
            pConfig.ndivisions =divs[idiv];
            int divx = divs[idiv];
            nDivs = {divx,divx};
            
            gmesh = CreateGeoMesh<tshape>(nDivs, EDomain, EBoundary, distortion);
            
            pConfig.gmesh=gmesh;
            
            
            SolveFEMProblemNew<pzshape::TPZShapeQuad>(xdiv,pConfig.porder,hdivfam, pConfig);
        }
    }
        
    }
    
}


template<class tshape>
void RunSmoothProblemTrapMesh(ProblemConfig &pConfig){
    
    pConfig.geometry = ProblemConfig::EGeometry::ETrap;
    pConfig.exactElast->fProblemType = TElasticity2DAnalytic::EHarmonic;
    pConfig.lambda= 123.;
    pConfig.mu= 79.3;
    pConfig.problemname="EHarmonic-ETrap";
    pConfig.dir_name = "SmoothProb-Trap";
    
    const int xdiv = 10; //Number of elements in each direction
    const int pOrder = 2;

    pConfig.ndivisions = xdiv;
    pConfig.hdivmais = 1;// internal order
    pConfig.isAdaptivity = false;
    pConfig.adaptivityStep = 2;//numero de steps no refinamento
    HDivFamily hdivfam = HDivFamily::EHDivStandard;
    TPZGeoMesh *gmesh;
    REAL distortion = 1./3;
    int DIM = tshape::Dimension;
    TPZVec<int> nDivs = {4,4};
   
    
    TPZVec<int> divs = {16};//4,8,16,32};//4,8,16,32,64,128,256,512};
    
    for (int64_t iorder=1; iorder< pOrder;iorder++) {
        pConfig.porder = iorder;
        
        for (int idiv = 0; idiv < divs.size(); idiv++){
            pConfig.ndivisions =divs[idiv];
            int divx = divs[idiv];
            nDivs = {divx,divx};
            
            gmesh = CreateGeoMesh<tshape>(nDivs, EDomain, EBoundary, distortion);
            
            pConfig.gmesh=gmesh;
            
            
            SolveFEMProblemNew<pzshape::TPZShapeQuad>(xdiv,pConfig.porder,hdivfam, pConfig);
        }
        
    }
    
}


template<class tshape>
void RunLShapeProblem(ProblemConfig &pConfig){
    pConfig.geometry = ProblemConfig::EGeometry::ELShape;
    pConfig.problemname = "ELShape";
    pConfig.exactElast->fProblemType = TElasticity2DAnalytic::ELShape;
    pConfig.mu = 1.;
    pConfig.lambda = 5.0;
    
    
    const int xdiv = 10; //Number of elements in each direction
    const int pOrder = 1;

    pConfig.ndivisions = xdiv;
    pConfig.hdivmais = 1;// internal order
    pConfig.isAdaptivity = false;
    pConfig.adaptivityStep = 1;//numero de steps no refinamento
    HDivFamily hdivfam = HDivFamily::EHDivStandard;
    TPZGeoMesh *gmesh;

    int DIM = tshape::Dimension;
    TPZVec<int> divs;
    TPZVec<int> nDivs = {4,4};
    if (pConfig.isAdaptivity) {
         divs = {1};
        pConfig.dir_name = "LShapeProblem-Adapt";
    }
    else{
        divs = {1,2,3,4,5,6};
        pConfig.dir_name = "LShapeProblem-Uni";
    }

    for (int64_t iorder = 1; iorder <= pOrder; iorder++) {
        pConfig.porder = iorder;

        for (int idiv = 0; idiv < divs.size(); idiv++) {
            pConfig.ndivisions = divs[idiv];
            int divx = divs[idiv];

            nDivs = {divx, divx};

            TPZVec<int> bcids(8, EBoundary);
            gmesh = Tools::CreateQuadLShapeMesh(bcids);

            int uniref = divx;
            Tools::UniformRefinement(uniref, gmesh);

            pConfig.gmesh = gmesh;

            SolveFEMProblemNew<pzshape::TPZShapeQuad>(xdiv, pConfig.porder, hdivfam, pConfig);
        }
    }
}


template<class tshape>
void SolveFEMProblem(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily,ProblemConfig &config)
{

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    int DIM = tshape::Dimension;
    TPZVec<int> nDivs = {4,4};
    TPZVec<int> divs = {4};//,8,16,32,64,128,256,512};
    
    int pend = 2;
    //config.ndivisions=divs[0];
    
    TPZGeoMesh *gmesh;

   
    
    for (int iorder = config.porder; iorder < config.porder+1; iorder++){
        
        printerrors <<  " porder " << " h " <<   " error stress "<< " error diplacement"<<std::endl;
        
        for (int idiv = 0; idiv < divs.size(); idiv++){
            config.ndivisions =divs[idiv];
                        config.porder= iorder;
                        int divx = divs[idiv];


                            if (DIM == 2) nDivs = {divx,divx};
                            if (DIM == 3) nDivs = {divx,divx,divx};

                        switch(config.geometry) {
                            case ProblemConfig::EGeometry::ECrack:
                            {
                                gmesh = ReadMeshFromGmsh("../../../Crack.msh");
                                break;
                            }
                            case ProblemConfig::EGeometry::ELShape:
                            {
                                TPZVec<int> bcids(8, EBoundary);
                                gmesh = Tools::CreateQuadLShapeMesh(bcids);

                                int uniref = divx;
                                Tools::UniformRefinement(uniref, gmesh);
                                break;
                            }
                            case ProblemConfig::EGeometry::ETrap:
                            {
                                REAL distortion = 1. / 3;
                                gmesh = CreateGeoMesh<tshape>(nDivs, EDomain, EBoundary, distortion);
                                break;
                            }
                            case ProblemConfig::EGeometry::EQuad:
                            {
                                REAL distortion = 0;
                                gmesh = CreateGeoMesh<tshape>(nDivs, EDomain, EBoundary, distortion);
                                break;
                            }
                        }

                        config.gmesh = gmesh;
                    
        for(int refsteps = 1; refsteps< config.adaptivityStep; refsteps ++){
            #ifdef ERRORESTIMATION_DEBUG
            {
                // Prints gmesh mesh properties
                std::string vtk_name = "geoMeshToSolveProbem.vtk";
                std::ofstream vtkfile(vtk_name.c_str());
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
            }
            #endif
            

                
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
                hdivCreator.SetExtraInternalOrder(config.hdivmais);
                //Sets if the resulting problem should or not be condensed
                hdivCreator.SetShouldCondense(true);
                // hdivCreator.SetShouldCondense(false);
                
                //Sets the type of hybridizantion desired.
                //The current options are HybridizationType::ENone, HybridizationType::EStandard
                //and HybridizationType::ESemi (the last only works with H(div)-constant spaces)
                hdivCreator.HybridType() = HybridizationType::ENone;
                // hdivCreator.HybridType() = HybridizationType::EStandard;
                
//                // Prints gmesh mesh properties
//                std::string vtk_name2 = "geoMeshLshape.vtk";
//                std::ofstream vtkfile2(vtk_name2.c_str());
//                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile2, true);
                
              //  config.gmesh=gmesh;
                
                //Creates an analytical solution to test
                TPZAnalyticSolution *gAnalytic = 0;
                if (hdivCreator.ProbType() == ProblemType::EDarcy){
                    TLaplaceExample1 *lap = new TLaplaceExample1;
                    lap->fExact = TLaplaceExample1::EHarmonic;
                    gAnalytic = lap;
                } else if (hdivCreator.ProbType() == ProblemType::EElastic){
                    if (DIM == 2){
                        TElasticity2DAnalytic *elas = new TElasticity2DAnalytic;
                        elas->gE = config.mu*(3*config.lambda+2*config.mu)/(config.lambda+config.mu);
                        elas->gPoisson = 0.5*config.lambda/(config.lambda+config.mu);
                        elas->fProblemType = config.exactElast->fProblemType;
                        elas->fPlaneStress=0;
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
         
                
                //Insert Materials
                InsertMaterials(DIM,hdivCreator,gAnalytic);
                
                //Gets the Multiphysics mesh from the HdivApproxCreator
                TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();

                
                //Create the analysis environment
                TPZLinearAnalysis an(cmesh,RenumType::ESloan);
                // TPZLinearAnalysis an(cmesh,RenumType::EMetis);
                an.SetExact(gAnalytic->ExactSolution(),4);
                
#ifdef PZ_USING_MKL
                TPZSSpStructMatrix<> strmat(cmesh);
                strmat.SetNumThreads(1);
#else
                TPZSkylineStructMatrix<STATE> strmat(cmesh);
                strmat.SetNumThreads(4);
#endif
                
                an.SetStructuralMatrix(strmat);
                TPZStepSolver<STATE> step;
                step.SetDirect(ELDLt);
                an.SetSolver(step);
                std::cout << "Assembling\n";
                an.Assemble();
                
                std::cout << "Solving\n";
                an.Solve();
                std::cout << "Finished\n";
                an.LoadSolution(); // compute internal dofs
                
            /*Compute error
                   std::ofstream anPostProcessFile("PostprocessFem.txt");
                   TPZManVector<REAL,7> error;
                   an.LoadSolution();
                    cmesh->LoadSolution(cmesh->Solution());
                    cmesh->ExpandSolution();
                   int64_t nelem = an.Mesh()->NElements();
                       
               //TODO: The element solution now has 7 columns, the error of the antisymmetric part has increased
                   
                   an.Mesh()->ElementSolution().Redim(nelem, 7);
                   an.SetExact(gAnalytic->ExactSolution(),5);
                   an.PostProcessError(error,true,anPostProcessFile);
                   
                   double haux=pow(2, -idiv);
                   printerrors << iorder << ", " << std::fixed << std::setprecision(5) <<  1./haux << ", "
                   << std::scientific << std::setprecision(15) << error[0] << ", "
                   << error[3] << std::endl;
             */
                
                EstimateErrorElasticity(config, cmesh, refsteps);
                
                // Prints gmesh mesh properties
//                std::string vtk_name = "geoMeshAfterAdapt.vtk";
//                std::ofstream vtkfile(vtk_name.c_str());
//                TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, vtkfile, true);
            
            
                
            }
            
        }
        
    }
    
   
}


//Create
template <class tshape>
TPZGeoMesh*
CreateGeoMesh(TPZVec<int> &nDivs, EMatid volId, EMatid bcId, REAL distortion)
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
    matIds[1] = bcId;
    matIds[2] = bcId;
    matIds[3] = bcId;
    matIds[4] = bcId;
    
    TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
                        matIds, nDivs, meshType,createBoundEls, distortion);
    // TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshSingleEl(meshType,
    //                     volId,createBoundEls, bcId);
    
    return gmesh;
    
}

void InsertMaterials(int &dim, TPZHDivApproxCreator& hdivCreator,TPZAnalyticSolution *fAn){

    if (hdivCreator.ProbType() == ProblemType::EDarcy){
        TPZMixedDarcyFlow* matdarcy = new TPZMixedDarcyFlow(EDomain,dim);
        // matdarcy->SetConstantPermeability(1.);
        TLaplaceExample1* lapl = dynamic_cast<TLaplaceExample1*> (fAn) ;
        matdarcy->SetExactSol(lapl->ExactSolution(),4);
        matdarcy->SetForcingFunction(lapl->ForceFunc(),4);

        hdivCreator.InsertMaterialObject(matdarcy);

        TPZFMatrix<STATE> val1(1,1,0.);
        TPZManVector<STATE> val2(1,0.);
        TPZBndCondT<STATE> *BCond1 = matdarcy->CreateBC(matdarcy, EBoundary, 0, val1, val2);
        BCond1->SetForcingFunctionBC(lapl->ExactSolution(),4);
        hdivCreator.InsertMaterialObject(BCond1);
    } else if (hdivCreator.ProbType() == ProblemType::EElastic){
        TElasticity2DAnalytic *elas2D;
        TElasticity3DAnalytic *elas3D;
        TPZMixedElasticityND* matelas;
        if (dim == 2) {
            elas2D = dynamic_cast<TElasticity2DAnalytic*> (fAn) ;
            matelas = new TPZMixedElasticityND(EDomain, elas2D->gE, elas2D->gPoisson, 0, 0, elas2D->fPlaneStress, dim);
            matelas->SetExactSol(elas2D->ExactSolution(),4);
            matelas->SetForcingFunction(elas2D->ForceFunc(),4);
            hdivCreator.InsertMaterialObject(matelas);

            TPZFMatrix<STATE> val1(dim,dim,0.);
            TPZManVector<STATE> val2(dim,0.);
            TPZBndCondT<STATE> *BCond1 = matelas->CreateBC(matelas, EBoundary, 0, val1, val2);
            BCond1->SetForcingFunctionBC(elas2D->ExactSolution(),4);
            hdivCreator.InsertMaterialObject(BCond1);
        }
        if (dim == 3) {
            elas3D = dynamic_cast<TElasticity3DAnalytic*> (fAn) ;
            matelas = new TPZMixedElasticityND(EDomain, elas3D->fE, elas3D->fPoisson, 0, 0, 0, dim);
            matelas->SetExactSol(elas3D->ExactSolution(),4);
            matelas->SetForcingFunction(elas3D->ForceFunc(),4);
            hdivCreator.InsertMaterialObject(matelas);

            TPZFMatrix<STATE> val1(dim,dim,0.);
            TPZManVector<STATE> val2(dim,0.);
            TPZBndCondT<STATE> *BCond1 = matelas->CreateBC(matelas, EBoundary, 0, val1, val2);
            BCond1->SetForcingFunctionBC(elas3D->ExactSolution(),4);
            hdivCreator.InsertMaterialObject(BCond1);
        }
        
    } else {
        DebugStop();//Material Not Implemented
    }
}

void InsertH1Materials(int &dim, TPZH1ApproxCreator& h1Creator,TPZAnalyticSolution *fAn){
    
    
    TElasticity2DAnalytic *elas2D;
    elas2D = dynamic_cast<TElasticity2DAnalytic*> (fAn) ;
    TPZElasticity2D *matelas= new TPZElasticity2D(EDomain, elas2D->gE, elas2D->gPoisson, 0, 0, elas2D->fPlaneStress);
    
    matelas->SetExactSol(elas2D->ExactSolution(),4);
    matelas->SetForcingFunction(elas2D->ForceFunc(),4);
    h1Creator.InsertMaterialObject(matelas);
    
    TPZFMatrix<STATE> val1(dim,dim,0.);
    TPZManVector<STATE> val2(dim,0.);
    TPZBndCondT<STATE> *BCond1 = matelas->CreateBC(matelas, EBoundary, 0, val1, val2);
    BCond1->SetForcingFunctionBC(elas2D->ExactSolution(),4);
    h1Creator.InsertMaterialObject(BCond1);
    
}

void EstimateErrorElasticity(ProblemConfig &config, TPZMultiphysicsCompMesh *originalMesh, int step) {

    std::cout << "\nError Estimation processing for Elasticity problem " << std::endl;

    // Error estimation
    if (!originalMesh) DebugStop();
    

    bool postProcWithHDiv = false;
    TPZElasticityErrorEstimator ErrorEstimator(config, *originalMesh, postProcWithHDiv);
    ErrorEstimator.SetAnalyticSolution(*config.exactElast);
    
    ErrorEstimator.PrimalReconstruction();
    
    TPZMultiphysicsCompMesh *PostProcMesh =ErrorEstimator.TPZHDivErrorEstimator<TPZMixedElasticityND>::PostProcMesh();
    
    TPZCompMesh *H1Mesh= PostProcMesh->MeshVector()[5];
    #ifdef ERRORESTIMATION_DEBUG
    {
        std::ofstream outTXT("H1MeshVector5.txt");
        H1Mesh->Print(outTXT);
    }
    #endif
    SolvingH1Displacement<pzshape::TPZShapeQuad>(H1Mesh,config);
    
    //PostProcMesh->MeshVector()[5]->ExpandSolution();
    #ifdef ERRORESTIMATION_DEBUG
    {
//        std::ofstream outTXT("SolutionH1.txt");
//        PostProcMesh->MeshVector()[5]->Solution().Print(outTXT);
        std::ofstream out1("SolutionH1.nb");
        PostProcMesh->MeshVector()[5]->ElementSolution().Print("SolutionH1 = ", out1, EMathematicaInput);
     //   std::ostream *out;
//        PostProcMesh->MeshVector()[5]->ElementSolution().Print("SolutionH1",std::cout);
        
    }
    #endif
    #ifdef ERRORESTIMATION_DEBUG
    {
        std::ofstream outTXT("H1MeshVector5AfterSol.txt");
        H1Mesh->Print(outTXT);
    }
    #endif
    
    #ifdef ERRORESTIMATION_DEBUG
    {
        std::ofstream outTXT("MultiMesh.txt");
        PostProcMesh->Print(outTXT);
    }
    #endif
    

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    TPZManVector<REAL, 6> errors;
    TPZManVector<REAL, 6> elementerrors;
    std::stringstream outVTK;
    outVTK << config.dir_name << "/" << config.problemname << "-" << config.ndivisions << "-k-" << config.porder<<"-lambda-"<<config.lambda
           << "-step "<<step << "-Errors.vtk";
    ErrorEstimator.ComputeErrors(errors, elementerrors, outVTK.str());

    #ifdef ERRORESTIMATION_DEBUG
    {
        std::string vtk_name = "geoMeshBeforeAdapt.vtk";
        std::ofstream vtkfile(vtk_name.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, vtkfile, true);
    }
    #endif
    
    if(config.isAdaptivity){
        Tools::hAdaptivity(ErrorEstimator.PostProcMesh(), config.gmesh, originalMesh, config);
    }
   
    int nel=config.gmesh->NElements();
    for (int iel=0; iel<nel; iel++) {
        TPZGeoEl * geo =config.gmesh->ElementVec()[iel];
        if (!geo) {
            continue;
        }
        int matId = geo->MaterialId();
        if ((matId==config.fSkeletonMatId || matId == config.fHangingNodeMatId) && geo->FatherIndex() == -1) {
            config.gmesh-> DeleteElement(geo);
        }
    }

   {
       std::string fileName = config.dir_name + "/" + config.problemname + "-GlobalErrors.txt";
       std::ofstream file(fileName, std::ios::app);
       Tools::PrintElasticityErrors(file, config, errors);
   }
    
    // Prints gmesh mesh properties
    std::string vtk_name = "geoMeshAfterAdapt_1.vtk";
    std::ofstream vtkfile(vtk_name.c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, vtkfile, true);
}

TPZGeoMesh*
ReadMeshFromGmsh(std::string file_name)
{
    //read mesh from gmsh
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    {
        TPZGmshReader reader;
        // essa interface permite voce mapear os nomes dos physical groups para
        // o matid que voce mesmo escolher
        TPZManVector<std::map<std::string,int>,4> stringtoint(4);
        stringtoint[2]["Domain"] = 1;
        stringtoint[1]["Boundaries"] = 2;
        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh(file_name,gmesh);
    }
    return gmesh;
}

template<class tshape>
void SolveFEMProblemNew(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily,ProblemConfig &config)
{

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    int DIM = tshape::Dimension;

        for(int refsteps = 1; refsteps <= config.adaptivityStep; refsteps ++){
            config.refStepCounter = refsteps;
            #ifdef ERRORESTIMATION_DEBUG
            {
                // Prints gmesh mesh properties
                std::string vtk_name = "geoMeshToSolveProblem.vtk";
                std::ofstream vtkfile(vtk_name.c_str());
                TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, vtkfile, true);
            }
            #endif
            
                // Creates an hdivApproxCreator object. It is an environment developped to
                // help creating H(div)-family possible approximation spaces.
                TPZHDivApproxCreator hdivCreator(config.gmesh);
                //Set the family of H(div) functions: Standard, Constant or Kernel
                hdivCreator.HdivFamily() = hdivfamily;
                //Set the problem type to be solved: Only EDarcy and EElastic are currently available
                hdivCreator.ProbType() = ProblemType::EElastic;
                //Includes the rigid body spaces (constant flux and pressure) if set as true
                hdivCreator.IsRigidBodySpaces() = false;
                //Set the default polynomial order
                hdivCreator.SetDefaultOrder(config.porder);
                //Set the extra polynomial order for the bubble functions. If zero, the polynomial degree
                //of the internal functions are the same as the default order
                hdivCreator.SetExtraInternalOrder(config.hdivmais);
                //Sets if the resulting problem should or not be condensed
              //  hdivCreator.SetShouldCondense(true);
                hdivCreator.SetShouldCondense(false);
                
                //Sets the type of hybridizantion desired.
                //The current options are HybridizationType::ENone, HybridizationType::EStandard
                //and HybridizationType::ESemi (the last only works with H(div)-constant spaces)
                hdivCreator.HybridType() = HybridizationType::ENone;
                //hdivCreator.HybridType() = HybridizationType::EStandard;
                
                //Creates an analytical solution to test
                TPZAnalyticSolution *gAnalytic = 0;
                if (hdivCreator.ProbType() == ProblemType::EDarcy){
                    TLaplaceExample1 *lap = new TLaplaceExample1;
                    lap->fExact = TLaplaceExample1::EHarmonic;
                    gAnalytic = lap;
                } else if (hdivCreator.ProbType() == ProblemType::EElastic){
                    if (DIM == 2){
                        TElasticity2DAnalytic *elas = new TElasticity2DAnalytic;
                        elas->gE = config.mu*(3*config.lambda+2*config.mu)/(config.lambda+config.mu);
                        elas->gPoisson = 0.5*config.lambda/(config.lambda+config.mu);
                        elas->fProblemType = config.exactElast->fProblemType;
                        elas->fPlaneStress=0;
                        gAnalytic = elas;
                        //std::cout<<"-----planestress no misto "<<elas->fPlaneStress<<std::endl;
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
         
                
                //Insert Materials
                InsertMaterials(DIM,hdivCreator,gAnalytic);
                
                //Gets the Multiphysics mesh from the HdivApproxCreator
            
                TPZMultiphysicsCompMesh *cmesh = nullptr;// = hdivCreator.CreateApproximationSpace();
                Tools::PRefinementNew(cmesh, config, hdivCreator);
            
            // TPZMultiphysicsCompMesh *cmesh =hdivCreator.CreateApproximationSpace();
                
                //Create the analysis environment
                TPZLinearAnalysis an(cmesh,RenumType::ESloan);
                // TPZLinearAnalysis an(cmesh,RenumType::EMetis);
                an.SetExact(gAnalytic->ExactSolution(),4);
                
#ifdef PZ_USING_MKL
                TPZSSpStructMatrix<> strmat(cmesh);
                strmat.SetNumThreads(1);
#else
                TPZSkylineStructMatrix<STATE> strmat(cmesh);
                strmat.SetNumThreads(4);
#endif
                
                an.SetStructuralMatrix(strmat);
                TPZStepSolver<STATE> step;
                step.SetDirect(ELDLt);
                an.SetSolver(step);
                std::cout << "Assembling\n";
                an.Assemble();
                
                std::cout << "Solving\n";
                an.Solve();
                std::cout << "Finished\n";
                an.LoadSolution(); // compute internal dofs
            
                #ifdef ERRORESTIMATION_DEBUG
                {
                    TPZStack<std::string> vecnames,scalnames;
                    vecnames.Push("Displacement");
                    scalnames.Push("Rotation");
                    scalnames.Push("SigmaX");
                    scalnames.Push("SigmaY");
                    scalnames.Push("SigmaXY");
                    scalnames.Push("SigmaYX");
                    scalnames.Push("TauXY");
                    vecnames.Push("ExactDisplacement");
                    vecnames.Push("ExactStrain");
                    vecnames.Push("ExactStress");
                // scalnames.Push("POrder");
                
                    int dim = 2;

                    std::string refsteps_str = std::to_string(refsteps);
                    std::string filename = "SolutionFEM_step" + refsteps_str + ".vtk";

                    an.DefineGraphMesh(dim, scalnames, vecnames, filename);
                    an.PostProcess(4, dim);
                }
                #endif
                
                std::ofstream anPostProcessFile("PostprocessFem.txt");
                TPZManVector<REAL,7> error(7,0);
                an.LoadSolution();
                cmesh->LoadSolution(cmesh->Solution());
                cmesh->ExpandSolution();
                int64_t nelem = an.Mesh()->NElements();
                
                
                an.Mesh()->ElementSolution().Redim(nelem, 7);
                // an.SetExact(elas->ExactSolution(),5);
                an.PostProcessError(error,true,anPostProcessFile);

                std::cout << "Error FEM: " << error << std::endl;



                EstimateErrorElasticity(config, cmesh, refsteps);
                
            }
    
   
}

template<class tshape>
void SolveH1Problem(const int &xdiv, const int &pOrder, H1Family &h1family,ProblemConfig &config)
{

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    int DIM = tshape::Dimension;
        
        TElasticity2DAnalytic *elas = new TElasticity2DAnalytic;
        elas->gE = config.mu*(3*config.lambda+2*config.mu)/(config.lambda+config.mu);
        elas->gPoisson = 0.5*config.lambda/(config.lambda+config.mu);
        elas->fProblemType = config.exactElast->fProblemType;
        elas->fPlaneStress=0;
        //gAnalytic = elas;
        
       // std::cout<<"fPlaneStress para H1 "<<elas->fPlaneStress<< std::endl;
    
        
        TPZCompMesh *cmesh=CreateH1CMesh(config.gmesh, config.porder,elas,config);
        {
            std::ofstream outTXT("H1MeshSolveH1Problem.txt");
            cmesh->Print(outTXT);
        }
       
        
        //Create the analysis environment
        TPZLinearAnalysis an(cmesh,RenumType::ESloan);
        // TPZLinearAnalysis an(cmesh,RenumType::EMetis);
        an.SetExact(elas->ExactSolution(),4);
        
#ifdef PZ_USING_MKL
        TPZSSpStructMatrix<> strmat(cmesh);
        strmat.SetNumThreads(1);
#else
        TPZSkylineStructMatrix<STATE> strmat(cmesh);
        strmat.SetNumThreads(4);
#endif
        
        an.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        an.SetSolver(step);
        std::cout << "Assembling\n";
        an.Assemble();
        
        std::cout << "Solving\n";
        an.Solve();
        std::cout << "Finished\n";
        an.LoadSolution(); // compute internal dofs
        
//        {
//            std::ofstream outTXT("H1MeshSolveH1ProblemPostSol.txt");
//            cmesh->Print(outTXT);
//        }
        
        
        
        //Compute error
        std::ofstream anPostProcessFile("PostprocessFem.txt");
        TPZManVector<REAL,7> error(6,0);
        an.LoadSolution();
        cmesh->LoadSolution(cmesh->Solution());
        cmesh->ExpandSolution();
        int64_t nelem = an.Mesh()->NElements();
        
        
        an.Mesh()->ElementSolution().Redim(nelem, 6);
        an.SetExact(elas->ExactSolution(),5);
        an.PostProcessError(error,true,anPostProcessFile);
        
       
        //std::cout<<"L2 error for displacement= "<<error[1]<<std::endl;
        
        
        std::ofstream outFile("TaxaAproxH1.txt", std::ios::app);
        
        outFile << "Erro Norma L2 disp = "<<error[1]<<"\n";
    
        
        //PosProcess Graph
        
//        TPZStack<std::string> vecnames,scalnames;
//        vecnames.Push("Displacement");
//        vecnames.Push("DisplacementExact");
//        std::stringstream out;
//        out << config.dir_name << "/" << "Fem_Solution_k" << config.porder << "Nref_" << config.ndivisions
//                << "NAdapStep_" << config.adaptivityStep << ".vtk";
//        an.DefineGraphMesh(DIM, scalnames, vecnames, out.str());
//        an.PostProcess(0, DIM);

}

TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity2DAnalytic *elas2D,ProblemConfig &config) {
    
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    const int dim = gmesh->Dimension();
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pord);
    cmesh->SetAllCreateFunctionsContinuous();
    
    // Domain elas mat
    const STATE E = elas2D->gE, nu = elas2D->gPoisson;
   

    TPZElasticity2D *mat= new TPZElasticity2D(EDomain, elas2D->gE, elas2D->gPoisson, 0, 0, elas2D->fPlaneStress);
    
    mat->SetExactSol(elas2D->ExactSolution(), 4);
    mat->SetForcingFunction(elas2D->ForceFunc(), 4);
    cmesh->InsertMaterialObject(mat);

    
    // BC
    
    TPZFMatrix<STATE> val1(dim,dim,0.);
    TPZManVector<STATE> val2(dim,0.);
    TPZBndCondT<STATE> *BCond1 = mat->CreateBC(mat, EBoundary, 0, val1, val2);
    BCond1->SetForcingFunctionBC(elas2D->ExactSolution(),4);
    cmesh->InsertMaterialObject(BCond1);

    // Constructs mesh
    cmesh->AutoBuild();
    
    return cmesh;
    
}
template<class tshape>
void SolvingH1Displacement(TPZCompMesh *cH1Mesh,ProblemConfig &config){
    
//    {
//        std::ofstream outTXT("H1MeshInitial-2.txt");
//        cH1Mesh->Print(outTXT);
//    }
    
            //Create the analysis environment
            TPZLinearAnalysis an(cH1Mesh,RenumType::ESloan);
            an.SetExact(config.exactElast->ExactSolution(),4);
    
    #ifdef PZ_USING_MKL
            TPZSSpStructMatrix<> strmat(cH1Mesh);
            strmat.SetNumThreads(1);
    #else
            TPZSkylineStructMatrix<STATE> strmat(cH1Mesh);
            strmat.SetNumThreads(4);
    #endif
    
            an.SetStructuralMatrix(strmat);
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt);
            an.SetSolver(step);
            std::cout<< "Assembling\n";
            an.Assemble();
            std::cout<< "Solving\n";
            an.Solve();
            std::cout << "Finished\n";
            an.LoadSolution(); // compute internal dofs
    
    {
        std::ofstream outTXT("H1MeshAfterSolution.txt");
        cH1Mesh->Print(outTXT);
    }
    
    
    
    TPZStack<std::string> vecnames,scalnames;
        vecnames.Push("Displacement");
    vecnames.Push("DisplacementExact");
    
        int dim = 2;

        an.DefineGraphMesh(dim, scalnames, vecnames, "SolutionH1.vtk");
        an.PostProcess(0, dim);
    
    //Compute error
  
    std::ofstream anPostProcessFile("PostprocessFemH1.txt");
    TPZManVector<REAL,7> error(6,0);
    an.LoadSolution();
    cH1Mesh->LoadSolution(cH1Mesh->Solution());
    cH1Mesh->ExpandSolution();
    int64_t nelem = an.Mesh()->NElements();
    
    
    an.Mesh()->ElementSolution().Redim(nelem, 6);
    an.PostProcessError(error,true,anPostProcessFile);
    
//    printerrors << config.porder <<  "L2 Error for Disp= "<<error[1]<<std::endl;
//    std::cout<<"L2 Error for Disp= "<<error[1]<<std::endl;
   // std::cout<<"Vector Error= "<<error<<std::endl;
    
    

    
}
