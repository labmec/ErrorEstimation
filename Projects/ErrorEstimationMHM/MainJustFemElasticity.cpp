//
//  MainJustFemElasticity.cpp
//  MHM_ElastPhil
//
//  Created by Denise De Siqueira on 12/12/23.
//

#include <stdio.h>


//#include "InputTreatment.h"
//#include "Solver.h"
//#include "Output.h"
//#include "tpzgeoelrefpattern.h"
//#include "DataStructure.h"
//#include "Tools.h"
//#include <tuple>

//%%%%%%

//#include "common_files.h"

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
std::ofstream printerrors("results_errors2.txt",std::ios::app);

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

void InsertMaterials(int &dim, TPZHDivApproxCreator& hdivc, TPZAnalyticSolution *fAn);
void EstimateError(ProblemConfig &config, TPZMultiphysicsCompMesh *multimesh);
void EstimateErrorElasticity(const ProblemConfig &config, TPZMultiphysicsCompMesh *originalMesh);

//%%%%


int main() {
    

    ProblemConfig pConfig;
    pConfig.geometry = ProblemConfig::EGeometry::ECrack;
    pConfig.exactElast = new TElasticity2DAnalytic;
    //pConfig.exactElast.operator*().fProblemType = TElasticity2DAnalytic::EDispy;
    switch (pConfig.geometry){
        case ProblemConfig::EGeometry::ECrack:
            pConfig.exactElast->fProblemType = TElasticity2DAnalytic::ECrack;
            break;
        case ProblemConfig::EGeometry::EQuad:
        case ProblemConfig::EGeometry::ETrap:
        default:
            pConfig.exactElast->fProblemType = TElasticity2DAnalytic::EThiago;
            break;
    }
    
    const int xdiv = 3; //Number of elements in each direction
    const int pOrder = 1; // Polynomial degree
    pConfig.porder = pOrder;
    pConfig.ndivisions = xdiv;
    pConfig.hdivmais = 2;// internal order
   
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
    
    int DIM = tshape::Dimension;
    TPZVec<int> nDivs;
    TPZVec<int> divs = {4};//{2,4,8,16,32,64};//{2,5,10,20,50,100};
    
    int pend = 2;
    
    TPZGeoMesh *gmesh;

    switch(config.geometry) {
        case ProblemConfig::EGeometry::ECrack:
        {
            gmesh = ReadMeshFromGmsh("../../../Crack.msh");
            // TPZVec<int> bcids(8, EBoundary);
            // gmesh = Tools::CreateQuadLShapeMesh(bcids);

            // int uniref = 3;
            // Tools::UniformRefinement(uniref, gmesh);
            break;
        }
        case ProblemConfig::EGeometry::ELShape:
        {
            TPZVec<int> bcids(8, EBoundary);
            gmesh = Tools::CreateQuadLShapeMesh(bcids);

            int uniref = 3;
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
   
    // Prints gmesh mesh properties
    std::string vtk_name = "geoMesBeforeAdapt.vtk";
    std::ofstream vtkfile(vtk_name.c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
    
    int nsteps =4;
    config.gmesh = gmesh;
    
    for (int iorder = 1; iorder < pend; iorder++){
        
        printerrors <<  " porder " << " h " <<   " error stress "<< " error diplacement"<<std::endl;
                    
        for(int refsteps = 1; refsteps< nsteps; refsteps ++){
        
//            for (int idiv = 0; idiv < divs.size(); idiv++){
//                config.ndivisions =divs[idiv];
//                config.porder= iorder;
//                int divx = divs[idiv];
                
                
                //    if (DIM == 2) nDivs = {divx,divx};
                //    if (DIM == 3) nDivs = {divx,divx,divx};
                
                // Creates/import a geometric mesh
                //  auto gmesh = CreateGeoMesh<tshape>(nDivs, EDomain, EBoundary);
                
                
                //  int uniref= divx/2;
                //  Tools::UniformRefinement(uniref, gmesh);
                
                // Prints gmesh mesh properties
                std::string vtk_name3 = "geoMeshLshape0.vtk";
                std::ofstream vtkfile3(vtk_name3.c_str());
                TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, vtkfile3, true);
                
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
                // hdivCreator.HybridType() = HybridizationType::ENone;
                hdivCreator.HybridType() = HybridizationType::EStandard;
                
                // Prints gmesh mesh properties
                std::string vtk_name2 = "geoMeshLshape.vtk";
                std::ofstream vtkfile2(vtk_name2.c_str());
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile2, true);
                
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
                        double lambda = 5;//123.;
                        double mu = 1;//79.3;
                        elas->gE = mu*(3*lambda+2*mu)/(lambda+mu);
                        elas->gPoisson = 0.5*lambda/(lambda+mu);
                        //Crack
                        elas->gE = 100.;
                        elas->gPoisson = 0.3;
                        // elas->fProblemType = TElasticity2DAnalytic::EDispx;
                        elas->fProblemType = config.exactElast->fProblemType;
                        elas->fPlaneStress = 0;
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

                config.fWrapMaterialId = hdivCreator.HybridData().fWrapMatId;
                config.fInterfaceMaterialId = hdivCreator.HybridData().fInterfaceMatId;
                config.fLagMultiplierMaterialId = hdivCreator.HybridData().fLagrangeMatId;
                
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
                
                //----
                
#ifdef PZ_USING_MKL
                TPZSSpStructMatrix<> strmat(cmesh);
                strmat.SetNumThreads(8);
#else
                TPZSkylineStructMatrix<STATE> strmat(cmesh);
                strmat.SetNumThreads(0);
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
                
                
                
                
                //Printing results in a vtk file
                TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(cmesh->MeshVector(), an.Mesh());
                cmesh->LoadSolution(cmesh->Solution());
                std::string plotfile = "res_h"+std::to_string(refsteps)+"p"+std::to_string(iorder);//sem o .vtk no final
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
                 cmesh->LoadSolution(cmesh->Solution());
                 cmesh->ExpandSolution();
                int64_t nelem = an.Mesh()->NElements();//0;
//                for (int64_t i = 0; i < an.Mesh()->NElements(); i++)
//                {
//                    auto el = an.Mesh()->ElementVec()[i];
//                    if (el->Dimension() == DIM) nelem++;
//                }
                
                an.Mesh()->ElementSolution().Redim(nelem, 7);
                an.SetExact(gAnalytic->ExactSolution(),5);
                an.PostProcessError(error,true,anPostProcessFile);
                
//                std::ofstream postVTK(plotfile+".0.vtk",std::ios::app);
//                TPZFMatrix<STATE> SolMatrix = an.Mesh()->ElementSolution();
//                postVTK << "CELL_DATA " << nelem << std::endl;
//                postVTK << "SCALARS ErrorStress float \nLOOKUP_TABLE default" << std::endl;
//                for (int64_t iel = 0; iel < nelem; iel++){
//                    postVTK << SolMatrix.GetVal(iel,0) << " ";
//                }
//                postVTK << std::endl;
//                postVTK << "SCALARS ErrorEnergy float \nLOOKUP_TABLE default" << std::endl;
//                for (int64_t iel = 0; iel < nelem; iel++){
//                    postVTK << SolMatrix.GetVal(iel,1) << " ";
//                }
//                postVTK << std::endl;
//                postVTK << "SCALARS ErrorDivStress float \nLOOKUP_TABLE default" << std::endl;
//                for (int64_t iel = 0; iel < nelem; iel++){
//                    postVTK << SolMatrix.GetVal(iel,2) << " ";
//                }
//                postVTK << std::endl;
//                postVTK << "SCALARS ErrorDisplacement float \nLOOKUP_TABLE default" << std::endl;
//                for (int64_t iel = 0; iel < nelem; iel++){
//                    postVTK << SolMatrix.GetVal(iel,3) << " ";
//                }
//                postVTK << std::endl;
//                postVTK << "SCALARS ErrorRotation float \nLOOKUP_TABLE default" << std::endl;
//                for (int64_t iel = 0; iel < nelem; iel++){
//                    postVTK << SolMatrix.GetVal(iel,4) << " ";
//                }
//                postVTK << std::endl;
//                postVTK << "SCALARS ErrorSymmetry float \nLOOKUP_TABLE default" << std::endl;
//                for (int64_t iel = 0; iel < nelem; iel++){
//                    postVTK << SolMatrix.GetVal(iel,5) << " ";
//                }
//                postVTK << std::endl;
                
                
                // //Print Errors
                //     std::cout << "POrder = " << iorder << std::endl;
                //     std::cout << "h = " << std::fixed <<  1./double(divx) << std::endl;
                //     std::cout << "L2 Stress = " << std::scientific << std::setprecision(15) << error[0] << std::endl;
                //     std::cout << "En Stress = " << error[1] << std::endl;
                //     std::cout << "L2 DivStr = " << error[2] << std::endl;
                //     std::cout << "L2 Displa = " << error[3] << std::endl;
                //     std::cout << "L2 Rotati = " << error[4] << std::endl;
                //     std::cout << "L2 Symmet = " << error[5] << std::endl;
                //     std::cout << "En Displa = " << error[6] << std::endl;
                //     printerrors << iorder << " " << std::fixed << std::setprecision(5) <<  1./double(divx) << " "
                //                 << std::scientific << std::setprecision(15) << error[0] << " "
                //                 << error[1] << " " << error[2] << " " << error[3] << " " << error[4]
                //                 << " " << error[5] << " " << error[6] << std::endl;
                
                printerrors << iorder << ", " << std::fixed << std::setprecision(5) <<  1./double(refsteps) << ", "
                << std::scientific << std::setprecision(15) << error[0] << ", "
                << error[3] << std::endl;
                
                EstimateErrorElasticity(config, cmesh);
                
                // Prints gmesh mesh properties
                std::string vtk_name = "geoMeshAfterAdapt.vtk";
                std::ofstream vtkfile(vtk_name.c_str());
                TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, vtkfile, true);
                
 //           }
            
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
    matIds[2] = EBoundary;
    matIds[3] = EBoundary;
    matIds[4] = EBoundary;
    
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
         //   matelas->SetForcingFunction(elas2D->ForceFunc(),4);
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



void EstimateErrorElasticity(const ProblemConfig &config, TPZMultiphysicsCompMesh *originalMesh) {

    std::cout << "\nError Estimation processing for MHM-Hdiv problem " << std::endl;

    // Error estimation
    if (!originalMesh) DebugStop();
    

    bool postProcWithHDiv = false;
    TPZElasticityErrorEstimator ErrorEstimator(config, *originalMesh, postProcWithHDiv);
    ErrorEstimator.SetAnalyticSolution(*config.exactElast);
    
    ErrorEstimator.PrimalReconstruction();

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    TPZManVector<REAL, 6> errors;
    TPZManVector<REAL, 6> elementerrors;
    std::stringstream outVTK;
    outVTK << config.dir_name << "/" << config.problemname << "-" << config.ndivisions << "-" << config.ninternalref
           << "-Errors.vtk";
    std::string outVTKstring = outVTK.str();
    ErrorEstimator.ComputeErrors(errors, elementerrors, outVTKstring);
    
    Tools::hAdaptivity(ErrorEstimator.PostProcMesh(), config.gmesh, config);
    //TODO
    int nel=config.gmesh->NElements();
    for (int iel=0; iel<nel; iel++) {
        TPZGeoEl * geo =config.gmesh->ElementVec()[iel];
        if (!geo) {
            continue;
        }
        int matId = geo->MaterialId();
        if (matId==3) {
            config.gmesh-> DeleteElement(geo);
        }
    }

//    {
//        std::string fileName = config.dir_name + "/" + config.problemname + "-GlobalErrors.txt";
//        std::ofstream file(fileName, std::ios::app);
//        Tools::PrintElasticityErrors(file, config, errors);
//    }
    
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