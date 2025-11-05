/**
 * @file This file implements an error estimator in space H1.
 */

#include "pzlog.h"

#include "TPZMultiphysicsCompMesh.h"
#include "TPZMultiphysicsInterfaceEl.h"

#include "TPZH1ApproxCreator.h"

#include "TPZLinearAnalysis.h"
#include "pzstepsolver.h"
#include "TPZBndCondT.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "TPZH1ErrorHybridH1EstimateMaterial.h"

#include "DarcyFlow/TPZHybridDarcyFlow.h"
#include "Elasticity/TPZElasticity2D.h"
#include "TPZHybridElasticity2D.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#include "TPZLagrangeMultiplierCS.h"
#include "pzgeoelbc.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZRefPatternDataBase.h"
#include "TPZVTKGeoMesh.h"
#include "TPZVTKGenerator.h"
#include "TPZGmshReader.h"
#include "TPZAnalyticSolution.h"
#include "TPZFrontSym.h"

#include "TPZBuildSBFemHybrid.h"
#include "TPZSBFemElementGroup.h"

#include "TPZPostProcessErrorSBFem.h"

#include <fstream>
#include <ctime>
#include <cstdio>
#include <cmath>

#include "TPZCreateHybridH1Space.h"

#include <iostream>
// Global variables

#include <nlohmann/json.hpp>

using json = nlohmann::json;

/// @brief  read the geometric mesh for the file name defined in the json input
/// @param input json object with the input data
/// @return pointer to the geometric mesh
TPZGeoMesh *ReadGmshFile(json &input);

/// @brief Refine the geometric mesh
/// @param geometric mesh object
/// @param input json data structure
/// @param number of refinements
void UniformRefine(TPZGeoMesh *gmesh, json &input, int iref);

/// @brief Create collapsed elements towards the singular element
std::set<int64_t> CreateCollapsedElements(TPZGeoMesh *gmesh, int64_t singular_el, TPZBuildSBFem &builder);

/// @brief divide randomly elements of the geometric mesh
/// @param gmesh pointer to the geometric mesh
/// @param nref number of refinements
void RefineRandomElement(TPZGeoMesh *gmesh, int nref);
/// @brief configure the SBFem builder according to the json input
/// @param builder reference to the SBFem builder
/// @param input json object with the input data
void ConfigureSBFemBuilder (TPZBuildSBFemHybrid &builder, json &input);

/// @brief configure the SBFem builder according to the json input
/// @param builder reference to the SBFem builder
/// @param input json object with the input data
void ConfigureSBFemBuilder (TPZBuildSBFem &builder, json &input);

/// @brief insert material objects into the H1 computational mesh according to the json input
/// @param cmesh reference to the computational mesh
/// @param input json object with the input data
void InsertMaterialObjectsH1(TPZCompMesh &cmesh, json &input);
/// @brief insert material objects into the Hybrid H1 computational mesh according to the json input
/// @param cmesh reference to the computational mesh
/// @param input json object with the input data
void InsertMaterialObjectsHybridH1(TPZCompMesh &cmesh, json &input);
/// @brief insert interface material objects into the computational mesh according to the json input
/// @param cmesh reference to the computational mesh
/// @param input json object with the input data
void InsertInterfaceMaterialObjects(TPZCompMesh &cmesh, json &input);

/// @brief Delete the condensed and group elements and restore the original elements
/// @param cmesh computational mesh with group elements
void UnwrapMesh(TPZCompMesh &cmesh);

/// @brief Create a computational mesh based on the SBFem hybrid H1 space
/// @param gmesh autopointer to the geometric mesh
/// @param input json object with the input data
/// @return pointer to the computational mesh
TPZCompMesh *CreateHybridSBFemSpace( TPZBuildSBFemHybrid &builder, json &input);

/// @brief Create a computational mesh based on the SBFem hybrid H1 space
/// @param gmesh autopointer to the geometric mesh
/// @param input json object with the input data
/// @return pointer to the computational mesh
TPZCompMesh *CreateHybridSBFemSpaceFromSBFem(TPZBuildSBFemHybrid &buildSBFem, json &input);

/// @brief Create a computational mesh based on the SBFem hybrid H1 space
/// @param gmesh autopointer to the geometric mesh
/// @param input json object with the input data
/// @return pointer to the computational mesh
TPZCompMesh *CreateSBFemSpace(TPZBuildSBFem &builder, json &input);

/// @brief Solve the linear system defined by the json input and the computational mesh
/// @param input json object with the input data
/// @param cmesh reference to the computational mesh
void SolveSystem(json &input, TPZCompMesh &cmesh);

/// @brief Post-process the solution of the computational mesh (compute the approximation error)
/// @param cmesh reference to the computational mesh
/// @param output json object with the output data
void PostProcessSolution(TPZCompMesh &cmesh, int iref, json &output, TPZVec<REAL> &errors);

/// @brief Plot the solution of the computational mesh
/// @param cmesh reference to the computational mesh
/// @param output json object with the output data
void PlotSolution(const std::string &rootname, TPZCompMesh &cmesh, int refstep, json &output);

/// @brief Perform a convergence study applying uniform refinements
/// @param gmesh autopointer to the geometric mesh
/// @param input json object with the input data
void ConvergenceStudy(TPZAutoPointer<TPZGeoMesh> gmesh, json &input, bool append = true);

/// @brief Compute the difference between two meshes
/// @param cmesh1 reference to the first computational mesh
/// @param cmesh2 reference to the second computational mesh
void CompareMeshes(TPZCompMesh &cmesh1, TPZCompMesh &cmesh2, TPZVec<REAL> &errors);
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.refine"));
#endif

std::string gProblem;
std::map<std::string,TLaplaceExample1> gExact1;
std::map<std::string,TElasticity2DAnalytic> gExact2;


/// @brief Add an exact solution to the map
void AddExact(const std::string &problemtype, const std::string &exactname) {
    if(problemtype == "Scalar") {
        if(gExact1.find(exactname) != gExact1.end()) return;
        TLaplaceExample1::EExactSol enumfunc = TLaplaceExample1::StringToExactSol(exactname);
        if(enumfunc != TLaplaceExample1::ENone) {
            gExact1[exactname].fExact = enumfunc;
        } else {
            DebugStop();
        }
    } else if (problemtype == "Elastic") {
        if(gExact2.find(exactname) != gExact2.end()) return;
        TElasticity2DAnalytic::EDefState enumfunc = TElasticity2DAnalytic::StringToExactSol(exactname);
        if(enumfunc != TElasticity2DAnalytic::ENone) {
            gExact2[exactname].fProblemType = enumfunc;
        } else {
            DebugStop();
        }
    }

}

std::complex<STATE> integrateF = 0;
//bool Print = false;

int main(int argc, char *argv[]) {
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
  	using json = nlohmann::json;
    std::cout << "Reading input file\n";
    std::string filename = "hybrid.json";
//    std::string filename = "SqrtSingular.json";
//    std::string filename = "elasticsmooth.json";
//    std::string filename = "SqrtSingularElastic.json";
#ifdef MACOSX
    filename = "../" + filename;
#endif
	// Read file
	std::ifstream file(filename);
	if (!file){
		std::cout << "\nCouldn't find file \"" << filename << "\""<< std::endl;
		DebugStop();
	}

	// Parse json
	json input;
	// file >> input;
	input = json::parse(file,nullptr,true,true); // to ignore comments in json file
    std::string problemtype = input["problem_type"];
    std::vector<int> refsteps = input["UniformRefinement"].get<std::vector<int>>();
;
    if(refsteps.size() != 2) DebugStop();
    TPZAutoPointer<TPZGeoMesh> gmeshorig = ReadGmshFile(input);
    for(int iref = refsteps[0]; iref <= refsteps[1]; iref++) {
        TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh(gmeshorig);
        UniformRefine(gmesh.operator->(), input, iref);
        //RefineRandomElement(gmesh.operator->(), 1);
        {
            std::ofstream out("GeoMesh.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
            std::ofstream out2("GeoMesh.txt");
            gmesh->Print(out2);
        }
        
        TPZBuildSBFem builder(gmesh);
        
        TPZCompMesh *cmeshH1ptr = nullptr;
        TPZCompMesh *cmeshHybridptr = nullptr;
        cmeshH1ptr = CreateSBFemSpace(builder, input);
        {
            TPZFMatrix<STATE> &sol = cmeshH1ptr->Solution();
            cmeshH1ptr->LoadSolution(sol);
        }
        {
            std::ofstream out("GeoMesh.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
            std::ofstream out2("GeoMesh.txt");
            gmesh->Print(out2);
            std::ofstream out3("CompMesh.txt");
            cmeshH1ptr->Print(out3);
        }
        TPZBuildSBFemHybrid builderHybrid(builder);
        gmesh->ResetReference();
        // cmeshHybridptr = CreateHybridSBFemSpace(builderHybrid, input);
        cmeshHybridptr = CreateHybridSBFemSpaceFromSBFem(builderHybrid, input);
        gmesh->ResetReference();
        {
            TPZFMatrix<STATE> &sol = cmeshH1ptr->Solution();
            cmeshH1ptr->LoadSolution(sol);
        }
        
        TPZCompMesh &cmeshH1 = *cmeshH1ptr;
        TPZCompMesh &cmeshHybrid = *cmeshHybridptr;
        {
            std::ofstream out("CompMeshHybrid.txt");
            cmeshHybrid.Print(out);
        }
        SolveSystem(input, cmeshH1);
        if(1)
        {
            // Plot the solution
            std::cout << "Plotting the solution\n";
            PlotSolution("H1",cmeshH1, iref, input);
            std::cout << "Plotting the solution done\n";
            int nerrors = 3;
            if(problemtype == "Scalar") {
                nerrors = 3;
            } else if(problemtype == "Elastic") {
                nerrors = 6;
            }
            TPZManVector<REAL> errors(nerrors,0.);
            cmeshH1.ExpandElementSolution(nerrors);
            bool store_errors = true;
            cmeshH1.EvaluateError(store_errors, errors, builder.GetMaterialIds());
            std::cout << "H1 Errors = " << errors << std::endl;
            PostProcessSolution(cmeshH1, iref, input, errors);
        }
        bool computehybrid = true;
        if(computehybrid) {
            SolveSystem(input, cmeshHybrid);
            int nerrors = 3;
            if(problemtype == "Scalar") {
                nerrors = 3;
            } else if(problemtype == "Elastic") {
                nerrors = 6;
            }
            TPZManVector<REAL> errors(nerrors,0.);
            cmeshHybrid.ExpandElementSolution(nerrors);
            bool store_errors = true;
            cmeshHybrid.EvaluateError(store_errors, errors, builder.GetMaterialIds());
            std::cout << "Hybrid H1 Errors = " << errors << std::endl;
            std::cout << "Plotting the solution\n";
            PlotSolution("Hybrid",cmeshHybrid, iref, input);
            std::cout << "Plotting the solution done\n";
            // Post-process the solution
            PostProcessSolution(cmeshHybrid, iref, input, errors);
            //        TPZVec<REAL> errors(5,0.);
            //        errors.Resize(5, 0.);
            //        errors.Fill(0.);
            //        CompareMeshes(cmeshH1, cmeshHybrid, errors);
            //        std::cout << "Errors between H1 and Hybrid H1\n";
            //        std::cout << "H1 = " << errors[0] << " Hybrid H1 = " << errors[1] << " Estimate = " << errors[2] <<
            //        std::endl;
            //        std::cout << "H1 sq = " << errors[0]*errors[0] << " Hybrid H1 sq = " << errors[1]*errors[1] << " Estimate sq = " << errors[2]*errors[2] <<
            //        std::endl;
        }
        
        
        int refstep = 0;
        // builder.Print();
        if(1)
        {
            std::ofstream out("SBFemGeoMesh.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
            std::ofstream out2("SBFemGeoMesh.txt");
            gmesh->Print(out2);
            std::ofstream out3("SBFemCompMeshH1.txt");
            cmeshH1.Print(out3);
            std::ofstream out4("SBFemCompMeshHybridH1.txt");
            cmeshHybrid.Print(out4);
        }
        // this is necessary to be able to include the hybrid sbfem mesh as an atomic mesh
        UnwrapMesh(*cmeshHybridptr);
        TPZPostProcessErrorSBFem postproc(&cmeshH1, &cmeshHybrid, builderHybrid);
        postproc.BuildPatchStructures2();
        bool plotpatches = true;
        if(plotpatches) {
            postproc.CreateMultiphysicsMesh();
            postproc.PlotPatches("Patches");
            std::ofstream out("Patches.txt");
            postproc.PrintPatchInformation(out);
        }
        if(computehybrid) {
            postproc.CreateMultiphysicsMesh();
            TPZVec<REAL> errors(5,0.);
            postproc.ComputeElementErrors(errors);
            PostProcessSolution(*postproc.MultiPhysicsMesh(), iref, input, errors);
            std::cout << "Errors between H1 and compute Hybrid H1\n";
            std::cout << "H1 = " << errors[0] << " Hybrid H1 = " << errors[1] << " Estimate = " << errors[2] <<
            std::endl;
            std::cout << "H1 sq = " << errors[0]*errors[0] << " Hybrid H1 sq = " << errors[1]*errors[1] << " Estimate sq = " << errors[2]*errors[2] <<
            std::endl;
            PlotSolution("ErrorFulHybrid", *postproc.MultiPhysicsMesh(), iref, input);
        }
        
        bool estimateerror = true;
        if(estimateerror) {
            // compute a locally conservative hybrid H1 approximation
            postproc.ReconstructHybridH1();
            //        TPZVec<REAL> elementerrors;
            postproc.CreateMultiphysicsMesh();
            TPZVec<REAL> errors(5,0.);
            postproc.ComputeElementErrors(errors);
            PostProcessSolution(*postproc.MultiPhysicsMesh(), iref, input, errors);
            //        CompareMeshes(cmeshH1, cmeshHybrid, errors);
            std::cout << "Errors between H1 and Hybrid H1 reconstructed locally\n";
            std::cout << "H1 = " << errors[0] << " Hybrid H1 = " << errors[1] << " Estimate = " << errors[2] <<
            std::endl;
            std::cout << "H1 sq = " << errors[0]*errors[0] << " Hybrid H1 sq = " << errors[1]*errors[1] << " Estimate sq = " << errors[2]*errors[2] <<
            std::endl;
            PlotSolution("ReconstructHybrid", *postproc.MultiPhysicsMesh(), iref, input);
        }
        
        /// Set the Lagrange multipliers
        /// The connects of SBFemGroups level 0, one connect level 2
        /// The connects of fluxes and boundary elements connect level 1
        /// Group sbfem elements with skeleton and interface
        /// Condense the elements keeping one connect external
        /// Compute a global stiffness matrix
        UnwrapMesh(*cmeshH1ptr);
        delete cmeshH1ptr;
        UnwrapMesh(*cmeshHybridptr);
        delete cmeshHybridptr;
    }
    return 0;
}

TPZGeoMesh *ReadGmshFile(json &input) {

    std::string filename;
    if(input.find("geometry") != input.end()){
		filename = input["geometry"];
    } else {
        DebugStop();
    }
    TPZGmshReader gmsh;
    int skeletonmatid = input["skeletonMatid"];
        std::set<int> boundmatids;
    for (const auto &item : input["boundary_conditions"]) {
        int matid = item["matid"];
        std::string typeName = item["name"];
        gmsh.GetDimNamePhysical()[1][typeName] = matid;
    }
    int pointmatid = 0;
    if(input.find("point_singularity") != input.end()) {
        pointmatid = input["point_singularity"];
        gmsh.GetDimNamePhysical()[0]["point_singularity"] = pointmatid;
    }
    gmsh.GetDimNamePhysical()[2]["domain"] = 1;
    gmsh.SetCreateRefPatterns(false);
    
#ifdef MACOSX
    filename = "../" + filename;
#endif

    auto gmesh = gmsh.GeometricGmshMesh(filename);
    gmesh->BuildConnectivity();
    

    return gmesh;
}

/// @brief Refine the geometric mesh
/// @param geometric mesh object
/// @param input json data structure
/// @param number of refinements
void UniformRefine(TPZGeoMesh *gmesh, json &input, int iref) {
    int pointmatid = 0;
    if(input.find("point_singularity") != input.end()) {
        pointmatid = input["point_singularity"];
    }
    std::set<int64_t> norefine;
    if(pointmatid != 0) {
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            int ncorner = gel->NCornerNodes();
            for(int in=0; in<ncorner; in++) {
                TPZGeoElSide gelside(gel,in);
                if(gelside.HasNeighbour(pointmatid)) {
                    norefine.insert(el);
                }
            }
        }
    }
    int uniform = iref;
    int64_t nel = gmesh->NElements();
    for(int el = 0; el<nel; el++) {
        if(norefine.find(el) != norefine.end()) continue;
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->Dimension() == 0) continue;
        std::list<TPZGeoEl *> refset = {gel};
        for(int iref = 0; iref<uniform; iref++) {
            std::list<TPZGeoEl *> subrefset;
            for(auto sub : refset) {
                TPZManVector<TPZGeoEl *> subels;
                sub->Divide(subels);
                if(subels.size()) {
                    subrefset.insert(subrefset.end(),&subels[0],&subels[0]+subels.size());
                }
            }
            refset = subrefset;
        }
    }

}

void ConfigureSBFemBuilder (TPZBuildSBFemHybrid &builder, json &input) {
    ConfigureSBFemBuilder(static_cast<TPZBuildSBFem &>(builder), input);
    int porder = input["porder_hybrid_bubblefunctions"];
    builder.SetPOrderBubbleFunctions(porder);

    int skeletonporder = input["skeleton porder_hybrid"];
    builder.SetSkeletonPOrder(skeletonporder);

    int fluxmatid = input["flux matid"];
    builder.SetFluxMaterialId(fluxmatid);

    std::pair<int,int> interfacematids;
    if (input.find("interface matids") != input.end()) {
        interfacematids.first = input["interface matids"][0];
        interfacematids.second = input["interface matids"][1];
    } else {
        DebugStop();
    }
    
    builder.SetInterfaceMaterialIds(interfacematids.first, interfacematids.second);
    int lagrangeporder = input["lagrange porder"];
    builder.SetLagrangePOrder(lagrangeporder);

}

/// @brief configure the SBFem builder according to the json input
/// @param builder reference to the SBFem builder
/// @param input json object with the input data
void ConfigureSBFemBuilder (TPZBuildSBFem &builder, json &input) {
    int skeletonmatid = input["skeleton matid"];
    builder.SetSkeletonMatid(skeletonmatid);

    std::pair<int,int> interfacematids;
    if (input.find("interface matids") != input.end()) {
        interfacematids.first = input["interface matids"][0];
        interfacematids.second = input["interface matids"][1];
    } else {
        DebugStop();
    }
    

    std::map<int,int> matidtranslation;
    std::set<int> sbfemids;
    for (const auto &item : input["matid_translation"]) {
        int from = item["from"];
        int to = item["to"];
        matidtranslation[from] = to;
        sbfemids.insert(to);
    }
    builder.SetMatIdTranslation(matidtranslation);
    std::set<int> boundmatids;
    for (const auto &item : input["boundary_conditions"]) {
        int matid = item["matid"];
        if(sbfemids.find(matid) != sbfemids.end()) continue;
        boundmatids.insert(matid);
    }
    builder.SetBoundaryMatIds(boundmatids);
    int porder = input["porder_bubblefunctions"];
    builder.SetPOrderBubbleFunctions(porder);
    int skeletonporder = input["skeleton porder_h1"];
    builder.SetSkeletonPOrder(skeletonporder);
    int nref = input["nref_skeleton"];
    builder.SetNRefSkeleton(nref);

}



void RefineRandomElement(TPZGeoMesh *gmesh, int nref) {
    int iref = 0;
    int dim = gmesh->Dimension();
    while(iref < nref) {
        int64_t nel = gmesh->NElements();
        if (nel == 0) break;
        int64_t el = rand() % nel;
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || gel->Dimension() != dim || gel->HasSubElement()) continue;
        TPZManVector<TPZGeoEl *> subels;
        gel->Divide(subels);
        int firstside = gel->FirstSide(dim-1);
        int lastside = gel->FirstSide(dim);
        for (int is = firstside; is < lastside; is++) {
            TPZGeoElSide thisside(gel, is);
            for(TPZGeoElSide neighbour = thisside.Neighbour(); neighbour != thisside; neighbour = neighbour.Neighbour()) {
                TPZGeoEl *neigh = neighbour.Element();
                if(neigh->Dimension() != dim-1) continue;
                if(!neigh->HasSubElement()) {
                    TPZManVector<TPZGeoEl *> subels;
                    neigh->Divide(subels);
                }
            }
        }
        iref++;
    }
}

void InsertMaterialObjectsHybridH1(TPZCompMesh &cmesh, json &input)
{
    int dim = cmesh.Dimension();
    TPZMaterialT<STATE> *mat = nullptr;
    TPZDarcyFlow *darcy = nullptr;
    TPZElasticity2D *matelast = nullptr;
    int nstate = 0;
    std::string problemtype = input["problem_type"];

    for (const auto &item : input["materials"]) {
        int matid = item["matid"];
        if(problemtype == "Scalar") {
            nstate = 1;
            REAL perm = item["permeability"];
            darcy = new TPZHybridDarcyFlow(matid, dim);
            //        auto darcy = new TPZDarcyFlow(matid, dim);
            darcy->SetConstantPermeability(perm);
            if(item.find("exactsolution") != item.end()) {
                std::string name = item["exactsolution"];
                AddExact(problemtype,name);
                darcy->SetForcingFunction(gExact1[name].ForceFunc(),5);
                darcy->SetExactSol(gExact1[name].ExactSolution(),5);
            }
            cmesh.InsertMaterialObject(darcy);
            mat = darcy;
        } else if (problemtype == "Elastic") {
            nstate = 2;
            REAL young = item["young"];
            REAL poisson = item["poisson"];
            REAL fx(0.), fy(0.);
            int planestress = 1;
            matelast = new TPZHybridElasticity2D(matid,young,poisson,fx,fy,planestress);
            if(item.find("exactsolution") != item.end()) {
                std::string name = item["exactsolution"];
                AddExact(problemtype,name);
                matelast->SetForcingFunction(gExact2[name].ForceFunc(),5);
                matelast->SetExactSol(gExact2[name].ExactSolution(),5);
                planestress = gExact2[name].fPlaneStress;
                if(planestress) matelast->SetPlaneStress();
                else matelast->SetPlaneStrain();
                gExact2[name].gE = young;
                gExact2[name].gPoisson = poisson;
                
            }
            cmesh.InsertMaterialObject(matelast);
            mat = matelast;
            
        }
    }
    for (const auto &item : input["boundary_conditions"]) {
        int matid = item["matid"];
        std::string typeName = item["name"];
        std::vector<double> valvec = item["value"];
        if(valvec.size() != nstate) DebugStop();

        int type = item["type"];
        TPZFNMatrix<4,STATE> val1(nstate,nstate,0.);
        TPZManVector<STATE,3> val2(nstate,0.);
        for(int i = 0; i<nstate; i++) val2[i] = valvec[i];
        if (type == 0) {
        } else if (type == 1) {
        } else if (type == 2) {
            DebugStop();
        }
        TPZBndCondT<STATE> *bc;
        if(darcy) {
            bc = darcy->TPZDarcyFlow::CreateBC(mat, matid, type, val1, val2);
        } else if(matelast) {
            bc = matelast->TPZElasticity2D::CreateBC(mat, matid, type, val1, val2);
        }
        if(item.find("exactsolution") != item.end()) {
            std::string name = item["exactsolution"];
            AddExact(problemtype,name);
            if(nstate == 1) {
                bc->SetForcingFunctionBC(gExact1[name].ExactSolution(),5);
            } else if(nstate == 2) {
                bc->SetForcingFunctionBC(gExact2[name].ExactSolution(),5);
            }
        }
        cmesh.InsertMaterialObject(bc);
    }

    int skeletonMatId = input["skeleton matid"];
    TPZNullMaterial<STATE> *skelmat = new TPZNullMaterial<STATE>(skeletonMatId, dim, nstate);
    cmesh.InsertMaterialObject(skelmat);
    int fluxmaterialid = input["flux matid"];
    TPZNullMaterial<STATE> *fluxmat = new TPZNullMaterial<STATE>(fluxmaterialid, dim, nstate);
    cmesh.InsertMaterialObject(fluxmat);
}

void InsertMaterialObjectsH1(TPZCompMesh &cmesh, json &input)
{
    int dim = cmesh.Dimension();
    TPZMaterialT<STATE> *mat = nullptr;
    int nstate = 0;
    std::string problemtype = input["problem_type"];

    for (const auto &item : input["materials"]) {
        int matid = item["matid"];
        if(problemtype == "Scalar") {
            nstate = 1;
            REAL perm = item["permeability"];
            auto darcy = new TPZDarcyFlow(matid, dim);
            //        auto darcy = new TPZDarcyFlow(matid, dim);
            darcy->SetConstantPermeability(perm);
            if(item.find("exactsolution") != item.end()) {
                std::string name = item["exactsolution"];
                AddExact(problemtype,name);
                darcy->SetForcingFunction(gExact1[name].ForceFunc(),5);
                darcy->SetExactSol(gExact1[name].ExactSolution(),5);
            }
            cmesh.InsertMaterialObject(darcy);
            if(!mat) mat = darcy;
        } else if (problemtype == "Elastic") {
            nstate = 2;
            REAL young = item["young"];
            REAL poisson = item["poisson"];
            REAL fx(0.), fy(0.);
            int planestress = 1;
            TPZElasticity2D *matelast = new TPZElasticity2D(matid,young,poisson,fx,fy,planestress);
            if(item.find("exactsolution") != item.end()) {
                std::string name = item["exactsolution"];
                AddExact(problemtype,name);
                matelast->SetForcingFunction(gExact2[name].ForceFunc(),5);
                matelast->SetExactSol(gExact2[name].ExactSolution(),5);
                planestress = gExact2[name].fPlaneStress;
                if(planestress) matelast->SetPlaneStress();
                else matelast->SetPlaneStrain();
                gExact2[name].gE = young;
                gExact2[name].gPoisson = poisson;
            }
            if(!mat) mat = matelast;
            cmesh.InsertMaterialObject(matelast);
            nstate = 2;
            
        }
    }
    for (const auto &item : input["boundary_conditions"]) {
        int matid = item["matid"];
        std::string typeName = item["name"];
        std::vector<double> valvec = item["value"];
        
        if(valvec.size() != nstate) DebugStop();

        int type = item["type"];
        TPZFNMatrix<4,STATE> val1(nstate,nstate,0.);
        TPZManVector<STATE,3> val2(nstate,0.);
        for(int i = 0; i<nstate; i++) val2[i] = valvec[i];
        if (type == 0) {
        } else if (type == 1) {
        } else if (type == 2) {
            DebugStop();
        }
        TPZBndCondT<STATE> *bc = mat->CreateBC(mat, matid, type, val1, val2);
        if(item.find("exactsolution") != item.end()) {
            std::string name = item["exactsolution"];
            AddExact(problemtype,name);
            if(nstate == 1) {
                bc->SetForcingFunctionBC(gExact1[name].ExactSolution(),5);
            } else if(nstate == 2) {
                bc->SetForcingFunctionBC(gExact2[name].ExactSolution(),5);
            }
        }
        cmesh.InsertMaterialObject(bc);
    }

    int skeletonMatId = input["skeleton matid"];
    TPZNullMaterial<STATE> *skelmat = new TPZNullMaterial<STATE>(skeletonMatId, dim,nstate);
    cmesh.InsertMaterialObject(skelmat);
}

void InsertInterfaceMaterialObjects(TPZCompMesh &cmesh, json &input) {
    int dim = cmesh.Dimension();
    std::pair<int,int> interfacematids;
    if (input.find("interface matids") != input.end()) {
        interfacematids.first = input["interface matids"][0];
        interfacematids.second = input["interface matids"][1];
    } else {
        DebugStop();
    }
    std::string problemtype = input["problem_type"];
    int nstate = 0;
    if(problemtype == "Scalar") {
        nstate = 1;
    } else if (problemtype == "Elastic") {
        nstate = 2;
    } else {
        DebugStop();
    }
    TPZLagrangeMultiplier<STATE> *lagrange = new TPZLagrangeMultiplier<STATE>(interfacematids.first, dim-1, nstate);
    cmesh.InsertMaterialObject(lagrange);
    lagrange = new TPZLagrangeMultiplier<STATE>(interfacematids.second, dim-1, nstate);
    lagrange->SetMultiplier(-1.);
    cmesh.InsertMaterialObject(lagrange);
}

#include "pzcondensedcompel.h"
#include "pzelementgroup.h"
void UnwrapMesh(TPZCompMesh &cmesh) {
    int64_t nel = cmesh.NElements();
    int npass = 3;
    for (int pass = 0; pass < npass; pass++) {
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = cmesh.Element(el);
            if (!cel) continue;
            TPZCondensedCompEl *condensed = dynamic_cast<TPZCondensedCompEl *>(cel);
            if (condensed) {
                condensed->Unwrap();
                continue;
            }
            TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
            TPZSBFemElementGroup *sbfemgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
            if (elgr && !sbfemgr) {
                elgr->Unwrap(false);
            }
        }
    }
}

#include "pzfstrmatrix.h"
#include "TPZLinearAnalysis.h"
#include "pzstepsolver.h"
#include "TPZVTKGenerator.h"
void SolveSystem(json &input, TPZCompMesh &cmesh) {
    TPZFStructMatrix<STATE> strmat(&cmesh);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    // Create the linear system
    TPZLinearAnalysis an(&cmesh);
    an.SetStructuralMatrix(strmat);
    an.SetSolver(step);

    an.Assemble();
    // Solve the system
    an.Solve();

}

/// @brief Create a computational mesh based on the SBFem hybrid H1 space
/// @param gmesh autopointer to the geometric mesh
/// @param input json object with the input data
/// @return pointer to the computational mesh
TPZCompMesh *CreateHybridSBFemSpaceFromSBFem(TPZBuildSBFemHybrid &builder, json &input) {

    TPZAutoPointer<TPZGeoMesh> gmesh = builder.GetGMesh();
    ConfigureSBFemBuilder(builder, input);
    if(0)
    {
        std::ofstream out("gmesh.txt");
        gmesh->Print(out);
        std::ofstream out2("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
    }
    builder.DuplicateSkeletonElements();

    TPZCompMesh *cmeshptr = new TPZCompMesh(gmesh);
    TPZCompMesh &cmesh = *cmeshptr;
    InsertMaterialObjectsHybridH1(cmesh, input);
    builder.CreateSkeletonApproximationSpace(cmesh);
    // this will create the SBFem group elements
    builder.CreateVolumetricElements(cmesh);
    InsertInterfaceMaterialObjects(cmesh, input);
    builder.CreateInterfaceElements(cmesh);
    builder.InitializeLagrangeLevels(cmesh);
    // this method groups the sbfem element group with the skeleton and interfaces
    builder.GroupAndCondenseElements(cmesh);
    {
        std::ofstream out("SBFemConfig.txt");
        builder.Print(out);
    }
    return cmeshptr;
}

TPZCompMesh *CreateHybridSBFemSpace(TPZBuildSBFemHybrid &builder, json &input){
    ConfigureSBFemBuilder(builder, input);
    TPZAutoPointer<TPZGeoMesh> gmesh = builder.GetGMesh();
    builder.StandardConfiguration();
    // builder.SetNRefSkeleton(0);
    int nref = builder.GetNRefSkeleton();
    builder.DivideSkeleton(nref);
    auto matidtranslation = builder.GetMatIdTranslation();
    std::set<int> matids;
    for(auto it : matidtranslation) {
        matids.insert(it.first);
    }
    builder.GenerateCollapsedGeometricElements(matids);
    builder.DuplicateSkeletonElements();

    TPZCompMesh *cmeshptr = new TPZCompMesh(gmesh);
    TPZCompMesh &cmesh = *cmeshptr;
    InsertMaterialObjectsHybridH1(cmesh, input);
    builder.CreateSkeletonApproximationSpace(cmesh);
    builder.CreateVolumetricElements(cmesh);
    InsertInterfaceMaterialObjects(cmesh, input);
    builder.CreateInterfaceElements(cmesh);
    builder.InitializeLagrangeLevels(cmesh);
    builder.GroupAndCondenseElements(cmesh);
    {
        std::ofstream out("SBFemConfig.txt");
        builder.Print(out);
    }
    return cmeshptr;
}

/// @brief Create a computational mesh based on the SBFem hybrid H1 space
/// @param gmesh autopointer to the geometric mesh
/// @param input json object with the input data
/// @return pointer to the computational mesh
TPZCompMesh *CreateSBFemSpace(TPZBuildSBFem &builder, json &input) {
    TPZAutoPointer<TPZGeoMesh> gmesh = builder.GetGMesh();
    int dim = gmesh->Dimension();
    ConfigureSBFemBuilder(builder, input);
    int64_t singular_partition = -1;
    int64_t singular_element = -1;
    if(input.find("point_singularity") != input.end()) {
        // two steps :
        // 1 create standard elements for all elements that do not touch the singularity
        int singular_matid = input["point_singularity"];
        TPZStack<int64_t> volels;
        std::map<int,int> matidtranslation = builder.GetMatIdTranslation();
        
        int64_t singular_node = -1;
        TPZGeoEl *singular_element = 0;
        int64_t nel = gmesh->NElements();
        std::set<int64_t> singular_group;
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel || gel->MaterialId() != singular_matid) continue;
            singular_node = gel->NodeIndex(0);
            singular_element = gel;
            break;
        }
        if(singular_node == -1) DebugStop();
        TPZGeoElSide singularside (singular_element);
        std::set<int64_t> to_invert;
        std::set<int64_t> to_delete;
        for (auto gelside = singularside.Neighbour(); gelside != singularside; gelside++) {
            TPZGeoEl *gel = gelside.Element();
            if(!gel || gel->HasSubElement()) continue;
            if(gel->Dimension() == dim) {
                singular_group.insert(gel->Index());
            } else if(gel->Dimension() == dim-1) {
                // boundary element touching the singularity
                to_delete.insert(gel->Index());
            }
        }
        for(auto it : to_delete) {
            TPZGeoEl *gel = gmesh->Element(it);
            gel->RemoveConnectivities();
            delete gel;
        }

        nel = gmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel || gel->Dimension() != dim) continue;
            if(singular_group.find(el) != singular_group.end()) continue;
            volels.Push(el);
        }
        builder.CreateElementCenterNodes(volels);
        singular_partition = builder.AddPartition(singular_group, singular_node);
        builder.AddSkeletonElements();
        // 2 create collapsed elements towards the singular point
    } else {
        builder.StandardConfiguration();
    }
    // builder.SetNRefSkeleton(0);
    int nref = builder.GetNRefSkeleton();
    builder.DivideSkeleton(nref);
    // auto matidtranslation = builder.GetMatIdTranslation();
    // std::set<int> matids;
    // for(auto it : matidtranslation) {
    //     matids.insert(it.first);
    // }
    // builder.GenerateCollapsedGeometricElements(matids);
    TPZCompMesh *cmeshptr = new TPZCompMesh(gmesh);
    TPZCompMesh &cmesh = *cmeshptr;
    InsertMaterialObjectsH1(cmesh, input);
    builder.BuildComputationMesh(cmesh);
    
    // std::set<int> matidsskel = builder.GetBoundaryMatIds();
    // matidsskel.insert(builder.GetSkeletonMatid());
    // cmesh.SetDefaultOrder(builder.GetSkeletonPOrder());
    // cmesh.SetAllCreateFunctionsContinuous();
    // cmesh.AutoBuild(matidsskel);
    // builder.CreateVolumetricElements(cmesh);
    // builder.CreateElementGroups(cmesh);

    return cmeshptr;
}


/// @brief Plot the solution of the computational mesh
/// @param cmesh reference to the computational mesh
/// @param output json object with the output data
void PlotSolution(const std::string &rootname, TPZCompMesh &cmesh, int refstep, json &output) {
    std::set<int> matiderror;
    if (output.find("matid_translation") != output.end()) {
        for (const auto &item : output["matid_translation"]) {
            int matid = item["to"];
            matiderror.insert(matid);
        }
    } else {
        DebugStop();
    }    // make a directory based on the exact solution
    
    std::string dirname = "PostProcess";
    if(output.find("directory") != output.end()) {
        std::string exactname = output["directory"];
        dirname = exactname;
    }
    int matid = *matiderror.rbegin();
    TPZMaterial *mat = cmesh.FindMaterial(matid);
    if (!mat) {
        std::cout << "Material id " << matid << " not found in the computational mesh\n";
        DebugStop();
    }
    TPZDarcyFlow *darcy = dynamic_cast<TPZDarcyFlow *>(mat);
    TPZElasticity2D *elast = dynamic_cast<TPZElasticity2D *>(mat);
    if (!darcy && !elast) {
        std::cout << "Material id " << matid << " is neither a Darcy flow or elasticity material\n";
        DebugStop();
    }
    std::system((std::string("mkdir -p ") + dirname).c_str());
    // std::string plotfile = dirname + std::string("/Solution.vtk");
    std::string plotfile = dirname + "/" + rootname;

    int porder = output["porder_bubblefunctions"];
    int skeletonporder = output["skeleton porder_h1"];
    int lagrangeporder = output["lagrange porder"];
    int nref_skeleton = output["nref_skeleton"];

    plotfile += "_r" + std::to_string(refstep) + "_rS" + std::to_string(nref_skeleton) + "_pB" + std::to_string(porder) + "_pS" + std::to_string(skeletonporder) + "_pL" + std::to_string(lagrangeporder);

    TPZStack<std::string> postprocess;
    if (output.find("PostProcess") != output.end()) {
        for (const auto &item : output["PostProcess"]) {
            std::string field = item;
            
            if(darcy && darcy->VariableIndex(field) != -1) {
                postprocess.Push(field);
            }
            if(elast && elast->VariableIndex(field) != -1) {
                postprocess.Push(field);
            }
        }
    }
    int dim = 2;
    int res = 2;
    std::cout << "Post processing in file " << plotfile << std::endl;
    TPZVTKGenerator vtk(&cmesh, postprocess, plotfile, res, dim);
    vtk.SetStep(refstep);
    vtk.Do();
}

/// @brief Post-process the solution of the computational mesh (compute the approximation error)
/// @param cmesh reference to the computational mesh
/// @param input json object with the output data
void PostProcessSolution(TPZCompMesh &cmesh, int iref, json &input, TPZVec<REAL> &errors) {
    bool store_errors = true;
    std::set<int> matiderror;
    if (input.find("matid_translation") != input.end()) {
        for (const auto &item : input["matid_translation"]) {
            int matid = item["to"];
            matiderror.insert(matid);
        }
    } else {
        DebugStop();
    }
    int matid = *matiderror.rbegin();
    TPZMaterial *mat = cmesh.FindMaterial(matid);
    if (!mat) {
        std::cout << "Material id " << matid << " not found in the computational mesh\n";
        DebugStop();
    }
    TPZManVector<std::string> errornames;
    std::string problemtype = input["problem_type"];
    if(problemtype == "Scalar") {
        TPZDarcyFlow *darcy = dynamic_cast<TPZDarcyFlow *>(mat);
        if (!darcy) {
            std::cout << "Material id " << matid << " is not a Darcy flow material\n";
            DebugStop();
        }
        int nerrors = darcy->NEvalErrors();
        errornames.resize(nerrors);
        darcy->ErrorNames(errornames);
    } else if(problemtype == "Elastic") {
        TPZElasticity2D *elast = dynamic_cast<TPZElasticity2D *>(mat);
        if (!elast) {
            std::cout << "Material id " << matid << " is not an elasticity material\n";
            DebugStop();
        }
        int nerrors = elast->NEvalErrors();
        errornames.resize(nerrors);
        elast->ErrorNames(errornames);

    }

    // make a directory based on the exact solution
    std::string dirname = "PostProcess";
    
    if(input.find("directory") != input.end()) {
        std::string exactname = input["directory"];
        dirname = exactname;
    }
    std::system((std::string("mkdir -p ") + dirname).c_str());
    
    std::string errorfilename = dirname + "/Errors.txt";
    std::ofstream errorfile(errorfilename, std::ios::app);

    int neq_condensed = cmesh.NEquations();
    int neq_total = cmesh.Solution().Rows();
    int porder = input["porder_bubblefunctions"];
    int skeletonporderH1 = input["skeleton porder_h1"];
    int skeletonporderHybrid = input["skeleton porder_hybrid"];
    int lagrangeporder = input["lagrange porder"];
    int nref_skeleton = input["nref_skeleton"];
    int nref_uniform = iref;
    errorfile << "Material_name " << mat->Name() << " ";
    for (int i = 0; i < 3; i++) {
        errorfile << errornames[i] << " " << errors[i] << " ";
    }
    errorfile << "nref " << nref_uniform << " ";
    errorfile << "neq_condensed " << neq_condensed << " neq_total " << neq_total << " ";
    errorfile << " nref_skeleton " << nref_skeleton << " bubble_order " << porder << " skeleton_order_h1 " << skeletonporderH1 << " skeleton_order_hybrid " << skeletonporderHybrid << " lagrange_order " << lagrangeporder << " ";
    errorfile << std::endl;
}

#include "pzcheckgeom.h"
/// @brief Perform a convergence study applying uniform refinements
/// @param gmesh autopointer to the geometric mesh
/// @param input json object with the input data
void ConvergenceStudy(TPZAutoPointer<TPZGeoMesh> gmeshin, json &input, bool append) {

    // make a directory based on the exact solution
    std::string dirname = "PostProcess";
    if(input.find("directory") != input.end()) {
        std::string exactname = input["directory"];
        dirname = exactname;
    }
    std::system((std::string("mkdir -p ") + dirname).c_str());

    if (!append)
    {
        std::string errorfilename = dirname + "/Errors.txt";
        std::ofstream errorfile(errorfilename);
    }

    int matid = input["matid_translation"][0]["to"];
    std::cout << "Material id for error computation " << matid << std::endl;

    std::vector<int> refsteps = input["UniformRefinement"];
    if(refsteps.size() != 2) DebugStop();
    for (int iref = refsteps[0]; iref <= refsteps[1]; iref++) {
        TPZGeoMesh *gmeshcopy = new TPZGeoMesh(*gmeshin);
        TPZCheckGeom check(gmeshcopy);
        check.UniformRefine(iref);
        TPZAutoPointer<TPZGeoMesh> gmesh(gmeshcopy);
        TPZBuildSBFemHybrid builder(gmesh);
        TPZCompMesh *cmeshptr = CreateHybridSBFemSpace(builder, input);
        TPZCompMesh &cmesh = *cmeshptr;
        SolveSystem(input, cmesh);
        // Plot the solution
        std::cout << "Plotting the solution\n";
        PlotSolution("hybrid", cmesh, iref, input);
        std::cout << "Plotting the solution done\n";
        // Post-process the solution
        std::cout << "Post-processing the solution\n";
        TPZManVector<REAL> errors;
        PostProcessSolution(cmesh, iref, input,errors);
        std::cout << "Post-processing the solution done\n";
        UnwrapMesh(cmesh);
        delete cmeshptr;
    }
}
#include "TPZSBFemVolume.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"
void AddSBFemElement(TPZCompEl *cel, std::map<int64_t, TPZSBFemVolume *> &sbfem) {
    if (!cel) return;
    TPZSBFemVolume *sbfemvol = dynamic_cast<TPZSBFemVolume *>(cel);
    if (sbfemvol) {
        int64_t index = sbfemvol->Reference()->Index();
        sbfem[index] = sbfemvol;
    }
    TPZCondensedCompEl *condensed = dynamic_cast<TPZCondensedCompEl *>(cel);
    if (condensed) {
        TPZCompEl *father = condensed->ReferenceCompEl();
        AddSBFemElement(father, sbfem);
    }
    TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
    if (elgr) {
        const TPZVec<TPZCompEl *> &elvec = elgr->GetElGroup();
        int nels = elvec.NElements();
        for (int iel = 0; iel < nels; iel++) {
            TPZCompEl *subel = elvec[iel];
            AddSBFemElement(subel, sbfem);
        }
    }

}

/// @brief Compute the difference between two meshes
/// @param cmesh1 reference to the first computational mesh
/// @param cmesh2 reference to the second computational mesh
void CompareMeshes(TPZCompMesh &cmesh1, TPZCompMesh &cmesh2, TPZVec<REAL> &errors)
{
    errors.Fill(0.);
    std::map<int64_t, TPZSBFemVolume *> sbfem1, sbfem2;
    int64_t nel1 = cmesh1.NElements();
    for (int64_t el = 0; el < nel1; el++) {
        TPZCompEl *cel = cmesh1.Element(el);
        AddSBFemElement(cel, sbfem1);
    }
    int64_t nel2 = cmesh2.NElements();
    for (int64_t el = 0; el < nel2; el++) {
        TPZCompEl *cel = cmesh2.Element(el);
        AddSBFemElement(cel, sbfem2);
    }
    std::map<int,TPZH1ErrorHybridH1EstimateMaterial *> matmap;
    for(auto itmesh1 : sbfem1) {
        int64_t index = itmesh1.first;
        auto itmesh2 = sbfem2.find(index);
        if (itmesh2 == sbfem2.end()) {
            std::cout << "Element index " << index << " not found in mesh 2\n";
            continue;
        }
        TPZSBFemVolume *sbfemvol1 = itmesh1.second;
        TPZSBFemVolume *sbfemvol2 = itmesh2->second;
        TPZMaterial *mat1 = sbfemvol1->Material();
        TPZMaterial *mat2 = sbfemvol2->Material();
        if (!mat1 || !mat2) {
            std::cout << "Element index " << index << " has no material\n";
            DebugStop();
        }
        if (mat1->Id() != mat2->Id()) {
            std::cout << "Element index " << index << " has different material ids " << mat1->Id() << " and " << mat2->Id() << std::endl;
            DebugStop();
        }
        TPZDarcyFlow *Darcy = dynamic_cast<TPZDarcyFlow *>(mat1);
        if (!Darcy) {
            std::cout << "Element index " << index << " material id " << mat1->Id() << " is not a DarcyFlow material\n";
            DebugStop();
        }
        if(matmap.find(mat1->Id()) == matmap.end()) {
            matmap[mat1->Id()] = new TPZH1ErrorHybridH1EstimateMaterial(*Darcy);
        }
        TPZH1ErrorHybridH1EstimateMaterial *material = matmap[mat1->Id()];
        TPZManVector<TPZMaterialDataT<STATE>,4> datavec(4);
        int dim = cmesh1.Dimension();
        TPZGeoEl *ref1 = sbfemvol1->Reference();
        TPZAutoPointer<TPZIntPoints> intrule = ref1->CreateSideIntegrationRule(ref1->NSides() - 1, 5);
        int maxIntOrder = intrule->GetMaxOrder();
        TPZManVector<int,3> intorders(dim,maxIntOrder);
        intrule->SetOrder(intorders);
        int problemdimension = cmesh1.Dimension();
        int NErrors = material->NEvalErrors();
        if(NErrors != errors.NElements()) {
            std::cout << "Number of error components " << NErrors << " different from allocated space " << errors.NElements() << std::endl;
            DebugStop();
        }
        TPZManVector<REAL> elerrors(NErrors,0.);
        int ndof = material->NStateVariables();
        TPZMaterialDataT<STATE> &data0 = datavec[0];
        data0.x.Resize(3);
        data0.sol.Resize(1);
        data0.sol[0].Resize(ndof);
        int nintpoints = intrule->NPoints();

        REAL weight;
        TPZManVector<REAL,3> intpoint(problemdimension,0.);

        for (int nint = 0; nint < nintpoints; nint++) {

            intrule->Point(nint, intpoint, weight);
            ref1->X(intpoint, data0.x);
            ref1->Jacobian(intpoint, data0.jacobian, data0.axes, data0.detjac, data0.jacinv);
            datavec[3].fNeedsSol = true;
            sbfemvol1->ComputeRequiredData(datavec[3],intpoint);
            weight *= fabs(data0.detjac);
            datavec[1].fNeedsSol = true;
            sbfemvol2->ComputeRequiredData(datavec[1],intpoint);
            TPZManVector<REAL> locerrors(NErrors,0.);
            material->Errors(datavec, locerrors);
            for (int i = 0; i < NErrors; i++) {
                elerrors[i] += locerrors[i] * weight;
            }
        }
        for (int i = 0; i < NErrors; i++) {
            errors[i] += elerrors[i];
        }
    }
    for (int i = 0; i < errors.NElements(); i++) {
        errors[i] = sqrt(errors[i]);
    }
}

/// @brief Create collapsed elements towards the singular element
std::set<int64_t> CreateCollapsedElements(TPZGeoMesh *gmesh, int64_t singular_el, TPZBuildSBFem &builder) {
    int dim = gmesh->Dimension();
    if(dim != 2) DebugStop();
    auto matidtranslation = builder.GetMatIdTranslation();
    std::set<int64_t> elgroup;
    TPZGeoEl *sing = gmesh->Element(singular_el);
    int64_t singnodeindex = sing->NodeIndex(0);
    TPZGeoElSide singside(sing);
    for(auto neighbour = singside.Neighbour(); neighbour != singside; neighbour++) {
        TPZStack<TPZGeoElSide> highelsides;
        TPZGeoEl *neighgel = neighbour.Element();
        int neighdim = neighgel->Dimension();
        int neighmatid = neighgel->MaterialId();
        if(neighdim == 1) {
            if(neighbour.Side() == 1) {
                elgroup.insert(neighgel->Index());
            } else {
                // flip the element
                DebugStop();
            }
        } else if(neighdim == 2) {
            if(matidtranslation.find(neighmatid) == matidtranslation.end()) DebugStop();
            int sbfemmatid = matidtranslation[neighmatid];
            neighgel->AllHigherDimensionSides(neighbour.Side(), 1, highelsides);
            std::set<int> highsides;
            for(auto &it : highelsides) highsides.insert(it.Side());
            int firstside = neighgel->FirstSide(1);
            int lastside = neighgel->FirstSide(2);
            for(int side = firstside; side<lastside; side++) {
                if(highsides.find(side) != highsides.end()) continue;
                TPZManVector<int64_t> nodeindexes(4,singnodeindex);
                nodeindexes[0] = neighgel->SideNodeIndex(side, 0);
                nodeindexes[1] = neighgel->SideNodeIndex(side, 1);
                int64_t index;
                gmesh->CreateGeoElement(EQuadrilateral, nodeindexes, sbfemmatid, index);
                elgroup.insert(index);
            }
        }
    }
}
