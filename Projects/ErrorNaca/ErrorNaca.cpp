/**
 * @file This file implements an error estimator in space H1.
 */

#include "pzlog.h"
#include "pzinterpolationspace.h"
#include "pzcheckgeom.h"
#include "pzstepsolver.h"
#include <pzbuildmultiphysicsmesh.h>
#include "pzgeoelbc.h"
#include "pzskylstrmatrix.h"
#include "pzskylmat.h"
#include "TPZLinearAnalysis.h"
#include "TPZBndCondT.h"
#include <TPZNullMaterial.h>
#include <pzmultiphysicscompel.h>
#include "TPZSBFemElementGroup.h"
#include "TPZSBFemVolume.h"
#include <TPZMultiphysicsCompMesh.h>
#include "TPZSSpStructMatrix.h"
#include "TPZSYSMPMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZRefPatternDataBase.h"
#include "TPZVTKGeoMesh.h"
#include "TPZGmshReader.h"
#include "TPZAnalyticSolution.h"
#include "TPZFrontSym.h"
#include <TPZSimpleTimer.h>
#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "Projection/TPZL2ProjectionCS.h"      
#include "TPZNullMaterialCS.h"             

#include <fstream>
#include <ctime>
#include <cstdio>
#include <cmath>

#include "TPZGenGrid2D.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "TPZGeoLinear.h"

#include "tpzblendnaca.h"
#include "tpznacaprofile.h"

#include "tpzchangeel.h"

#include <iostream>

#define USING_MKL

using namespace pzgeom;
/// @brief verify is the derivative of the NACA coordinate is correct
/// @param naca profile object
/// @param point parametric coordinate around which the derivative will be verified
void VerifyDerivative(TPZBlendNACA &naca, REAL point);

/// @brief Read a gmsh file and return a geometric mesh
/// @param filename 
/// @return a geometric mesh
TPZGeoMesh *ReadGmsh(const std::string &filename, const TPZNacaProfile &naca);

/// @brief Create a computational mesh with H1 elements
TPZCompMesh *CreateH1CompMesh(TPZGeoMesh *gmesh, TPZVec<int> &porders, int64_t &newcon);

/// @brief Adjust the integration rule for Equarterpoints H1 elements
void AdjustH1Equarterpointintrule(TPZGeoMesh *gmesh);
// void AdjustH1Equarterpointintrule(TPZGeoMesh *gmesh);

/// @brief Create a computational mesh with L2 elements
TPZCompMesh *CreateL2CompMesh(TPZGeoMesh *gmesh);

/// @brief Create the computational mesh with HDiv elements
TPZCompMesh *CreateHDivCompMesh(TPZGeoMesh *gmesh, TPZVec<int> &porders, int64_t &newcon);

/// @brief Adjust the integration rule for Equarterpoints HDiv elements
void AdjustHDivEquarterpointintrule(TPZGeoMesh *gmesh);
// void AdjustHDivEquarterpointintrule(TPZGeoMesh *gmesh);

/// @brief Create the computational "multiphysics" mesh with only HDiv elements
TPZMultiphysicsCompMesh *CreateMultiphysicsMesh(TPZCompMesh *cmeshHDiv, TPZCompMesh *cmeshL2, TPZGeoMesh *gmesh);

/// @brief Simulate the NACA profile using H1 approximation and Joukowski condition to find Beta
TPZCompMesh *SimulateNacaProfileH1(TPZGeoMesh *gmesh, TPZVec<int> &porders, REAL &circulation_H1);

/// @brief Simulate the NACA profile using H1 approximation and minimization to find Beta
TPZCompMesh *SimulateNacaProfileH1_Minimization(TPZGeoMesh *gmesh, TPZVec<int> &porders, REAL &circulation_H1);

/// @brief Simulate the NACA profile using H(div) approximation and Joukowski condition to find Beta
TPZMultiphysicsCompMesh *SimulateNacaProfileHDiv(TPZGeoMesh *gmesh, TPZVec<int> &porders, REAL &circulation_HDiv);

/// @brief Simulate the NACA profile using H(div) approximation and minimization to find Beta
TPZMultiphysicsCompMesh *SimulateNacaProfileHDiv_Minimization(TPZGeoMesh *gmesh, TPZVec<int> &porders, REAL &circulation_HDiv);

/// @brief Compute the error as an average error per element and save as an elemental solution 
void ComputeErrorEstimator(TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *cmesh_m,TPZVec<REAL> &ErrorEstimator);

/// @brief Compute the global error the norm of the vector ErrorEstimator 
void ComputeGlobalError(TPZVec<REAL> &ErrorEstimator, REAL &GlobalError);

/// @brief Compute the error as an average error per element and save as an elemental solution 
void Hrefinement(TPZMultiphysicsCompMesh *cmesh_m,TPZVec<REAL> &ErrorEstimator,TPZVec<REAL> &RefinementIndicator,TPZVec<int> &porders);

/// @brief Performe h-refinements at the computational meshes "cmesh" and "cmesh_m" based on the ErrorEstimator 
// @param cmesh_m multiphysics computational mesh
// @param ErrorEstimator vector with the error estimator for each computational element
// @param minh minimum level of refinement
// @param RefinementIndicator vector with the refinement indicator for each computational element
// @param porders vector with the polynomial order of each geometric element 
void HPrefinement(TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator, int minh, TPZVec<REAL> &RefinementIndicator, TPZVec<int> &porders);

/// @brief Smoothen the volumetric elements around the trailingedge by adopting the same level of refinement  
void Smoothentrailingedgeelements(TPZMultiphysicsCompMesh *cmesh_m,TPZVec<REAL> &RefinementIndicator);

/// @brief Substitute the trailing edge quadrilateral elements with colapsed quadrilateral elements with or without quarterpoint elements 
void Changetrailingedgeelements(TPZGeoMesh *gmesh);

/// @brief If the element is not linear mapping, change it to linear mapping
/// @param gel element that will be substituted
void ChangeToLinearQuad(TPZGeoEl *gel);
void ChangeToLinearOned(TPZGeoEl *gel);

/// @brief change the element to a blend element
void ChangeToBlend(TPZGeoEl *gel);

/// @brief Adjust the geometric mesh to create blend elements only for neighbours of the profile
void AdjustGeometry(TPZGeoMesh *gmesh);

/// @brief adjust the elements neighbouring the trailing edge to accomodate SBFem simulation
void AdjustToSBFemGeometry(TPZGeoMesh *gmesh);

/// @brief remove computational elements generated for SBFem simulation
void RemoveSBFemElements(TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *cmesh_m, TPZGeoMesh *gmesh);

/// @brief print the geometry of the trailing edge elements
void PrintTrailingEdgeElements(TPZGeoMesh *gmesh);

/// @brief Change the elements the touch the trailing edge to quarterpoint elements
void CreateQuarterPointElements(TPZGeoMesh *gmesh);

/// @brief Change the elements the touch the trailing edge to H1 SBFEM elements
void CreateH1SBFEMelements(TPZCompMesh *cmesh);

/// @brief Change the elements the touch the trailing edge to Hdiv SBFEM elements
void CreateHdivSBFEMelements(TPZMultiphysicsCompMesh *cmesh, int64_t newcon);

/// @brief create a refinement patter that cuts the element in longitudinal direction
/// @return refinement pattern created
TPZAutoPointer<TPZRefPattern> CreateRefPattern();

/// @brief divide a geometric element. If it is a trailing edge element, set the refinement pattern first
void DivideGeoEl(TPZGeoEl *gel, TPZVec<TPZGeoEl *> &subels);

/// @brief returns true if the element is counter clockwise
bool IsCounterClockwise(TPZGeoEl *gel);

/// @brief make the elements counter clockwise
void MakeCounterClockwise(TPZGeoEl *gel);

/// @brief Adjust the colapsed connects of the atomic computational meshes   
void Adjustcompmesh(TPZCompMesh *cmesh);

/// @brief Divide Trailing Edge Neighbours   
void DivideTrailingEdgeNeighbours(TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator, REAL tol, TPZVec<REAL> &RefinementIndicator, TPZVec<int> &porders);

/// @brief Divide Profile Neighbours   
void DivideProfileNeighbours(TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator, REAL tol, int minlevel, TPZVec<REAL> &RefinementIndicator);

/// @brief SmoothenGeometry   
void SmoothenGeometry(TPZGeoMesh *gmesh);

/// @brief Check the NACA profile, generate coordinates and write to a file
void CheckNacaProfile(const TPZNacaProfile &naca);

/// @brief print the results of the analysis
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh);

/// @brief Adjust...
void AdjustSystemofEquations (TPZLinearAnalysis &an, TPZCompMesh *cmesh, int64_t &newcon);
// void AdjustSystemofEquations (TPZLinearAnalysis &an, TPZCompMesh *cmesh, int64_t &newcon);

/// @brief Get the trailing edge elements
void GetTrailingEdgeElements(TPZGeoMesh *gmesh, TPZGeoElSide &Geosideminus, TPZGeoElSide &Geosideplus);
// void GetTrailingEdgeElements(TPZGeoMesh *gmesh, TPZGeoElSide Geosideplus, TPZGeoElSide Geosideminus);

/// @brief Evaluate...
void EvaluateSolutionGradientsH1(TPZGeoMesh *gmesh, TPZManVector<REAL,3> &gradplus, TPZManVector<REAL,3> &gradminus);
// void EvaluateSolutionGradientsH1(TPZGeoMesh *gmesh, TPZManVector<REAL,3> &gradplus, TPZManVector<REAL,3> &gradminus);

/// @brief Evaluate...
void EvaluateSolutionHDiv(TPZGeoMesh *gmesh, TPZMultiphysicsCompMesh *cmesh_m, TPZManVector<REAL,3> &u_plus, TPZManVector<REAL,3> &u_minus);
// void EvaluateSolutionHDiv(TPZGeoMesh *gmesh, TPZManVector<REAL,3> &u_plus, TPZManVector<REAL,3> &u_minus);

/// @brief Evaluate...
void EvaluateCirculationH1(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, int matid, REAL &circulation);
// void EvaluateSolutionHDiv(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, int matid, REAL &circulation);

/// @brief Evaluate...
void EvaluateCirculationHDiv(TPZGeoMesh *gmesh, TPZMultiphysicsCompMesh *cmesh_m, int matid, REAL &circulation);
// void EvaluateSolutionHDiv(TPZGeoMesh *gmesh, TPZMultiphysicsCompMesh *cmesh_m, int matid, REAL &circulation);

/// @brief Compute...
void ComputeBetaH1(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, TPZFMatrix<STATE> &phi_0, TPZFMatrix<STATE> &phi_1, REAL &Beta);
// void ComputeBetaH1(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, TPZFMatrix<STATE> &phi_0, TPZFMatrix<STATE> &phi_1, REAL &Beta);

/// @brief Compute...
void ComputeBetaHDiv(TPZGeoMesh *gmesh, TPZMultiphysicsCompMesh *cmesh_m, TPZFMatrix<STATE> &u_0, TPZFMatrix<STATE> &u_1, REAL &Beta);
// void ComputeBetaHDiv(TPZGeoMesh *gmesh, TPZMultiphysicsCompMesh *cmesh_m, TPZFMatrix<STATE> &u_0, TPZFMatrix<STATE> &u_1, REAL &Beta);

/// Add the SBFemVolume elements to the computational mesh
void AddSBFemVolumeElements();

/// Hide the SBFemVolume elements from the computational mesh
void HideSBFemVolumeElements();

int volmat = 1;
int cutmat = 2;
int profilemat = 3;
int boundmat = 4;
int pointmat = 5;
int trailingedgemat = 6;
int nacamat = 7;

int sbfem_skeleton = 8;
int sbfem_domain = 9;

int sbfem_highperm = 10;
REAL shift_distance = 1.e-2;
TPZSBFemElementGroup *sbfem_groupH1 = 0;
TPZSBFemElementGroup *sbfem_groupHdiv = 0;
enum MMeshStyle {ETraditional, ECollapsed, EQuarterPoint, ESBFem};
MMeshStyle meshstyle = ECollapsed;


enum BBetaDetermination {Joukowski, Minimization};
BBetaDetermination betadetermination = Minimization;

enum RRefinementStyle {h, hp};
RRefinementStyle refinementstyle = hp;

std::map<int,TPZAutoPointer<TPZRefPattern>> refpattern;
int64_t trailingedge_element_index = -1;

auto f_profile = [](const TPZVec<REAL> &loc, TPZVec<STATE> &rhsVal, TPZFMatrix<STATE> &matVal)
{
    matVal(0,0) = 0.0;
    matVal(1,0) = 0.0;
    rhsVal[0] = -loc[1]*matVal(0,0)+loc[0]*matVal(1,0);

};

auto f_infinity = [](const TPZVec<REAL> &loc, TPZVec<STATE> &rhsVal, TPZFMatrix<STATE> &matVal)
{
    REAL signal;
    if(betadetermination == Minimization) {
        signal = -1.0;
    } else {
        signal = 1.0;
    }
    matVal(0,0) = 1.0;
    matVal(1,0) = 0.5;
    rhsVal[0] = -signal*loc[1]*matVal(0,0)+signal*loc[0]*matVal(1,0);

};

int main() {

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    if(meshstyle != ETraditional) {
        gRefDBase.InitializeRefPatterns(2);
    }
    TPZManVector<REAL, 3> x0(3, 0.);
    REAL length = 10.;
    REAL angle = 0.0;
    TPZNacaProfile naca(length, 12, angle, x0);

    TPZGeoMesh *gmesh = ReadGmsh("naca.msh", naca);
    {
        std::ofstream out("gmesh.txt");
        gmesh->Print(out);

        std::ofstream out2("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
    }
    {
        int64_t nel = gmesh->NElements();
        for(int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel && gel->MaterialId() == trailingedgemat) {
                trailingedge_element_index = el;
                break;
            }
        }
        if(trailingedge_element_index == -1) DebugStop();
    }
    // look for a quadrilateral element and keep the result in a global variable
    if(meshstyle != ETraditional) {
        auto manual = CreateRefPattern();
        refpattern[4] = manual;
        refpattern[6] = manual;
    }
    TPZCheckGeom check(gmesh);
    int uniform = 3;
    if (uniform)
    {
        check.UniformRefine(uniform);
        if(0)
        {
            std::ofstream out4("gmeshfine.txt");
            gmesh->Print(out4);

            std::ofstream out5("gmeshfine.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out5);
        }
    }
    // switch the mapped elements to linear or blend elements
    // make the elements neighbouring the trailing edge linear
    AdjustGeometry(gmesh);
    if(0)
    {
        std::ofstream out("gmeshadjusted.txt");
        gmesh->Print(out);

        std::ofstream out2("gmeshadjusted.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
    }
#ifdef PZDEBUG
    // check the node coordinates
    // this broke before correcting the element type switch
    {
        const int64_t nel = gmesh->NElements();
        for(int64_t el = 0; el < nel; el++)
        {
            TPZGeoEl * gel = gmesh->ElementVec()[el];
            if(!gel) continue;
            gel->VerifyNodeCoordinates();
        }//for el
    }
#endif

    // change the elements touching the trailing edge to collapsed elements
    if(meshstyle != ETraditional)
    {
        // false : create the midside node in the middle of the trailing edge
        Changetrailingedgeelements(gmesh);
        if(0) {
            std::ofstream out("gmeshchanged.txt");
            gmesh->Print(out);

            std::ofstream out2("gmeshchanged.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);

            PrintTrailingEdgeElements(gmesh);
        }
    }
    int nrefinements = 13;
    int minh = uniform + 1;
    // indicating the flux order
    int64_t nel = gmesh->NElements();
    TPZVec<int> porders(nel, defaultporder);
    REAL GlobalError;
    REAL circulation_H1;
    REAL circulation_HDiv;
    std::ofstream outGE("GlobalError.txt");
    for (int64_t i = 0; i < nrefinements; i++) {
        TPZCompMesh *cmesh = 0;
        TPZMultiphysicsCompMesh *cmesh_m = 0;
        TPZGeoMesh *gmeshcopy = new TPZGeoMesh(*gmesh);
        // change the elements of gmeshcopy to quadratic elements
        // observe that is applied to the copy. It does not affect the original mesh
        if(meshstyle == EQuarterPoint) {
            CreateQuarterPointElements(gmeshcopy);
            PrintTrailingEdgeElements(gmeshcopy);
        }
        // this is where skeleton and high permeability elements are created
        // this is applied to the temporary copy of the mesh
        // these elements will be deleted in the method RemoveSBFemElements
        if(meshstyle == ESBFem) {
            AdjustToSBFemGeometry(gmeshcopy);
            porders.Resize(gmeshcopy->NElements(),defaultporder);
        }

        if(betadetermination == Minimization) {
            cmesh = SimulateNacaProfileH1_Minimization(gmeshcopy, porders,circulation_H1);
        } else {
            cmesh = SimulateNacaProfileH1(gmeshcopy, porders,circulation_H1);
        }
        if(0)
        {
            std::ofstream out("cmeshH1.txt");
            cmesh->Print(out);
        }
        if(betadetermination == Minimization) {
            cmesh_m = SimulateNacaProfileHDiv_Minimization(gmeshcopy, porders,circulation_HDiv);
        } else {
            cmesh_m = SimulateNacaProfileHDiv(gmeshcopy, porders,circulation_HDiv);
        }
        {
            std::ofstream out("cmeshHdiv.txt");
            cmesh_m->Print(out);
        }
        {
            const std::string plotfile = "postprocess_H1";
            constexpr int vtkRes{2};
            TPZVec<std::string> fields = {"Pressure", "Flux"};
            auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
            vtk.SetNThreads(0);
            vtk.Do();
        }
        {
            const std::string plotfile = "postprocess_Hdiv";
            constexpr int vtkRes{2};
            TPZVec<std::string> fields = {"DivFlux", "Flux"};
            auto vtk = TPZVTKGenerator(cmesh_m, fields, plotfile, vtkRes);
            vtk.SetNThreads(0);
            vtk.Do();
        }

        TPZVec<REAL> Error;
        ComputeErrorEstimator(cmesh, cmesh_m, Error);
        ComputeGlobalError(Error, GlobalError);

        // this step is essential to regain compatibility with the original mesh
        // otherwise, computational elements would not have corresponding geometric elements
        RemoveSBFemElements(cmesh, cmesh_m, gmesh);
        /// reset the reference to the original geometric mesh
        cmesh->SetReference(gmesh);
        cmesh_m->SetReference(gmesh);
        for(int i = 0; i<2; i++) {
            cmesh_m->MeshVector()[i]->SetReference(gmesh);
        }
        delete gmeshcopy;
        gmeshcopy = 0;

        int64_t nDOF = cmesh->NEquations();
        // int64_t nDOF_m = cmesh_m->NEquations();
        outGE << i << "  " << nDOF << "  " << GlobalError << "  " << circulation_H1 << "  " << circulation_HDiv << std::endl;

        TPZVec<REAL> RefinementIndicator;
        if(refinementstyle == hp) {
            HPrefinement(cmesh_m,Error,minh,RefinementIndicator,porders);
        } else {
            Hrefinement(cmesh_m, Error, RefinementIndicator,porders);
        }
        {
            std::ofstream out2("gmeshrefined.txt");
            gmesh->Print(out2);
            std::ofstream out6("ErrorEstimator.vtk");
            TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m, out6, Error, "Error");
        }

        {
            int64_t nnod = cmesh->NConnects();
            for (int64_t i = 0; i < nnod; i++) {
                TPZConnect &c = cmesh->ConnectVec()[i];
                if (c.HasDependency()) {
                    c.RemoveDepend();
                }
            }
        }
        delete cmesh;
        TPZVec<TPZCompMesh *> meshvec = cmesh_m->MeshVector();
        {
            int64_t nnod = meshvec[0]->NConnects();
            for (int64_t i = 0; i < nnod; i++) {
                TPZConnect &c = meshvec[0]->ConnectVec()[i];
                if (c.HasDependency()) {
                    c.RemoveDepend();
                }
            }
        }
        delete cmesh_m;
        delete meshvec[0];
        delete meshvec[1];
    }
    delete gmesh;
    return 0;
}

/// @brief verify is the derivative of the NACA coordinate is correct
/// @param naca profile object
/// @param point parametric coordinate around which the derivative will be verified
void VerifyDerivative(TPZBlendNACA &naca, REAL point)
{
    TPZFNMatrix<20,REAL> pointvalues(10,2,0.), expected(10,2,0.);
    TPZManVector<REAL,10> error(10,0.), distance(10,0.);
    TPZManVector<REAL,2> derivative = {naca.dxla(point),naca.dyla(point)};
    std::cout << "derivative " << derivative << std::endl;
    REAL delta = 0.001;
    for(int i=0; i<10; i++) {
        REAL pos = point + delta * i;
        pointvalues(i,0) = naca.xla(pos);
        pointvalues(i,1) = naca.yla(pos);
        distance[i] = delta*i;
        expected(i,0) = pointvalues(0,0)+derivative[0]*distance[i];
        expected(i,1) = pointvalues(0,1)+derivative[1]*distance[i];
        error[i] = sqrt((expected(i,0)-pointvalues(i,0))*(expected(i,0)-pointvalues(i,0))+
                        (expected(i,1)-pointvalues(i,1))*(expected(i,1)-pointvalues(i,1)));
    }
    std::cout << "error " << error << std::endl;
    TPZManVector<REAL,9> rate(9,0.);
    for(int i=1; i<9; i++) {
        rate[i] = log(error[i+1]/error[i])/log(distance[i+1]/distance[i]);
    }
    std::cout << "rate " << rate << std::endl;
}

void CheckNacaProfile(const TPZNacaProfile &nacaorig) {
    TPZNacaProfile naca(nacaorig);
    naca.lastpos = 10.;
    TPZManVector<REAL,3> x0(3,0.);
    for(int i=0; i<3; i++) x0[i] = nacaorig.fX0[i];
    REAL length = naca.fCord;
    naca.ParametricDomainNodeCoord(0, x0);
    
    if(0)
    {
        REAL par = 1.5;
        int uplow = 0;
        int maxpt = 1000;
        TPZManVector<REAL,2> coord = {naca.xla(par),naca.yla(par)};
        naca.NearestParameter(coord, uplow, maxpt, par);

        std::cout << "Computed par = " << par << std::endl;

        VerifyDerivative(naca, 1.);
    }
    int count = 6;
    TPZFMatrix<REAL> coup(count,2,0.),colow(count,2,0.);
    for (size_t i = 0; i < count; i++)
    {
        /* code */
        REAL par = length*pow(i*1./(count-1),1);
        colow(i,0) = naca.xla(par);
        colow(i,1) = naca.yla(par);
        coup(i,0) = naca.xua(par);
        coup(i,1) = naca.yua(par);
        std::cout << "par = " << par << " xl " << naca.xla(par) << " yl " << naca.yla(par) << std::endl;
        std::cout << "par = " << par << " xu " << naca.xua(par) << " yu " << naca.yua(par) << std::endl;
    }
    std::ofstream out("naca.nb");

    colow.Print("colow = ",out,EMathematicaInput);
    coup.Print("coup = ",out,EMathematicaInput);
    out << "ListPlot[{colow, coup}]\n";
}

/// @brief Read a gmsh file and return a geometric mesh
/// @param filename 
/// @return a geometric mesh
TPZGeoMesh *ReadGmsh(const std::string &meshfilename, const TPZNacaProfile &naca)
{
    TPZGmshReader gmsh;
    gmsh.SetVerbose(1);
    gmsh.GetDimNamePhysical()[0]["FixPoint"] = pointmat;
    gmsh.GetDimNamePhysical()[0]["Trailingedge"] = trailingedgemat;
    gmsh.GetDimNamePhysical()[1]["Profile"] = profilemat;
    gmsh.GetDimNamePhysical()[1]["OuterBoundary"] = boundmat;  
    gmsh.GetDimNamePhysical()[1]["Cut"] = cutmat;
    gmsh.GetDimNamePhysical()[2]["Domain"] = volmat;
    auto gmesh = gmsh.GeometricGmshMesh(meshfilename);

    for (int64_t el = 0; el < gmesh->NElements(); el++) {
        TPZGeoEl* gel = gmesh->Element(el);
        // Your code here
        if(gel->MaterialId() == profilemat) {
            // std::cout << "gel " << el << " index " << gel->Index() << std::endl;
            // gel->Print(std::cout);
            TPZGeoElSide gelside(gel);
            TPZManVector<int64_t,2> nodes(2);
            nodes[0] = gel->NodeIndex(0);
            nodes[1] = gel->NodeIndex(1);
            TPZGeoElRefPattern< TPZNacaProfile> *nacael = new TPZGeoElRefPattern< TPZNacaProfile>(nodes,nacamat,*gmesh);
            //std::cout << "nacael index " << nacael->Index() << std::endl;
            nacael->Geom() = naca;
            nacael->Geom().fNodeIndexes[0] = nodes[0];
            nacael->Geom().fNodeIndexes[1] = nodes[1];
            nacael->Initialize();
            TPZManVector<REAL,3> coord(3,0.);
            TPZManVector<REAL,1> ksi(1,-1.);
            nacael->X(ksi,coord);
            TPZFNMatrix<3,REAL> gradx(3,1);
            nacael->GradX(ksi,gradx);
            std::cout << "gradx " << gradx;
            //std::cout << "ksi " << ksi << " coord " << coord << std::endl;
            ksi[0] = 1.;
            nacael->X(ksi,coord);
            nacael->GradX(ksi,gradx);
            std::cout << "gradx " << gradx;
            //std::cout << "ksi " << ksi << " coord " << coord << std::endl;
            TPZGeoElSide nacaside(nacael);
            nacaside.CenterX(coord);
            REAL par = 0.;
            int uplow = 0;
            int maxpt = 1000;
            //std::cout << "center coord " << coord << std::endl;
            // naca.NearestParameter(coord, uplow, maxpt, par);
            // std::cout << "par = " << par << " coord " << coord << std::endl;

            nacaside.InsertConnectivity(gelside);

            // nacael->Print(std::cout);

            TPZStack<TPZGeoElSide> neighbours;
            TPZGeoElSide neighbour = nacaside.Neighbour();
            while(neighbour != nacaside)
            {
                neighbours.Push(neighbour);
                neighbour = neighbour.Neighbour();
            }
            for (int i = 0; i < neighbours.size(); i++)
            {
                TPZGeoElSide neighbour = neighbours[i];
                TPZGeoEl *neighgel = neighbour.Element();
                if(neighgel->IsGeoBlendEl()) continue;
                pzgeom::SwitchToBlend(neighgel);
            }
        }
    }
    return gmesh;

}

void BuildBlueRedElements(TPZGeoMesh *gmesh, std::set<int64_t> &blue, std::set<int64_t> &red) {
    int64_t nel = gmesh->NElements();
    std::set<TPZGeoEl *> cutelements;
    for(int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->MaterialId() == cutmat) {
            cutelements.insert(gel);
            TPZGeoElSide gelside(gel);
            TPZGeoElSide neighbour = gelside.Neighbour();
            int nfound = 0;
            // look at the neighbours of the cut elements
            // compute the parameter transformation. If the det is positive
            while(neighbour != gelside)
            {
                int64_t neighindex = neighbour.Element()->Index();
                if(neighbour.Element()->MaterialId() == volmat)
                {
                    TPZTransform<REAL> tr(1);
                    gelside.SideTransform3(neighbour, tr);
                    if(tr.Mult()(0,0) > 0.) {
                        blue.insert(neighindex);
                    } else {
                        red.insert(neighindex);
                    }
                    nfound++;
                }
                neighbour = neighbour.Neighbour();
            }
            // if(nfound != 2) DebugStop();
        }
    }
    // loop over the dim-1 sides of the cut elements
    for(auto gelcut : cutelements) {
        for(int side = gelcut->FirstSide(1); side < gelcut->NSides()-1; side++) {
            bool change = true;
            while(change) {
                TPZGeoElSide gelside(gelcut,side);
                // loop over the neighbours of the dim-1 sides of the gelcut element
                for(auto neighbour = gelside.Neighbour(); neighbour != gelside; neighbour++) {
                    int64_t neighindex = neighbour.Element()->Index();
                    if(blue.find(neighindex) != blue.end() || red.find(neighindex) != red.end()) {
                        continue;
                    }
                    // loop over the face sides of the neighbour element. 
                    // if there is a neighbour of a colour, adopt that colour
                    bool found = false;
                    for(int faceside = neighbour.Element()->FirstSide(1); faceside < neighbour.Element()->NSides()-1; faceside++) {
                        TPZGeoElSide face(neighbour.Element(),faceside);
                        for(auto neighbourface = face.Neighbour(); neighbourface != face; neighbourface++) {
                            int64_t neighfaceindex = neighbourface.Element()->Index();
                            if(blue.find(neighfaceindex) != blue.end()) {
                                blue.insert(neighindex);
                                found = true;
                                change = true;
                            }
                            if(red.find(neighfaceindex) != red.end()) {
                                red.insert(neighindex);
                                found = true;
                                change=true;
                            }
                            if(found) break;
                        }
                        if(found) break;
                    }
                }
                // if no change was made, the all neighours should be either blue or red
                if(change == false) {
                    for(auto neighbour = gelside.Neighbour(); neighbour != gelside; neighbour++) {
                        int64_t neighindex = neighbour.Element()->Index();
                        if(blue.find(neighindex) == blue.end() && red.find(neighindex) == red.end()) {
                            DebugStop();
                        }
                    }
                }
            }
        }
    }
}

#include "pzintel.h"
/// @brief Create a computational mesh with H1 elements
TPZCompMesh *CreateH1CompMesh(TPZGeoMesh *gmesh, TPZVec<int> &porders, int64_t &newcon)
{
    if(porders.size() != gmesh->NElements()) DebugStop();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    int dim = 2;
    TPZDarcyFlow *material = new TPZDarcyFlow(volmat,dim);
    cmesh->InsertMaterialObject(material);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(2);

    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);

    auto bnd1 = material->CreateBC(material, profilemat, 1, val1, val2);
    TPZBndCondT<STATE> *bndcond = dynamic_cast<TPZBndCondT<STATE> *>(bnd1);
    cmesh->InsertMaterialObject(bnd1);
    auto bnd2 = material->CreateBC(material, boundmat, 1, val1, val2);
    bnd2->SetForcingFunctionBC(f_infinity,1);
    cmesh->InsertMaterialObject(bnd2);
    auto bnd3 = material->CreateBC(material, pointmat, 0, val1, val2);
    cmesh->InsertMaterialObject(bnd3);
    if(meshstyle == ESBFem) {
        auto bnd4 = material->CreateBC(material, sbfem_skeleton, 1, val1, val2);
        cmesh->InsertMaterialObject(bnd4);
    }
    cmesh->SetAllCreateFunctionsContinuous();
    std::set<int> matidsh1 = {volmat,profilemat,boundmat};
    if(meshstyle == ESBFem) {
        matidsh1.insert(sbfem_skeleton);
    }
    {
        int64_t nel = gmesh->NElements();
        TPZStack<int64_t > gelstack;
        TPZStack<int> elp;
        for(int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->HasSubElement()) continue;
            int matid = gel->MaterialId();
            if(matidsh1.find(matid) != matidsh1.end()) {
                gelstack.Push(gel->Index());
                if(porders[el] < 0) DebugStop();
                elp.Push(porders[el]+1);
            }
        }
        // cmesh->ApproxSpace().CreateCompEl(gmesh->Element(507),*cmesh);
        cmesh->ApproxSpace().BuildMesh(*cmesh, gelstack, elp);
    }
    cmesh->ExpandSolution();
    cmesh->ComputeNodElCon();
    {
        std::fstream out("cmesh.txt");
        cmesh->Print(out);
    }

    // cmesh->AutoBuild(matidsh1);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    std::set<int> matidpoint = {pointmat};
    cmesh->AutoBuild(matidpoint);
    // // insert the cut boundary condition
    newcon = cmesh->AllocateNewConnect(1,1,1);
    std::cout << "newcon " << newcon << std::endl;
    int64_t nel = gmesh->NElements();
    // // loop over the geometric cut elements to duplicate the connects
    // now we are debugging!
    // Create a map of the connect indexes to the new connect indexes
    // create a set of computational element indexes which are one the side of the cut boundary
    std::map<int64_t,int64_t> connectmap, connectmapinverse;
    std::set<int64_t> blueelements, redelements;
    BuildBlueRedElements(gmesh, blueelements, redelements);

    {
        TPZVec<int> elcolor(nel,0);
        for(auto el : blueelements) elcolor[el] = -1;
        for(auto el : redelements) elcolor[el] = 1;
        std::ofstream out("cutcolor.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, elcolor);
    }
    for(int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->MaterialId() != cutmat) continue;
        TPZGeoElSide gelside(gel);
        TPZGeoElSide neighbour = gelside.Neighbour();
        // look at the neighbours of the cut elements
        // compute the parameter transformation. If the det is positive then the element is blue
        while(neighbour != gelside)
        {
            int64_t neighindex = neighbour.Element()->Index();
            if(blueelements.find(neighindex) != blueelements.end()) {
                break;
            }
            neighbour = neighbour.Neighbour();
        }
        // if(neighbour == gelside) DebugStop();
        TPZCompEl *cel = neighbour.Element()->Reference();
        if(!cel) continue;
        int nsidenodes = gel->NSideNodes(gelside.Side());
        // put the connect indexes in the map
        for(int is=0; is<nsidenodes; is++)
        {
            int nodelocindex = neighbour.Element()->SideNodeLocIndex(neighbour.Side(),is);
            int64_t cindex = cel->ConnectIndex(nodelocindex);
            TPZConnect &c = cel->Connect(nodelocindex);
            if(c.HasDependency()) continue;
            connectmap[cindex] = -1;
        }
    }
    // duplicate the connects and initialize the dependency
    for(auto it : connectmap)
    {
        int64_t cindex = it.first;
        int64_t cindex2 = cmesh->AllocateNewConnect(1,1,1);
        TPZConnect &c2 = cmesh->ConnectVec()[cindex2];
        TPZFNMatrix<1,STATE> val(1,1,1.);
        c2.AddDependency(cindex2,cindex,val,0,0,1,1);
        c2.AddDependency(cindex2,newcon,val,0,0,1,1);
        connectmap[cindex] = cindex2;
        connectmapinverse[cindex2] = cindex;
    }
    // change the connect indexes of the blue elements
    for(auto blue : blueelements) {
        TPZGeoEl *gel = gmesh->Element(blue);
        TPZCompEl *cel = gel->Reference();
        if(!cel) continue;
        int nc = gel->NCornerNodes();
        for(int ic=0; ic<nc; ic++)
        {
            int64_t cindex = cel->ConnectIndex(ic);
            TPZConnect &c = cmesh->ConnectVec()[cindex];
            if(c.HasDependency()) {
                auto dep = c.FirstDepend();
                while(dep) {
                    int64_t cindex2 = dep->fDepConnectIndex;
                    if(connectmap.find(cindex2) != connectmap.end()) {
                        dep->fDepConnectIndex = connectmap[cindex2];
                    }
                    dep = dep->fNext;
                }
            }
            if(connectmap.find(cindex) == connectmap.end()) continue;
            cel->SetConnectIndex(ic,connectmap[cindex]);
        }
    }

    int64_t nelem = cmesh->NElements();
    for(int64_t el = 0; el<nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        if(gel->Type() != EQuadrilateral) continue;
        for(int ic = 0; ic<4; ic++) {
            int ic1 = (ic+1)%4;
            if(cel->ConnectIndex(ic) == cel->ConnectIndex(ic1)) {
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                if(!intel) DebugStop();
                intel->SetSideOrder(ic+4,1);
            }
        }
    }

    cmesh->ExpandSolution();
    cmesh->ComputeNodElCon();
    {
        std::ofstream out("cmesh.txt");
        cmesh->Print(out);
    }
    cmesh->CleanUpUnconnectedNodes();

    if(meshstyle == ESBFem) {
        CreateH1SBFEMelements(cmesh);
    }
    {
        std::fstream out("cutsol.vtk");
        TPZConnect &c = cmesh->ConnectVec()[newcon];
        int64_t pos = cmesh->Block().Position(c.SequenceNumber());
        TPZFMatrix<STATE> &meshsol = cmesh->Solution();
        meshsol(pos,0) = 1;
        cmesh->LoadSolution(meshsol);
        TPZVec<std::string> fields = {
            "Pressure",
            "Flux"
        };
        auto vtk = TPZVTKGenerator(cmesh, fields, "cutsolution", 0);
        vtk.SetNThreads(0);
        vtk.Do();
    }
    return cmesh;
}

/// @brief Adjust the integration rule for Equarterpoints H1 elements
void AdjustH1Equarterpointintrule(TPZGeoMesh *gmesh)

{
    int64_t nel = gmesh->NElements();
    // std::set<TPZCompEl *> TrailingedgeH1Elements;
    TPZCompEl *TrailingedgeH1Element;
    TPZGeoElSide TrailingSide;

    // Finding and build a Stack of the volumetric Neighbors of the Trailingedge;
    TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
    TrailingSide = TPZGeoElSide(trailingedge_element);
    int64_t trailingnode = trailingedge_element->NodeIndex(0);
    auto Neighbor = TrailingSide.Neighbour();
    for(; Neighbor != TrailingSide; Neighbor++)
    {
        if(!Neighbor.Element()) DebugStop();
        if(Neighbor.Element()->HasSubElement()) continue;
        if (Neighbor.Element()->MaterialId() == volmat || Neighbor.Element()->MaterialId() == profilemat)  
        {
            TrailingedgeH1Element = Neighbor.Element()->Reference();
            TPZGeoElSide NeighborEl(Neighbor.Element());
            TPZIntPoints* intrule = NeighborEl.CreateIntegrationRule(10);
            TrailingedgeH1Element->SetIntegrationRule(intrule);
        }  
    }
}

/// @brief Create a computational mesh with L2 elements
 TPZCompMesh *CreateL2CompMesh(TPZGeoMesh *gmesh)
 {
     TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
     int dim = 2;

     cmesh->SetDimModel(dim);
     cmesh->SetDefaultOrder(1);

     cmesh->SetAllCreateFunctionsContinuous();
     cmesh->ApproxSpace().CreateDisconnectedElements(true);
     //o método ApproxSpace().CreateDisconnectedElements(true) é chamado para transformar os elementos em elementos discontinuos, criando assim um espaço L2 de aproximação.

     cmesh->AutoBuild();
    
     return cmesh;
 }

/// @brief Create a computational mesh with HDiv elements
TPZCompMesh *CreateHDivCompMesh(TPZGeoMesh *gmesh, TPZVec<int> &porders, int64_t &newcon)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    int dim = 2;
    TPZDarcyFlow *material = new TPZDarcyFlow(volmat,dim);
    cmesh->InsertMaterialObject(material);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(1);

    HDivFamily hdiv = HDivFamily::EHDivKernel;
    cmesh->ApproxSpace().SetHDivFamily(hdiv);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);

    TPZFMatrix<STATE> val1(2,2,0.);
    TPZManVector<STATE> val2(2,0.);

    // Neumann condition on the profile
    auto bnd1 = material->CreateBC(material, profilemat, 5, val1, val2);
    bnd1->SetForcingFunctionBC(f_profile,1);
    cmesh->InsertMaterialObject(bnd1);

    // Neumann condition on the infinity boundary
    auto bnd2 = material->CreateBC(material, boundmat, 5, val1, val2);
    bnd2->SetForcingFunctionBC(f_infinity,1);
    cmesh->InsertMaterialObject(bnd2);

    if(meshstyle == ESBFem) {
        auto bnd3 = material->CreateBC(material, sbfem_skeleton, 0, val1, val2);
        cmesh->InsertMaterialObject(bnd3);
    }

    std::set<int> matidshdiv = {volmat,profilemat,boundmat};
    // build the computational mesh. HDivKernel elements for matidshdiv
    // H1 elements for the skeleton
    {
        int64_t nel = gmesh->NElements();
        TPZStack<int64_t > gelstack, gelstack_skel;
        TPZStack<int> elp, elp_skel;
        for(int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->HasSubElement()) continue;
            int matid = gel->MaterialId();
            if(matidshdiv.find(matid) != matidshdiv.end()) {
                gelstack.Push(gel->Index());
                if(porders[el] < 0) DebugStop();
                elp.Push(porders[el]);
            } else if(matid == sbfem_skeleton && meshstyle == ESBFem) {
                gelstack_skel.Push(gel->Index());
                if(porders[el] < 0) DebugStop();
                elp_skel.Push(SBFemOrder);
            }
        }
        cmesh->ApproxSpace().BuildMesh(*cmesh, gelstack, elp);
        if(meshstyle == ESBFem) {
            // insert the skeleton elements as H1 elements
            cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
            cmesh->ApproxSpace().BuildMesh(*cmesh, gelstack_skel, elp_skel);
        }
        // adjust the connect order of the skeleton elements
        for(int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->Dimension() != 1) continue;
            if(gel->HasSubElement()) continue;
            int matid = gel->MaterialId();
            if(matid == sbfem_skeleton && meshstyle == ESBFem) {
                TPZCompEl *cel = gel->Reference();
                if(!cel) DebugStop();
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                if(!intel) DebugStop();
                intel->SetSideOrder(2,SBFemOrder);
            }
        }
        cmesh->ExpandSolution();
    }

    // cmesh->AutoBuild(matidshdiv);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();

    // Fixed point condition   
    // auto bnd3 = new TPZNullMaterial(pointmat);
    // bnd3->SetDimension(0);
    // cmesh->InsertMaterialObject(bnd3);
    // std::set<int> Pointmat = {pointmat};
    // cmesh->AutoBuild(Pointmat);
// insert the circulation condition
    {
        std::ofstream out("cmesh_hdiv.txt");
        cmesh->Print(out);
    }
    newcon = -1;
    int64_t nel = gmesh->NElements();
    std::set<int64_t> profileconnects;
    // // loop over the geometric cut elements
    for(int64_t el = 0; el<nel; el++)
    {
       TPZGeoEl *gel = gmesh->Element(el);
       if(gel->HasSubElement()) continue;
       if(gel->MaterialId() != profilemat) continue;
       TPZCompEl *cel = gel->Reference();
       if(!cel) DebugStop();
       int nsidenodes = gel->NCornerNodes();
       // loop over the connects of the computational element
       // if the connect has no dependency, make it dependent on the new connect
       for(int is=0; is<nsidenodes; is++)
       {
            int64_t cindex = cel->ConnectIndex(is);
            profileconnects.insert(cindex);
            TPZConnect &c = cel->Connect(is);
            if(c.HasDependency()) DebugStop();
        }

        // make sure the connect on the profile is of order 1
        {
            int64_t cindex = cel->ConnectIndex(2);
            TPZConnect &c = cel->Connect(2);
            if(c.HasDependency()) DebugStop();
            int64_t seq = c.SequenceNumber();
            c.SetOrder(0,cindex);
            c.SetNState(1);
            c.SetNShape(0);
            // reset the size of the block of the connect
            cmesh->Block().Set(seq,0);
        }
    }
    newcon = *profileconnects.begin();
    std::cout << "newcon " << newcon << std::endl;

    int64_t nelem = cmesh->NElements();
    for(int64_t el = 0; el<nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        int ncorner = gel->NCornerNodes();
        for(int ic = 0; ic<ncorner; ic++) {
            int64_t cindex = cel->ConnectIndex(ic);
            if(profileconnects.find(cindex) != profileconnects.end()) {
                cel->SetConnectIndex(ic,newcon);
            }
        }
    }
    // force the side corresponding to a collapsed collect to have order 1
    for(int64_t el = 0; el<nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        if(gel->Type() != EQuadrilateral) continue;
        for(int ic = 0; ic<4; ic++) {
            int ic1 = (ic+1)%4;
            if(cel->ConnectIndex(ic) == cel->ConnectIndex(ic1)) {
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                if(!intel) DebugStop();
                intel->SetSideOrder(ic+4,1);
            }
        }
    }
    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    {
        std::ofstream out("cmeshHdiv.txt");
        cmesh->Print(out);
        std::ofstream out2("cmeshHDiv.vtk"); 
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out2);
    }

    return cmesh;
}

/// @brief Adjust the integration rule for Equarterpoints HDiv elements
void AdjustHDivEquarterpointintrule(TPZGeoMesh *gmesh)

{
    int64_t nel = gmesh->NElements();
    // std::set<TPZMultiphysicsElement *> TrailingedgeHDivElements;
    TPZMultiphysicsElement *TrailingedgeHDivElement;
    TPZGeoElSide TrailingSide;

    // Finding and build a Stack of the volumetric Neighbors of the Trailingedge;
    TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
    TrailingSide = TPZGeoElSide(trailingedge_element);
    int64_t trailingnode = trailingedge_element->NodeIndex(0);
    auto Neighbor = TrailingSide.Neighbour();
    for(; Neighbor != TrailingSide; Neighbor++)
    {
        if(!Neighbor.Element()) DebugStop();
        if(Neighbor.Element()->HasSubElement()) continue;
        if (Neighbor.Element()->MaterialId() == volmat || Neighbor.Element()->MaterialId() == profilemat)
        {
            auto NeighborCompel = Neighbor.Element()->Reference();
            TrailingedgeHDivElement = dynamic_cast<TPZMultiphysicsElement*>(NeighborCompel);
            if (!TrailingedgeHDivElement) DebugStop();
            TPZGeoElSide NeighborEl(Neighbor.Element());
            TPZIntPoints* intrule = NeighborEl.CreateIntegrationRule(10);
            TrailingedgeHDivElement->SetIntegrationRule(intrule);
        }  
    }
}

/// @brief Create a computational "multiphysics" mesh with only HDiv elements
 TPZMultiphysicsCompMesh *CreateMultiphysicsMesh(TPZCompMesh* cmeshHDiv,TPZCompMesh* cmeshL2, TPZGeoMesh *gmesh){
//
     TPZMultiphysicsCompMesh *cmesh_m = new TPZMultiphysicsCompMesh(gmesh);

     const int dim = cmeshHDiv->Dimension();
     const int pord = cmeshHDiv->GetDefaultOrder();
     cmesh_m->SetDimModel(dim);
     cmesh_m->SetDefaultOrder(pord);
         
     cmesh_m->SetAllCreateFunctionsMultiphysicElem();

     // 1. Materials
     std::set<int> materialIDs;
     TPZMixedDarcyFlow *material = new TPZMixedDarcyFlow(volmat,dim);
    //  material->SetBigNumber(10e8);
     cmesh_m->InsertMaterialObject(material);
     materialIDs.insert(volmat);

     // 2. Boundary Conditions
    TPZFMatrix<STATE> val1(2,2,0.);
    TPZManVector<STATE> val2(2,0.);

    // 2.1 Neumann condition on the profile
    auto bnd1 = material->CreateBC(material, profilemat, 5, val1, val2);
    bnd1->SetForcingFunctionBC(f_profile,1);
    cmesh_m->InsertMaterialObject(bnd1);

    // // 2.2 Neumann condition on the infity boundary
    auto bnd2 = material->CreateBC(material, boundmat, 5, val1, val2);
    bnd2->SetForcingFunctionBC(f_infinity,1);
    cmesh_m->InsertMaterialObject(bnd2);

    TPZNullMaterialCS<STATE> *bnd3 = new TPZNullMaterialCS<STATE>(sbfem_skeleton);
    cmesh_m->InsertMaterialObject(bnd3);

    // 2.3 Fixed point condition
    // auto bnd3 = new TPZL2ProjectionCS(pointmat,0);
    // cmesh_m->InsertMaterialObject(bnd3);

    // 3. VECTOR OF COMPUTATIONAL MESHES (datavec)
     TPZManVector<int, 2> active_approx_spaces(2, 1);
     TPZManVector<TPZCompMesh *, 2> mesh_vec(2);
     mesh_vec[0] = cmeshHDiv;
     mesh_vec[1] = cmeshL2;

     cmesh_m->BuildMultiphysicsSpace(active_approx_spaces, mesh_vec);

     cmesh_m->ExpandSolution();
     cmesh_m->LoadReferences();
     return cmesh_m;
 }

/// @brief Simulate the NACA profile using H1 approximation and Joukowski condition to find Beta
TPZCompMesh *SimulateNacaProfileH1(TPZGeoMesh *gmesh, TPZVec<int> &porders, REAL &circulation_H1)
{
    int64_t newcon;
    auto cmeshH1 = CreateH1CompMesh(gmesh, porders, newcon);
    if(meshstyle == EQuarterPoint) {
        AdjustH1Equarterpointintrule(gmesh);
    }

    TPZLinearAnalysis an(cmeshH1,RenumType::EMetis);
    if (1)
    {
        std::ofstream out("cmeshH1.txt");
        cmeshH1->Print(out);
    }
#ifdef USING_MKL
    TPZSSpStructMatrix<STATE> strmat(cmeshH1);
#else
    TPZSkylineStructMatrix<STATE> strmat(cmeshH1);
#endif
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Assemble();

    int nvar = an.Solution().Rows();

//  Adjust the System of Equation to compute phi_0 and phi_1 directly 
    AdjustSystemofEquations(an,cmeshH1,newcon);

//  Solve the system of equations and save the solutions: phi_0 e phi_1 
    an.Solve();
    int64_t numeq = cmeshH1->NEquations();
    TPZFMatrix<STATE> fsol(numeq,2,0.),phi_0(numeq,1,0.),phi_1(numeq,1,0.),phi(numeq,1,0.);
    fsol = an.Solution();
    // fsol.Print("Solution",std::cout);
    for (int64_t eq = 0; eq < numeq; eq++) 
    {
        phi_0(eq) =  -fsol(eq,0);
        phi_1(eq) =  fsol(eq,1);
    }

//  Compute the "Beta" constant
    REAL Beta;
    ComputeBetaH1(gmesh, cmeshH1, phi_0, phi_1, Beta);

    phi = phi_0+Beta*phi_1;
    cmeshH1->LoadSolution(phi);

    {
        int nc = sbfem_groupH1->NConnects();
        for(int ic=0; ic<nc; ic++) {
            TPZConnect &c = sbfem_groupH1->Connect(ic);
            c.Print(*cmeshH1, std::cout);
        }
    }
    an.LoadSolution(phi);

    circulation_H1 = -Beta;
    // EvaluateCirculationH1(gmesh, cmeshH1, profilemat, circulation_H1);

    std::cout << "--------- PostProcess H1 ---------" << std::endl;
    std::cout << "Circulation = " << circulation_H1 << std::endl;
    //printa na tela "--------- PostProcess ---------", indicando que a simulação está em processamento.
    PrintResults(an,cmeshH1);
    //chama a função PrintResults para realizar o pós-processamento dos resultados. Essa função provavelmente gera saídas com os resultados da simulação.
    return cmeshH1;
}

/// @brief Simulate the NACA profile using H1 approximation and minimization to find Beta
TPZCompMesh *SimulateNacaProfileH1_Minimization(TPZGeoMesh *gmesh, TPZVec<int> &porders, REAL &circulation_H1)
{
    int64_t newcon;
    auto cmeshH1 = CreateH1CompMesh(gmesh, porders, newcon);
    AdjustH1Equarterpointintrule(gmesh);
    if (1)
    {
        std::ofstream out("cmeshH1.txt");
        cmeshH1->Print(out);
    }

    TPZLinearAnalysis an(cmeshH1,RenumType::EMetis);
#ifdef USING_MKL
    TPZSSpStructMatrix<STATE> strmat(cmeshH1);
#else
    TPZSkylineStructMatrix<STATE> strmat(cmeshH1);
#endif
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Assemble();

    int nvar = an.Solution().Rows();

//  Solve the system of equations and save the solutions: phi_0 e phi_1 
    an.Solve();
    int64_t numeq = cmeshH1->NEquations();
    TPZFMatrix<STATE> phi(numeq,1,0.);
    phi = an.Solution();
    cmeshH1->LoadSolution(phi);
    an.LoadSolution(phi);

    EvaluateCirculationH1(gmesh, cmeshH1, profilemat, circulation_H1);

    std::cout << "--------- PostProcess H1 ---------" << std::endl;
    std::cout << "Circulation = " << circulation_H1 << std::endl;
    //printa na tela "--------- PostProcess ---------", indicando que a simulação está em processamento.
    PrintResults(an,cmeshH1);
    //chama a função PrintResults para realizar o pós-processamento dos resultados. Essa função provavelmente gera saídas com os resultados da simulação.
    return cmeshH1;
}

/// @brief Simulate the NACA profile using H(div) approximation and Joukowski condition to find Beta
 TPZMultiphysicsCompMesh *SimulateNacaProfileHDiv(TPZGeoMesh *gmesh, TPZVec<int> &porders, REAL &circulation_HDiv)
 {
    int64_t newcon;
    gmesh->ResetReference();
     auto cmeshHDiv = CreateHDivCompMesh(gmesh, porders, newcon);
     auto cmeshL2 = CreateL2CompMesh(gmesh);

     TPZMultiphysicsCompMesh* cmesh_m = CreateMultiphysicsMesh(cmeshHDiv,cmeshL2, gmesh);
    if(meshstyle == EQuarterPoint) {
        AdjustHDivEquarterpointintrule(gmesh);
    }
    if(meshstyle == ESBFem) {
        CreateHdivSBFEMelements(cmesh_m,newcon);
    }
     cmesh_m->CleanUpUnconnectedNodes();
     // Define o pointer chamado cmesh_m relacionado à classe TPZMultiphysicsCompMesh, associando-o a função CreateMultiphysicsMesh, cujos parâmetros (já antes declarados) são: cmeshHdiv, gmesh.
    
    TPZLinearAnalysis an(cmesh_m,RenumType::EMetis);
    {
        std::ofstream out("cmesh_m.txt");
        cmesh_m->Print(out);
    }


#ifdef USING_MKL
    TPZSSpStructMatrix<STATE> strmat(cmesh_m);
#else
    TPZSkylineStructMatrix<STATE> strmat(cmesh_m);
#endif
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Assemble();

    int nvar = an.Solution().Rows();

    //  Adjust the System of Equation to compute phi_0 and phi_1 directly 
    AdjustSystemofEquations(an,cmesh_m,newcon);

    //  Solve the system of equations and save the solutions: phi_0 e phi_1 
    an.Solve();
    int64_t numeq = cmesh_m->NEquations();
    TPZFMatrix<STATE> fsol(numeq,2,0.),u_0(numeq,1,0.),u_1(numeq,1,0.),u(numeq,1,0.);
    fsol = an.Solution();

    // fsol.Print("HDiv solution",std::cout);

    for (int64_t eq = 0; eq < numeq; eq++) 
    {
        u_0(eq) =  fsol(eq,0);
        u_1(eq) =  fsol(eq,1);
    }

    //  Compute the "Kappa" constant
    REAL Beta;
    ComputeBetaHDiv(gmesh, cmesh_m, u_0, u_1, Beta);

    u = u_0+Beta*u_1;
    cmesh_m->LoadSolution(u);
    cmesh_m->TransferMultiphysicsSolution();
    an.LoadSolution(u);

    EvaluateCirculationHDiv(gmesh, cmesh_m, profilemat, circulation_HDiv);

    std::cout << "--------- PostProcess HDiv---------" << std::endl;
    std::cout << "Circulation = " << circulation_HDiv << std::endl;
     PrintResults(an,cmesh_m);

     return cmesh_m;
 }

 /// @brief Simulate the NACA profile using H(div) approximation and minimization to find Beta
 TPZMultiphysicsCompMesh *SimulateNacaProfileHDiv_Minimization(TPZGeoMesh *gmesh, TPZVec<int> &porders, REAL &circulation_HDiv)
 {
    int64_t newcon;
    gmesh->ResetReference();
     auto cmeshHDiv = CreateHDivCompMesh(gmesh, porders, newcon);
     auto cmeshL2 = CreateL2CompMesh(gmesh);

     TPZMultiphysicsCompMesh* cmesh_m = CreateMultiphysicsMesh(cmeshHDiv,cmeshL2, gmesh);
    if(meshstyle == EQuarterPoint) {
        AdjustHDivEquarterpointintrule(gmesh);
    }
    if(meshstyle == ESBFem) {
        CreateHdivSBFEMelements(cmesh_m,newcon);
    }
     cmesh_m->CleanUpUnconnectedNodes();
     // Define o pointer chamado cmesh_m relacionado à classe TPZMultiphysicsCompMesh, associando-o a função CreateMultiphysicsMesh, cujos parâmetros (já antes declarados) são: cmeshHdiv, gmesh.
    
    TPZLinearAnalysis an(cmesh_m,RenumType::EMetis);
    {
        std::ofstream out("cmesh_m.txt");
        cmesh_m->Print(out);
    }


#ifdef USING_MKL
    TPZSSpStructMatrix<STATE> strmat(cmesh_m);
#else
    TPZSkylineStructMatrix<STATE> strmat(cmesh_m);
#endif
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Assemble();

    int nvar = an.Solution().Rows();

    //  Solve the system of equations and save the solutions: phi_0 e phi_1 
    an.Solve();
    int64_t numeq = cmesh_m->NEquations();
    TPZFMatrix<STATE> u(numeq,1,0.);
    u = an.Solution();
    cmesh_m->LoadSolution(u);
    cmesh_m->TransferMultiphysicsSolution();
    an.LoadSolution(u);

    EvaluateCirculationHDiv(gmesh, cmesh_m, profilemat, circulation_HDiv);

    std::cout << "--------- PostProcess HDiv---------" << std::endl;
    std::cout << "Circulation = " << circulation_HDiv << std::endl;
     PrintResults(an,cmesh_m);

     return cmesh_m;
 }

/// @brief Compute the error as an average error per element and save as an elemental solution in "Error"  
void ComputeErrorEstimator(TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator)
{
    if(meshstyle == ESBFem) {
        AddSBFemVolumeElements();
    }
    cmesh->LoadReferences();
    int64_t nel_m = cmesh_m->NElements();
    ErrorEstimator.Resize(nel_m,0.);
    ErrorEstimator = 0.;
    // TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh->ElementVec();
    TPZAdmChunkVector<TPZCompEl *> &elementvec_m = cmesh_m->ElementVec();

    for (int64_t iel = 0; iel < nel_m; iel++) 
    {
        TPZCompEl *cell_m = elementvec_m[iel];
        if(!cell_m) continue;
        TPZGeoEl *gel_m = cell_m->Reference();
        if(!gel_m) continue;
        if(gel_m->Dimension() != 2) continue;

        TPZCompEl *cell = gel_m->Reference();
        int matid = gel_m->MaterialId();

        TPZGeoElSide gelside(gel_m);
        auto intrule = gelside.CreateIntegrationRule(5);
        int intrulepoints = intrule->NPoints();
        TPZManVector<REAL,4> intpoint(2,0.);
        REAL weight = 0.;
        TPZFMatrix<REAL> jac, axe, jacInv;
        TPZManVector<REAL,3> x(3,0.);
        REAL detJac;

        for(int int_ind = 0; int_ind < intrulepoints; ++int_ind)
        {
            intrule->Point(int_ind,intpoint,weight);
            gelside.Jacobian(intpoint, jac, axe, detJac , jacInv);
            weight *= fabs(detJac);
            TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(cell);
            TPZManVector<STATE,3> flux(3,0.);
            msp ->Solution(intpoint,7,flux);

                TPZManVector<STATE,3> sol(3,0.);
            if(matid == volmat) {
                TPZMultiphysicsElement *msp_m  = dynamic_cast <TPZMultiphysicsElement *>(cell_m);
                if(!msp_m) DebugStop();
                msp_m ->Solution(intpoint,1,sol);
            } else if(matid == sbfem_domain) {
                TPZSBFemVolume *msp_m  = dynamic_cast <TPZSBFemVolume *>(cell_m);
                if(!msp_m) DebugStop();
                msp_m ->Solution(intpoint,1,sol);
                // std::cout << "flux " << flux << std::endl;
                // std::cout << "sol " << sol << std::endl;
            }

            ErrorEstimator[iel] += ((flux[0]-sol[0])*(flux[0]-sol[0])+(flux[1]-sol[1])*(flux[1]-sol[1]))*weight;
        }//loop over integratin points
    }//loop over cemsh_m elements

    {
        std::ofstream out("ErrorEstimator.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m,out,ErrorEstimator,"ErrorEstimator");
    }
    HideSBFemVolumeElements();
}


/// @brief Compute the global error the norm of the vector ErrorEstimator  
void ComputeGlobalError(TPZVec<REAL> &ErrorEstimator, REAL &GlobalError)
{
    GlobalError = 0.;
    int64_t nel_m = ErrorEstimator.size();
    for (int64_t iel = 0; iel < nel_m; iel++) 
    {
        GlobalError += ErrorEstimator[iel]*ErrorEstimator[iel];
    }//loop over the elements of the vector ErrorEstimator
    GlobalError = sqrt(GlobalError);
}


/// @brief Performe h-refinements at the computational meshes "cmesh" and "cmesh_m" based on the ErrorEstimator  
void Hrefinement(TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator, TPZVec<REAL> &RefinementIndicator, TPZVec<int> &porders)
{
    int64_t nel_m = cmesh_m->NElements();
    auto gmesh = cmesh_m->Reference();
    int64_t nel = gmesh->NElements();
    RefinementIndicator.Resize(nel,0.);
    RefinementIndicator = 0.;
    auto tolel = std::max_element(ErrorEstimator.begin(), ErrorEstimator.end());
    REAL tol = *tolel/5.0;
    // divide all two dimensional elements with error larger than tol
    for (int64_t iel = 0; iel < nel_m; iel++) 
    {
        TPZCompEl *cel = cmesh_m->Element(iel);
        TPZGeoEl *gel = cel->Reference();
        // grouped elements have no reference (SBFem simulation)
        if(!gel) continue;
        // int64_t index = gel->Index();
        if(ErrorEstimator[iel]> tol)
        {
            RefinementIndicator[iel] = 1.0;
            TPZManVector<TPZGeoEl*> subel; 
            DivideGeoEl(gel,subel);
        } 
    }//loop over cmesh_m elements
    Smoothentrailingedgeelements(cmesh_m,RefinementIndicator);

    //Change the refined gmesh to blend
    AdjustGeometry(gmesh);

    // refinar os elementos de contorno
    nel = gmesh->NElements();
    for (int64_t iel = 0; iel < nel; iel++) 
    {
        TPZGeoEl *gel = gmesh->Element(iel);
        if(gel->HasSubElement()) continue;
        if(gel->Dimension() != 1) continue;
        if(gel->MaterialId() == nacamat) continue;
        if(gel->MaterialId() == sbfem_highperm) continue;
        if(gel->MaterialId() ==  boundmat) continue;
        if(gel->MaterialId() ==  cutmat) continue;
        //Percorrer os Neighbors de gel e encontar os volumetricos que tem subelementos
        TPZGeoElSide gelside(gel); 
        auto Neighbor = gelside.Neighbour();
        for(; Neighbor != gelside; Neighbor++)
        {
            if(!Neighbor.Element()) DebugStop();
            if(!Neighbor.HasSubElement()) continue;
            TPZManVector<TPZGeoEl*> subel; 
            DivideGeoEl(gel,subel);
            DebugStop();
            break;
        }
    }//loop over cemsh_m elements
    std::ofstream out1("RefinementIndicator.vtk"); 
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m,out1,RefinementIndicator,"RefinementIndicator");

    nel = gmesh->NElements();
    RefinementIndicator.Resize(nel,0.);
    std::ofstream out1("HRefinementIndicator.vtk"); 
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m,out1,RefinementIndicator,"RefinementIndicator");

    porders.Resize(nel);
    porders.Fill(defaultporder);

    {
        std::ofstream out("gmeshrefined.txt");
        gmesh->Print(out);
    }

    {
        std::ofstream out2("gmesh_Hrefinement.vtk"); 
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
    }

}

/// @brief Performe hp-refinements at the computational meshes "cmesh" and "cmesh_m" based on the ErrorEstimator  
// @param cmesh_m multiphysics computational mesh
// @param ErrorEstimator vector with the error estimator for each computational element
// @param minh minimum level of refinement
// @param RefinementIndicator vector with the refinement indicator for each computational element
// @param porders vector with the polynomial order of each geometric element 
void HPrefinement(TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator, int minh, TPZVec<REAL> &RefinementIndicator, TPZVec<int> &porders)
{

    int64_t nel_m = cmesh_m->NElements();

    auto gmesh = cmesh_m->Reference();
    gmesh->ResetReference();
    cmesh_m->LoadReferences();
    int64_t nel = gmesh->NElements();
    RefinementIndicator.Resize(nel,0.);
    RefinementIndicator = 0.;
    auto tolel = std::max_element(ErrorEstimator.begin(), ErrorEstimator.end());
    REAL tol = *tolel/5.0;
    // if the error of a trailing edge is larger than the tolerance, divide it
    DivideTrailingEdgeNeighbours(cmesh_m, ErrorEstimator, tol, RefinementIndicator, porders);
    // make sure all the trailing edge elements have the same level of refinement
    Smoothentrailingedgeelements(cmesh_m,RefinementIndicator);
    // if the error of an element neighbouring the profile is larger than the tolerance AND
    // the number of divisions is less than minh, divide it
    // DivideProfileNeighbours(cmesh_m, ErrorEstimator, tol, minh, RefinementIndicator);
    // enforce the two on one constraint
    // divide the one dimensional elements
    SmoothenGeometry(gmesh);
    nel = gmesh->NElements();
    RefinementIndicator.Resize(nel,0.);
    int64_t psize = porders.size();
    porders.Resize(nel);
    for(int64_t i=psize; i<nel; i++) {
        porders[i] = defaultporder;
    }
    int sbfem_maxorder = 0;
    for (int64_t iel = 0; iel < nel_m; iel++) 
    {
        TPZCompEl *cel = cmesh_m->Element(iel);
        if(!cel) continue;
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        TPZSBFemVolume *sbvol = dynamic_cast<TPZSBFemVolume *>(cel);
        TPZSBFemElementGroup *sbfem_group = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if(!mfcel && !sbvol && !sbfem_group) {
            DebugStop();
        }
        if(mfcel) {
            TPZInterpolationSpace *msp = dynamic_cast<TPZInterpolationSpace *>(mfcel->Element(0));
            if(!msp) DebugStop();
            TPZGeoEl *gel = cel->Reference();
            int64_t index = gel->Index();
            int porder = porders[index];
            if(ErrorEstimator[iel] > tol)
            {
                RefinementIndicator[iel] = 1.0;
                if(gel->HasSubElement()) {
                    // only trailing edge neighbours are divided
                    TPZManVector<TPZGeoEl*> subel; 
                    DivideGeoEl(gel,subel);
                    for(auto sub : subel) {
                        porders[sub->Index()] = porder;
                    }
                } else {
                    // apply p refinement
                    porders[index] = porder+1;
                }
            } else {
                if(gel->HasSubElement()) {
                    // the element was divided to satisfy the two on one constraint
                    TPZManVector<TPZGeoEl*> subel; 
                    DivideGeoEl(gel,subel);
                    for(auto sub : subel) {
                        porders[sub->Index()] = porder;
                    }
                } else {
                    porders[index] = porder;
                }
            }
        } else if(sbvol) {
            REAL error = ErrorEstimator[iel];
            TPZGeoEl *gel = cel->Reference();
            int64_t index = gel->Index();
            if(error > tol)
            {
                RefinementIndicator[iel] = 1.0;
                porders[index]++;
                if(porders[index] > sbfem_maxorder) sbfem_maxorder = porders[index];
            } else {
                porders[index] = SBFemOrder;
            } 
        }
    }//loop over cemsh_m elements
    // set all sbfem volume elements to the same order
    if(meshstyle == ESBFem) 
    {
        if(sbfem_maxorder > SBFemOrder) {
            std::cout << "Increasing the order of the SBFem elements from " << SBFemOrder << " to " << sbfem_maxorder << std::endl;
            SBFemOrder = sbfem_maxorder;
        }
        // update the porders data structure
        for (int64_t iel = 0; iel < nel_m; iel++) 
        {
            TPZCompEl *cel = cmesh_m->Element(iel);
            if(!cel) continue;
            TPZSBFemVolume *sbvol = dynamic_cast<TPZSBFemVolume *>(cel);
            if(!sbvol) continue;
            TPZGeoEl *gel = cel->Reference();
            int64_t index = gel->Index();
            porders[index] = SBFemOrder;
        }
    }
    {
        std::ofstream out("HRefinementIndicator.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m,out,RefinementIndicator,"RefinementIndicator");
        TPZVec<double> porderscopy(cmesh_m->NElements(),0);
        for(int64_t i=0; i<cmesh_m->NElements(); i++) {
            TPZCompEl *cel = cmesh_m->Element(i);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if(!gel) continue;
            porderscopy[i] = porders[gel->Index()];
        }
        std::ofstream out2("porders.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m, out2, porderscopy, "porders");
    }
    {
        std::ofstream out("porders.txt");
        for(auto p : porders) out << p << std::endl;
    }
    {
        std::ofstream out("gmeshrefined.txt");
        gmesh->Print(out);
    }
    {
        std::ofstream out2("gmesh_HPrefinement.vtk"); 
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2, porders, "Porders",true);
    }
    if(gmesh->NElements() != nel) DebugStop();
    if(porders.size() != gmesh->NElements()) DebugStop();

    {
        std::ofstream out("gmeshrefined.txt");
        gmesh->Print(out);
    }
    {
        std::ofstream out1("RefinementIndicator.vtk"); 
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m,out1,RefinementIndicator,"RefinementIndicator");
    }
    {
        std::ofstream out2("gmesh_HPrefinement.vtk"); 
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2, porders, "Porders",true);
    }

}

/// @brief Smoothen the volumetric elements around the trailingedge by adopting the same level of refinement  
void Smoothentrailingedgeelements(TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &RefinementIndicator)
{
    auto gmesh = cmesh_m->Reference();
    int64_t nel = gmesh->NElements();
    TPZGeoEl *gel_TrailingEdge = gmesh->Element(trailingedge_element_index);
    // if(gel_TrailingEdge->HasSubElement()) continue;
    if (gel_TrailingEdge->MaterialId() != trailingedgemat) DebugStop(); 

    TPZGeoElSide TrailingSide = TPZGeoElSide(gel_TrailingEdge);
    TPZGeoElSide Neighbor = TrailingSide.Neighbour();
    TPZGeoEl *GeoNeighbor;
    int maxlevel = 0;
    for(; Neighbor != gel_TrailingEdge; Neighbor++)
    {
        if(!Neighbor.Element()) DebugStop();
        if(Neighbor.Element()->HasSubElement()) continue;
        if (Neighbor.Element()->MaterialId() != volmat) continue;  
        GeoNeighbor = Neighbor.Element();
        if(maxlevel <= GeoNeighbor->Level())
        {
            maxlevel = GeoNeighbor->Level();    
        }  
    }
    Neighbor = gel_TrailingEdge->Neighbour(0);
    for(; Neighbor != gel_TrailingEdge; Neighbor++)
    {
        if(!Neighbor.Element()) DebugStop();
        if(Neighbor.Element()->HasSubElement()) continue;
        if (Neighbor.Element()->MaterialId() != volmat) continue;  
        GeoNeighbor = Neighbor.Element();
        if(GeoNeighbor->Level() < maxlevel)
        {
            auto iel_GeoNeighbor = GeoNeighbor->Index();
            RefinementIndicator[iel_GeoNeighbor] = 1.0;
            TPZManVector<TPZGeoEl*> subel_GeoNeighbor;
            DivideGeoEl(GeoNeighbor,subel_GeoNeighbor);
        }
    }
}

#include "tpzquadraticquad.h"
#include "tpzquadraticline.h"
#include "tpzgeoblend.h"
/// @brief create a collapsed quadratic quad element
/// @param sn the index of the singular node
/// @param edge1 the index of the first edge node
/// @param edge2 the index of the second edge node
/// @param gmesh the geometric mesh
static int64_t CreateCollapsedQuad(int64_t sn, int64_t edge1, int64_t edge2, bool blend, TPZGeoMesh *gmesh)
{
    TPZManVector<int64_t,4> nodeindices(4);
    nodeindices[0] = edge1;
    nodeindices[1] = edge2;
    nodeindices[2] = sn;
    nodeindices[3] = sn;
    int64_t index;
    int matid = volmat;
    if(blend) {
        new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>>(nodeindices,matid,*gmesh,index);
    } else {
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(nodeindices,matid,*gmesh,index);
    }
    return index;
}

/// @brief change the element to a blend element
void ChangeToBlend(TPZGeoEl *gel) {
        TPZGeoMesh *gmesh = gel->Mesh();
    if(gel->Type() != EQuadrilateral) DebugStop();
    TPZManVector<int64_t,4> nodeindices(4);
    for(int i=0; i<4; i++) nodeindices[i] = gel->NodeIndex(i);
    TPZManVector<TPZGeoElSide,9> neighbours(9);
    for(int i=0; i<8; i++) {
        neighbours[i] = gel->Neighbour(i);
        while(neighbours[i].Element() == gel) {
            neighbours[i] = neighbours[i].Neighbour();
        }
    }
    int nsubel = gel->NSubElements();
    TPZManVector<TPZGeoEl *,4> subels(nsubel);
    for(int is=0; is<nsubel; is++) {
        subels[is] = gel->SubElement(is);
    }
    auto refpat = gel->GetRefPattern();
    int whichsubel = -1;
    TPZGeoEl *father = gel->Father();
    if(father) {
        whichsubel = gel->WhichSubel();
    }
    gel->RemoveConnectivities();
    int64_t gelindex = gel->Index();
    delete gel;
    int64_t index;
    auto newgel = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>>(nodeindices,volmat,*gmesh,index);    
    if(index != gelindex) DebugStop();
    for(int i=0; i<8; i++) {
        if(neighbours[i].Element()) {
            TPZGeoElSide gelside(newgel,i);
            gelside.SetConnectivity(neighbours[i]);
        }
    }
    newgel->SetNeighbour(8,newgel);
    if(father) {
        father->SetSubElement(whichsubel,newgel);
    }
    newgel->SetFather(father);
    newgel->SetRefPattern(refpat);
    for(int is=0; is<nsubel; is++) {
        newgel->SetSubElement(is,subels[is]);
    }
    newgel->BuildBlendConnectivity();
}

/// @brief If the element is not linear mapping, change it to linear mapping
/// @param gel element that will be substituted
void ChangeToLinearQuad(TPZGeoEl *gel) {
    TPZGeoMesh *gmesh = gel->Mesh();
    if(gel->Type() != EQuadrilateral) DebugStop();
    TPZManVector<int64_t,4> nodeindices(4);
    for(int i=0; i<4; i++) nodeindices[i] = gel->NodeIndex(i);
    TPZManVector<TPZGeoElSide,9> neighbours(9);
    for(int i=0; i<8; i++) {
        neighbours[i] = gel->Neighbour(i);
        while(neighbours[i].Element() == gel) {
            neighbours[i] = neighbours[i].Neighbour();
        }
    }
    int nsubel = gel->NSubElements();
    TPZManVector<TPZGeoEl *,4> subels(nsubel);
    for(int is=0; is<nsubel; is++) {
        subels[is] = gel->SubElement(is);
    }
    auto refpat = gel->GetRefPattern();
    int whichsubel = -1;
    TPZGeoEl *father = gel->Father();
    if(father) {
        whichsubel = gel->WhichSubel();
    }
    gel->RemoveConnectivities();
    int64_t gelindex = gel->Index();
    int matid = gel->MaterialId();
    delete gel;
    int64_t index;
    auto newgel = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(nodeindices,matid,*gmesh,index);
    if(index != gelindex) DebugStop();
    for(int i=0; i<8; i++) {
        if(neighbours[i].Element()) {
            TPZGeoElSide gelside(newgel,i);
            gelside.SetConnectivity(neighbours[i]);
        }
    }
    newgel->SetNeighbour(8,newgel);
    if(father) {
        father->SetSubElement(whichsubel,newgel);
    }
    newgel->SetFather(father);
    newgel->SetRefPattern(refpat);
    for(int is=0; is<nsubel; is++) {
        newgel->SetSubElement(is,subels[is]);
    }
}

/// @brief If the element is not linear mapping, change it to linear mapping
/// @param gel element that will be substituted
void ChangeToLinearOned(TPZGeoEl *gel) {
    TPZGeoMesh *gmesh = gel->Mesh();
    if(gel->Type() != EOned) DebugStop();
    TPZManVector<int64_t,4> nodeindices(2);
    for(int i=0; i<2; i++) nodeindices[i] = gel->NodeIndex(i);
    TPZManVector<TPZGeoElSide,3> neighbours(3);
    for(int i=0; i<3; i++) {
        neighbours[i] = gel->Neighbour(i);
        while(neighbours[i].Element() == gel) {
            neighbours[i] = neighbours[i].Neighbour();
        }
    }
    int nsubel = gel->NSubElements();
    TPZManVector<TPZGeoEl *,4> subels(nsubel);
    for(int is=0; is<nsubel; is++) {
        subels[is] = gel->SubElement(is);
    }
    auto refpat = gel->GetRefPattern();
    int whichsubel = -1;
    TPZGeoEl *father = gel->Father();
    if(father) {
        whichsubel = gel->WhichSubel();
    }
    gel->RemoveConnectivities();
    int64_t gelindex = gel->Index();
    int matid = gel->MaterialId();
    delete gel;
    int64_t index;
    auto newgel = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear>(nodeindices,matid,*gmesh,index);
    if(index != gelindex) DebugStop();
    for(int i=0; i<3; i++) {
        if(neighbours[i].Element()) {
            TPZGeoElSide gelside(newgel,i);
            gelside.SetConnectivity(neighbours[i]);
        }
    }
    if(father) {
        father->SetSubElement(whichsubel,newgel);
    }
    newgel->SetFather(father);
    newgel->SetRefPattern(refpat);
    for(int is=0; is<nsubel; is++) {
        newgel->SetSubElement(is,subels[is]);
    }
}

/// @brief divide a geometric element. If it is a trailing edge element, set the refinement pattern first
void DivideGeoEl(TPZGeoEl *gel, TPZVec<TPZGeoEl *> &subels)
{
    subels.resize(0);
    if(!gel) return;
    if(gel->MaterialId() < 0) {
        DebugStop();
    }
    int64_t gelindex = gel->Index();
    if(gel->Dimension() == 0) return;
    if(gel->Type() != EQuadrilateral) {
        gel->Divide(subels);
        return;
    }
    bool iscollapsed = false;
    for(int side = 4; side<8; side++) {
        int n1 = side-4;
        int n2 = (n1+1)%4;
        if(gel->NodeIndex(n1) == gel->NodeIndex(n2)) {
            if(side != 6) DebugStop();
            iscollapsed = true;
            gel->SetRefPattern(refpattern[side]);
        }
    }
    gel->Divide(subels);
    if(0 && iscollapsed && !(gel->IsLinearMapping())) {
        ChangeToBlend(subels[0]);
        ChangeToBlend(subels[1]);
    }
    int ns = subels.size();
    for(int is = 0; is<ns; is++) {
        if(!IsCounterClockwise(subels[is])) DebugStop();
    }
}

/// @brief Substitute the trailing edge quadrilateral elements with colapsed quadrilateral elements with or without quarterpoint elements  
void Changetrailingedgeelements(TPZGeoMesh *gmesh)
{
    int64_t nel = gmesh->NElements();
    // identify the point element and trailing edge node
    int64_t singular = -1;
    TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
    TPZGeoElSide trailingedge(trailingedge_element);
    singular = trailingedge_element->NodeIndex(0);
    if(singular == -1) DebugStop();
    TPZStack<TPZGeoElSide> neighbours;
    // create nodes along the line from the singular node to the edge nodes
    TPZManVector<REAL,3> x0(3);
    gmesh->NodeVec()[singular].GetCoordinates(x0);
    // for all neighbours of the trailing edge element
    for(TPZGeoElSide neigh = trailingedge.Neighbour(); neigh != trailingedge; neigh = neigh.Neighbour()) {
        TPZGeoEl *gel = neigh.Element();
        if(!gel) DebugStop();
        if(gel->HasSubElement()) continue;
        neighbours.Push(neigh);
    }
    // create the new elements
    std::set<int64_t> sbfem_created;
    int numblend = 0;
    for(int64_t el = 0; el < neighbours.size(); el++) {
        TPZGeoEl *gel = neighbours[el].Element();
        if(!gel) DebugStop();
        if(gel->Dimension() != 2) continue;
        int nn = gel->NNodes();
        int64_t singular = gel->NodeIndex(neighbours[el].Side());

        // create collapsed elements splitting each quadrilateral element into 2 quads
        for(int in = 1; in < nn-1; in++) {
            int side1d = -1;
            if(in == 1) side1d = neighbours[el].Side()+4;
            else if(in == 2) side1d = neighbours[el].Side()+3;
            else DebugStop();
            if(side1d < 4) side1d += 4;
            TPZGeoElSide gelside(neighbours[el].Element(),side1d);
            bool blend = false;
            if(gelside.HasNeighbour(profilemat)) blend = true;
            // if the element is a blend element, create a linear element representing a high permeability element in H(div) simulations
            if(blend) {
                numblend++;
            }
            int locindex = (neighbours[el].Side()+in)%nn;
            int64_t edge1 = gel->NodeIndex(locindex);
            int64_t edge2 = gel->NodeIndex((locindex+1)%nn);
            int64_t newindex;
            if(nn == 4) {
                // create a geometrically consistent quad element
                newindex = CreateCollapsedQuad(singular,edge1,edge2,blend,gmesh);
                sbfem_created.insert(newindex);
            } else {
                DebugStop();
            }
        }
        if(nn > 1) gel->SetMaterialId(-gel->MaterialId());
    }
    if(numblend != 2) DebugStop();
    gmesh->BuildConnectivity();
    gmesh->BuildConnectivity();
    if(1)
    {
        std::ofstream out("gmesh_Changetrailingedgeelements.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmesh_Changetrailingedgeelements.txt");
        gmesh->Print(out2);
    }
}

void DivideTrailingEdgeNeighbours(TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator, REAL tol, TPZVec<REAL> &RefinementIndicator, TPZVec<int> &porders) 
{
    auto gmesh = cmesh_m->Reference();
    int64_t nel = gmesh->NElements();
    // compute the set of elements neighbouring the trailing edge
    std::set<int64_t> trailel;
    for(int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->MaterialId() == trailingedgemat) {
            TPZGeoElSide gelside(gel);
            for(TPZGeoElSide neighbour = gelside.Neighbour(); neighbour != gelside; neighbour = neighbour.Neighbour()) {
                TPZGeoEl *neighgel = neighbour.Element();
                if(neighgel->Reference() && neighgel->MaterialId() == volmat) {
                    trailel.insert(neighgel->Index());
                }
            }
        }
    }
    // divide the trailing edge element neighbours if necessary
    for(int64_t el : trailel) {
        TPZGeoEl *gel = gmesh->Element(el);
        TPZCompEl *cel = gel->Reference();
        int64_t index = cel->Index();
        if(ErrorEstimator[index] > tol && meshstyle != ESBFem) {
            RefinementIndicator[index] = 1.;
            TPZManVector<TPZGeoEl*> subel;
            DivideGeoEl(gel,subel);
        } else if(ErrorEstimator[index] > tol && meshstyle == ESBFem) {
            RefinementIndicator[index] = 1.;
            porders[index]++;
        }
    }
}

void DivideProfileNeighbours(TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator, REAL tol, int minlevel, TPZVec<REAL> &RefinementIndicator) 
{
    auto gmesh = cmesh_m->Reference();
    int64_t nel = gmesh->NElements();
    // compute the set of elements neighbouring the trailing edge
    std::set<int64_t> profileel;
    for(int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->MaterialId() == profilemat) {
            int nsides = gel->NSides();
            for(int is = 0; is < nsides; is++) {
                TPZGeoElSide gelside(gel,is);
                for(TPZGeoElSide neighbour = gelside.Neighbour(); neighbour != gelside; neighbour = neighbour.Neighbour()) {
                    TPZGeoEl *neighgel = neighbour.Element();
                    if(neighgel->Reference() && neighgel->MaterialId() == volmat) {
                        profileel.insert(neighgel->Index());
                    }
                }
            }
        }
    }
    // divide the trailing edge element neighbours if necessary
    for(int64_t el : profileel) {
        TPZGeoEl *gel = gmesh->Element(el);
        int level = gel->Level();
        if(level >= minlevel) continue;
        TPZCompEl *cel = gel->Reference();
        int64_t index = cel->Index();
        if(ErrorEstimator[index] > tol) {
            RefinementIndicator[index] = 1.;
            TPZManVector<TPZGeoEl*> subel;
            DivideGeoEl(gel,subel);
        }
    }
    Smoothentrailingedgeelements(cmesh_m,RefinementIndicator);
    // smoothen the neighbours of the cut elements
    // first divide the cut elements if necessary
    bool changed = true;
    while(changed) {
        changed = false;
        nel = gmesh->NElements();
        for(int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->MaterialId() != cutmat) continue;
            if(gel->HasSubElement()) continue;
            TPZGeoElSide gelside(gel);
            bool shoulddivide = false;
            for(TPZGeoElSide neighbour = gelside.Neighbour(); neighbour != gelside; neighbour = neighbour.Neighbour()) {
                TPZGeoEl *neighgel = neighbour.Element();
                if(neighgel->HasSubElement()) shoulddivide = true;
            }
            if(shoulddivide) {
                TPZManVector<TPZGeoEl*> subel;
                DivideGeoEl(gel,subel);
                changed = true;
            }
        }
    }
    // if the cut elements are divided, make sure the neighbours are divided
    changed = true;
    while(changed) {
        changed = false;
        nel = gmesh->NElements();
        for(int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->MaterialId() != cutmat) continue;
            // if the element is not divided, nothing to check
            if(!gel->HasSubElement()) continue;
            TPZGeoElSide gelside(gel);
            for(TPZGeoElSide neighbour = gelside.Neighbour(); neighbour != gelside; neighbour = neighbour.Neighbour()) {
                TPZGeoEl *neighgel = neighbour.Element();
                if(!neighgel->HasSubElement()) {
                    TPZManVector<TPZGeoEl*> subel;
                    DivideGeoEl(neighgel,subel);
                    changed = true;
                }
            }
        }
    }

}

/// @brief Adjust the geometric mesh to create blend elements only for neighbours of the profile
void AdjustGeometry(TPZGeoMesh *gmesh) {
    int64_t nel = gmesh->NElements();
    // make all 2d elements linear
    for(int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->HasSubElement()) continue;
        if(gel->MaterialId() == volmat) {
            bool linear = gel->IsLinearMapping();
            if(linear == true) {
                continue;
            }
            ChangeToLinearQuad(gel);
        }
    }
    // divide the profile elements if they have a neighbour element that is divided
    bool changed = true;
    while(changed) {
        changed = false;
        for(int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->MaterialId() != profilemat) continue;
            if(gel->HasSubElement()) continue;
            TPZGeoElSide gelside(gel);
            for(auto neighbour = gelside.Neighbour(); neighbour != gelside; neighbour++) {
                TPZGeoEl *neighgel = neighbour.Element();
                if(!neighgel->HasSubElement()) continue;
                if(neighgel->NSideSubElements(neighbour.Side())== 1) continue;
                TPZManVector<TPZGeoEl*> subel;
                DivideGeoEl(gel,subel);
                changed = true;
                break;
            }
        }
    }
    // change the neighbouring elements of the profile elements to blend elements
    for(int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->MaterialId() != profilemat) continue;
        TPZGeoElSide gelside(gel);
        TPZStack<TPZGeoElSide> allneighs;
        for(auto neighbour = gelside.Neighbour(); neighbour != gelside; neighbour++) {
            allneighs.Push(neighbour);
        }
        int nneigh = allneighs.size();
        for(int in = 0; in < nneigh; in++) {
            TPZGeoEl *neighgel = allneighs[in].Element();
            if(neighgel->MaterialId() != volmat) continue;
            if(neighgel->IsLinearMapping() == false) {
                ChangeToBlend(neighgel);
            }
        }
    }
}
static void AllSubElements(TPZGeoElSide &gelside, TPZStack<TPZGeoElSide> &smallest, int dimension) {
    // we push in the smallest datastructure only element/sides that have no subelements
    TPZStack<TPZGeoElSide> locsmall;
    // locsmall will be empty if the element is not divided
    gelside.GetSubElements2(locsmall);
    for(int is = 0; is<locsmall.size(); is++) {
        if(locsmall[is].Dimension() != dimension) continue;
        if(locsmall[is].HasSubElement()) {
            // recursive call to add the subelements of the subelements
            AllSubElements(locsmall[is],smallest,dimension);
        } else {
            // we have a "leaf" subelement
            smallest.Push(locsmall[is]);
        }
    }
}

void SmoothenGeometry(TPZGeoMesh *gmesh) 
{
    bool changed = true;
    while(changed) {
        changed = false;
        int64_t nel = gmesh->NElements();
        // enforce the two on one constraint
        for(int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->MaterialId() < 0) continue;
            if(gel->Dimension() != 2) continue;
            if(gel->HasSubElement()) continue;
            bool shouldrefine = false;
            for(int is = gel->FirstSide(1); is < gel->NSides()-1; is++) {
                TPZGeoElSide gelside(gel,is);
                for(TPZGeoElSide neighbour = gelside.Neighbour(); neighbour != gelside; neighbour = neighbour.Neighbour()) {
                    TPZGeoEl *neighgel = neighbour.Element();
                    if(!neighgel->HasSubElement()) continue;
                    TPZStack<TPZGeoElSide> subel;
                    AllSubElements(neighbour,subel,1);
                    // if there are 2 subelements or less the current element satisfies the
                    // two on one rule
                    if(subel.size() <= 2) continue;
                    shouldrefine = true;
                    break;
                }
                if(shouldrefine) break;
            }
            if(shouldrefine) {
                TPZManVector<TPZGeoEl *> subel;
                DivideGeoEl(gel,subel);
                changed = true;
            }
        }
        // divide an element that has more than half of its neighbours divided
        nel = gmesh->NElements();
        for(int el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->MaterialId() < 0) continue;
            if(gel->HasSubElement()) continue;
            int nside1 = gel->NSides(1);
            int nneighdiv = 0;
            for(int is = gel->FirstSide(1); is < gel->FirstSide(2); is++) {
                TPZGeoElSide gelside(gel,is);
                bool found = false;
                for(TPZGeoElSide neighbour = gelside.Neighbour(); neighbour != gelside; neighbour = neighbour.Neighbour()) {
                    TPZGeoEl *neighgel = neighbour.Element();
                    if(neighgel->MaterialId() < 0) continue;
                    int neighside = neighbour.Side();
                    if(neighgel->HasSubElement() && neighgel->NSideSubElements(neighside) > 1) found = true;
                }
                if(found) nneighdiv++;
            }
            if(nneighdiv > nside1/2) {
                TPZManVector<TPZGeoEl *> subel;
                DivideGeoEl(gel,subel);
                changed = true;
            }
        }
    }
    // adjust the one dimensional elements
    // they should have no divided neighbours
    changed = true;
    while(changed) {
        changed = false;
        int64_t nel = gmesh->NElements();
        for(int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->HasSubElement()) continue;
            if(gel->Dimension() != 1) continue;
            //Percorrer os Neighbors de gel e encontar os volumetricos que tem subelementos
            TPZGeoElSide gelside(gel); 
            auto Neighbor = gelside.Neighbour();
            for(; Neighbor != gelside; Neighbor++) {
                if(!Neighbor.Element()) DebugStop();
                if(Neighbor.HasSubElement() && Neighbor.NSubElements() > 1) {
                    TPZManVector<TPZGeoEl *> subel;
                    DivideGeoEl(gel,subel);
                    changed = true;
                    break;
                }
            }
        }
    }
}
  
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
//Define a função do tipo void chamada PrintResults, que recebe como parâmetros TPZLinearAnalysis &an e  TPZCompMesh *cmesh
{
 
    std::cout << "--------- Post Process ---------" << std::endl;
    //printa para o usuário "--------- Post Process ---------" indicando o início da fase de pós-processamento. 
    TPZSimpleTimer postProc("Post processing time");
    //declara uma variável chamada postProc, do tipo TPZSimpleTimer, chamando um construtor com uma string como argumento, igual a "Post processing time".
    //inicializa um temporizador chamado postProc que será usado para medir o tempo gasto no pós-processamento.
    const std::string plotfile = "postprocess";
    //define o nome base do arquivo de saída para o pós-processamento. O nome base é "postprocess".
    constexpr int vtkRes{3};
    //define a variável do tipo inteiro denominada vtkRes, do tipo constexpr, que significa que é uma expressão constante, ou seja,  vtkRes é um valor constante e não pode ser alterado. Ainda, {0} indica o valor associado a essa constante, e portanto não será alterado, com valor determinado na hora de compilação.
    //define a resolução para o formato de arquivo VTK. Neste caso, a resolução é definida como 0, o que geralmente significa que a resolução será automática.
    TPZVec<std::string> fields = {
        "Pressure",
        "Flux"
    };
    //nesse conjunto de linhas de código, temos que TPZVec é uma estrutura do tipo vetor que contém como argumento uma variável chamda "fields" que é uma lista de strings, que, pelo que se chamam, são relacionadas à pressão e ao fluxo.
    //cria um vetor de strings chamado fields que contém os nomes dos campos que serão pós-processados. Neste caso, os campos incluem "Pressure" (pressão) e "Flux" (fluxo). Esses campos representam propriedades do problema que desejamos visualizar após a simulação.
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    //essa linha de código declara uma variável chamada vtk do tipo auto, o que significa que o compilador irá deduzir o tipo que ela terá a depender do que ela é igual. No caso, ela é igual a função TPZVTKGenerator, de parâmetros cmesh, fields, plotfile, vtkRes.
    //cria um objeto vtk da classe TPZVTKGenerator, que é usado para gerar arquivos VTK a partir dos dados da malha computacional cmesh. Os argumentos passados para o construtor incluem a malha computacional, os campos a serem pós-processados, o nome base do arquivo de saída (plotfile) e a resolução VTK (vtkRes).
    vtk.SetNThreads(0);
    //define o número de threads a serem usadas durante o pós-processamento. A variável global_nthread provavelmente contém o número desejado de threads.
    vtk.Do();
    //inicia o processo de geração dos arquivos VTK. Esta função gera arquivos de saída contendo informações sobre os campos especificados na malha computacional.
    std::cout << "Total time = " << postProc.ReturnTimeDouble()/1000. << " s" << std::endl;
    //imprime o tempo gasto no pós-processamento, convertido para segundos.
    
    return;
    //a função é concluída e retorna.
}

/// @brief Adjust the system of equations: build the second collum of the rhs for Beta=1 and remove the Beta equations from the structural matrix 
/// @return Adjusted rhs and structural matrices: TPZFMatrix<STATE> &rhs, TPZSYsmpMatrix<STATE> &spmat.
void AdjustSystemofEquations (TPZLinearAnalysis &an, TPZCompMesh *cmesh, int64_t &newcon)
{
    ///Geting the structural matrix and the rhs
    cmesh->ConnectVec()[newcon].Print(*cmesh);
    TPZConnect &c = cmesh->ConnectVec()[newcon];
    auto seq = c.SequenceNumber();
    auto eqbeta = cmesh->Block().Position(seq);
    std::cout << "eqbeta " << eqbeta << std::endl;
    auto solver = an.Solver();
    TPZMatrixSolver<STATE>* matsolv = dynamic_cast<TPZMatrixSolver<STATE>*>(solver);
    if(!matsolv) DebugStop();
    auto mat = matsolv->Matrix();
#if defined(USING_MKL) || defined(USING_EIGEN)
    TPZSYsmpMatrix<STATE>* spmat = dynamic_cast<TPZSYsmpMatrix<STATE>*> (mat.operator->());
    if(!spmat) DebugStop();
#else
    TPZSkylMatrix<STATE> *spmat = dynamic_cast<TPZSkylMatrix<STATE>*>(mat.operator->());
    if(!spmat) DebugStop();
#endif
    TPZFMatrix<STATE> &rhs = an.Rhs(); 

    ///Variables
    auto nrows =  spmat->Rows();
    auto ncols =  spmat->Cols();
    rhs.Resize(nrows,2);
    for(int64_t row = 0; row < nrows; row++)
    {
       rhs(row,1) = 0.0; 
    }
 
    for (int64_t row = 0; row < nrows; row++) 
    {
        if (row == eqbeta)
        {
            rhs(row,1) = 1.0;
            rhs(row,0) = 0.;
            spmat->PutVal(row,row,1.);
        } else
        {
            rhs(row,1) = -spmat->GetVal(row,eqbeta);
            spmat->PutVal(row,eqbeta,0.);    
        }
    }
}

/// @brief Get the trailing edge elements
/// @param gmeshfilename 
/// @return The GeoElSides of the trailing edge: TPZGeoElSide &Geosideminus (T.E. plus) and TPZGeoElSide &Geosideplus (T.E. minus).
void GetTrailingEdgeElements(TPZGeoMesh *gmesh, TPZGeoElSide &Geosideminus, TPZGeoElSide &Geosideplus)

{
    int64_t nel = gmesh->NElements();
    TPZStack<TPZGeoEl *, 2> NeighborNacaElements;
    TPZStack<int, 2> NeighborNacaElementsSide;
    TPZGeoElSide TrailingSide;

    // Put all neighbours of the trailing edge element with matid profilemat in the stack NeighborNacaElements
    TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
    TrailingSide = TPZGeoElSide(trailingedge_element);
    int64_t trailingnode = trailingedge_element->NodeIndex(0);
    auto Neighbor = TrailingSide.Neighbour();
    for(; Neighbor != TrailingSide; Neighbor++)
    {
        if(!Neighbor.Element()) DebugStop();
        if(Neighbor.Element()->HasSubElement()) continue;
        if (Neighbor.Element()->MaterialId() == profilemat)  
        {
            NeighborNacaElements.Push (Neighbor.Element());
        }  
    }
    // NeighborNacaElements are one dimensional elements with profilemat
    if(NeighborNacaElements.size() != 2) {
        // debugging code, ends with DebugStop
        auto Neighbor = TrailingSide.Neighbour();
        for(; Neighbor != TrailingSide; Neighbor++)
        {
            if(!Neighbor.Element()) DebugStop();
            if(Neighbor.Element()->HasSubElement()) continue;
            std::cout << "Neighbor.Element()->MaterialId() " << Neighbor.Element()->MaterialId() << std::endl;
            if (Neighbor.Element()->MaterialId() == profilemat)  
            {
                NeighborNacaElements.Push (Neighbor.Element());
            }  
        }

        DebugStop();
    }

    // make sure the NeighborNacaElements are ordered counter clockwise
    if(NeighborNacaElements[0]->NodeIndex(1) == trailingnode) {
        TPZGeoEl *tmp = NeighborNacaElements[0];
        NeighborNacaElements[0] = NeighborNacaElements[1];
        NeighborNacaElements[1] = tmp;
    }
    // the first element is on the top of the profile
    // the second element is on the bottom of the profile
    if(NeighborNacaElements[0]->NodeIndex(0) != trailingnode) DebugStop();
    if(NeighborNacaElements[1]->NodeIndex(1) != trailingnode) DebugStop();
    if(! TrailingSide) DebugStop();


    ////Fazer um loop para cada um dos filhos dos dois elementos do vetor NeighborNacaElements e 
    //verificar qual dos filhos tem o Trailing Edge como um neighbor

    // Find the volumetric elements neighbouring the trailing edge
    {
        // take the profile element
        auto gel = NeighborNacaElements[0];
        TPZGeoElSide gelside = TPZGeoElSide(gel);
        TPZGeoElSide Neighbor = gelside.Neighbour();
        for(; Neighbor != gelside; Neighbor++) 
        {
            if (Neighbor.Element()->MaterialId() != volmat && Neighbor.Element()->MaterialId() != sbfem_domain) continue;
            if(Neighbor.Element()->HasSubElement()) continue;
            // neighgel is neighbouring volumetric element
            Geosideminus = Neighbor;
        }
    }
    {
        auto gel = NeighborNacaElements[1];
        TPZGeoElSide gelside = TPZGeoElSide(gel);
        TPZGeoElSide Neighbor = gelside.Neighbour();
        for(; Neighbor != gelside; Neighbor++) 
        {
            if (Neighbor.Element()->MaterialId() != volmat && Neighbor.Element()->MaterialId() != sbfem_domain) continue;
            if(Neighbor.Element()->HasSubElement()) continue;
            Geosideplus = Neighbor;
        }
    }
    if(!Geosideminus || !Geosideplus) DebugStop();
    TPZCompMesh *cmesh = gmesh->Reference();
    TPZMultiphysicsCompMesh *cmesh_m = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);
    if(meshstyle == ESBFem && !cmesh_m) {
        TPZStack<TPZCompEl *> compels;
        sbfem_groupH1->GetCompElList(compels);
        for(auto compel : compels) {
            TPZGeoEl *gel = compel->Reference();
            if(gel == Geosideminus.Element()) Geosideminus.Element()->SetReference(compel);
            if(gel == Geosideplus.Element()) Geosideplus.Element()->SetReference(compel);
        }
    } else if(meshstyle == ESBFem) {
        TPZStack<TPZCompEl *> compels;
        sbfem_groupHdiv->GetCompElList(compels);
        for(auto compel : compels) {
            TPZGeoEl *gel = compel->Reference();
            if(gel == Geosideminus.Element()) Geosideminus.Element()->SetReference(compel);
            if(gel == Geosideplus.Element()) Geosideplus.Element()->SetReference(compel);
        }
    }
}

#include "pzvec_extras.h"
/// @brief Evaluate the gradient of the H1 potential solution at the trailing edge points: T.E. plus and T.E. minus
/// @return The gradient of the solution at the trailling edge points: TPZManVector<REAL,3> gradminus, TPZManVector<REAL,3> gradminus.
void EvaluateSolutionGradientsH1(TPZGeoMesh *gmesh, TPZManVector<REAL,3> &gradminus, TPZManVector<REAL,3> &gradplus)
{
    ///Variables
    const int dim = gmesh->Dimension();
    TPZManVector<REAL,3> qsi(dim,0.);

    TPZGeoElSide Geosideminus;
    TPZGeoElSide Geosideplus;
    GetTrailingEdgeElements(gmesh, Geosideminus, Geosideplus);
    TPZCompEl*Compelminus = Geosideminus.Element()->Reference();
    TPZCompEl*Compelplus = Geosideplus.Element()->Reference();

    //Compute the gradient of the potential solution for the Compelminus
    auto tr = Geosideminus.Element()->SideToSideTransform(Geosideminus.Side(),Geosideminus.Element()->NSides()-1);
    TPZManVector<REAL,1> qsiminus = {-1.};
    tr.Apply(qsiminus,qsi);
    {
        TPZManVector<REAL,3> x1(3),x2(3);
        Geosideminus.Element()->X(qsi,x1);
        TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
        auto trailnode = trailingedge_element->NodePtr(0);
        trailnode->GetCoordinates(x2);
        if(dist(x1,x2) > 1.e-10) {
            for(int i=0; i<4; i++) std::cout << Geosideminus.Element()->NodeIndex(i) << " ";
            std::cout << std::endl;
            DebugStop();
        }
    }
    if(meshstyle != ETraditional) {
        qsiminus[0] = -1.+1.e-2;
        tr.Apply(qsiminus,qsi);
    }

    Compelminus->Solution(qsi,2,gradminus); 

    //Compute the gradient of the potential solution for the Compelplus
    tr = Geosideplus.Element()->SideToSideTransform(Geosideplus.Side(),Geosideplus.Element()->NSides()-1);
    TPZManVector<REAL,1> qsiplus = {1.};
    tr.Apply(qsiplus,qsi);
    {
        TPZManVector<REAL,3> x1(3),x2(3);
        Geosideplus.Element()->X(qsi,x1);
        TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
        auto trailnode = trailingedge_element->NodePtr(0);
        trailnode->GetCoordinates(x2);
        if(dist(x1,x2) > 1.e-10) DebugStop();
    }
    if(meshstyle != ETraditional) {
        qsiplus[0] = 1.-1.e-2;
        tr.Apply(qsiplus,qsi);
    }

    Compelplus->Solution(qsi,2,gradplus); 
}

/// @brief Evaluate the H(div) Flux solution at the trailing edge points: T.E. plus and T.E. minus
/// @return The solution at the trailling edge points: TPZManVector<REAL,3> u_minus, TPZManVector<REAL,3> u_plus.
void EvaluateSolutionHDiv(TPZGeoMesh *gmesh, TPZMultiphysicsCompMesh *cmesh_m, TPZManVector<REAL,3> &u_minus, TPZManVector<REAL,3> &u_plus)

{
    cmesh_m->LoadReferences();
    ///Variables
    const int dim = gmesh->Dimension();
    TPZManVector<REAL,3> qsi(dim,0.);

    TPZGeoElSide Geosideminus;
    TPZGeoElSide Geosideplus;
    GetTrailingEdgeElements(gmesh, Geosideminus, Geosideplus);
    TPZCompEl*Compelminus = Geosideminus.Element()->Reference();
    TPZCompEl*Compelplus = Geosideplus.Element()->Reference();

    //Compute the Flus solution for the Compelminus
    auto tr = Geosideminus.Element()->SideToSideTransform(Geosideminus.Side(),Geosideminus.Element()->NSides()-1);
        TPZManVector<REAL,1> qsiminus = {-1.};
    tr.Apply(qsiminus,qsi);
    {
        TPZManVector<REAL,3> x1(3),x2(3);
        Geosideminus.Element()->X(qsi,x1);
        TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
        auto trailnode = trailingedge_element->NodePtr(0);
        trailnode->GetCoordinates(x2);
        if(dist(x1,x2) > 1.e-10) {
            for(int i=0; i<4; i++) std::cout << Geosideminus.Element()->NodeIndex(i) << " ";
            std::cout << std::endl;
            DebugStop();
        }
    }
    if(meshstyle != ETraditional) {
        qsiminus[0] = -1.+1.e-2;
        tr.Apply(qsiminus,qsi);
    }

    Compelminus->Solution(qsi,1,u_minus); 

    //Compute the Flux solution for the Compelplus
    tr = Geosideplus.Element()->SideToSideTransform(Geosideplus.Side(),Geosideplus.Element()->NSides()-1);
    TPZManVector<REAL,1> qsiplus = {1.};
    tr.Apply(qsiplus,qsi);
    {
        TPZManVector<REAL,3> x1(3),x2(3);
        Geosideplus.Element()->X(qsi,x1);
        TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
        auto trailnode = trailingedge_element->NodePtr(0);
        trailnode->GetCoordinates(x2);
        if(dist(x1,x2) > 1.e-10) DebugStop();
    }
    if(meshstyle != ETraditional) {
        qsiplus[0] = 1.-1.e-2;
        tr.Apply(qsiplus,qsi);
    }
    Compelplus->Solution(qsi,1,u_plus); 
}

/// @brief Evaluate the circulation constant "Kappa" for the flux potential solution loaded into the TPZCompMesh "cmesh"
/// @return The circulation constant: REAL &Kappa.
void EvaluateCirculationH1(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, int matid, REAL &circulation)

{
    circulation = 0.0;
    cmesh->LoadReferences();
    TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh->ElementVec();
    int64_t nelem = cmesh->NElements();
    for (int iel = 0; iel < nelem; iel++) {
        TPZCompEl *el = elementvec[iel];
        TPZGeoEl *gel = el->Reference();
        if (gel->MaterialId() != matid) continue;
        // NUMERICAL INTEGRATION DATA;
        TPZGeoElSide gelside(gel);
        auto intrule = gelside.CreateIntegrationRule(20);
        int intrulepoints = intrule->NPoints();
        TPZManVector<REAL,4> intpointtemp(1,0.);
        TPZManVector<REAL,4> intpointvol(2,0.);
        REAL weight = 0.;
        TPZFMatrix<REAL> jac, axe, jacInv;
        TPZManVector<REAL,3> x(3,0.);
        REAL detJac;

        TPZGeoElSide neighbor = gelside.Neighbour();
        for(; neighbor != gel; neighbor++)
        {
            int matid = neighbor.Element()->MaterialId();
            if (matid == volmat || matid == sbfem_domain) 
            {
                break;
            }   
        }

        auto neighel = neighbor.Element();
        auto tr = neighel->SideToSideTransform(neighbor.Side(),neighel->NSides()-1);  
        for(int int_ind = 0; int_ind < intrulepoints; ++int_ind)
        {
            intrule->Point(int_ind,intpointtemp,weight);
            neighbor.Jacobian(intpointtemp, jac, axe, detJac , jacInv);
            weight *= fabs(detJac);
            tr.Apply(intpointtemp,intpointvol);
            auto ref = neighel->Reference();
            TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(ref);
            TPZManVector<STATE,3> flux(3,0.);
            msp ->Solution(intpointvol,7,flux);

            REAL ut = flux[0]*axe(0,0) + flux[1]*axe(0,1);
            circulation += ut*weight;
        }//loop over integration points
    }//loop over elements
}

/// @brief Evaluate the circulation constant "Kappa" for the flux solution loaded into the TPZMultiphysicsCompMesh "cmesh_m"
/// @return The circulation constant: REAL &Kappa.
void EvaluateCirculationHDiv(TPZGeoMesh *gmesh, TPZMultiphysicsCompMesh *cmesh_m, int matid, REAL &circulation)

{
    circulation = 0.0;
    cmesh_m->LoadReferences();
    TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh_m->ElementVec();
    int64_t nelem = cmesh_m->NElements();
    for (int iel = 0; iel < nelem; iel++) {
        TPZCompEl *el = elementvec[iel];
        if(!el) continue;
        TPZGeoEl *gel = el->Reference();
        if(!gel) continue;
        if (gel->MaterialId() != matid) continue;

        //NUMERICAL INTEGRATION DATA;
        TPZGeoElSide gelside(gel);
        auto intrule = gelside.CreateIntegrationRule(20);
        int intrulepoints = intrule->NPoints();
        TPZManVector<REAL,4> intpointtemp(1,0.);
        TPZManVector<REAL,4> intpointvol(2,0.);
        REAL weight = 0.;
        TPZFMatrix<REAL> jac, axe, jacInv;
        REAL detJac;

        TPZGeoElSide neighbor = gelside.Neighbour();
        for(; neighbor != gel; neighbor++)
        {
            if (neighbor.Element()->MaterialId() == volmat || neighbor.Element()->MaterialId() == sbfem_domain)
            {
                break;
            }   
        }

        auto neighel = neighbor.Element();
        if(neighel->MaterialId() != sbfem_domain && neighel->MaterialId() != volmat) DebugStop();
        auto tr = neighel->SideToSideTransform(neighbor.Side(),neighel->NSides()-1);  
        for(int int_ind = 0; int_ind < intrulepoints; ++int_ind)
        {
            intrule->Point(int_ind,intpointtemp,weight);
            neighbor.Jacobian(intpointtemp, jac, axe, detJac , jacInv);
            weight *= fabs(detJac);
            tr.Apply(intpointtemp,intpointvol);
            auto ref = neighel->Reference();
            TPZMultiphysicsElement *msp  = dynamic_cast <TPZMultiphysicsElement *>(ref);
            TPZManVector<STATE,3> sol;
            ref->Solution(intpointvol,1,sol);
            REAL ut = sol[0]*axe(0,0) + sol[1]*axe(0,1);
            circulation += ut*weight;
        }//loop over integration points
    }//loop over elements
}

/// @brief Compute the Beta constant
/// @return The value of Beta: REAL Beta.
void ComputeBetaH1(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, TPZFMatrix<STATE> &phi_0, TPZFMatrix<STATE> &phi_1, REAL &Beta)

{
    AddSBFemVolumeElements();
    cmesh->LoadSolution(phi_0);
    TPZManVector<REAL,3> grad_phi0_plus(3,0.);
    TPZManVector<REAL,3> grad_phi0_minus(3,0.);
    EvaluateSolutionGradientsH1(gmesh, grad_phi0_plus, grad_phi0_minus);

    cmesh->LoadSolution(phi_1);
    TPZManVector<REAL,3> grad_phi1_plus(3,0.);
    TPZManVector<REAL,3> grad_phi1_minus(3,0.);
    EvaluateSolutionGradientsH1(gmesh, grad_phi1_plus, grad_phi1_minus);

    REAL A = grad_phi1_plus[0]*grad_phi1_plus[0]+grad_phi1_plus[1]*grad_phi1_plus[1]-grad_phi1_minus[0]*grad_phi1_minus[0]-grad_phi1_minus[1]*grad_phi1_minus[1];
    REAL B = 2.0*(grad_phi0_plus[0]*grad_phi1_plus[0]+grad_phi0_plus[1]*grad_phi1_plus[1]-grad_phi0_minus[0]*grad_phi1_minus[0]-grad_phi0_minus[1]*grad_phi1_minus[1]);
    REAL C = grad_phi0_plus[0]*grad_phi0_plus[0]+grad_phi0_plus[1]*grad_phi0_plus[1]-grad_phi0_minus[0]*grad_phi0_minus[0]-grad_phi0_minus[1]*grad_phi0_minus[1];
    REAL Delta = B*B-4.0*A*C;

    if (abs((-B+sqrt(Delta))/(2*A))<abs((-B-sqrt(Delta))/(2*A)))
    {
        Beta = (-B+sqrt(Delta))/(2*A);
    }
    else 
    {
        Beta = (-B-sqrt(Delta))/(2*A);
    }
    HideSBFemVolumeElements();
}

/// @brief Compute the Kappa constant
/// @return The value of Kappa: REAL Kappa.
void ComputeBetaHDiv(TPZGeoMesh *gmesh, TPZMultiphysicsCompMesh *cmesh_m, TPZFMatrix<STATE> &u_0, TPZFMatrix<STATE> &u_1, REAL &Beta)

{
    AddSBFemVolumeElements();
    cmesh_m->LoadSolution(u_0);
    cmesh_m->TransferMultiphysicsSolution();
    TPZManVector<REAL,3> u0_plus(3,0.);
    TPZManVector<REAL,3> u0_minus(3,0.);
    EvaluateSolutionHDiv(gmesh, cmesh_m, u0_plus, u0_minus);

    cmesh_m->LoadSolution(u_1);
    cmesh_m->TransferMultiphysicsSolution();
    TPZManVector<REAL,3> u1_plus(3,0.);
    TPZManVector<REAL,3> u1_minus(3,0.);
    EvaluateSolutionHDiv(gmesh, cmesh_m, u1_plus, u1_minus);

    REAL A = u1_plus[0]*u1_plus[0]+u1_plus[1]*u1_plus[1]-u1_minus[0]*u1_minus[0]-u1_minus[1]*u1_minus[1];
    REAL B = 2.0*(u0_plus[0]*u1_plus[0]+u0_plus[1]*u1_plus[1]-u0_minus[0]*u1_minus[0]-u0_minus[1]*u1_minus[1]);
    REAL C = u0_plus[0]*u0_plus[0]+u0_plus[1]*u0_plus[1]-u0_minus[0]*u0_minus[0]-u0_minus[1]*u0_minus[1];
    REAL Delta = B*B-4.0*A*C;

    if (abs((-B+sqrt(Delta))/(2*A))<abs((-B-sqrt(Delta))/(2*A)))
    {
        Beta = (-B+sqrt(Delta))/(2*A);
    }
    else 
    {
        Beta = (-B-sqrt(Delta))/(2*A);
    }
    HideSBFemVolumeElements();
}

/// @brief returns true if the element is counter clockwise
bool IsCounterClockwise(TPZGeoEl *gel) {
    if(gel->Dimension() != 2) DebugStop();
    TPZGeoElSide gelside(gel);
    TPZManVector<REAL,3> qsi(2),ax1(3),ax2(3),ax3(3);
    TPZFNMatrix<6,REAL> gradx(3,2);
    gelside.CenterPoint(qsi);
    gel->GradX(qsi,gradx);
    for(int i=0; i<3; i++) {
        ax1[i] = gradx(i,0);
        ax2[i] = gradx(i,1);
    }
    if(!IsZero(ax1[2]) || !IsZero(ax2[2])) DebugStop();
    Cross(ax1,ax2,ax3);
    if(ax3[2] > 0.) return true;
    return false;
}



/// @brief make the elements counter clockwise
void MakeCounterClockwise(TPZGeoEl *gel) {
    if(IsCounterClockwise(gel)) return;
    if(gel->Type() != EQuadrilateral) DebugStop();
    TPZManVector<TPZGeoElSide,9> gelsides(gel->NSides());
    for(int is = 0; is<gel->NSides(); is++) {
        TPZGeoElSide thisside(gel,is);
        TPZGeoElSide neigh = thisside.Neighbour();
        if(neigh!= thisside) gelsides[is] = neigh;
    }
    gel->RemoveConnectivities();
    // flip nodes 1 and 3
    {
        int64_t node1 = gel->NodeIndex(1);
        int64_t node3 = gel->NodeIndex(3);
        gel->SetNodeIndex(1,node3);
        gel->SetNodeIndex(3,node1);
    }
    // connectivities 0->0 1->3 2->2 3->1 4->7 5->6 6->5 7->4 8->8
    for(int is = 0; is<9; is++) {
        int mapside[] = {0,3,2,1,7,6,5,4,8};
        TPZGeoElSide gelside(gel,is);
        if(gelsides[is]) {
            gelside.SetConnectivity(gelsides[is]);
        } else {
            gelside.SetConnectivity(gelside);
        }
    }
#ifdef PZDEBUG
    if(IsCounterClockwise(gel) == false) DebugStop();
    for(int is = 0; is<gel->NSides(); is++) {
        std::set<int64_t> nodes;
        for(int i=0; i<gel->NSideNodes(is); i++) nodes.insert(gel->SideNodeIndex(is,i));
        TPZGeoElSide gelside(gel,is);
        for(auto neigh = gelside.Neighbour(); neigh != gelside; neigh++) {
            for(int in = 0; in<neigh.NSideNodes(); in++) {
                if(nodes.find(neigh.SideNodeIndex(in)) == nodes.end()) {
                    DebugStop();
                }
            }
        }
    }
#endif
}



TPZAutoPointer<TPZRefPattern> CreateRefPattern() {
    char buf[] =
    "6     3  "
    "1010       UnifQua	"
    "-1.    -1.     0. "
    " 1.    -1.     0. "
    " 1.     1.     0. "
    "-1.     1.     0.	"
    " -1.    0.     0. "
    " 1.     0.     0. "
    " 3     4     0     1     2     3 "
    " 3     4     0     1     5     4 "
    " 3     4     4     5     2     3 ";
    std::istringstream str(buf);
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(str);
    return refpat;
}


/// @brief Change the elements the touch the trailing edge to quarterpoint elements
void CreateQuarterPointElements(TPZGeoMesh *gmesh) {
    TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
    TPZGeoElSide trailingedge(trailingedge_element);
    std::map<TPZGeoEl *,int> connected_elements;
    for(TPZGeoElSide neigh = trailingedge.Neighbour(); neigh != trailingedge; neigh++) {
        TPZGeoEl *neighgel = neigh.Element();
        if(neighgel->HasSubElement()) continue;
        if(neighgel->Dimension() == 2) {
            if(neighgel->NodeIndex(2) == neighgel->NodeIndex(3)) {
                connected_elements[neighgel] = 6;
            }
        } else if(neighgel->Dimension() == 1) {
            connected_elements[neighgel] = neigh.Side();
        }

    }
    for(auto it : connected_elements) {
        TPZGeoEl *gel = it.first;
        int side = it.second;
        gel = TPZChangeEl::ChangeToQuadratic(gel->Mesh(), gel->Index());
        TPZChangeEl::ChangeToQuarterPoint(gel->Mesh(), gel->Index(),side);
    }
}

/// @brief remove computational elements generated for SBFem simulation
void RemoveSBFemElements(TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *cmesh_m, TPZGeoMesh *gmesh)
{
    if(sbfem_groupHdiv) {
        sbfem_groupHdiv->Unwrap();
        sbfem_groupHdiv = 0;
    }
    if(sbfem_groupH1) {
        sbfem_groupH1->Unwrap();
        sbfem_groupH1 = 0;
    }
    auto Remove= [gmesh] (TPZCompMesh *cmesh) {
            int64_t nel = cmesh->NElements();
            std::cout << "gmesh nelements " << gmesh->NElements() << std::endl;
            for(int64_t iel = 0; iel < nel; iel++) 
            {
                TPZCompEl *cel = cmesh->Element(iel);
                if(!cel) continue;
                int64_t index = cel->ReferenceIndex();
                if(index < 0) continue;
                TPZGeoEl *gel = cel->Reference();
                int matid = gel->MaterialId();
                if(index >= gmesh->NElements()) {
                    std::cout << "gel index " << index << " matid " << gel->MaterialId() << " is out of range " << " celindex " << cel->Index() << "\n";
                    int64_t celindex = cel->Index();
                    delete cel;
                    cel = cmesh->Element(celindex);
                    if(cel) DebugStop();
                }
            }

        };
    Remove(cmesh);
    Remove(cmesh_m);
    auto meshvec = cmesh_m->MeshVector();
    for(auto mesh : meshvec) {
        Remove(mesh);
    }
}



/// @brief print the geometry of the trailing edge elements
void PrintTrailingEdgeElements(TPZGeoMesh *gmesh) {
    TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
    TPZGeoElSide trailingedge(trailingedge_element);
    std::set<TPZGeoEl *> connected_elements;
    for(TPZGeoElSide neigh = trailingedge.Neighbour(); neigh != trailingedge; neigh++) {
        TPZGeoEl *neighgel = neigh.Element();
        if(neighgel->HasSubElement()) continue;
        if(neighgel->Dimension() == 2) {
            if(neighgel->NodeIndex(2) == neighgel->NodeIndex(3)) {
                connected_elements.insert(neighgel);
            }
        } else if(neighgel->Dimension() == 1) {
            connected_elements.insert(neighgel);
        }

    }
    for(auto it : connected_elements) {
        TPZGeoEl *gel = it;
        std::cout << "Element " << gel->Index() << std::endl;
        int nnodes = gel->NNodes();
        for(int i=0; i<nnodes; i++) {
            TPZManVector<REAL,3> x(3);
            gel->NodePtr(i)->GetCoordinates(x);
            std::cout << "Node " << i << " " << x << std::endl;
        }
    }
}


/// @brief Change the elements the touch the trailing edge to H1 SBFEM elements
void CreateH1SBFEMelements(TPZCompMesh *cmesh) {
    cmesh->LoadReferences();
    TPZMaterialT<STATE> *mat2d = dynamic_cast<TPZMaterialT<STATE> *>(cmesh->FindMaterial(volmat));
    if(!mat2d) DebugStop();
    TPZDarcyFlow *darcy = dynamic_cast<TPZDarcyFlow *>(mat2d);
    if(!darcy) DebugStop();
    TPZDarcyFlow *darcysbfem = new TPZDarcyFlow(*darcy);
    darcysbfem->SetId(sbfem_domain);
    cmesh->InsertMaterialObject(darcysbfem);
    TPZGeoMesh *gmesh = cmesh->Reference();
    TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
    TPZGeoElSide trailingedge(trailingedge_element);
    std::set<TPZGeoEl *> connected_elements;
    for(TPZGeoElSide neigh = trailingedge.Neighbour(); neigh != trailingedge; neigh++) {
        TPZGeoEl *neighgel = neigh.Element();
        if(neighgel->HasSubElement()) continue;
        if(neighgel->MaterialId() < 0) continue;
        if(connected_elements.find(neighgel) != connected_elements.end()) continue;
        if(neighgel->Dimension() == 2) {
            if(neighgel->NodeIndex(2) == neighgel->NodeIndex(3)) {
                TPZCompEl *cel = neighgel->Reference();
                if(cel) DebugStop();
                connected_elements.insert(neighgel);
                TPZGeoElSide gelside(neighgel,4);
                TPZGeoElSide neigh = gelside.HasNeighbour(sbfem_skeleton);
                if(!neigh) {
                    DebugStop();
                }
                if(neigh.Element()->Reference() == 0) {
                    DebugStop();
                }
            } else {
                DebugStop();
            }
        }
    }
    TPZSBFemElementGroup *elgr = new TPZSBFemElementGroup(*cmesh);
    sbfem_groupH1 = elgr;
    TPZStack<TPZSBFemVolume *> sbvols;
    for(auto it : connected_elements) {
        TPZGeoEl *gel = it;
        if(gel->MaterialId() != sbfem_domain) DebugStop();
        TPZSBFemVolume *sbvol = new TPZSBFemVolume(*cmesh,gel);
        sbvols.Push(sbvol);
        TPZGeoElSide gelside(gel,4);
        TPZGeoElSide neigh = gelside.HasNeighbour(sbfem_skeleton);
        TPZGeoEl *skel = neigh.Element();
        if(skel->MaterialId() != sbfem_skeleton) DebugStop();
        // create an H1 element
        TPZCompEl *cskel = skel->Reference();
        if(!cskel) DebugStop();
        TPZInterpolatedElement *icskel = dynamic_cast<TPZInterpolatedElement *>(cskel);
        if(!cskel) DebugStop();
        int sideorder = icskel->Connect(2).Order();
        if(sideorder < SBFemOrder) icskel->SetSideOrder(2,SBFemOrder);
        std::cout << "gel " << gel->Index() << " is linear " << gel->IsLinearMapping() << std::endl;
        std::cout << "skel is linear " << skel->IsLinearMapping() << std::endl;
        if(!gel->IsLinearMapping()) {
            ChangeToLinearQuad(gel);
            // gel->Print(std::cout);
            // std::ofstream out("gmesh.txt");
            // gmesh->Print(out);
            // DebugStop();
        }
        sbvol->SetSkeleton(cskel->Index());
        elgr->AddElement(sbvol);
    }
    for(auto sbvol : sbvols) {
        sbvol->SetElementGroupIndex(elgr->Index());
    }
    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    TPZElementMatrixT<STATE> ek,ef;
    elgr->CalcStiff(ek,ef);
    if(0) {
        ek.fMat.Print("ek = ",std::cout,EMathematicaInput);
    }
}

/// @brief Change the elements the touch the trailing edge to Hdiv SBFEM elements
void CreateHdivSBFEMelements(TPZMultiphysicsCompMesh *m_cmesh, int64_t newcon) {
    TPZCompMesh *cmesh = m_cmesh;
    cmesh->LoadReferences();
    // cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(2);
    // cmesh->ApproxSpace().SetHDivFamily(HDivFamily::EHDivKernel);
    TPZGeoMesh *gmesh = cmesh->Reference();

    // insert the skeleton material
    if(cmesh->FindMaterial(sbfem_skeleton) == 0) {
        DebugStop();
    }
    TPZMixedDarcyFlow *mat1d = new TPZMixedDarcyFlow(sbfem_highperm,1);
    cmesh->InsertMaterialObject(mat1d);
    mat1d->SetConstantPermeability(1.e-6);
    TPZMaterialT<STATE> *mat2d = dynamic_cast<TPZMaterialT<STATE> *>(cmesh->FindMaterial(volmat));
    if(!mat2d) DebugStop();
    TPZMixedDarcyFlow *darcy = dynamic_cast<TPZMixedDarcyFlow *>(mat2d);
    if(!darcy) DebugStop();
    TPZMixedDarcyFlow *darcysbfem = new TPZMixedDarcyFlow(*darcy);
    darcysbfem->SetId(sbfem_domain);
    cmesh->InsertMaterialObject(darcysbfem);

    TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
    TPZGeoElSide trailingedge(trailingedge_element);
    std::set<TPZGeoEl *> connected_elements;
    // create the skeleton elements
    for(TPZGeoElSide neigh = trailingedge.Neighbour(); neigh != trailingedge; neigh++) {
        TPZGeoEl *neighgel = neigh.Element();
        if(neighgel->HasSubElement()) continue;
        if(neighgel->MaterialId() < 0) continue;
        if(connected_elements.find(neighgel) != connected_elements.end()) continue;
        if(neighgel->Dimension() == 2) {
            if(neighgel->NodeIndex(2) == neighgel->NodeIndex(3)) {
                TPZCompEl *cel = neighgel->Reference();
                if(cel) DebugStop();
                if(!neighgel->IsLinearMapping()) {
                    int64_t index = neighgel->Index();
                    ChangeToLinearQuad(neighgel);
                    neighgel = gmesh->Element(index);
                }
                connected_elements.insert(neighgel);
                TPZGeoElSide gelside(neighgel,4);
                TPZGeoElSide neighbour = gelside.HasNeighbour(sbfem_skeleton);
                if(!neighbour) DebugStop();
                if(!neighbour.Element()->Reference()) DebugStop();
            } else {
                DebugStop();
            }
        } else if (neighgel->Dimension() == 1 && neighgel->MaterialId() == sbfem_highperm) {
            // create a zero dimensional element
            if(!neighgel->IsLinearMapping()) {
                int64_t index = neighgel->Index();
                ChangeToLinearOned(neighgel);
                neighgel = gmesh->Element(index);
            }
            connected_elements.insert(neighgel);
            TPZGeoElSide gelside(neighgel,0);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while(neighbour != gelside) {
                TPZGeoEl *neighgel = neighbour.Element();
                if(neighgel->MaterialId() == sbfem_skeleton && neighgel->Dimension() == 0) {
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
            if(neighbour == gelside) DebugStop();
            if(neighbour.Element()->Dimension() != 0) DebugStop();
            if(!neighbour.Element()->Reference()) DebugStop();
        }
    }
    TPZSBFemElementGroup *elgr = new TPZSBFemElementGroup(*cmesh);
    sbfem_groupHdiv = elgr;
    TPZStack<TPZSBFemVolume *> sbvols;
    for(auto it : connected_elements) {
        TPZGeoEl *gel = it;
        TPZSBFemVolume *sbvol = new TPZSBFemVolume(*cmesh,gel);
        sbvols.Push(sbvol);
        if(gel->Dimension() == 2) {
            TPZGeoElSide gelside(gel,4);
            TPZGeoElSide neigh = gelside.HasNeighbour(sbfem_skeleton);
            TPZGeoEl *skel = neigh.Element();
            if(skel->MaterialId() != sbfem_skeleton) DebugStop();
            TPZCompEl *cskel = skel->Reference();
            if(!cskel) DebugStop();
            TPZMultiphysicsElement *mcel = dynamic_cast<TPZMultiphysicsElement *>(cskel);
            if(!mcel) DebugStop();
            TPZCompEl *cskelat = mcel->ElementVec()[0].Element();
            TPZInterpolatedElement *icskel = dynamic_cast<TPZInterpolatedElement *>(cskelat);
            if(!icskel) DebugStop();
            int sideorder = icskel->Connect(2).Order();
            if(sideorder != SBFemOrder) DebugStop();
            sbvol->SetSkeleton(cskel->Index());
        } else if(gel->Dimension() == 1) {
            TPZGeoElSide gelside(gel,0);
            TPZGeoElSide neigh = gelside.Neighbour();
            while(neigh != gelside) {
                TPZGeoEl *neighgel = neigh.Element();
                if(neighgel->MaterialId() == sbfem_skeleton && neighgel->Dimension() == 0) {
                    break;
                }
                neigh = neigh.Neighbour();
            }
            if(neigh.Element()->MaterialId() != sbfem_skeleton) DebugStop();
                        // create an H1 element
            TPZCompEl *cskel = neigh.Element()->Reference();
            if(!cskel) DebugStop();
            sbvol->SetSkeleton(cskel->Index());
        }
        elgr->AddElement(sbvol);
    }
    for(auto sbvol : sbvols) {
        sbvol->SetElementGroupIndex(elgr->Index());
    }
    TPZElementMatrixT<STATE> ek,ef;
    elgr->CalcStiff(ek,ef);
    if(0) {
        ek.fMat.Print("ek = ",std::cout,EMathematicaInput);
    }
    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();

}

/// Add the SBFemVolume elements to the computational mesh
void AddSBFemVolumeElements()
{
    if(meshstyle == ESBFem) {

        TPZCompMesh *cmesh;
        if(sbfem_groupH1) cmesh = sbfem_groupH1->Mesh();
        TPZMultiphysicsCompMesh *cmesh_m;
        if(sbfem_groupHdiv) cmesh_m = dynamic_cast<TPZMultiphysicsCompMesh *>(sbfem_groupHdiv->Mesh());
        if(sbfem_groupH1 && !cmesh) DebugStop();
        if(sbfem_groupHdiv && !cmesh_m) DebugStop();
        if(sbfem_groupH1) {
            TPZStack<TPZCompEl *> compels;
            sbfem_groupH1->GetCompElList(compels);
            for(auto compel : compels) {
                int64_t index = compel->Index();
                cmesh->ElementVec()[index] = compel;
            } 
        }
        if(sbfem_groupHdiv) {
            TPZStack<TPZCompEl *> compels;
            sbfem_groupHdiv->GetCompElList(compels);
            for(auto compel : compels) {
                int64_t index = compel->Index();
                cmesh_m->ElementVec()[index] = compel;
            } 
        }
    }
}

/// Hide the SBFemVolume elements from the computational mesh
void HideSBFemVolumeElements() {
    if(meshstyle == ESBFem) {

        TPZCompMesh *cmesh;
        if(sbfem_groupH1) cmesh = sbfem_groupH1->Mesh();
        TPZMultiphysicsCompMesh *cmesh_m;
        if(sbfem_groupHdiv) cmesh_m = dynamic_cast<TPZMultiphysicsCompMesh *>(sbfem_groupHdiv->Mesh());
        if(sbfem_groupH1 && !cmesh) DebugStop();
        if(sbfem_groupHdiv && !cmesh_m) DebugStop();
        if(sbfem_groupH1) {
            TPZStack<TPZCompEl *> compels;
            sbfem_groupH1->GetCompElList(compels);
            for(auto compel : compels) {
                int64_t index = compel->Index();
                cmesh->ElementVec()[index] = 0;
            } 
        }
        if(sbfem_groupHdiv)
        {
            TPZStack<TPZCompEl *> compels;
            sbfem_groupHdiv->GetCompElList(compels);
            for(auto compel : compels) {
                int64_t index = compel->Index();
                cmesh_m->ElementVec()[index] = 0;
            } 
        }
    }
}

/// @brief adjust the elements neighbouring the trailing edge to accomodate SBFem simulation
void AdjustToSBFemGeometry(TPZGeoMesh *gmesh) {
    TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
    int64_t singular = trailingedge_element->NodeIndex(0);
    TPZGeoElSide trailingedge(trailingedge_element);
    std::set<TPZGeoEl *> connected_elements;
    for(TPZGeoElSide neigh = trailingedge.Neighbour(); neigh != trailingedge; neigh++) {
        TPZGeoEl *neighgel = neigh.Element();
        if(neighgel->HasSubElement()) continue;
        // the element with matid < 0 is the element subsituted by two collapsed elements
        if(neighgel->MaterialId() < 0) continue;
        if(connected_elements.find(neighgel) != connected_elements.end()) continue;
        // only insert collapsed elements
        if(neighgel->Dimension() == 2) {
            if(neighgel->NodeIndex(2) == neighgel->NodeIndex(3)) {
                connected_elements.insert(neighgel);
            }
        } else if(neighgel->Dimension() == 1) {
            connected_elements.insert(neighgel);
        }
    }
    // change the material id of the volumetric elements to sbfem_domain
    // create the skeleton elements
    for(auto it : connected_elements) {
        TPZGeoEl *gel = it;
        if(gel->Dimension() == 2) {
            gel->SetMaterialId(sbfem_domain);
            TPZGeoElBC gbc(gel,4,sbfem_skeleton);
            TPZGeoEl *gbcgel = gbc.CreatedElement();
            if(gbcgel->IsLinearMapping() == false) {
                gel->Print(std::cout);
                gbcgel->Print(std::cout);
                DebugStop();
            }
        }
    }
    // create sbfem_highperm elements
    for(auto it : connected_elements) {
        TPZGeoEl *gel = it;
        if(gel->MaterialId() == profilemat) {
            TPZGeoElSide gelside(gel);
            TPZManVector<int64_t,4> nodeindices(2);
            if(gelside.SideNodeIndex(0) == singular) {
                nodeindices[0] = gelside.SideNodeIndex(1);
                nodeindices[1] = gelside.SideNodeIndex(0);
            } else {
                nodeindices[0] = gelside.SideNodeIndex(0);
                nodeindices[1] = gelside.SideNodeIndex(1);
            }
            // creating a linear element will lead to an inconsistent geometry
            int64_t index;
            gmesh->CreateGeoElement(EOned, nodeindices, sbfem_highperm, index);
            // create the skeleton element
            gmesh->CreateGeoElement(EPoint, nodeindices, sbfem_skeleton, index);
        }
    }
    // insert the connectivity of the sbfem_highperm elements
    gmesh->BuildConnectivity();
    // make the elements linear
    std::set<TPZGeoEl *> removed;
    for(auto it : connected_elements) {
        TPZGeoEl *gel = it;
        if(!gel->IsLinearMapping()) {
            if(gel->Dimension() == 2) {
                removed.insert(gel);
                ChangeToLinearQuad(gel);
            }
        } else if(gel->Dimension() == 2) {
            removed.insert(gel);
        }
    }
    for(auto it : removed) {
        connected_elements.erase(it);
    }
    // Print the 1d elements that are not linear
    if(0) {
        for(auto it : connected_elements) {
            TPZGeoEl *gel = it;
            if(gel->Dimension() != 1) {
                DebugStop();
            }
            TPZGeoElSide gelside(gel);
            std::cout << "Nonlinear neighbouring element " << gel->Index() << " matid " << gel->MaterialId() << std::endl;
            // gelside.RemoveConnectivity();
            // gelside.SetConnectivity(gelside);
        }
        PrintTrailingEdgeElements(gmesh);
    }
#ifdef PZDEBUG
    {
        const int64_t nel = gmesh->NElements();
        for(int64_t el = 0; el < nel; el++)
        {
            TPZGeoEl * gel = gmesh->ElementVec()[el];
            if(!gel) continue;
            gel->VerifyNodeCoordinates();
        }//for el
    }
#endif

}

