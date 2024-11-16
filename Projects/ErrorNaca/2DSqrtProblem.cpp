/**
 * @file This file implements an error estimator in space H1.
 */

#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "Projection/TPZL2ProjectionCS.h"
#include "TPZAnalyticSolution.h"
#include "TPZBndCondT.h"
#include "TPZFrontSym.h"
#include "TPZGmshReader.h"
#include "TPZLinearAnalysis.h"
#include "TPZNullMaterialCS.h"
#include "TPZRefPatternDataBase.h"
#include "TPZSBFemElementGroup.h"
#include "TPZSBFemVolume.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSYSMPMatrix.h"
#include "TPZTensor.h"
#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"
#include "pzcheckgeom.h"
#include "pzgeoelbc.h"
#include "pzinterpolationspace.h"
#include "pzlog.h"
#include "pzskylmat.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include <TPZMultiphysicsCompMesh.h>
#include <TPZNullMaterial.h>
#include <TPZSimpleTimer.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzmultiphysicscompel.h>

#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>

#include "TPZGenGrid2D.h"

#include "TPZGeoLinear.h"
#include "tpzarc3d.h"
#include "tpzgeoblend.h"

#include "tpzblendnaca.h"

#include "tpzchangeel.h"

#include <iostream>

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.material.darcy");
#endif



//#define USING_MKL

using namespace pzgeom;
/// @brief verify is the derivative of the NACA coordinate is correct
/// @param naca profile object
/// @param point parametric coordinate around which the derivative will be verified
void VerifyDerivative(TPZBlendNACA &naca, REAL point);

/// @brief Read a gmsh file and return a geometric mesh
/// @param filename
/// @return a geometric mesh
TPZGeoMesh *ReadGmsh(const std::string &filename);

/// @brief Create a computational mesh with H1 elements
TPZCompMesh *CreateH1CompMesh(TPZGeoMesh *gmesh, TPZVec<int> &porders, TPZAnalyticSolution *analyticSol);

/// @brief Adjust the integration rule for Equarterpoints H1 elements
void AdjustH1Equarterpointintrule(TPZGeoMesh *gmesh);
// void AdjustH1Equarterpointintrule(TPZGeoMesh *gmesh);

/// @brief Create a computational mesh with L2 elements
TPZCompMesh *CreateL2CompMesh(TPZGeoMesh *gmesh);

/// @brief Create the computational mesh with HDiv elements
TPZCompMesh *CreateHDivCompMesh(TPZGeoMesh *gmesh, TPZVec<int> &porders, TPZAnalyticSolution *analyticSol);

/// @brief Adjust the integration rule for Equarterpoints HDiv elements
void AdjustHDivEquarterpointintrule(TPZGeoMesh *gmesh);
// void AdjustHDivEquarterpointintrule(TPZGeoMesh *gmesh);

/// @brief Create the computational "multiphysics" mesh with only HDiv elements
TPZMultiphysicsCompMesh *CreateMultiphysicsMesh(TPZCompMesh *cmeshHDiv, TPZCompMesh *cmeshL2, TPZGeoMesh *gmesh,
                                                TPZAnalyticSolution *analyticSol);

/// @brief Simulate the NACA profile using H1 approximation for Beta = 0
TPZCompMesh *SimulateH1(TPZGeoMesh *gmesh, TPZVec<int> &porders, TPZAnalyticSolution *analyticSol);

/// @brief Simulate the NACA profile using H(div) approximation for Beta = 0
TPZMultiphysicsCompMesh *SimulateHDiv(TPZGeoMesh *gmesh, TPZVec<int> &porders, TPZAnalyticSolution *analyticSol);

/// @brief Compute the error as an average error per element and save as an elemental solution
void ComputeErrorEstimator(TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator, REAL &GlobalError, REAL &ErrorH1, REAL &ErrorHDiv);

/// @brief Compute the global error the norm of the vector ErrorEstimator
void ComputeGlobalError(TPZVec<REAL> &ErrorEstimator, REAL &GlobalError);

/// @brief Compute the error as an average error per element and save as an elemental solution
void Hrefinement(TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator, TPZVec<REAL> &RefinementIndicator,
                 TPZVec<int> &porders);

/// @brief Performe h-refinements at the computational meshes "cmesh" and "cmesh_m" based on the ErrorEstimator
// @param cmesh_m multiphysics computational mesh
// @param ErrorEstimator vector with the error estimator for each computational element
// @param minh minimum level of refinement
// @param RefinementIndicator vector with the refinement indicator for each computational element
// @param porders vector with the polynomial order of each geometric element
void HPrefinement(TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator, int minh,
                  TPZVec<REAL> &RefinementIndicator, TPZVec<int> &porders);

/// @brief Smoothen the volumetric elements around the trailingedge by adopting the same level of refinement
void Smoothentrailingedgeelements(TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &RefinementIndicator);

/// @brief Substitute the trailing edge quadrilateral elements with colapsed quadrilateral elements with or without
/// quarterpoint elements
void Changetrailingedgeelements(TPZGeoMesh *gmesh);


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
void CreateHdivSBFEMelements(TPZMultiphysicsCompMesh *cmesh, TPZAnalyticSolution *analyticSol);

/// @brief create a refinement patter that cuts the element in horizontal direction
/// @return refinement pattern created
TPZAutoPointer<TPZRefPattern> CreateRefPattern_collapsed();

/// @brief create a refinement patter that cuts the element in vertical direction
/// @return refinement pattern created
TPZAutoPointer<TPZRefPattern> CreateRefPattern_sbfem();

/// @brief divide a geometric element. If it is a trailing edge element, set the refinement pattern first
void DivideGeoEl(TPZGeoEl *gel, TPZVec<TPZGeoEl *> &subels);

/// @brief returns true if the element is counter clockwise
bool IsCounterClockwise(TPZGeoEl *gel);

/// @brief make the elements counter clockwise
void MakeCounterClockwise(TPZGeoEl *gel);

/// @brief Divide Trailing Edge Neighbours
void DivideTrailingEdgeNeighbours(TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator, REAL tol,
                                  TPZVec<REAL> &RefinementIndicator, TPZVec<int> &porders);

/// @brief SmoothenGeometry
void SmoothenGeometry(TPZGeoMesh *gmesh);

/// @brief print the results of the analysis
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh);

/// Add the SBFemVolume elements to the computational mesh
void AddSBFemVolumeElements();

/// Hide the SBFemVolume elements from the computational mesh
void HideSBFemVolumeElements();

int volmat = 1;
int dirichletmat = 2;
int neumannmat = 3;
int trailingedgemat = 4;

int sbfem_skeleton = 8;
int sbfem_domain = 9;

int sbfem_highperm_hdiv = 10;
int sbfem_highperm_h1 = 11;
REAL shift_distance = 1.e-2;
TPZSBFemElementGroup *sbfem_groupH1 = 0;
TPZSBFemElementGroup *sbfem_groupHdiv = 0;
enum MMeshStyle { ETraditional, ECollapsed, EQuarterPoint, ESBFem };
MMeshStyle meshstyle = ESBFem;
int defaultporder = 1;
int SBFemOrder = 2;
int nuniform = 1;
int nrefinements = 13;
// set of geometric element indices that are SBFem elements
std::set<int64_t> sbfem_elements;

enum RRefinementStyle { h, hp , huniform, puniform};
RRefinementStyle refinementstyle = h;

TPZAutoPointer<TPZRefPattern> refpattern_collapsed;
TPZAutoPointer<TPZRefPattern> refpattern_sbfem;
int64_t trailingedge_element_index = -1;

int main() {

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    if (meshstyle != ETraditional) {
        gRefDBase.InitializeRefPatterns(2);
    }

    TPZGeoMesh *gmesh = ReadGmsh("2DSqrtProblem.msh");
    {
        std::ofstream out("gmesh.txt");
        gmesh->Print(out);

        std::ofstream out2("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
    }
    {
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (gel && gel->MaterialId() == trailingedgemat) {
                trailingedge_element_index = el;
                break;
            }
        }
        if (trailingedge_element_index == -1) DebugStop();
    }
    // look for a quadrilateral element and keep the result in a global variable
    refpattern_collapsed = CreateRefPattern_collapsed();
    refpattern_sbfem = CreateRefPattern_sbfem();
    {
        TPZCheckGeom check(gmesh);
        if (nuniform) {
            check.UniformRefine(nuniform);
            if (0) {
                std::ofstream out4("gmeshfine.txt");
                gmesh->Print(out4);

                std::ofstream out5("gmeshfine.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out5);
            }
        }
    }
    
    // change the elements touching the trailing edge to collapsed elements
    if (meshstyle != ETraditional) {
        // false : create the midside node in the middle of the trailing edge
        Changetrailingedgeelements(gmesh);
        if (0) {
            std::ofstream out("gmeshchanged.txt");
            gmesh->Print(out);

            std::ofstream out2("gmeshchanged.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);

            PrintTrailingEdgeElements(gmesh);
        }
    }
    // this is where the SBFem elements are created
    // this is where skeleton and high permeability elements are created
    // this is applied to the temporary copy of the mesh
    // these elements will be deleted in the method RemoveSBFemElements
    if (meshstyle == ESBFem) {
        // Adjust the geometry to accomodate the SBFem simulation
        AdjustToSBFemGeometry(gmesh);
        if (0) {
            std::ofstream out("gmeshsbfem.txt");
            gmesh->Print(out);

            std::ofstream out2("gmeshsbfem.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
        }
    } else if (meshstyle == EQuarterPoint) {
        // change the elements of gmeshcopy to quadratic elements
        // observe that is applied to the copy. It does not affect the original mesh
        CreateQuarterPointElements(gmesh);
        PrintTrailingEdgeElements(gmesh);
    }

    // Creating the analytical solution
    TLaplaceExample1 *analytic = new TLaplaceExample1();
    analytic->fExact = TLaplaceExample1::ESquareRoot;

    int minh = 1;
    // indicating the flux order
    int64_t nel = gmesh->NElements();
    TPZVec<int> porders(nel, defaultporder);
    REAL GlobalError;
    REAL ErrorH1;
    REAL ErrorHDiv;
    REAL circulation_H1;
    REAL circulation_HDiv;
    std::ofstream outGE0("GlobalError.txt");
    std::ofstream outGE1("ErrorEstimatorH1HDiv.txt");
    std::ofstream outGE2("ErrorH1.txt");
    std::ofstream outGE3("ErrorHDiv.txt");
    for (int64_t i = 0; i < nrefinements; i++) {
        TPZCompMesh *cmesh = 0;
        TPZMultiphysicsCompMesh *cmesh_m = 0;
        if (1) {
            std::ofstream out("gmesh.txt");
            gmesh->Print(out);
            std::ofstream out2("gmesh.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
        }

        cmesh_m = SimulateHDiv(gmesh, porders, analytic);
        {
            std::ofstream out("cmeshHdiv.txt");
            cmesh_m->Print(out);
        }

        cmesh = SimulateH1(gmesh, porders, analytic);
        {
            std::ofstream out("cmeshH1.txt");
            cmesh->Print(out);
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
        ComputeErrorEstimator(cmesh, cmesh_m, Error, GlobalError, ErrorH1, ErrorHDiv);
        // ComputeGlobalError(Error, GlobalError);
        std::cout << "GlobalError: " << GlobalError << std::endl;

#ifdef PZDEBUG
        {
            gmesh->ResetReference();
            cmesh_m->LoadReferences();
            std::cout << "Number of SBFem elements " << sbfem_elements.size() << std::endl;
            for (auto el : sbfem_elements) {
                TPZGeoEl *gel = gmesh->Element(el);
                if (gel->Reference() == 0) {
                    std::cout << "Element " << el << " has no reference\n";
                }
            }
        }
#endif
        int64_t nDOF = cmesh->NEquations();
        // int64_t nDOF_m = cmesh_m->NEquations();
        outGE0 << i << "  " << nDOF << "  " << GlobalError << "  " << circulation_H1 << "  " << circulation_HDiv
              << std::endl;
        outGE1 << nDOF << "  " << GlobalError << std::endl;
        outGE2 << nDOF << "  " << ErrorH1 << std::endl;
        outGE3 << nDOF << "  " << ErrorHDiv << std::endl;

        TPZVec<REAL> RefinementIndicator;
        if (refinementstyle == hp) {
            HPrefinement(cmesh_m, Error, minh, RefinementIndicator, porders);
        } else if (refinementstyle == h) {
            Hrefinement(cmesh_m, Error, RefinementIndicator, porders);
        } else if (refinementstyle == huniform) {
            int64_t nel = gmesh->NElements();
            for (int64_t i = 0; i < nel; i++) {
                TPZGeoEl *gel = gmesh->Element(i);
                if(!gel) continue;
                int matid = gel->MaterialId();
                if(gel->Dimension() == 0) continue;
                if(gel->HasSubElement()) continue;
                if(matid < 0) continue;
                if(matid == sbfem_highperm_h1) continue;
                if(matid == sbfem_highperm_hdiv) continue;
                TPZGeoElSide gelside(gel);
                if(gelside.HasNeighbour(sbfem_highperm_h1)) continue;
                if(gelside.HasNeighbour(sbfem_highperm_hdiv)) continue;
                TPZManVector<TPZGeoEl *, 4> subels;
                DivideGeoEl(gel, subels);
            }
            porders.Resize(gmesh->NElements(), defaultporder);
        } else if (refinementstyle == puniform) {
            porders.Fill(++defaultporder);
            SBFemOrder++;
        }
        if(meshstyle == ESBFem)
        {
            sbfem_elements.clear();
            TPZGeoEl *gel = gmesh->Element(trailingedge_element_index);
            TPZGeoElSide gelside(gel);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                TPZGeoEl *neighgel = neighbour.Element();
                int dim = neighgel->Dimension();
                int matid = neighgel->MaterialId();

                if(dim == 2 && !neighgel->HasSubElement() && matid >= 0) {
                    sbfem_elements.insert(neighgel->Index());
                }
                neighbour = neighbour.Neighbour();
            }
        }
        {
            std::ofstream out2("gmeshrefined.txt");
            gmesh->Print(out2);
            std::ofstream out3("gmeshrefined.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out3);
            std::ofstream out6("ErrorEstimator.vtk");
            AddSBFemVolumeElements();
            TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m, out6, Error, "Error");
            HideSBFemVolumeElements();
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
        sbfem_groupH1 = 0;
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
        sbfem_groupHdiv = 0;
    }
    delete gmesh;
    return 0;
}

/// @brief verify is the derivative of the NACA coordinate is correct
/// @param naca profile object
/// @param point parametric coordinate around which the derivative will be verified
void VerifyDerivative(TPZBlendNACA &naca, REAL point) {
    TPZFNMatrix<20, REAL> pointvalues(10, 2, 0.), expected(10, 2, 0.);
    TPZManVector<REAL, 10> error(10, 0.), distance(10, 0.);
    TPZManVector<REAL, 2> derivative = {naca.dxla(point), naca.dyla(point)};
    std::cout << "derivative " << derivative << std::endl;
    REAL delta = 0.001;
    for (int i = 0; i < 10; i++) {
        REAL pos = point + delta * i;
        pointvalues(i, 0) = naca.xla(pos);
        pointvalues(i, 1) = naca.yla(pos);
        distance[i] = delta * i;
        expected(i, 0) = pointvalues(0, 0) + derivative[0] * distance[i];
        expected(i, 1) = pointvalues(0, 1) + derivative[1] * distance[i];
        error[i] = sqrt((expected(i, 0) - pointvalues(i, 0)) * (expected(i, 0) - pointvalues(i, 0)) +
                        (expected(i, 1) - pointvalues(i, 1)) * (expected(i, 1) - pointvalues(i, 1)));
    }
    std::cout << "error " << error << std::endl;
    TPZManVector<REAL, 9> rate(9, 0.);
    for (int i = 1; i < 9; i++) {
        rate[i] = log(error[i + 1] / error[i]) / log(distance[i + 1] / distance[i]);
    }
    std::cout << "rate " << rate << std::endl;
}

/// @brief Read a gmsh file and return a geometric mesh
/// @param filename
/// @return a geometric mesh
TPZGeoMesh *ReadGmsh(const std::string &meshfilename) {
    TPZGmshReader gmsh;
    gmsh.SetVerbose(1);
    gmsh.GetDimNamePhysical()[0]["Trailingedge"] = trailingedgemat;
    gmsh.GetDimNamePhysical()[1]["Dirichlet"] = dirichletmat;
    gmsh.GetDimNamePhysical()[1]["Neumman"] = neumannmat;
    gmsh.GetDimNamePhysical()[2]["Domain"] = volmat;
    auto gmesh = gmsh.GeometricGmshMesh(meshfilename);

    return gmesh;
}

#include "pzintel.h"
/// @brief Create a computational mesh with H1 elements
TPZCompMesh *CreateH1CompMesh(TPZGeoMesh *gmesh, TPZVec<int> &porders, TPZAnalyticSolution *analyticSol) {
    if (porders.size() != gmesh->NElements()) DebugStop();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    int dim = 2;
    TPZDarcyFlow *material = new TPZDarcyFlow(volmat, dim);
    material->SetExactSol(analyticSol->ExactSolution(), 4);
    cmesh->InsertMaterialObject(material);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(defaultporder);

    TPZFMatrix<STATE> val1(1, 1, 0.);
    TPZManVector<STATE> val2(1, 0.);

    auto bnd1 = material->CreateBC(material, dirichletmat, 0, val1, val2);
    bnd1->SetForcingFunctionBC(analyticSol->ExactSolution(), 4);
    cmesh->InsertMaterialObject(bnd1);
    auto bnd2 = material->CreateBC(material, neumannmat, 1, val1, val2);
    // bnd2->SetForcingFunctionBC(analyticSol->ExactSolution(),4);
    cmesh->InsertMaterialObject(bnd2);
    if (meshstyle == ESBFem) {
        auto bnd4 = material->CreateBC(material, sbfem_skeleton, 1, val1, val2);
        cmesh->InsertMaterialObject(bnd4);
    }
    cmesh->SetAllCreateFunctionsContinuous();
    std::set<int> matidsh1 = {volmat, dirichletmat, neumannmat};
    if (meshstyle == ESBFem) {
        matidsh1.insert(sbfem_skeleton);
    }
    {
        int64_t nel = gmesh->NElements();
        TPZStack<int64_t> gelstack;
        TPZStack<int> elp;
        for (int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (gel->HasSubElement()) continue;
            int matid = gel->MaterialId();
            if (matidsh1.find(matid) != matidsh1.end()) {
                gelstack.Push(gel->Index());
                if (porders[el] < 0) DebugStop();
                elp.Push(porders[el] + 1);
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

    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild(matidsh1);
    cmesh->ExpandSolution();
    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();

    if (meshstyle == ESBFem) {
        CreateH1SBFEMelements(cmesh);
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
    for (; Neighbor != TrailingSide; Neighbor++) {
        if (!Neighbor.Element()) DebugStop();
        if (Neighbor.Element()->HasSubElement()) continue;
        if (Neighbor.Element()->MaterialId() == volmat) {
            TrailingedgeH1Element = Neighbor.Element()->Reference();
            TPZGeoElSide NeighborEl(Neighbor.Element());
            TPZIntPoints *intrule = NeighborEl.CreateIntegrationRule(10);
            TrailingedgeH1Element->SetIntegrationRule(intrule);
        }
    }
}

/// @brief Create a computational mesh with L2 elements
TPZCompMesh *CreateL2CompMesh(TPZGeoMesh *gmesh) {
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    int dim = 2;

    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(1);

    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    // o método ApproxSpace().CreateDisconnectedElements(true) é chamado para transformar os elementos em elementos
    // discontinuos, criando assim um espaço L2 de aproximação.

    cmesh->AutoBuild();

    return cmesh;
}

/// @brief Create a computational mesh with HDiv elements
TPZCompMesh *CreateHDivCompMesh(TPZGeoMesh *gmesh, TPZVec<int> &porders, TPZAnalyticSolution *analyticSol) {
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    int dim = 2;
    TPZDarcyFlow *material = new TPZDarcyFlow(volmat, dim);
    material->SetExactSol(analyticSol->ExactSolution(), 4);
    cmesh->InsertMaterialObject(material);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(defaultporder);

    HDivFamily hdiv = HDivFamily::EHDivKernel;
    cmesh->ApproxSpace().SetHDivFamily(hdiv);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);

    TPZFMatrix<STATE> val1(2, 2, 0.);
    TPZManVector<STATE> val2(2, 0.);

    // Dirichlet condition
    auto bnd1 = material->CreateBC(material, dirichletmat, 0, val1, val2);
    bnd1->SetForcingFunctionBC(analyticSol->ExactSolution(), 4);
    cmesh->InsertMaterialObject(bnd1);

    // Neumann condition
    auto bnd2 = material->CreateBC(material, neumannmat, 5, val1, val2);
    // bnd1->SetForcingFunctionBC(analyticSol->ExactSolution(),4);
    cmesh->InsertMaterialObject(bnd2);

    if (meshstyle == ESBFem) {
        auto bnd3 = material->CreateBC(material, sbfem_skeleton, 0, val1, val2);
        cmesh->InsertMaterialObject(bnd3);
    }

    std::set<int> matidshdiv = {volmat, dirichletmat, neumannmat};
    // build the computational mesh. HDivKernel elements for matidshdiv
    // H1 elements for the skeleton
    {
        int64_t nel = gmesh->NElements();
        TPZStack<int64_t> gelstack, gelstack_skel;
        TPZStack<int> elp, elp_skel;
        for (int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel) continue;
            if (gel->HasSubElement()) continue;
            int matid = gel->MaterialId();
            if (matidshdiv.find(matid) != matidshdiv.end()) {
                gelstack.Push(gel->Index());
                if (porders[el] < 0) DebugStop();
                elp.Push(porders[el] + 1);
            } else if (matid == sbfem_skeleton && meshstyle == ESBFem) {
                gelstack_skel.Push(gel->Index());
                if (porders[el] < 0) DebugStop();
                elp_skel.Push(SBFemOrder);
            }
        }
        cmesh->ApproxSpace().BuildMesh(*cmesh, gelstack, elp);
        if (meshstyle == ESBFem) {
            // insert the skeleton elements as H1 elements
            cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
            cmesh->SetDefaultOrder(SBFemOrder);
            cmesh->ApproxSpace().BuildMesh(*cmesh, gelstack_skel, elp_skel);
        }
        // adjust the connect order of the skeleton elements
        for (int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (gel->Dimension() != 1) continue;
            if (gel->HasSubElement()) continue;
            int matid = gel->MaterialId();
            if (matid == sbfem_skeleton && meshstyle == ESBFem) {
                TPZCompEl *cel = gel->Reference();
                if (!cel) DebugStop();
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                if (!intel) DebugStop();
                intel->SetSideOrder(2, SBFemOrder);
            }
        }
        cmesh->ExpandSolution();
    }

    // cmesh->AutoBuild(matidshdiv);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();

    // force the side corresponding to a collapsed connect to have order 1
    int64_t nelem = cmesh->NElements();
    for (int64_t el = 0; el < nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (!gel) continue;
        if (gel->Type() != EQuadrilateral) continue;
        for (int ic = 0; ic < 4; ic++) {
            int ic1 = (ic + 1) % 4;
            if (cel->ConnectIndex(ic) == cel->ConnectIndex(ic1)) {
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                if (!intel) DebugStop();
                intel->SetSideOrder(ic + 4, 1);
            }
        }
    }
    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();

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
    for (; Neighbor != TrailingSide; Neighbor++) {
        if (!Neighbor.Element()) DebugStop();
        if (Neighbor.Element()->HasSubElement()) continue;
        if (Neighbor.Element()->MaterialId() == volmat) {
            auto NeighborCompel = Neighbor.Element()->Reference();
            TrailingedgeHDivElement = dynamic_cast<TPZMultiphysicsElement *>(NeighborCompel);
            if (!TrailingedgeHDivElement) DebugStop();
            TPZGeoElSide NeighborEl(Neighbor.Element());
            TPZIntPoints *intrule = NeighborEl.CreateIntegrationRule(10);
            TrailingedgeHDivElement->SetIntegrationRule(intrule);
        }
    }
}

/// @brief Create a computational "multiphysics" mesh with only HDiv elements
TPZMultiphysicsCompMesh *CreateMultiphysicsMesh(TPZCompMesh *cmeshHDiv, TPZCompMesh *cmeshL2, TPZGeoMesh *gmesh,
                                                TPZAnalyticSolution *analyticSol) {
    //
    TPZMultiphysicsCompMesh *cmesh_m = new TPZMultiphysicsCompMesh(gmesh);

    const int dim = cmeshHDiv->Dimension();
    const int pord = cmeshHDiv->GetDefaultOrder();
    cmesh_m->SetDimModel(dim);
    cmesh_m->SetDefaultOrder(pord);

    cmesh_m->SetAllCreateFunctionsMultiphysicElem();

    // 1. Materials
    std::set<int> materialIDs;
    TPZMixedDarcyFlow *material = new TPZMixedDarcyFlow(volmat, dim);
    material->SetExactSol(analyticSol->ExactSolution(), 4);
    cmesh_m->InsertMaterialObject(material);
    materialIDs.insert(volmat);

    // 2. Boundary Conditions
    TPZFMatrix<STATE> val1(2, 2, 0.);
    TPZManVector<STATE> val2(2, 0.);

    // 2.1 Dirichlet condition
    auto bnd1 = material->CreateBC(material, dirichletmat, 0, val1, val2);
    bnd1->SetForcingFunctionBC(analyticSol->ExactSolution(), 4);
    cmesh_m->InsertMaterialObject(bnd1);

    // // 2.2 Neumann condition
    auto bnd2 = material->CreateBC(material, neumannmat, 5, val1, val2);
    // bnd2->SetForcingFunctionBC(analyticSol->ExactSolution(),4);
    cmesh_m->InsertMaterialObject(bnd2);

    TPZNullMaterialCS<STATE> *bnd3 = new TPZNullMaterialCS<STATE>(sbfem_skeleton);
    cmesh_m->InsertMaterialObject(bnd3);

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

/// @brief Simulate the NACA profile using H1 approximation for Beta = 0
TPZCompMesh *SimulateH1(TPZGeoMesh *gmesh, TPZVec<int> &porders, TPZAnalyticSolution *analyticSol) {
    auto cmeshH1 = CreateH1CompMesh(gmesh, porders, analyticSol);
    if (meshstyle == EQuarterPoint) {
        AdjustH1Equarterpointintrule(gmesh);
    }

    TPZLinearAnalysis an(cmeshH1, RenumType::EMetis);
    if (1) {
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

    an.Solve();
    TPZFMatrix<STATE> &fsol = an.Solution();
    cmeshH1->LoadSolution(fsol);

    {
        int64_t nel = cmeshH1->NElements();
        TPZFMatrix<STATE> &elsol = cmeshH1->ElementSolution();
        elsol.Redim(nel, 6);
    }
    TPZVec<REAL> Errors;
    an.PostProcessError(Errors, true);

    std::cout << "The H1 error is " << Errors[2] << std::endl;

    AddSBFemVolumeElements();
    {
        REAL trad_error = 0.;
        REAL sbfem_error = 0.;
        int64_t nel = cmeshH1->NElements();
        TPZFMatrix<STATE> &elsol = cmeshH1->ElementSolution();
       for(int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = cmeshH1->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if(!gel) continue;
            if(gel->Dimension() != 2) continue;
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
            if(intel) {
                trad_error += elsol(el,2)*elsol(el,2);
            } else if(sbfem) {
                sbfem_error += elsol(el,2)*elsol(el,2);
            }
        }
        std::cout << "The interior H1 error is " << sqrt(trad_error) << std::endl;
        std::cout << "The SBFem H1 error is " << sqrt(sbfem_error) << std::endl;
    }
    HideSBFemVolumeElements();

    return cmeshH1;
}

/// @brief Simulate the NACA profile using H(div) approximation for Beta = 0
TPZMultiphysicsCompMesh *SimulateHDiv(TPZGeoMesh *gmesh, TPZVec<int> &porders, TPZAnalyticSolution *analyticSol) {
    gmesh->ResetReference();
    auto cmeshHDiv = CreateHDivCompMesh(gmesh, porders, analyticSol);
    auto cmeshL2 = CreateL2CompMesh(gmesh);

    TPZMultiphysicsCompMesh *cmesh_m = CreateMultiphysicsMesh(cmeshHDiv, cmeshL2, gmesh, analyticSol);
    if (meshstyle == EQuarterPoint) {
        AdjustHDivEquarterpointintrule(gmesh);
    }
    if (meshstyle == ESBFem) {
        CreateHdivSBFEMelements(cmesh_m, analyticSol);
    }
    cmesh_m->CleanUpUnconnectedNodes();
    // Define o pointer chamado cmesh_m relacionado à classe TPZMultiphysicsCompMesh, associando-o a função
    // CreateMultiphysicsMesh, cujos parâmetros (já antes declarados) são: cmeshHdiv, gmesh.

    TPZLinearAnalysis an(cmesh_m, RenumType::EMetis);
    {
        std::ofstream out("cmesh_m.txt");
        cmesh_m->Print(out);
    }

#ifdef USING_MKL
    TPZSSpStructMatrix<STATE> strmat(cmesh_m);
    strmat.SetNumThreads(0);xx
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
    TPZFMatrix<STATE> &fsol = an.Solution();
    cmesh_m->LoadSolution(fsol);
    cmesh_m->TransferMultiphysicsSolution();

    {
        int64_t nel = cmesh_m->NElements();
        TPZFMatrix<STATE> &elsol = cmesh_m->ElementSolution();
        elsol.Redim(nel, 5);
    }
    TPZVec<REAL> Errors;
    an.PostProcessError(Errors, true);
    std::cout << "The HDiv error is " << Errors[1] << std::endl;


    AddSBFemVolumeElements();
    {
        REAL trad_error = 0.;
        REAL sbfem_error = 0.;
        int64_t nel = cmesh_m->NElements();
        TPZFMatrix<STATE> &elsol = cmesh_m->ElementSolution();
       for(int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = cmesh_m->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if(!gel) continue;
            if(gel->Dimension() != 2) continue;
            TPZMultiphysicsElement *intel = dynamic_cast<TPZMultiphysicsElement *>(cel);
            TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
            if(intel) {
                trad_error += elsol(el,1)*elsol(el,1);
            } else if(sbfem) {
                sbfem_error += elsol(el,1)*elsol(el,1);
            }
        }
        std::cout << "The interior Hdiv error is " << sqrt(trad_error) << std::endl;
        std::cout << "The SBFem Hdiv error is " << sqrt(sbfem_error) << std::endl;
    }
    HideSBFemVolumeElements();
    return cmesh_m;
}

/// @brief Compute the error as an average error per element and save as an elemental solution in "Error"
void ComputeErrorEstimator(TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator, REAL &GlobalError, REAL &ErrorH1, REAL &ErrorHDiv) {
    if (meshstyle == ESBFem) {
        AddSBFemVolumeElements();
    }
    cmesh->LoadReferences();
    int64_t nel_m = cmesh_m->NElements();
    ErrorEstimator.Resize(nel_m, 0.);
    ErrorEstimator = 0.;
    // TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh->ElementVec();
    TPZAdmChunkVector<TPZCompEl *> &elementvec_m = cmesh_m->ElementVec();
    TLaplaceExample1 *analytic = new TLaplaceExample1();
    analytic->fExact = TLaplaceExample1::ESquareRoot;
    TPZFMatrix<STATE> &MFelementsol = cmesh_m->ElementSolution();
    TPZFMatrix<STATE> &Elementsol = cmesh->ElementSolution();
    REAL total_est_trad = 0.;
    REAL total_est_sbfem = 0.;

    REAL H1_trad = 0.;
    REAL H1_sbfem = 0.;
    REAL HDiv_trad = 0.;
    REAL HDiv_sbfem = 0.;
    REAL Orthogonal_sbfem = 0;
    REAL Orthogonal_trad = 0;

    for (int64_t iel = 0; iel < nel_m; iel++) {
#ifdef PZ_LOG
    if(logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "index " << iel;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
        REAL H1_trad_el = 0.;
        REAL H1_sbfem_el = 0.;
        REAL HDiv_trad_el = 0.;
        REAL HDiv_sbfem_el = 0.;
        REAL Orthogonal_sbfem_el = 0;
        REAL Orthogonal_trad_el = 0;
        TPZCompEl *cell_m = elementvec_m[iel];
        if (!cell_m) continue;
        TPZGeoEl *gel_m = cell_m->Reference();
        if (!gel_m) continue;
        if (gel_m->Dimension() != 2) continue;

        TPZCompEl *cell = gel_m->Reference();
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cell);
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cell);
        if (!intel && !sbfem) DebugStop();
        TPZManVector<REAL,9> errorvec(5,0.);
        cell->EvaluateError(errorvec, false);
        int porder = 0;
        if (intel) {
            porder = intel->PreferredSideOrder(gel_m->NSides() - 1);
        } else if (sbfem) {
            porder = 8+porder;//18;
            //std::cout << "sbfem iel = " << iel << std::endl;
        }
        int matid = gel_m->MaterialId();

        TPZGeoElSide gelside(gel_m);
        auto intrule = gelside.CreateIntegrationRule(2 * 18);
        int intrulepoints = intrule->NPoints();
        TPZManVector<REAL, 4> intpoint(2, 0.);
        REAL weight = 0.;
        TPZFMatrix<REAL> jac, axe, jacInv;
        TPZManVector<REAL, 3> x(3, 0.);
        REAL detJac;
        TPZManVector<STATE, 3> solEx(1, 0.);
        TPZFNMatrix<9, STATE> dsolEx(2, 1, 0.);

        for (int int_ind = 0; int_ind < intrulepoints; ++int_ind) {
            intrule->Point(int_ind, intpoint, weight);
            gel_m->Jacobian(intpoint, jac, axe, detJac, jacInv);
            gel_m->X(intpoint, x);
            analytic->Solution(x, solEx, dsolEx);
            weight *= fabs(detJac);
            TPZInterpolationSpace *msp = dynamic_cast<TPZInterpolationSpace *>(cell);
            TPZManVector<STATE, 3> flux(3, 0.);
            msp->Solution(intpoint, 7, flux);

            TPZManVector<STATE, 3> sol(3, 0.);
            if (matid == volmat) {
                TPZMultiphysicsElement *msp_m = dynamic_cast<TPZMultiphysicsElement *>(cell_m);
                if (!msp_m) DebugStop();
                msp_m->Solution(intpoint, 1, sol);
            } else if (matid == sbfem_domain) {
                TPZSBFemVolume *msp_m = dynamic_cast<TPZSBFemVolume *>(cell_m);
                if (!msp_m) DebugStop();
                msp_m->Solution(intpoint, 1, sol);
                // std::cout << "flux " << flux << std::endl;
                // std::cout << "sol " << sol << std::endl;
            }

            REAL contr =
                ((flux[0] - sol[0]) * (flux[0] - sol[0]) + (flux[1] - sol[1]) * (flux[1] - sol[1])) * weight;
            ErrorEstimator[iel] += contr;
            REAL h1_loc = ((dsolEx(0,0)+flux[0])*(dsolEx(0,0)+flux[0])+(dsolEx(1,0)+flux[1])*(dsolEx(1,0)+flux[1]));
            REAL hdiv_loc = ((dsolEx(0,0)+sol[0])*(dsolEx(0,0)+sol[0])+(dsolEx(1,0)+sol[1])*(dsolEx(1,0)+sol[1]));
            if(sbfem) {
                H1_sbfem_el += h1_loc*weight;
                HDiv_sbfem_el += hdiv_loc*weight;
                Orthogonal_sbfem_el += ((dsolEx(0,0)+flux[0])*(dsolEx(0,0)+sol[0])+(dsolEx(1,0)+flux[1])*(dsolEx(1,0)+sol[1]))*weight;
            } else {
                H1_trad_el += h1_loc*weight;
                HDiv_trad_el += hdiv_loc*weight;
                Orthogonal_trad_el += ((dsolEx(0,0)+flux[0])*(dsolEx(0,0)+sol[0])+(dsolEx(1,0)+flux[1])*(dsolEx(1,0)+sol[1]))*weight;
            }
#ifdef PZ_LOG
            if(logger.isDebugEnabled()) {
                std::stringstream sout;
                sout << "i " << int_ind << " x " << x << " exact " << dsolEx(0,0) << " " << dsolEx(1,0) << std::endl;
                sout << "flux " << flux << " err loc H1 " << h1_loc  << " acc " << (H1_trad_el+H1_sbfem_el) << std::endl;
                sout << "dsol " << sol  << " err loc HDiv " << hdiv_loc  << " acc " << (HDiv_trad_el+HDiv_sbfem_el) << std::endl;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif

        } // loop over integration points
        
        H1_trad += H1_trad_el;
        H1_sbfem += H1_sbfem_el;
        HDiv_trad += HDiv_trad_el;
        HDiv_sbfem += HDiv_sbfem_el;
        Orthogonal_sbfem += Orthogonal_sbfem_el;
        Orthogonal_trad += Orthogonal_trad_el;
        
        H1_trad_el = sqrt(H1_trad_el);
        H1_sbfem_el = sqrt(H1_sbfem_el);
        HDiv_trad_el = sqrt(HDiv_trad_el);
        HDiv_sbfem_el = sqrt(HDiv_sbfem_el);

        if(sbfem) {
            total_est_sbfem += ErrorEstimator[iel];
        } else {
            total_est_trad += ErrorEstimator[iel];
        }
        REAL elest = sqrt(ErrorEstimator[iel]);
        ErrorEstimator[iel] = elest;
        {
            REAL diff = HDiv_trad_el+HDiv_sbfem_el - MFelementsol(iel,1);
            if(abs(diff) > 1.e-9) {
                std::cout << "HDiv error iel = " << iel << " diff = " << diff << std::endl;
            }
        }
        {
            REAL diff = H1_trad_el+H1_sbfem_el - Elementsol(iel,2);
            if(abs(diff) > 1.e-9) {
                std::cout << "H1 error iel = " << iel << " diff = " << diff << std::endl;
                cell->EvaluateError(errorvec, false);
                REAL test = sqrt(errorvec[2]);
            }
        }
    } // loop over cemsh_m elements
    GlobalError = sqrt(total_est_trad+total_est_sbfem);
    ErrorH1 = sqrt(H1_trad+H1_sbfem);
    ErrorHDiv = sqrt(HDiv_trad+HDiv_sbfem);

    {
        std::ofstream out("ErrorEstimator.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m, out, ErrorEstimator, "ErrorEstimator");
        std::cout << "Total estimated traditional error " << sqrt(total_est_trad) << std::endl;
        std::cout << "Total real traditional H1 error " << sqrt(H1_trad) << std::endl;
        std::cout << "Total real traditional HDiv error " << sqrt(HDiv_trad) << std::endl;
        std::cout << "***\n";
        std::cout << "Total estimated sbfem error " << sqrt(total_est_sbfem) << std::endl;
        std::cout << "Total real sbfem H1 error " << sqrt(H1_sbfem) << std::endl;
        std::cout << "Total real sbfem HDiv error " << sqrt(HDiv_sbfem) << std::endl;
        std::cout << "Orthogonal sbfem error " << Orthogonal_sbfem << std::endl;
        std::cout << "***\n";
        std::cout << "Total estimated error " << sqrt(total_est_trad+total_est_sbfem) << std::endl;
        std::cout << "Total real H1 error " << sqrt(H1_trad+H1_sbfem) << std::endl;
        std::cout << "Total real HDiv error " << sqrt(HDiv_trad+HDiv_sbfem) << std::endl;
        std::cout << "Orthogonal trad error " << Orthogonal_trad << std::endl;
        std::cout << "Teste " << H1_trad + HDiv_trad +2*Orthogonal_trad - total_est_trad << std::endl;
        std::cout << "Prager Synger = " << H1_trad + HDiv_trad - total_est_trad << std::endl;
    }
    HideSBFemVolumeElements();
}

/// @brief Compute the global error the norm of the vector ErrorEstimator
void ComputeGlobalError(TPZVec<REAL> &ErrorEstimator, REAL &GlobalError) {
    GlobalError = 0.;
    int64_t nel_m = ErrorEstimator.size();
    for (int64_t iel = 0; iel < nel_m; iel++) {
        GlobalError += ErrorEstimator[iel] * ErrorEstimator[iel];
    } // loop over the elements of the vector ErrorEstimator
    GlobalError = sqrt(GlobalError);
}

/// @brief Performe h-refinements at the computational meshes "cmesh" and "cmesh_m" based on the ErrorEstimator
void Hrefinement(TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator, TPZVec<REAL> &RefinementIndicator,
                 TPZVec<int> &porders) {
    int64_t nel_m = cmesh_m->NElements();
    auto gmesh = cmesh_m->Reference();
    int64_t nel = gmesh->NElements();
    RefinementIndicator.Resize(nel, 0.);
    RefinementIndicator = 0.;
    auto tolel = std::max_element(ErrorEstimator.begin(), ErrorEstimator.end());
    REAL tol = *tolel / 5.0;
    // divide all two dimensional elements with error larger than tol
    for (int64_t iel = 0; iel < nel_m; iel++) {
        TPZCompEl *cel = cmesh_m->Element(iel);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        // group elements have no reference (SBFem simulation)
        if (!gel) continue;
        // int64_t index = gel->Index();
        if (ErrorEstimator[iel] > tol) {
            RefinementIndicator[iel] = 1.0;
            TPZManVector<TPZGeoEl *> subel;
            DivideGeoEl(gel, subel);
        }
    } // loop over cmesh_m elements
    Smoothentrailingedgeelements(cmesh_m, RefinementIndicator);

    // refinar os elementos de contorno
    nel = gmesh->NElements();
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->Element(iel);
        if (!gel) continue;
        if (gel->HasSubElement()) continue;
        if (gel->Dimension() != 1) continue;
        int matid = gel->MaterialId();
        if (matid < 0) continue;
        if (matid == sbfem_highperm_hdiv || matid == sbfem_highperm_h1) continue;
        // if(gel->MaterialId() ==  boundmat) continue;
        // if(gel->MaterialId() ==  cutmat) continue;
        // Percorrer os Neighbors de gel e encontar os volumetricos que tem subelementos
        TPZGeoElSide gelside(gel);
        auto Neighbor = gelside.Neighbour();
        for (; Neighbor != gelside; Neighbor++) {
            if (!Neighbor.Element()) DebugStop();
            if (!Neighbor.HasSubElement()) continue;
            TPZManVector<TPZGeoEl *> subel;
            DivideGeoEl(gel, subel);
            if(matid == sbfem_skeleton){
                auto sbfem = gelside.HasNeighbour(sbfem_domain);
                if(!sbfem) DebugStop();
                if(!sbfem.HasSubElement()) {
                    DivideGeoEl(sbfem.Element(), subel);
                }
            }
            break;
        }
    } // loop over cmesh_m elements
    nel = gmesh->NElements();
    RefinementIndicator.Resize(nel, 0.);
    std::ofstream out1("HRefinementIndicator.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m, out1, RefinementIndicator, "RefinementIndicator");

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
void HPrefinement(TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator, int minh,
                  TPZVec<REAL> &RefinementIndicator, TPZVec<int> &porders) {

    int64_t nel_m = cmesh_m->NElements();

    auto gmesh = cmesh_m->Reference();
    gmesh->ResetReference();
    cmesh_m->LoadReferences();
    int64_t nel = gmesh->NElements();
    RefinementIndicator.Resize(nel, 0.);
    RefinementIndicator = 0.;
    auto tolel = std::max_element(ErrorEstimator.begin(), ErrorEstimator.end());
    std::cout << "max error " << *tolel << std::endl;
    REAL tol = *tolel / 5.0;
    // if the error of a trailing edge is larger than the tolerance, divide it
    DivideTrailingEdgeNeighbours(cmesh_m, ErrorEstimator, tol, RefinementIndicator, porders);
    // make sure all the trailing edge elements have the same level of refinement
    Smoothentrailingedgeelements(cmesh_m, RefinementIndicator);
    SmoothenGeometry(gmesh);
    nel = gmesh->NElements();
    RefinementIndicator.Resize(nel, 0.);
    int64_t psize = porders.size();
    porders.Resize(nel);
    for (int64_t i = psize; i < nel; i++) {
        porders[i] = defaultporder;
    }
    int sbfem_maxorder = 0;
    for (int64_t iel = 0; iel < nel_m; iel++) {
        TPZCompEl *cel = cmesh_m->Element(iel);
        if (!cel) continue;
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        TPZSBFemVolume *sbvol = dynamic_cast<TPZSBFemVolume *>(cel);
        TPZSBFemElementGroup *sbfem_group = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (!mfcel && !sbvol && !sbfem_group) {
            DebugStop();
        }
        if (mfcel) {
            TPZInterpolationSpace *msp = dynamic_cast<TPZInterpolationSpace *>(mfcel->Element(0));
            if (!msp) DebugStop();
            TPZGeoEl *gel = cel->Reference();
            int64_t index = gel->Index();
            int porder = porders[index];
            REAL Elerror = ErrorEstimator[iel];
            if (Elerror > tol) {
                RefinementIndicator[iel] = 1.0;
                if (gel->HasSubElement()) {
                    // only trailing edge neighbours are divided
                    TPZManVector<TPZGeoEl *> subel;
                    DivideGeoEl(gel, subel);
                    for (auto sub : subel) {
                        porders[sub->Index()] = porder;
                    }
                } else {
                    // apply p refinement
                    porders[index] = porder + 1;
                }
            } else {
                if (gel->HasSubElement()) {
                    // the element was divided to satisfy the two on one constraint
                    TPZManVector<TPZGeoEl *> subel;
                    DivideGeoEl(gel, subel);
                    for (auto sub : subel) {
                        porders[sub->Index()] = porder;
                    }
                } else {
                    porders[index] = porder;
                }
            }
        } else if (sbvol) {
            REAL error = ErrorEstimator[iel];
            TPZGeoEl *gel = cel->Reference();
            int64_t index = gel->Index();
            if (error > tol) {
                RefinementIndicator[iel] = 1.0;
                porders[index] = SBFemOrder + 1;
                sbfem_maxorder = SBFemOrder + 1;
            } else {
                porders[index] = SBFemOrder;
            }
        } else {
            std::cout << "I don't understand\n";
        }
    } // loop over cemsh_m elements
    // set all sbfem volume elements to the same order
    if (meshstyle == ESBFem) {
        if (sbfem_maxorder > SBFemOrder) {
            std::cout << "Increasing the order of the SBFem elements from " << SBFemOrder << " to " << sbfem_maxorder
                      << std::endl;
            SBFemOrder = sbfem_maxorder;
        }
        // update the porders data structure
        for (int64_t iel : sbfem_elements) {
            porders[iel] = SBFemOrder;
        }
    }
    {
        std::ofstream out("HRefinementIndicator.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m, out, RefinementIndicator, "RefinementIndicator");
        TPZVec<double> porderscopy(cmesh_m->NElements(), 0);
        for (int64_t i = 0; i < cmesh_m->NElements(); i++) {
            TPZCompEl *cel = cmesh_m->Element(i);
            if (!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if (!gel) continue;
            porderscopy[i] = porders[gel->Index()];
        }
        std::ofstream out2("porders.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m, out2, porderscopy, "porders");
    }
    {
        std::ofstream out("porders.txt");
        for (auto p : porders)
            out << p << std::endl;
    }
    {
        std::ofstream out("gmeshrefined.txt");
        gmesh->Print(out);
    }
    {
        std::ofstream out2("gmesh_HPrefinement.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2, porders, "Porders", true);
    }
    if (gmesh->NElements() != nel) DebugStop();
    if (porders.size() != gmesh->NElements()) DebugStop();

    {
        std::ofstream out("gmeshrefined.txt");
        gmesh->Print(out);
    }
    {
        std::ofstream out1("RefinementIndicator.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m, out1, RefinementIndicator, "RefinementIndicator");
    }
    {
        std::ofstream out2("gmesh_HPrefinement.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2, porders, "Porders", true);
    }
}

/// @brief Smoothen the volumetric elements around the trailingedge by adopting the same level of refinement
void Smoothentrailingedgeelements(TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &RefinementIndicator) {
    auto gmesh = cmesh_m->Reference();
    if (0) {
        std::ofstream out("gmeshbefore.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }

    int64_t nel = gmesh->NElements();
    TPZGeoEl *gel_TrailingEdge = gmesh->Element(trailingedge_element_index);
    // if(gel_TrailingEdge->HasSubElement()) continue;
    if (gel_TrailingEdge->MaterialId() != trailingedgemat) DebugStop();

    TPZGeoElSide TrailingSide = TPZGeoElSide(gel_TrailingEdge);
    TPZGeoElSide Neighbor = TrailingSide.Neighbour();
    TPZGeoEl *GeoNeighbor;
    int maxlevel = 0;
    for (; Neighbor != gel_TrailingEdge; Neighbor++) {
        if (!Neighbor.Element()) DebugStop();
        if (Neighbor.Element()->HasSubElement()) continue;
        if (Neighbor.Element()->MaterialId() != volmat) continue;
        GeoNeighbor = Neighbor.Element();
        if (maxlevel <= GeoNeighbor->Level()) {
            maxlevel = GeoNeighbor->Level();
        }
    }
    Neighbor = gel_TrailingEdge->Neighbour(0);
    for (; Neighbor != gel_TrailingEdge; Neighbor++) {
        if (!Neighbor.Element()) DebugStop();
        if (Neighbor.Element()->HasSubElement()) continue;
        if (Neighbor.Element()->MaterialId() != volmat) continue;
        GeoNeighbor = Neighbor.Element();
        if (GeoNeighbor->Level() < maxlevel) {
            auto iel_GeoNeighbor = GeoNeighbor->Index();
            RefinementIndicator[iel_GeoNeighbor] = 1.0;
            TPZManVector<TPZGeoEl *> subel_GeoNeighbor;
            DivideGeoEl(GeoNeighbor, subel_GeoNeighbor);
            RefinementIndicator.Resize(gmesh->NElements(), 0.);
        }
    }
    if (0) {
        std::ofstream out("gmeshafter.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
}

#include "tpzgeoblend.h"
#include "tpzquadraticline.h"
#include "tpzquadraticquad.h"
/// @brief create a collapsed quadratic quad element
/// @param sn the index of the singular node
/// @param edge1 the index of the first edge node
/// @param edge2 the index of the second edge node
/// @param gmesh the geometric mesh
static int64_t CreateCollapsedQuad(int64_t sn, int64_t edge1, int64_t edge2, bool blend, TPZGeoMesh *gmesh) {
    TPZManVector<int64_t, 4> nodeindices(4);
    nodeindices[0] = edge1;
    nodeindices[1] = edge2;
    nodeindices[2] = sn;
    nodeindices[3] = sn;
    int64_t index;
    int matid = volmat;
    if (blend) {
        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>>(nodeindices, matid, *gmesh, index);
    } else {
        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeindices, matid, *gmesh, index);
    }
    return index;
}


/// @brief divide a geometric element. If it is a trailing edge element, set the refinement pattern first
void DivideGeoEl(TPZGeoEl *gel, TPZVec<TPZGeoEl *> &subels) {
    subels.resize(0);
    if (!gel) return;
    if (gel->MaterialId() < 0) {
        DebugStop();
    }
    int64_t gelindex = gel->Index();
    if (gel->Dimension() == 0) return;
    if (gel->Type() != EQuadrilateral) {
        gel->Divide(subels);
        return;
    }
    bool iscollapsed = false;
    for (int side = 4; side < 8; side++) {
        int n1 = side - 4;
        int n2 = (n1 + 1) % 4;
        if (gel->NodeIndex(n1) == gel->NodeIndex(n2)) {
            if (side != 6) DebugStop();
            iscollapsed = true;
            if(meshstyle == ESBFem) {
                gel->SetRefPattern(refpattern_sbfem);
            } else {
                //gel->SetRefPattern(refpattern_collapsed);
            }
        }
    }
    gel->Divide(subels);
    
    if(iscollapsed) {
        // set all nodes of the collapsed element to the singular node
        // update the neighbouring information of the elements
        int64_t singular = -1;
        TPZGeoMesh *gmesh = gel->Mesh();
        TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
        TPZGeoElSide trailingedge(trailingedge_element);
        singular = trailingedge_element->NodeIndex(0);
        for (auto sub : subels) {
            int64_t node1 = sub->NodeIndex(2);
            int64_t node2 = sub->NodeIndex(3);
            if(node1 == singular || node2 == singular) {
                if(gel->NodeIndex(2) != singular) DebugStop();
                if(gel->NodeIndex(3) != singular) DebugStop();
                sub->SetNodeIndex(2,singular);
                sub->SetNodeIndex(3,singular);
                {
                    TPZGeoElSide gelside(sub,2);
                    gelside.RemoveConnectivity();
                    gelside.SetConnectivity(trailingedge);
                }
                {
                    TPZGeoElSide gelside(sub,3);
                    gelside.RemoveConnectivity();
                    gelside.SetConnectivity(trailingedge);
                }
                // this breaks the code for LowerLevelCompElementList
                if(0)
                {
                    TPZGeoElSide gelside(sub,6);
                    gelside.RemoveConnectivity();
                    TPZGeoElSide fatherside(gel,6);
                    gelside.SetConnectivity(fatherside);
                }
            }


        }
    }
    int ns = subels.size();
    for (int is = 0; is < ns; is++) {
        if (!IsCounterClockwise(subels[is])) DebugStop();
    }
}

/// @brief Substitute the trailing edge quadrilateral elements with colapsed quadrilateral elements with or without
/// quarterpoint elements
void Changetrailingedgeelements(TPZGeoMesh *gmesh) {
    int64_t nel = gmesh->NElements();
    // identify the point element and trailing edge node
    int64_t singular = -1;
    TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
    TPZGeoElSide trailingedge(trailingedge_element);
    singular = trailingedge_element->NodeIndex(0);
    if (singular == -1) DebugStop();
    TPZStack<TPZGeoElSide> neighbours;
    // create nodes along the line from the singular node to the edge nodes
    TPZManVector<REAL, 3> x0(3);
    gmesh->NodeVec()[singular].GetCoordinates(x0);
    // for all neighbours of the trailing edge element
    for (TPZGeoElSide neigh = trailingedge.Neighbour(); neigh != trailingedge; neigh = neigh.Neighbour()) {
        TPZGeoEl *gel = neigh.Element();
        if (!gel) DebugStop();
        if (gel->HasSubElement()) continue;
        neighbours.Push(neigh);
    }
    // create the new elements
    std::set<int64_t> sbfem_created;
    for (int64_t el = 0; el < neighbours.size(); el++) {
        TPZGeoEl *gel = neighbours[el].Element();
        if (!gel) DebugStop();
        if (gel->Dimension() != 2) continue;
        int nn = gel->NNodes();
        int64_t singular = gel->NodeIndex(neighbours[el].Side());

        // create collapsed elements splitting each quadrilateral element into 2 quads
        for (int in = 1; in < nn - 1; in++) {
            int side1d = -1;
            if (in == 1)
                side1d = neighbours[el].Side() + 4;
            else if (in == 2)
                side1d = neighbours[el].Side() + 3;
            else
                DebugStop();
            if (side1d < 4) side1d += 4;
            TPZGeoElSide gelside(neighbours[el].Element(), side1d);
            int locindex = (neighbours[el].Side() + in) % nn;
            int64_t edge1 = gel->NodeIndex(locindex);
            int64_t edge2 = gel->NodeIndex((locindex + 1) % nn);
            int64_t newindex;
            if (nn == 4) {
                // create a geometrically consistent quad element
                newindex = CreateCollapsedQuad(singular, edge1, edge2, false, gmesh);
                sbfem_created.insert(newindex);
            } else {
                DebugStop();
            }
        }
        if (nn > 1) gel->SetMaterialId(-gel->MaterialId());
    }
    gmesh->BuildConnectivity();
    if (1) {
        std::ofstream out("gmesh_Changetrailingedgeelements.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmesh_Changetrailingedgeelements.txt");
        gmesh->Print(out2);
    }
}

void DivideTrailingEdgeNeighbours(TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator, REAL tol,
                                  TPZVec<REAL> &RefinementIndicator, TPZVec<int> &porders) {
    auto gmesh = cmesh_m->Reference();
    int64_t nel = gmesh->NElements();
    // compute the set of elements neighbouring the trailing edge
    std::set<int64_t> trailel;
    {
        TPZGeoEl *gel = gmesh->Element(trailingedge_element_index);
        TPZGeoElSide gelside(gel);
        for (TPZGeoElSide neighbour = gelside.Neighbour(); neighbour != gelside; neighbour = neighbour.Neighbour()) {
            TPZGeoEl *neighgel = neighbour.Element();
            if (neighgel->Reference() && neighgel->MaterialId() == volmat) {
                trailel.insert(neighgel->Index());
            }
        }
    }
    // divide the trailing edge element neighbours if necessary
    for (int64_t el : trailel) {
        TPZGeoEl *gel = gmesh->Element(el);
        TPZCompEl *cel = gel->Reference();
        int64_t index = cel->Index();
        if (ErrorEstimator[index] > tol && meshstyle != ESBFem) {
            RefinementIndicator[index] = 1.;
            TPZManVector<TPZGeoEl *> subel;
            DivideGeoEl(gel, subel);
        } else if (ErrorEstimator[index] > tol && meshstyle == ESBFem) {
            RefinementIndicator[index] = 1.;
            porders[index]++;
        }
    }
}

static void AllSubElements(TPZGeoElSide &gelside, TPZStack<TPZGeoElSide> &smallest, int dimension) {
    // we push in the smallest datastructure only element/sides that have no subelements
    TPZStack<TPZGeoElSide> locsmall;
    // locsmall will be empty if the element is not divided
    gelside.GetSubElements2(locsmall);
    for (int is = 0; is < locsmall.size(); is++) {
        if (locsmall[is].Dimension() != dimension) continue;
        if (locsmall[is].HasSubElement()) {
            // recursive call to add the subelements of the subelements
            AllSubElements(locsmall[is], smallest, dimension);
        } else {
            // we have a "leaf" subelement
            smallest.Push(locsmall[is]);
        }
    }
}

void SmoothenGeometry(TPZGeoMesh *gmesh) {
    bool changed = true;
    while (changed) {
        changed = false;
        int64_t nel = gmesh->NElements();
        // enforce the two on one constraint
        for (int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (gel->MaterialId() < 0) continue;
            if (gel->Dimension() != 2) continue;
            if (gel->HasSubElement()) continue;
            bool shouldrefine = false;
            for (int is = gel->FirstSide(1); is < gel->NSides() - 1; is++) {
                TPZGeoElSide gelside(gel, is);
                for (TPZGeoElSide neighbour = gelside.Neighbour(); neighbour != gelside;
                     neighbour = neighbour.Neighbour()) {
                    TPZGeoEl *neighgel = neighbour.Element();
                    if (!neighgel->HasSubElement()) continue;
                    TPZStack<TPZGeoElSide> subel;
                    AllSubElements(neighbour, subel, 1);
                    // if there are 2 subelements or less the current element satisfies the
                    // two on one rule
                    if (subel.size() <= 2) continue;
                    shouldrefine = true;
                    break;
                }
                if (shouldrefine) break;
            }
            if (shouldrefine) {
                TPZManVector<TPZGeoEl *> subel;
                DivideGeoEl(gel, subel);
                changed = true;
            }
        }
        // divide an element that has more than half of its neighbours divided
        nel = gmesh->NElements();
        for (int el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (gel->MaterialId() < 0) continue;
            if (gel->HasSubElement()) continue;
            int nside1 = gel->NSides(1);
            int nneighdiv = 0;
            for (int is = gel->FirstSide(1); is < gel->FirstSide(2); is++) {
                TPZGeoElSide gelside(gel, is);
                bool found = false;
                for (TPZGeoElSide neighbour = gelside.Neighbour(); neighbour != gelside;
                     neighbour = neighbour.Neighbour()) {
                    TPZGeoEl *neighgel = neighbour.Element();
                    if (neighgel->MaterialId() < 0) continue;
                    int neighside = neighbour.Side();
                    if (neighgel->HasSubElement() && neighgel->NSideSubElements(neighside) > 1) found = true;
                }
                if (found) nneighdiv++;
            }
            if (nneighdiv > nside1 / 2) {
                TPZManVector<TPZGeoEl *> subel;
                DivideGeoEl(gel, subel);
                changed = true;
            }
        }
    }
    // adjust the one dimensional elements
    // they should have no divided neighbours
    changed = true;
    while (changed) {
        changed = false;
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (gel->HasSubElement()) continue;
            if (gel->Dimension() != 1) continue;
            // Percorrer os Neighbors de gel e encontar os volumetricos que tem subelementos
            TPZGeoElSide gelside(gel);
            auto Neighbor = gelside.Neighbour();
            for (; Neighbor != gelside; Neighbor++) {
                if (!Neighbor.Element()) DebugStop();
                if (Neighbor.HasSubElement() && Neighbor.NSubElements() > 1) {
                    TPZManVector<TPZGeoEl *> subel;
                    DivideGeoEl(gel, subel);
                    changed = true;
                    break;
                }
            }
        }
    }
}

void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
// Define a função do tipo void chamada PrintResults, que recebe como parâmetros TPZLinearAnalysis &an e  TPZCompMesh
// *cmesh
{

    std::cout << "--------- Post Process ---------" << std::endl;
    // printa para o usuário "--------- Post Process ---------" indicando o início da fase de pós-processamento.
    TPZSimpleTimer postProc("Post processing time");
    // declara uma variável chamada postProc, do tipo TPZSimpleTimer, chamando um construtor com uma string como
    // argumento, igual a "Post processing time". inicializa um temporizador chamado postProc que será usado para medir
    // o tempo gasto no pós-processamento.
    const std::string plotfile = "postprocess";
    // define o nome base do arquivo de saída para o pós-processamento. O nome base é "postprocess".
    constexpr int vtkRes{3};
    // define a variável do tipo inteiro denominada vtkRes, do tipo constexpr, que significa que é uma expressão
    // constante, ou seja,  vtkRes é um valor constante e não pode ser alterado. Ainda, {0} indica o valor associado a
    // essa constante, e portanto não será alterado, com valor determinado na hora de compilação. define a resolução
    // para o formato de arquivo VTK. Neste caso, a resolução é definida como 0, o que geralmente significa que a
    // resolução será automática.
    TPZVec<std::string> fields = {"Pressure", "Flux"};
    // nesse conjunto de linhas de código, temos que TPZVec é uma estrutura do tipo vetor que contém como argumento uma
    // variável chamda "fields" que é uma lista de strings, que, pelo que se chamam, são relacionadas à pressão e ao
    // fluxo. cria um vetor de strings chamado fields que contém os nomes dos campos que serão pós-processados. Neste
    // caso, os campos incluem "Pressure" (pressão) e "Flux" (fluxo). Esses campos representam propriedades do problema
    // que desejamos visualizar após a simulação.
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    // essa linha de código declara uma variável chamada vtk do tipo auto, o que significa que o compilador irá deduzir
    // o tipo que ela terá a depender do que ela é igual. No caso, ela é igual a função TPZVTKGenerator, de parâmetros
    // cmesh, fields, plotfile, vtkRes. cria um objeto vtk da classe TPZVTKGenerator, que é usado para gerar arquivos
    // VTK a partir dos dados da malha computacional cmesh. Os argumentos passados para o construtor incluem a malha
    // computacional, os campos a serem pós-processados, o nome base do arquivo de saída (plotfile) e a resolução VTK
    // (vtkRes).
    vtk.SetNThreads(0);
    // define o número de threads a serem usadas durante o pós-processamento. A variável global_nthread provavelmente
    // contém o número desejado de threads.
    vtk.Do();
    // inicia o processo de geração dos arquivos VTK. Esta função gera arquivos de saída contendo informações sobre os
    // campos especificados na malha computacional.
    std::cout << "Total time = " << postProc.ReturnTimeDouble() / 1000. << " s" << std::endl;
    // imprime o tempo gasto no pós-processamento, convertido para segundos.

    return;
    // a função é concluída e retorna.
}

/// @brief Adjust the system of equations: build the second collum of the rhs for Beta=1 and remove the Beta equations
/// from the structural matrix
/// @return Adjusted rhs and structural matrices: TPZFMatrix<STATE> &rhs, TPZSYsmpMatrix<STATE> &spmat.
void AdjustSystemofEquations(TPZLinearAnalysis &an, TPZCompMesh *cmesh, int64_t &newcon) {
    /// Geting the structural matrix and the rhs
    cmesh->ConnectVec()[newcon].Print(*cmesh);
    TPZConnect &c = cmesh->ConnectVec()[newcon];
    auto seq = c.SequenceNumber();
    auto eqbeta = cmesh->Block().Position(seq);
    std::cout << "eqbeta " << eqbeta << std::endl;
    auto solver = an.Solver();
    TPZMatrixSolver<STATE> *matsolv = dynamic_cast<TPZMatrixSolver<STATE> *>(solver);
    if (!matsolv) DebugStop();
    auto mat = matsolv->Matrix();
#if defined(USING_MKL) || defined(USING_EIGEN)
    TPZSYsmpMatrix<STATE> *spmat = dynamic_cast<TPZSYsmpMatrix<STATE> *>(mat.operator->());
    if (!spmat) DebugStop();
#else
    TPZSkylMatrix<STATE> *spmat = dynamic_cast<TPZSkylMatrix<STATE> *>(mat.operator->());
    if (!spmat) DebugStop();
#endif
    TPZFMatrix<STATE> &rhs = an.Rhs();

    /// Variables
    auto nrows = spmat->Rows();
    auto ncols = spmat->Cols();
    rhs.Resize(nrows, 2);
    for (int64_t row = 0; row < nrows; row++) {
        rhs(row, 1) = 0.0;
    }

    for (int64_t row = 0; row < nrows; row++) {
        if (row == eqbeta) {
            rhs(row, 1) = 1.0;
            rhs(row, 0) = 0.;
            spmat->PutVal(row, row, 1.);
        } else {
            rhs(row, 1) = -spmat->GetVal(row, eqbeta);
            spmat->PutVal(row, eqbeta, 0.);
        }
    }
}

/// @brief returns true if the element is counter clockwise
bool IsCounterClockwise(TPZGeoEl *gel) {
    if (gel->Dimension() != 2) DebugStop();
    TPZGeoElSide gelside(gel);
    TPZManVector<REAL, 3> qsi(2), ax1(3), ax2(3), ax3(3);
    TPZFNMatrix<6, REAL> gradx(3, 2);
    gelside.CenterPoint(qsi);
    gel->GradX(qsi, gradx);
    for (int i = 0; i < 3; i++) {
        ax1[i] = gradx(i, 0);
        ax2[i] = gradx(i, 1);
    }
    if (!IsZero(ax1[2]) || !IsZero(ax2[2])) DebugStop();
    Cross(ax1, ax2, ax3);
    if (ax3[2] > 0.) return true;
    return false;
}

/// @brief make the elements counter clockwise
void MakeCounterClockwise(TPZGeoEl *gel) {
    if (IsCounterClockwise(gel)) return;
    if (gel->Type() != EQuadrilateral) DebugStop();
    TPZManVector<TPZGeoElSide, 9> gelsides(gel->NSides());
    for (int is = 0; is < gel->NSides(); is++) {
        TPZGeoElSide thisside(gel, is);
        TPZGeoElSide neigh = thisside.Neighbour();
        if (neigh != thisside) gelsides[is] = neigh;
    }
    gel->RemoveConnectivities();
    // flip nodes 1 and 3
    {
        int64_t node1 = gel->NodeIndex(1);
        int64_t node3 = gel->NodeIndex(3);
        gel->SetNodeIndex(1, node3);
        gel->SetNodeIndex(3, node1);
    }
    // connectivities 0->0 1->3 2->2 3->1 4->7 5->6 6->5 7->4 8->8
    for (int is = 0; is < 9; is++) {
        int mapside[] = {0, 3, 2, 1, 7, 6, 5, 4, 8};
        TPZGeoElSide gelside(gel, is);
        if (gelsides[is]) {
            gelside.SetConnectivity(gelsides[is]);
        } else {
            gelside.SetConnectivity(gelside);
        }
    }
#ifdef PZDEBUG
    if (IsCounterClockwise(gel) == false) DebugStop();
    for (int is = 0; is < gel->NSides(); is++) {
        std::set<int64_t> nodes;
        for (int i = 0; i < gel->NSideNodes(is); i++)
            nodes.insert(gel->SideNodeIndex(is, i));
        TPZGeoElSide gelside(gel, is);
        for (auto neigh = gelside.Neighbour(); neigh != gelside; neigh++) {
            for (int in = 0; in < neigh.NSideNodes(); in++) {
                if (nodes.find(neigh.SideNodeIndex(in)) == nodes.end()) {
                    DebugStop();
                }
            }
        }
    }
#endif
}

TPZAutoPointer<TPZRefPattern> CreateRefPattern_collapsed() {
    char buf[] = "6     3  "
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

TPZAutoPointer<TPZRefPattern> CreateRefPattern_sbfem() {
    char buf[] = "6     3  "
                 "1010       UnifQua	"
                 "-1.    -1.     0. "
                 " 1.    -1.     0. "
                 " 1.     1.     0. "
                 "-1.     1.     0.	"
                 " 0.    -1.     0. "
                 " 0.     1.     0. "
                 " 3     4     0     1     2     3 "
                 " 3     4     0     4     5     3 "
                 " 3     4     4     1     2     5 ";
    std::istringstream str(buf);
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(str);
    return refpat;
}

/// @brief Change the elements the touch the trailing edge to quarterpoint elements
void CreateQuarterPointElements(TPZGeoMesh *gmesh) {
    TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
    TPZGeoElSide trailingedge(trailingedge_element);
    std::map<TPZGeoEl *, int> connected_elements;
    for (TPZGeoElSide neigh = trailingedge.Neighbour(); neigh != trailingedge; neigh++) {
        TPZGeoEl *neighgel = neigh.Element();
        if (neighgel->HasSubElement()) continue;
        if (neighgel->Dimension() == 2) {
            if (neighgel->NodeIndex(2) == neighgel->NodeIndex(3)) {
                connected_elements[neighgel] = 6;
            }
        } else if (neighgel->Dimension() == 1) {
            connected_elements[neighgel] = neigh.Side();
        }
    }
    for (auto it : connected_elements) {
        TPZGeoEl *gel = it.first;
        int side = it.second;
        gel = TPZChangeEl::ChangeToQuadratic(gel->Mesh(), gel->Index());
        TPZChangeEl::ChangeToQuarterPoint(gel->Mesh(), gel->Index(), side);
    }
}

/// @brief remove computational elements generated for SBFem simulation
void RemoveSBFemElements(TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *cmesh_m, TPZGeoMesh *gmesh) {
    if (sbfem_groupHdiv) {
        sbfem_groupHdiv->Unwrap();
        sbfem_groupHdiv = 0;
    }
    if (sbfem_groupH1) {
        sbfem_groupH1->Unwrap();
        sbfem_groupH1 = 0;
    }
    auto Remove = [gmesh](TPZCompMesh *cmesh) {
        int64_t nel = cmesh->NElements();
        // std::cout << "gmesh nelements " << gmesh->NElements() << std::endl;
        for (int64_t iel = 0; iel < nel; iel++) {
            TPZCompEl *cel = cmesh->Element(iel);
            if (!cel) continue;
            int64_t index = cel->ReferenceIndex();
            if (index < 0) continue;
            TPZGeoEl *gel = cel->Reference();
            int matid = gel->MaterialId();
            if (index >= gmesh->NElements()) {
                // std::cout << "gel index " << index << " matid " << gel->MaterialId() << " is out of range " << "
                // celindex " << cel->Index() << "\n";
                int64_t celindex = cel->Index();
                delete cel;
                cel = cmesh->Element(celindex);
                if (cel) DebugStop();
            }
        }
    };
    Remove(cmesh);
    Remove(cmesh_m);
    auto meshvec = cmesh_m->MeshVector();
    for (auto mesh : meshvec) {
        Remove(mesh);
    }
}

/// @brief print the geometry of the trailing edge elements
void PrintTrailingEdgeElements(TPZGeoMesh *gmesh) {
    TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
    TPZGeoElSide trailingedge(trailingedge_element);
    std::set<TPZGeoEl *> connected_elements;
    for (TPZGeoElSide neigh = trailingedge.Neighbour(); neigh != trailingedge; neigh++) {
        TPZGeoEl *neighgel = neigh.Element();
        if (neighgel->HasSubElement()) continue;
        if (neighgel->Dimension() == 2) {
            if (neighgel->NodeIndex(2) == neighgel->NodeIndex(3)) {
                connected_elements.insert(neighgel);
            }
        } else if (neighgel->Dimension() == 1) {
            connected_elements.insert(neighgel);
        }
    }
    for (auto it : connected_elements) {
        TPZGeoEl *gel = it;
        std::cout << "Element " << gel->Index() << std::endl;
        int nnodes = gel->NNodes();
        for (int i = 0; i < nnodes; i++) {
            TPZManVector<REAL, 3> x(3);
            gel->NodePtr(i)->GetCoordinates(x);
            std::cout << "Node " << i << " " << x << std::endl;
        }
    }
}

/// @brief Change the elements the touch the trailing edge to H1 SBFEM elements
void CreateH1SBFEMelements(TPZCompMesh *cmesh) {
    cmesh->LoadReferences();
    TPZMaterialT<STATE> *mat2d = dynamic_cast<TPZMaterialT<STATE> *>(cmesh->FindMaterial(volmat));
    if (!mat2d) DebugStop();
    TPZDarcyFlow *darcy = dynamic_cast<TPZDarcyFlow *>(mat2d);
    if (!darcy) DebugStop();
    TPZDarcyFlow *darcysbfem = new TPZDarcyFlow(*darcy);
    darcysbfem->SetId(sbfem_domain);
    cmesh->InsertMaterialObject(darcysbfem);
    TPZDarcyFlow *mat1d = new TPZDarcyFlow(sbfem_highperm_h1, 1);
    mat1d->SetConstantPermeability(1.e9);
    cmesh->InsertMaterialObject(mat1d);
    TPZGeoMesh *gmesh = cmesh->Reference();
    TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
    TPZGeoElSide trailingedge(trailingedge_element);
    std::set<TPZGeoEl *> connected_elements;
    for (TPZGeoElSide neigh = trailingedge.Neighbour(); neigh != trailingedge; neigh++) {
        TPZGeoEl *neighgel = neigh.Element();
        if (neighgel->HasSubElement()) continue;
        if (neighgel->MaterialId() < 0) continue;
        if (connected_elements.find(neighgel) != connected_elements.end()) continue;
        if (neighgel->Dimension() == 2) {
            if (neighgel->NodeIndex(2) == neighgel->NodeIndex(3)) {
                TPZCompEl *cel = neighgel->Reference();
                if (cel) DebugStop();
                connected_elements.insert(neighgel);
                TPZGeoElSide gelside(neighgel, 4);
                TPZGeoElSide neigh = gelside.HasNeighbour(sbfem_skeleton);
                if (!neigh) {
                    DebugStop();
                }
                if (neigh.Element()->Reference() == 0) {
                    DebugStop();
                }
            } else {
                DebugStop();
            }
        } else if (neighgel->Dimension() == 1 && neighgel->MaterialId() == sbfem_highperm_h1) {
            // create a zero dimensional element
            connected_elements.insert(neighgel);
            TPZGeoElSide gelside(neighgel, 0);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                TPZGeoEl *neighgel = neighbour.Element();
                if (neighgel->MaterialId() == sbfem_skeleton && neighgel->Dimension() == 0) {
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
            if (neighbour == gelside) DebugStop();
            if (neighbour.Element()->Dimension() != 0) DebugStop();
            if (!neighbour.Element()->Reference()) DebugStop();
        }
    }
    TPZSBFemElementGroup *elgr = new TPZSBFemElementGroup(*cmesh);
    sbfem_groupH1 = elgr;
    TPZStack<TPZSBFemVolume *> sbvols;
    for (auto it : connected_elements) {
        TPZGeoEl *gel = it;
        int matid = gel->MaterialId();
        if (matid != sbfem_domain && matid != sbfem_highperm_h1) DebugStop();
        TPZSBFemVolume *sbvol = new TPZSBFemVolume(*cmesh, gel);
        sbvols.Push(sbvol);
        if(gel->Dimension() == 2) {
            TPZGeoElSide gelside(gel, 4);
            TPZGeoElSide neigh = gelside.HasNeighbour(sbfem_skeleton);
            TPZGeoEl *skel = neigh.Element();
            if (skel->MaterialId() != sbfem_skeleton) DebugStop();
            // create an H1 element
            TPZCompEl *cskel = skel->Reference();
            if (!cskel) DebugStop();
            TPZInterpolatedElement *icskel = dynamic_cast<TPZInterpolatedElement *>(cskel);
            if (!cskel) DebugStop();
            int sideorder = icskel->Connect(2).Order();
            if (sideorder < SBFemOrder) icskel->SetSideOrder(2, SBFemOrder);
            // std::cout << "gel " << gel->Index() << " is linear " << gel->IsLinearMapping() << std::endl;
            // std::cout << "skel is linear " << skel->IsLinearMapping() << std::endl;
            if (!gel->IsLinearMapping()) {
                // gel->Print(std::cout);
                // std::ofstream out("gmesh.txt");
                // gmesh->Print(out);
                DebugStop();
            }
            sbvol->SetSkeleton(cskel->Index());
        } else if(gel->Dimension() == 1) {
            TPZGeoElSide gelside(gel, 0);
            TPZGeoElSide neigh = gelside.Neighbour();
            while (neigh != gelside) {
                TPZGeoEl *neighgel = neigh.Element();
                if (neighgel->MaterialId() == sbfem_skeleton && neighgel->Dimension() == 0) {
                    break;
                }
                neigh = neigh.Neighbour();
            }
            if (neigh.Element()->MaterialId() != sbfem_skeleton) DebugStop();
            // create an H1 element
            TPZCompEl *cskel = neigh.Element()->Reference();
            if (!cskel) DebugStop();
            sbvol->SetSkeleton(cskel->Index());
        }
        elgr->AddElement(sbvol);
    }
    for (auto sbvol : sbvols) {
        sbvol->SetElementGroupIndex(elgr->Index());
    }
    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    TPZElementMatrixT<STATE> ek, ef;
    elgr->CalcStiff(ek, ef);
    if (0) {
        ek.fMat.Print("ek = ", std::cout, EMathematicaInput);
    }
}

/// @brief Change the elements the touch the trailing edge to Hdiv SBFEM elements
void CreateHdivSBFEMelements(TPZMultiphysicsCompMesh *m_cmesh, TPZAnalyticSolution *analyticSol) {
    TPZCompMesh *cmesh = m_cmesh;
    cmesh->LoadReferences();
    // cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(2);
    // cmesh->ApproxSpace().SetHDivFamily(HDivFamily::EHDivKernel);
    TPZGeoMesh *gmesh = cmesh->Reference();

    // insert the skeleton material
    if (cmesh->FindMaterial(sbfem_skeleton) == 0) {
        DebugStop();
    }
    TPZMixedDarcyFlow *mat1d = new TPZMixedDarcyFlow(sbfem_highperm_hdiv, 1);
    cmesh->InsertMaterialObject(mat1d);
    mat1d->SetConstantPermeability(1.e-9);
    TPZMaterialT<STATE> *mat2d = dynamic_cast<TPZMaterialT<STATE> *>(cmesh->FindMaterial(volmat));
    if (!mat2d) DebugStop();
    TPZMixedDarcyFlow *darcy = dynamic_cast<TPZMixedDarcyFlow *>(mat2d);
    if (!darcy) DebugStop();
    TPZMixedDarcyFlow *darcysbfem = new TPZMixedDarcyFlow(*darcy);
    darcysbfem->SetId(sbfem_domain);
    cmesh->InsertMaterialObject(darcysbfem);

    TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
    TPZGeoElSide trailingedge(trailingedge_element);
    std::set<TPZGeoEl *> connected_elements;
    // create the skeleton elements
    for (TPZGeoElSide neigh = trailingedge.Neighbour(); neigh != trailingedge; neigh++) {
        TPZGeoEl *neighgel = neigh.Element();
        if (neighgel->HasSubElement()) continue;
        if (neighgel->MaterialId() < 0) continue;
        if (connected_elements.find(neighgel) != connected_elements.end()) continue;
        if (neighgel->Dimension() == 2) {
            if (neighgel->NodeIndex(2) == neighgel->NodeIndex(3)) {
                TPZCompEl *cel = neighgel->Reference();
                if (cel) DebugStop();
                if (!neighgel->IsLinearMapping()) {
                    int64_t index = neighgel->Index();
                    DebugStop();
                }
                connected_elements.insert(neighgel);
                TPZGeoElSide gelside(neighgel, 4);
                TPZGeoElSide neighbour = gelside.HasNeighbour(sbfem_skeleton);
                if (!neighbour) DebugStop();
                if (!neighbour.Element()->Reference()) DebugStop();
            } else {
                DebugStop();
            }
        } else if (neighgel->Dimension() == 1 && neighgel->MaterialId() == sbfem_highperm_hdiv) {
            // create a zero dimensional element
            if (!neighgel->IsLinearMapping()) {
                int64_t index = neighgel->Index();
                DebugStop();
                neighgel = gmesh->Element(index);
            }
            connected_elements.insert(neighgel);
            TPZGeoElSide gelside(neighgel, 0);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                TPZGeoEl *neighgel = neighbour.Element();
                if (neighgel->MaterialId() == sbfem_skeleton && neighgel->Dimension() == 0) {
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
            if (neighbour == gelside) DebugStop();
            if (neighbour.Element()->Dimension() != 0) DebugStop();
            if (!neighbour.Element()->Reference()) DebugStop();
        }
    }
    TPZSBFemElementGroup *elgr = new TPZSBFemElementGroup(*cmesh);
    sbfem_groupHdiv = elgr;
    TPZStack<TPZSBFemVolume *> sbvols;
    for (auto it : connected_elements) {
        TPZGeoEl *gel = it;
        TPZSBFemVolume *sbvol = new TPZSBFemVolume(*cmesh, gel);
        sbvols.Push(sbvol);
        if (gel->Dimension() == 2) {
            TPZGeoElSide gelside(gel, 4);
            TPZGeoElSide neigh = gelside.HasNeighbour(sbfem_skeleton);
            TPZGeoEl *skel = neigh.Element();
            if (skel->MaterialId() != sbfem_skeleton) DebugStop();
            TPZCompEl *cskel = skel->Reference();
            if (!cskel) DebugStop();
            TPZMultiphysicsElement *mcel = dynamic_cast<TPZMultiphysicsElement *>(cskel);
            if (!mcel) DebugStop();
            TPZCompEl *cskelat = mcel->ElementVec()[0].Element();
            TPZInterpolatedElement *icskel = dynamic_cast<TPZInterpolatedElement *>(cskelat);
            if (!icskel) DebugStop();
            int sideorder = icskel->Connect(2).Order();
            if (sideorder != SBFemOrder) DebugStop();
            sbvol->SetSkeleton(cskel->Index());
        } else if (gel->Dimension() == 1) {
            TPZGeoElSide gelside(gel, 0);
            TPZGeoElSide neigh = gelside.Neighbour();
            while (neigh != gelside) {
                TPZGeoEl *neighgel = neigh.Element();
                if (neighgel->MaterialId() == sbfem_skeleton && neighgel->Dimension() == 0) {
                    break;
                }
                neigh = neigh.Neighbour();
            }
            if (neigh.Element()->MaterialId() != sbfem_skeleton) DebugStop();
            // create an H1 element
            TPZCompEl *cskel = neigh.Element()->Reference();
            if (!cskel) DebugStop();
            sbvol->SetSkeleton(cskel->Index());
        }
        elgr->AddElement(sbvol);
    }
    for (auto sbvol : sbvols) {
        sbvol->SetElementGroupIndex(elgr->Index());
    }
    TPZElementMatrixT<STATE> ek, ef;
    elgr->CalcStiff(ek, ef);
    if (0) {
        ek.fMat.Print("ek = ", std::cout, EMathematicaInput);
    }
    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();
}

/// Add the SBFemVolume elements to the computational mesh
void AddSBFemVolumeElements() {
    if (meshstyle == ESBFem) {

        TPZCompMesh *cmesh = 0;
        if (sbfem_groupH1) cmesh = sbfem_groupH1->Mesh();
        TPZMultiphysicsCompMesh *cmesh_m = 0;
        if (sbfem_groupHdiv) cmesh_m = dynamic_cast<TPZMultiphysicsCompMesh *>(sbfem_groupHdiv->Mesh());
        if (sbfem_groupH1 && !cmesh) DebugStop();
        if (sbfem_groupHdiv && !cmesh_m) DebugStop();
        if (sbfem_groupH1) {
            TPZStack<TPZCompEl *> compels;
            sbfem_groupH1->GetCompElList(compels);
            for (auto compel : compels) {
                int64_t index = compel->Index();
                cmesh->ElementVec()[index] = compel;
            }
        }
        if (sbfem_groupHdiv) {
            TPZStack<TPZCompEl *> compels;
            sbfem_groupHdiv->GetCompElList(compels);
            for (auto compel : compels) {
                int64_t index = compel->Index();
                cmesh_m->ElementVec()[index] = compel;
            }
        }
    }
}

/// Hide the SBFemVolume elements from the computational mesh
void HideSBFemVolumeElements() {
    if (meshstyle == ESBFem) {

        TPZCompMesh *cmesh;
        if (sbfem_groupH1) cmesh = sbfem_groupH1->Mesh();
        TPZMultiphysicsCompMesh *cmesh_m;
        if (sbfem_groupHdiv) cmesh_m = dynamic_cast<TPZMultiphysicsCompMesh *>(sbfem_groupHdiv->Mesh());
        if (sbfem_groupH1 && !cmesh) DebugStop();
        if (sbfem_groupHdiv && !cmesh_m) DebugStop();
        if (sbfem_groupH1) {
            TPZStack<TPZCompEl *> compels;
            sbfem_groupH1->GetCompElList(compels);
            for (auto compel : compels) {
                int64_t index = compel->Index();
                cmesh->ElementVec()[index] = 0;
            }
        }
        if (sbfem_groupHdiv) {
            TPZStack<TPZCompEl *> compels;
            sbfem_groupHdiv->GetCompElList(compels);
            for (auto compel : compels) {
                int64_t index = compel->Index();
                cmesh_m->ElementVec()[index] = 0;
            }
        }
    }
}

/// @brief adjust the elements neighbouring the trailing edge to accomodate SBFem simulation
void AdjustToSBFemGeometry(TPZGeoMesh *gmesh) {
    if(meshstyle != ESBFem) DebugStop();
    TPZGeoEl *trailingedge_element = gmesh->Element(trailingedge_element_index);
    int64_t singular = trailingedge_element->NodeIndex(0);
    TPZGeoElSide trailingedge(trailingedge_element);
    std::set<TPZGeoEl *> connected_elements;
    for (TPZGeoElSide neigh = trailingedge.Neighbour(); neigh != trailingedge; neigh++) {
        TPZGeoEl *neighgel = neigh.Element();
        if (neighgel->HasSubElement()) continue;
        // the element with matid < 0 is the element subsituted by two collapsed elements
        if (neighgel->MaterialId() < 0) continue;
        if (connected_elements.find(neighgel) != connected_elements.end()) continue;
        // only insert collapsed elements
        if (neighgel->Dimension() == 2) {
            if (neighgel->NodeIndex(2) == neighgel->NodeIndex(3)) {
                connected_elements.insert(neighgel);
            }
        } else if (neighgel->Dimension() == 1) {
            connected_elements.insert(neighgel);
        }
    }
    // change the material id of the volumetric elements to sbfem_domain
    // create the skeleton elements
    for (auto it : connected_elements) {
        TPZGeoEl *gel = it;
        if (gel->Dimension() == 2) {
            gel->SetMaterialId(sbfem_domain);
            sbfem_elements.insert(gel->Index());
            TPZGeoElBC gbc(gel, 4, sbfem_skeleton);
            TPZGeoEl *gbcgel = gbc.CreatedElement();
        }
    }
    // create sbfem_highperm elements
    for (auto it : connected_elements) {
        TPZGeoEl *gel = it;
        if(gel->Dimension() != 1) continue;
        int matid = gel->MaterialId();
        TPZGeoElSide gelside(gel);
        TPZManVector<int64_t, 4> nodeindices(2);
        if (gelside.SideNodeIndex(0) == singular) {
            nodeindices[0] = gelside.SideNodeIndex(1);
            nodeindices[1] = gelside.SideNodeIndex(0);
        } else if (gelside.SideNodeIndex(1) == singular) {
            nodeindices[0] = gelside.SideNodeIndex(0);
            nodeindices[1] = gelside.SideNodeIndex(1);
        } else {
            DebugStop();
        }
        if(matid == dirichletmat) {
            // creating a linear element will lead to an inconsistent geometry
            int64_t index;
            gel->SetMaterialId(-1);
            gmesh->CreateGeoElement(EOned, nodeindices, sbfem_highperm_h1, index);
            // create the skeleton element
            gmesh->CreateGeoElement(EPoint, nodeindices, sbfem_skeleton, index);
        } else if(matid == neumannmat) {
            int64_t index;
            gel->SetMaterialId(-1);
            gmesh->CreateGeoElement(EOned, nodeindices, sbfem_highperm_hdiv, index);
            // create the skeleton element
            gmesh->CreateGeoElement(EPoint, nodeindices, sbfem_skeleton, index);
        } else {
            DebugStop();
        }
    }

    // insert the connectivity of the sbfem_highperm elements
    gmesh->BuildConnectivity();
    PrintTrailingEdgeElements(gmesh);
#ifdef PZDEBUG
    {
        const int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->ElementVec()[el];
            if (!gel) continue;
            gel->VerifyNodeCoordinates();
        } // for el
    }
#endif
}
