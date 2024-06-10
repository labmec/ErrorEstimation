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
#include "TPZLinearAnalysis.h"
#include "TPZBndCondT.h"
#include <TPZNullMaterial.h>
#include <TPZMultiphysicsCompMesh.h>
#include "TPZSSpStructMatrix.h"
#include "TPZSYSMPMatrix.h"
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

#include <iostream>

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
TPZCompMesh *CreateH1CompMesh(TPZGeoMesh *gmesh, int64_t &newcon);

/// @brief Create a computational mesh with L2 elements
TPZCompMesh *CreateL2CompMesh(TPZGeoMesh *gmesh);

/// @brief Create the computational mesh with HDiv elements
TPZCompMesh *CreateHDivCompMesh(TPZGeoMesh *gmesh, int64_t &newcon);

/// @brief Create the computational "multiphysics" mesh with only HDiv elements
TPZMultiphysicsCompMesh *CreateMultiphysicsMesh(TPZCompMesh *cmeshHDiv, TPZCompMesh *cmeshL2, TPZGeoMesh *gmesh);

/// @brief Simulate the NACA profile using H1 approximation
TPZCompMesh *SimulateNacaProfileH1(TPZGeoMesh *gmesh);

/// @brief Simulate the NACA profile using H(div) approximation
TPZMultiphysicsCompMesh *SimulateNacaProfileHDiv(TPZGeoMesh *gmesh);

/// @brief Compute the error as an average error per element and save as an elemental solution 
void ComputeErrorEstimator(TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *cmesh_m,TPZVec<REAL> &ErrorEstimator);

/// @brief Compute the error as an average error per element and save as an elemental solution 
void Hrefinement(TPZMultiphysicsCompMesh *cmesh_m,TPZVec<REAL> &ErrorEstimator,TPZVec<REAL> &RefinementIndicator);

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


int volmat = 1;
int cutmat = 2;
int profilemat = 3;
int boundmat = 4;
int pointmat = 5;
int trailingedgemat = 6;
int blendmat = 7;

auto f_profile = [](const TPZVec<REAL> &loc, TPZVec<STATE> &rhsVal, TPZFMatrix<STATE> &matVal)
{
    matVal(0,0) = 0.0;
    matVal(1,0) = 0.0;
    rhsVal[0] = -loc[1]*matVal(0,0)+loc[0]*matVal(1,0);

};

auto f_infinity = [](const TPZVec<REAL> &loc, TPZVec<STATE> &rhsVal, TPZFMatrix<STATE> &matVal)
{
    matVal(0,0) = 1.0;
    matVal(1,0) = 0.5;
    rhsVal[0] = -loc[1]*matVal(0,0)+loc[0]*matVal(1,0);

};

int main() {

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    TPZManVector<REAL,3> x0(3,0.);
    REAL length = 10.;
    REAL angle = 0.0;
    TPZNacaProfile naca(length, 12, angle, x0);

    TPZGeoMesh *gmesh = ReadGmsh("naca.msh",naca);

    {
        std::ofstream out("gmesh.txt");
        gmesh->Print(out);

        std::ofstream out2("gmesh.vtk"); 
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
    }

    TPZCheckGeom check(gmesh);
    check.UniformRefine(4);

    {
        std::ofstream out4("gmeshfine.txt");
        gmesh->Print(out4);

        std::ofstream out5("gmeshfine.vtk"); 
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out5);
    }

    int nrefinements = 5;
    for(int64_t i=0; i<nrefinements; i++)
    {
        auto cmesh = SimulateNacaProfileH1(gmesh);
        auto cmesh_m = SimulateNacaProfileHDiv(gmesh);
        TPZVec<REAL> Error;
        ComputeErrorEstimator(cmesh,cmesh_m,Error);
        TPZVec<REAL> RefinementIndicator;
        Hrefinement(cmesh_m,Error,RefinementIndicator);
        std::ofstream out6("ErrorEstimator.vtk"); 
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m,out6,Error,"Error");  
        {
            const std::string plotfile = "postprocess_H1";
            constexpr int vtkRes{0};
            TPZVec<std::string> fields = 
            {
            "Pressure",
            "Flux"
            };
            auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
            vtk.SetNThreads(0);
            vtk.Do();
        }

        {
            std::ofstream out("cmesh.txt");
            cmesh->Print(out);
        }

        {
            int64_t nnod = cmesh->NConnects();
            for(int64_t i=0; i<nnod; i++)
            {
                TPZConnect &c = cmesh->ConnectVec()[i];
                if(c.HasDependency())
                {
                    std::cout << "Connect " << i << " has dependency\n";
                    c.RemoveDepend();
                }
            }    
        }
        delete cmesh;
        TPZVec<TPZCompMesh*> meshvec = cmesh_m->MeshVector();
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
            TPZGeoElRefPattern< TPZNacaProfile> *nacael = new TPZGeoElRefPattern< TPZNacaProfile>(nodes,blendmat,*gmesh);
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

/// @brief Create a computational mesh with H1 elements
TPZCompMesh *CreateH1CompMesh(TPZGeoMesh *gmesh, int64_t &newcon)
{
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

    std::set<int> matidsh1 = {volmat,profilemat,boundmat};
    cmesh->AutoBuild(matidsh1);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    std::set<int> matidpoint = {pointmat};
    cmesh->AutoBuild(matidpoint);
    // // insert the cut boundary condition
    newcon = cmesh->AllocateNewConnect(1,1,1);
    std::cout << "newcon " << newcon << std::endl;
    int64_t nel = gmesh->NElements();
    // // loop over the geometric cut elements
    for(int64_t el = 0; el<nel; el++)
    {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->MaterialId() != cutmat) continue;
        TPZGeoElSide gelside(gel);
        TPZGeoElSide neighbour = gelside.Neighbour();
        // look at the neighbours of the cut elements
        // compute the parameter transformation. If the det is positive
        while(neighbour != gelside)
        {
            if(neighbour.Element()->MaterialId() == volmat)
            {
                TPZTransform<REAL> tr(1);
                gelside.SideTransform3(neighbour, tr);
                if(tr.Mult()(0,0) > 0.) break;
            }
            neighbour = neighbour.Neighbour();
        }
        if(neighbour == gelside) DebugStop();
        TPZCompEl *cel = neighbour.Element()->Reference();
        if(!cel) continue;
        int nsidenodes = gel->NSideNodes(gelside.Side());
        // loop over the connects of the computational element
        // if the connect has no dependency, create a new connect
        for(int is=0; is<nsidenodes; is++)
        {
            int nodelocindex = neighbour.Element()->SideNodeLocIndex(neighbour.Side(),is);
            int64_t cindex = cel->ConnectIndex(nodelocindex);
            TPZConnect &c = cel->Connect(nodelocindex);
            if(c.HasDependency()) continue;
            // insert the dependency
            int64_t cindex2 = cmesh->AllocateNewConnect(c.NShape(),c.NState(),c.Order());
            TPZConnect &c2 = cmesh->ConnectVec()[cindex2];
            TPZFNMatrix<1,STATE> val(1,1,1.);
            c2.AddDependency(cindex2,cindex,val,0,0,1,1);
            c2.AddDependency(cindex2,newcon,val,0,0,1,1);
            cel->SetConnectIndex(nodelocindex,cindex2);
        }
        auto sidecut = neighbour.Side();
        TPZGeoEl *neighel = neighbour.Element();
        TPZCompEl*neighcel = neighel->Reference();
        int ncorner = neighel->NCornerNodes();
        for(int side = ncorner; side<neighel->NSides()-1; side++)
        {
            if(side == sidecut) continue;
            TPZGeoElSide neighneigh = neighel->Neighbour(side);
            if(neighneigh.NSideNodes()!=2) DebugStop();
            {
                int indexneighneigh = neighneigh.SideNodeLocIndex(0);
                int indexneigh = neighel->SideNodeLocIndex(side,1);
                int64_t cindex = neighcel->ConnectIndex(indexneigh); 
                TPZCompEl *neighneighcel = neighneigh.Element()->Reference(); 
                neighneighcel->SetConnectIndex(indexneighneigh,cindex);  
            }
            {
                int indexneighneigh = neighneigh.SideNodeLocIndex(1);
                int indexneigh = neighel->SideNodeLocIndex(side,0);
                int64_t cindex = neighcel->ConnectIndex(indexneigh); 
                TPZCompEl *neighneighcel = neighneigh.Element()->Reference(); 
                neighneighcel->SetConnectIndex(indexneighneigh,cindex);
            }
        }
    }
    cmesh->ExpandSolution();
    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
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
TPZCompMesh *CreateHDivCompMesh(TPZGeoMesh *gmesh, int64_t &newcon)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    int dim = 2;
    auto material = new TPZNullMaterial(volmat);
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

    std::set<int> matidshdiv = {volmat,profilemat,boundmat};
    cmesh->AutoBuild(matidshdiv);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();

    // Fixed point condition   
    // auto bnd3 = new TPZNullMaterial(pointmat);
    // bnd3->SetDimension(0);
    // cmesh->InsertMaterialObject(bnd3);
    // std::set<int> Pointmat = {pointmat};
    // cmesh->AutoBuild(Pointmat);
// insert the circulation condition
    newcon = cmesh->AllocateNewConnect(1,1,1);
    std::cout << "newcon " << newcon << std::endl;
    int64_t nel = gmesh->NElements();
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
            TPZConnect &c = cel->Connect(is);
            if(c.HasDependency()) continue;
            // insert the dependency
            TPZFNMatrix<1,STATE> val(1,1,1.);
            c.AddDependency(cindex,newcon,val,0,0,1,1);
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
    cmesh->ExpandSolution();
    std::ofstream out("cmeshHDiv.vtk"); 
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out);
    return cmesh;
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
     return cmesh_m;
 }

/// @brief Simulate the NACA profile using H1 approximation
TPZCompMesh *SimulateNacaProfileH1(TPZGeoMesh *gmesh)
{
    int64_t newcon;
    auto cmeshH1 = CreateH1CompMesh(gmesh, newcon);
    if (1)
    {
        std::ofstream out("cmeshH1.txt");
        cmeshH1->Print(out);
    }

    TPZLinearAnalysis an(cmeshH1,RenumType::EMetis);
    TPZSSpStructMatrix<STATE> strmat(cmeshH1);
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
    for (int64_t eq = 0; eq < numeq; eq++) 
    {
        phi_0(eq) =  fsol(eq,0);
        phi_1(eq) =  fsol(eq,1);
    }

//  Compute the "Beta" constant
    REAL Beta;
    ComputeBetaH1(gmesh, cmeshH1, phi_0, phi_1, Beta);

    phi = phi_0+Beta*phi_1;
    cmeshH1->LoadSolution(phi);
    an.LoadSolution(phi);

    REAL circulation;
    EvaluateCirculationH1(gmesh, cmeshH1, profilemat, circulation);

    std::cout << "--------- PostProcess ---------" << std::endl;
    //printa na tela "--------- PostProcess ---------", indicando que a simulação está em processamento.
    PrintResults(an,cmeshH1);
    //chama a função PrintResults para realizar o pós-processamento dos resultados. Essa função provavelmente gera saídas com os resultados da simulação.
    return cmeshH1;
}

/// @brief Simulate the NACA profile using H(div) approximation
 TPZMultiphysicsCompMesh *SimulateNacaProfileHDiv(TPZGeoMesh *gmesh)
 {
    int64_t newcon;
    gmesh->ResetReference();
     auto cmeshHDiv = CreateHDivCompMesh(gmesh, newcon);
     auto cmeshL2 = CreateL2CompMesh(gmesh);

     TPZMultiphysicsCompMesh* cmesh_m = CreateMultiphysicsMesh(cmeshHDiv,cmeshL2, gmesh);
     // Define o pointer chamado cmesh_m relacionado à classe TPZMultiphysicsCompMesh, associando-o a função CreateMultiphysicsMesh, cujos parâmetros (já antes declarados) são: cmeshHdiv, gmesh.
    
     {
         std::ofstream out("cmesh_m.txt");
         cmesh_m->Print(out);
     }

    TPZLinearAnalysis an(cmesh_m,RenumType::EMetis);
    TPZSSpStructMatrix<STATE> strmat(cmesh_m);
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

    REAL circulation;
    EvaluateCirculationHDiv(gmesh, cmesh_m, profilemat, circulation);

     std::cout << "--------- PostProcess ---------" << std::endl;

     PrintResults(an,cmesh_m);

     return cmesh_m;
 }

/// @brief Compute the error as an average error per element and save as an elemental solution in "Error"  
void ComputeErrorEstimator(TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator)
{
    cmesh->LoadReferences();
    int64_t nel_m = cmesh_m->NElements();
    ErrorEstimator.Resize(nel_m,0.);
    ErrorEstimator = 0.;
    // TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh->ElementVec();
    TPZAdmChunkVector<TPZCompEl *> &elementvec_m = cmesh_m->ElementVec();

    for (int64_t iel = 0; iel < nel_m; iel++) 
    {
        TPZCompEl *cell_m = elementvec_m[iel];
        TPZGeoEl *gel_m = cell_m->Reference();
        if(gel_m->Dimension() != 2) continue;

        TPZCompEl *cell = gel_m->Reference();

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
            TPZManVector<STATE,3> flux;
            msp ->Solution(intpoint,7,flux);

            TPZMultiphysicsElement *msp_m  = dynamic_cast <TPZMultiphysicsElement *>(cell_m);
            TPZManVector<STATE,3> sol;
            msp_m ->Solution(intpoint,1,sol);

            ErrorEstimator[iel] += ((flux[0]-sol[0])*(flux[0]-sol[0])+(flux[1]-sol[1])*(flux[1]-sol[1]))*weight;
        }//loop over integratin points
    }//loop over cemsh_m elements
}

/// @brief Performe h-refinements at the computational meshes "cmesh" and "cmesh_m" based on the ErrorEstimator  
void Hrefinement(TPZMultiphysicsCompMesh *cmesh_m, TPZVec<REAL> &ErrorEstimator, TPZVec<REAL> &RefinementIndicator)
{
    int64_t nel_m = cmesh_m->NElements();
    auto gmesh = cmesh_m->Reference();
    int64_t nel = gmesh->NElements();
    RefinementIndicator.Resize(nel,0.);
    RefinementIndicator = 0.;
    auto tol = std::max_element(ErrorEstimator.begin(), ErrorEstimator.end());
    *tol = *tol/5.0;

    for (int64_t iel = 0; iel < nel_m; iel++) 
    {
        TPZCompEl *cel = cmesh_m->Element(iel);
        TPZGeoEl *gel = cel->Reference();
        int64_t index = gel->Index();
        if(ErrorEstimator[iel]>*tol)
        {
            RefinementIndicator[iel] = 1.0;
            TPZManVector<TPZGeoEl*> subel; 
            gel->Divide(subel);
        } 
    }//loop over cemsh_m elements
    nel = gmesh->NElements();
    for (int64_t iel = 0; iel < nel; iel++) 
    {
        TPZGeoEl *gel = gmesh->Element(iel);
        if(gel->HasSubElement()) continue;
        if(gel->Dimension() != 1) continue;
        //Percorrer os Neighbors de gel e encontar os volumetricos que tem subelementos
        TPZGeoElSide gelside(gel); 
        auto Neighbor = gelside.Neighbour();
        for(; Neighbor != gelside; Neighbor++)
        {
            if(!Neighbor.Element()) DebugStop();
            if(!Neighbor.HasSubElement()) continue;
            TPZManVector<TPZGeoEl*> subel; 
            gel->Divide(subel);
            break;
        }
    }//loop over cemsh_m elements

    std::ofstream out1("RefinementIndicator.vtk"); 
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m,out1,RefinementIndicator,"RefinementIndicator");

    nel = gmesh->NElements();
    RefinementIndicator.Resize(nel,0.);

    {
        std::ofstream out2("gmesh_Hrefinement.vtk"); 
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
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
    constexpr int vtkRes{0};
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
    TPZSYsmpMatrix<STATE>* spmat = dynamic_cast<TPZSYsmpMatrix<STATE>*> (mat.operator->());
    if(!spmat) DebugStop();
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
            std::cout << "rhs before " << rhs(row,0) << std::endl;
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

    // First coment;
    for (int64_t el = 0; el < nel; el++) 
    {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->HasSubElement()) continue;

        if (gel->MaterialId() == trailingedgemat) 
        {
            TrailingSide = TPZGeoElSide(gel);
            auto Neighbor = gel->Neighbour(0);
            for(; Neighbor != gel; Neighbor++)
            {
                if(!Neighbor.Element()) DebugStop();
                if(Neighbor.Element()->HasSubElement()) continue;
                if (Neighbor.Element()->MaterialId() == profilemat)  
                {
                    NeighborNacaElements.Push (Neighbor.Element());
                }  
            }
        }
    }
    if(NeighborNacaElements.size() != 2) DebugStop();
    if(! TrailingSide) DebugStop();


    ////Fazer um loop para cada um dos filhos dos dois elementos do vetor NeighborNacaElements e 
    //verificar qual dos filhos tem o Trailing Edge como um neighbor

    // Second coment;
    TPZManVector<int,2> TE_elementsindex;
    TE_elementsindex.resize(2);
    int cont = -1;
    for(auto el : NeighborNacaElements)
    {
        TPZGeoElSide side = TPZGeoElSide(el, 2);
        TPZGeoElSide Neighbor = el->Neighbour(2);
        for(; Neighbor != side; Neighbor++) 
        {
            if (Neighbor.Element()->MaterialId() != volmat) continue;
            TE_elementsindex[++cont] = Neighbor.Element()->Index();
        }
    }

    // Third coment;
    TPZGeoElSide Neighbor = TrailingSide.Neighbour();
    for(; Neighbor !=TrailingSide; Neighbor++) 
    {
        if (Neighbor.Element()->MaterialId() != volmat) continue;  
        if(Neighbor.Element()->Index() == TE_elementsindex[0])
        {
            NeighborNacaElementsSide.Push (Neighbor.Side());
        }
        if(Neighbor.Element()->Index() == TE_elementsindex[1])
        {
            NeighborNacaElementsSide.Push (Neighbor.Side());
        }
    }
    
    TPZGeoEl *gel0 = gmesh->Element(TE_elementsindex[0]);
    Geosideminus = TPZGeoElSide (gel0,NeighborNacaElementsSide[0]);
    TPZGeoEl *gel1 = gmesh->Element(TE_elementsindex[1]);
    Geosideplus = TPZGeoElSide (gel1,NeighborNacaElementsSide[1]);
}

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
    qsi[0] = tr.Sum()(0,0);
    qsi[1] = tr.Sum()(1,0);

    Compelminus->Solution(qsi,2,gradminus); 

    //Compute the gradient of the potential solution for the Compelplus
    tr = Geosideplus.Element()->SideToSideTransform(Geosideplus.Side(),Geosideplus.Element()->NSides()-1);
    qsi[0] = tr.Sum()(0,0);
    qsi[1] = tr.Sum()(1,0);

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
    qsi[0] = tr.Sum()(0,0);
    qsi[1] = tr.Sum()(1,0);

    Compelminus->Solution(qsi,1,u_minus); 

    //Compute the Flux solution for the Compelplus
    tr = Geosideplus.Element()->SideToSideTransform(Geosideplus.Side(),Geosideplus.Element()->NSides()-1);
    qsi[0] = tr.Sum()(0,0);
    qsi[1] = tr.Sum()(1,0);

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
        auto intrule = gelside.CreateIntegrationRule(5);
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
            if (neighbor.Element()->MaterialId() == volmat) 
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
            TPZManVector<STATE,3> flux;
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
        TPZGeoEl *gel = el->Reference();
        if (gel->MaterialId() != matid) continue;

        //NUMERICAL INTEGRATION DATA;
        TPZGeoElSide gelside(gel);
        auto intrule = gelside.CreateIntegrationRule(5);
        int intrulepoints = intrule->NPoints();
        TPZManVector<REAL,4> intpointtemp(1,0.);
        TPZManVector<REAL,4> intpointvol(2,0.);
        REAL weight = 0.;
        TPZFMatrix<REAL> jac, axe, jacInv;
        REAL detJac;

        TPZGeoElSide neighbor = gelside.Neighbour();
        for(; neighbor != gel; neighbor++)
        {
            if (neighbor.Element()->MaterialId() == volmat) 
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
            TPZMultiphysicsElement *msp  = dynamic_cast <TPZMultiphysicsElement *>(ref);
            TPZManVector<STATE,3> sol;
            msp ->Solution(intpointvol,1,sol);
            REAL ut = sol[0]*axe(0,0) + sol[1]*axe(0,1);
            circulation += ut*weight;
        }//loop over integration points
    }//loop over elements
}

/// @brief Compute the Beta constant
/// @return The value of Beta: REAL Beta.
void ComputeBetaH1(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, TPZFMatrix<STATE> &phi_0, TPZFMatrix<STATE> &phi_1, REAL &Beta)

{
    cmesh->LoadSolution(phi_0);
    TPZManVector<REAL,3> grad_phi0_plus;
    TPZManVector<REAL,3> grad_phi0_minus;
    EvaluateSolutionGradientsH1(gmesh, grad_phi0_plus, grad_phi0_minus);

    cmesh->LoadSolution(phi_1);
    TPZManVector<REAL,3> grad_phi1_plus;
    TPZManVector<REAL,3> grad_phi1_minus;
    EvaluateSolutionGradientsH1(gmesh, grad_phi1_plus, grad_phi1_minus);

    REAL A = grad_phi1_plus[0]*grad_phi1_plus[0]+grad_phi1_plus[1]*grad_phi1_plus[1]-grad_phi1_minus[0]*grad_phi1_minus[0]-grad_phi1_minus[1]*grad_phi1_minus[1];
    REAL B = 2.0*(grad_phi0_plus[0]*grad_phi1_plus[0]+grad_phi0_plus[1]*grad_phi1_plus[1]-grad_phi0_minus[0]*grad_phi1_minus[0]-grad_phi0_minus[1]*grad_phi1_minus[1]);
    REAL C = grad_phi0_plus[0]*grad_phi0_plus[0]+grad_phi0_plus[1]*grad_phi0_plus[1]-grad_phi0_minus[0]*grad_phi0_minus[0]-grad_phi0_minus[1]*grad_phi0_minus[1];
    REAL Delta = B*B-4.0*A*C;

    if ((-B+sqrt(Delta))/(2*A)>(-B-sqrt(Delta))/(2*A))
    {
        Beta = (-B+sqrt(Delta))/(2*A);
    }
    else 
    {
        Beta = (-B-sqrt(Delta))/(2*A);
    }
}

/// @brief Compute the Kappa constant
/// @return The value of Kappa: REAL Kappa.
void ComputeBetaHDiv(TPZGeoMesh *gmesh, TPZMultiphysicsCompMesh *cmesh_m, TPZFMatrix<STATE> &u_0, TPZFMatrix<STATE> &u_1, REAL &Beta)

{
    cmesh_m->LoadSolution(u_0);
    cmesh_m->TransferMultiphysicsSolution();
    TPZManVector<REAL,3> u0_plus;
    TPZManVector<REAL,3> u0_minus;
    EvaluateSolutionHDiv(gmesh, cmesh_m, u0_plus, u0_minus);

    cmesh_m->LoadSolution(u_1);
    cmesh_m->TransferMultiphysicsSolution();
    TPZManVector<REAL,3> u1_plus;
    TPZManVector<REAL,3> u1_minus;
    EvaluateSolutionHDiv(gmesh, cmesh_m, u1_plus, u1_minus);

    REAL A = u1_plus[0]*u1_plus[0]+u1_plus[1]*u1_plus[1]-u1_minus[0]*u1_minus[0]-u1_minus[1]*u1_minus[1];
    REAL B = 2.0*(u0_plus[0]*u1_plus[0]+u0_plus[1]*u1_plus[1]-u0_minus[0]*u1_minus[0]-u0_minus[1]*u1_minus[1]);
    REAL C = u0_plus[0]*u0_plus[0]+u0_plus[1]*u0_plus[1]-u0_minus[0]*u0_minus[0]-u0_minus[1]*u0_minus[1];
    REAL Delta = B*B-4.0*A*C;

    if ((-B+sqrt(Delta))/(2*A)<(-B-sqrt(Delta))/(2*A))
    {
        Beta = (-B+sqrt(Delta))/(2*A);
    }
    else 
    {
        Beta = (-B-sqrt(Delta))/(2*A);
    }
}



