/**
 * @file This file implements an error estimator in space H1.
 */

#include "pzlog.h"
#include "TPZLinearAnalysis.h"
#include "pzstepsolver.h"
#include "TPZBndCondT.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "pzgeoelbc.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZRefPatternDataBase.h"
#include "TPZVTKGeoMesh.h"
#include "TPZGmshReader.h"
#include "TPZAnalyticSolution.h"
#include "TPZFrontSym.h"
#include "pzcheckgeom.h"

#include <fstream>
#include <ctime>
#include <cstdio>
#include <cmath>

//#include "pzgengrid.h"
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
TPZCompMesh *CreateH1CompMesh(TPZGeoMesh *gmesh);

/// @brief Create de computational mesh with HDiv elements
TPZCompMesh *CreateHDivCompMesh(TPZGeoMesh *gmesh);

/// @brief Simulate the NACA profile using H1 approximation
TPZCompMesh *SimulateNacaProfileH1(TPZGeoMesh *gmesh);

/// @brief Check the NACA profile, generate coordinates and write to a file
void CheckNacaProfile(const TPZNacaProfile &naca);

int volmat = 1;
int cutmat = 2;
int profilemat = 3;
int boundmat = 4;
int pointmat = 5;

int main() {

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    TPZManVector<REAL,3> x0(3,0.);
    REAL length = 10.;
    REAL angle = 0.;
    TPZNacaProfile naca(length, 12, angle, x0);

    TPZGeoMesh *gmesh = ReadGmsh("naca.msh",naca);

    {
        std::ofstream out("gmesh.txt");
        gmesh->Print(out);
    }
    auto cmeshH1 = SimulateNacaProfileH1(gmesh);


    std::ofstream out3("nacacoarse.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out3);

    TPZCheckGeom check(gmesh);
    check.UniformRefine(1);

    std::ofstream out4("gmeshfine.txt");
    gmesh->Print(out4);


    std::ofstream out2("naca.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);

    {
        int64_t nnod = cmeshH1->NConnects();
        for(int64_t i=0; i<nnod; i++)
        {
            TPZConnect &c = cmeshH1->ConnectVec()[i];
            if(c.HasDependency())
            {
                std::cout << "Connect " << i << " has dependency\n";
                c.RemoveDepend();
            }
        }    
    }
    delete cmeshH1;
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
    gmsh.GetDimNamePhysical()[1]["Profile"] = profilemat;
    gmsh.GetDimNamePhysical()[1]["OuterBoundary"] = boundmat;
    gmsh.GetDimNamePhysical()[1]["Cut"] = cutmat;
    gmsh.GetDimNamePhysical()[2]["Domain"] = volmat;
    auto gmesh = gmsh.GeometricGmshMesh(meshfilename);

    for (int64_t el = 0; el < gmesh->NElements(); el++) {
        TPZGeoEl* gel = gmesh->Element(el);
        // Your code here
        if(gel->MaterialId() == -1) {
            // std::cout << "gel " << el << " index " << gel->Index() << std::endl;
            // gel->Print(std::cout);
            TPZGeoElSide gelside(gel);
            TPZManVector<int64_t,2> nodes(2);
            nodes[0] = gel->NodeIndex(0);
            nodes[1] = gel->NodeIndex(1);
            TPZGeoElRefPattern< TPZNacaProfile> *nacael = new TPZGeoElRefPattern< TPZNacaProfile>(nodes,-4,*gmesh);
            //std::cout << "nacael index " << nacael->Index() << std::endl;
            nacael->Geom() = naca;
            nacael->Geom().fNodeIndexes[0] = nodes[0];
            nacael->Geom().fNodeIndexes[1] = nodes[1];
            nacael->Initialize();
            TPZManVector<REAL,3> coord(3,0.);
            TPZManVector<REAL,1> ksi(1,-1.);
            nacael->X(ksi,coord);
            //std::cout << "ksi " << ksi << " coord " << coord << std::endl;
            ksi[0] = 1.;
            nacael->X(ksi,coord);
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
TPZCompMesh *CreateH1CompMesh(TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    int dim = 2;
    TPZDarcyFlow *material = new TPZDarcyFlow(volmat,dim);
    cmesh->InsertMaterialObject(material);
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);
    auto bnd = material->CreateBC(material, profilemat, 1, val1, val2);
    TPZBndCondT<STATE> *bndcond = dynamic_cast<TPZBndCondT<STATE> *>(bnd);
    cmesh->InsertMaterialObject(bnd);
//        std::function<void (const TPZVec<REAL> &loc,
//                                                   TPZVec<TVar> &rhsVal,
//                                                   TPZFMatrix<TVar> &matVal)>;

    auto f = [](const TPZVec<REAL> &loc, TPZVec<STATE> &rhsVal, TPZFMatrix<STATE> &matVal)
    {
        rhsVal[0] = 0.;
        matVal(0,0) = 1.;
        matVal(1,0) = 1.;
    };
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(1);
    auto bnd2 = material->CreateBC(material, boundmat, 1, val1, val2);
    bnd2->SetForcingFunctionBC(f,1);
    cmesh->InsertMaterialObject(bnd2);
    auto bnd3 = material->CreateBC(material, pointmat, 0, val1, val2);
    cmesh->InsertMaterialObject(bnd3);
    std::set<int> matidshdiv = {volmat,profilemat,boundmat};
    cmesh->AutoBuild(matidshdiv);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    std::set<int> matidpoint = {pointmat};
    cmesh->AutoBuild(matidpoint);
    // insert the cut boundary condition
    int64_t newcon = cmesh->AllocateNewConnect(1,1,1);
    std::cout << "newcon " << newcon << std::endl;
    int64_t nel = gmesh->NElements();
    // loop over the geometric cut elements
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
        if(!cel) DebugStop();
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
    }
    cmesh->ExpandSolution();
    cmesh->ComputeNodElCon();
    return cmesh;
}

/// @brief Create a computational mesh with HDiv elements
TPZCompMesh *CreateHDivCompMesh(TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    int dim = 2;
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    HDivFamily hdiv = HDivFamily::EHDivKernel;
    cmesh->ApproxSpace().SetHDivFamily(hdiv);
    auto material = new TPZMixedDarcyFlow(volmat,dim);
    cmesh->InsertMaterialObject(material);
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);
    // we set a Neumann condition on the profile
    // this will be modified later
    auto bnd = material->CreateBC(material, profilemat, 0, val1, val2);
    TPZBndCondT<STATE> *bndcond = dynamic_cast<TPZBndCondT<STATE> *>(bnd);
    cmesh->InsertMaterialObject(bnd);
//        std::function<void (const TPZVec<REAL> &loc,
//                                                   TPZVec<TVar> &rhsVal,
//                                                   TPZFMatrix<TVar> &matVal)>;

    auto f = [](const TPZVec<REAL> &loc, TPZVec<STATE> &rhsVal, TPZFMatrix<STATE> &matVal)
    {
        rhsVal[0] = 0.;
        matVal(0,0) = 1.;
        matVal(1,0) = 1.;
    };
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(1);
    auto bnd2 = material->CreateBC(material, boundmat, 1, val1, val2);
    bnd2->SetForcingFunctionBC(f,1);
    cmesh->InsertMaterialObject(bnd2);
    auto bnd3 = material->CreateBC(material, pointmat, 0, val1, val2);
    cmesh->InsertMaterialObject(bnd3);
    cmesh->AutoBuild();
        // insert the cut boundary condition
    int64_t newcon = cmesh->AllocateNewConnect(1,1,1);
    std::cout << "newcon " << newcon << std::endl;
    int64_t nel = gmesh->NElements();
    // loop over the geometric cut elements
    for(int64_t el = 0; el<nel; el++)
    {
        TPZGeoEl *gel = gmesh->Element(el);
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
            c.SetOrder(1,cindex);
            c.SetNState(1);
            c.SetNShape(0);
            // reset the size of the block of the connect
            cmesh->Block().Set(seq,0);
        }
    }
    cmesh->ExpandSolution();
    return cmesh;
}

/// @brief Simulate the NACA profile using H1 approximation
TPZCompMesh *SimulateNacaProfileH1(TPZGeoMesh *gmesh)
{
    auto cmeshH1 = CreateH1CompMesh(gmesh);
    {
        std::ofstream out("cmeshH1.txt");
        cmeshH1->Print(out);
    }

    {
        std::ofstream out("cmesh.txt");
        cmeshH1->Print(out);
    }

    TPZLinearAnalysis an(cmeshH1);
    TPZSkylineStructMatrix<STATE> strmat(cmeshH1);
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Run();
    return cmeshH1;
}
