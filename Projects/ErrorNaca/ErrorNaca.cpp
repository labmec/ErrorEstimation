/**
 * @file This file implements an error estimator in space H1.
 */

#include "pzlog.h"
#include "TPZLinearAnalysis.h"
#include "pzstepsolver.h"
#include "TPZBndCondT.h"
#include "DarcyFlow/TPZDarcyFlow.h"
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

int main() {
    TPZManVector<REAL,3> x0(3,0.);
    REAL length = 10.;
    REAL angle = 0.;
    TPZNacaProfile naca(length, 12, angle, x0);
    naca.lastpos = 10.;
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

    TPZGmshReader gmsh;
    
    gmsh.GetDimNamePhysical()[1]["dirichlet"] = -1;
    gmsh.GetDimNamePhysical()[1]["neuman"] = -2;
    gmsh.GetDimNamePhysical()[2]["domain"] = 1;
    std::string meshfilename = "naca.msh";
    auto gmesh = gmsh.GeometricGmshMesh(meshfilename);

    for (int64_t el = 0; el < gmesh->NElements(); el++) {
        TPZGeoEl* gel = gmesh->Element(el);
        // Your code here
        if(gel->MaterialId() == 2) {
            std::cout << "gel " << el << " index " << gel->Index() << std::endl;
            // gel->Print(std::cout);
            TPZGeoElSide gelside(gel);
            TPZManVector<int64_t,2> nodes(2);
            nodes[0] = gel->NodeIndex(0);
            nodes[1] = gel->NodeIndex(1);
            naca.fNodeIndexes[0] = nodes[0];
            naca.fNodeIndexes[1] = nodes[1];
            TPZGeoElRefPattern< TPZNacaProfile> *nacael = new TPZGeoElRefPattern< TPZNacaProfile>(nodes,4,*gmesh);
            std::cout << "nacael index " << nacael->Index() << std::endl;
            nacael->Geom() = naca;
            nacael->Initialize();
            TPZManVector<REAL,3> coord(3,0.);
            TPZManVector<REAL,1> ksi(1,-1.);
            nacael->X(ksi,coord);
            std::cout << "ksi " << ksi << " coord " << coord << std::endl;
            ksi[0] = 1.;
            nacael->X(ksi,coord);
            std::cout << "ksi " << ksi << " coord " << coord << std::endl;
            TPZGeoElSide nacaside(nacael);
            nacaside.CenterX(coord);
            REAL par = 0.;
            int uplow = 0;
            int maxpt = 1000;
            std::cout << "center coord " << coord << std::endl;
            naca.NearestParameter(coord, uplow, maxpt, par);
            std::cout << "par = " << par << " coord " << coord << std::endl;

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
    std::ofstream out1("gmesh.txt");
    gmesh->Print(out1);

    TPZCheckGeom check(gmesh);
    check.UniformRefine(4);

    std::ofstream out2("naca.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
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

