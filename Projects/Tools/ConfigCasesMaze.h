/**
 * @file
 * @brief Contains the TPZMatLaplacian class.
 */

#ifndef ConfigCasesMaze_H
#define ConfigCasesMaze_H

#include "DarcyFlow/TPZDarcyFlow.h"
//#include "pzdiscgal.h"
#include "pzfmatrix.h"
#include <Pre/TPZAnalyticSolution.h>
#include <iostream>
#include <string>
#include <utility>

class TPZMultiphysicsCompMesh;

class ConfigCasesMaze {
    
private:
    std::string ImageName;
    double ImperviousPermeability = 1.0;
    double PermeablePermeability = 100000;
    int fluxOrder=1;
    int PressureOrder=1;
    double CCpressureIn= 1000;
    double CCpressureOut = 10;
    bool MhmOpenChannel = false;
    std::string VTKName = "Salida.vtk";
    TLaplaceExample1 * exact = nullptr;

public:
    void SetImageName(const std::string &name){
        ImageName = name;
    }

    void SetExactSolution(TLaplaceExample1 * exact_sol){
        exact = exact_sol;
    }

    void SetImperviousMatPermeability(double perm){
        ImperviousPermeability=perm;
    }
    void SetPermeableMatPermeability(double perm){
        PermeablePermeability=perm;
    }
    void SetFluxOrder(int order){
        fluxOrder=order;
    }
    void SetPressureOrder(int order){
        PressureOrder=order;
    }
    void SetCCPressureIn(double pressureval){
        CCpressureIn=pressureval;
    }
    void SetCCPressureOut(double pressureval){
        CCpressureOut=pressureval;
    }
    void SetMHMOpenChannel(bool value){
        MhmOpenChannel = value;
    }
    void SetVTKName(std::string name){
        VTKName = name;
    }
    
    //Get
    
    std::string GetImageName(){
        return ImageName;
    }
    
    double GetImperviousMatPermeability(){
        return ImperviousPermeability;
    }
    double GetPermeableMatPermeability(){
        return PermeablePermeability;
    }
    int GetFluxOrder(){
        return fluxOrder;
    }
    int GetPressureOrder(){
        return PressureOrder;
    }
    double GetCCPressureIn(){
       return  CCpressureIn;
    }
    double GetCCPressureOut(){
        return CCpressureOut;
    }
    bool GetMHMOpenChannel(){
       return MhmOpenChannel;
    }
    std::string GetVTKName(){
       return VTKName;
    }

    TLaplaceExample1* GetExactSolution(){
        return exact;
    }

    TPZGeoMesh *GenerateGeoMesh(const std::string &name, int nx, int ny);

    // this will create a geometric mesh of size L by H
    TPZGeoMesh *GeoMeshFromPng(const std::string &name, int &pix_x, int &pix_y);
    
    TPZCompMesh *CMeshPressure(TPZGeoMesh * gmesh, int pOrder);
    
    TPZCompMesh *CMeshFlux(TPZGeoMesh * gmesh,int pOrder);
    
    /// create a constant value mesh. The singleconnect parameter indicates the mesh is associated
    /// with a subdomain and should have a single connect
    TPZCompMesh *CMeshAverage(TPZGeoMesh *gmesh, bool singleconnect);

    TPZMultiphysicsCompMesh *CMeshMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> &meshvec);

    TPZMultiphysicsCompMesh *BuildMesh(int nx, int ny);

    /// set the lagrange multipliers so that the global stiffness can be inverted without pivoting
    void ConfigureLagrangeLevels(TPZMultiphysicsCompMesh *cmesh);
};


#endif

