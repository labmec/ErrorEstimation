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

class ConfigCasesMaze {
    
private:
    std::string ImageName = "Salida.png";
    double ImperviousPermeability = 1.0;
    double PermeablePermeability = 100000;
    int fluxOrder=1;
    int PressureOrder=1;
    double CCpressureIn= 1000;
    double CCpressureOut = 10;
    bool MhmOpenChannel = false;
    std::string VTKName = "Salida.vtk";
    TLaplaceExample1 * exact = nullptr;
    int nSubDomains = 1;
    int SkeletonDivision = 1;

public:
    void SetImageName(std::string name){
        ImageName=std::move(name);
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

    void SetNumberOfSubdomains(int nsub){
        nSubDomains = nsub;
    }

    int GetNumberOfSubdomains(){
        return nSubDomains;
    }

    void SetSkeletonDivision(int nsub){
        SkeletonDivision = nsub;
    }

    int GetSkeletonDivision(){
        return SkeletonDivision;
    }

};


#endif

