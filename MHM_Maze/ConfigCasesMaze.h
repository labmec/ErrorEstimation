/**
 * @file
 * @brief Contains the TPZMatLaplacian class.
 */

#ifndef ConfigCasesMaze_H
#define ConfigCasesMaze_H

#include <iostream>
#include "pzdiscgal.h"
#include "pzfmatrix.h"
#include "TPZMatLaplacian.h"
#include <string>

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
    
public:
    void SetImageName( std::string name){
        ImageName=name;
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
    
    
};


#endif

