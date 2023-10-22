//
//  TPZHDivErrorEstimateMaterial.cpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 11/06/19.
//

#include "TPZHDivErrorEstimateMaterial.h"
#include "pzaxestools.h"
#include "TPZAnalyticSolution.h"
#include "TPZMaterialDataT.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "Elasticity/TPZMixedElasticityND.h"



template<>
TPZHDivErrorEstimateMaterial<TPZMixedDarcyFlow>::TPZHDivErrorEstimateMaterial(int matid, int dim) : TPZMixedDarcyFlow(matid,dim)
{
    
}

template<>
TPZHDivErrorEstimateMaterial<TPZMixedElasticityND>::TPZHDivErrorEstimateMaterial(int matid, int dim)
{
    DebugStop();
}


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.errorestimation.hdiv"));
#endif

template <typename MixedMaterial>
int TPZHDivErrorEstimateMaterial<MixedMaterial>::FirstNonNullApproxSpaceIndex(const TPZVec<TPZMaterialDataT<STATE>> &datavec) const { 
    signed long int nvec = datavec.NElements();
    for (unsigned long int ivec = 0; ivec < nvec ; ivec++){
        if (datavec[ivec].fShapeType != TPZMaterialData::EEmpty){
            return ivec;
        }
    }
    return -1;
}

template <typename MixedMaterial>
void TPZHDivErrorEstimateMaterial<MixedMaterial>::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const {
    //fem solution for flux and potential
    datavec[2].SetAllRequirements(false);
    datavec[2].fNeedsSol = true;
    datavec[2].fNeedsNormal = true;

    datavec[3].SetAllRequirements(false);
    datavec[3].fNeedsSol = true;
}

template <typename MixedMaterial>
void TPZHDivErrorEstimateMaterial<MixedMaterial>::FillBoundaryConditionDataRequirements(int type,TPZVec<TPZMaterialDataT<STATE> > &datavec) const {
    datavec[2].SetAllRequirements(false);
    datavec[2].fNeedsSol = true;
    datavec[2].fNeedsNormal = true;
}

template <typename MixedMaterial>
int TPZHDivErrorEstimateMaterial<MixedMaterial>::VariableIndex(const std::string &name) const
{
    if(name == "FluxFem") return 40;
    if(name == "FluxReconstructed") return 41;
    if(name == "FluxExact") return 42;
    if(name == "PressureFEM") return 43;
    if(name == "PressureReconstructed") return 44;
    if(name == "PressureExact") return 45;
    if(name == "PressureErrorExact") return 100;
    if(name == "PressureErrorEstimate") return 101;
    if(name == "EnergyErrorExact") return 102;
    if(name == "EnergyErrorEstimate") return 103;
    if(name == "ResidualError") return 104;
    if(name == "PressureEffectivityIndex") return 105;
    if(name == "EnergyEffectivityIndex") return 106;
    if(name == "POrder") return 46;
     
    return -1;
}


template <typename MixedMaterial>
int TPZHDivErrorEstimateMaterial<MixedMaterial>::NSolutionVariables(int var) const
{
    switch (var) {
        case 40:
        case 41:
        case 42:
            return 3;
            break;
        case 43:
        case 44:
        case 45:
        case 46:
        case 47:
        case 100:
        case 101:
        case 102:
        case 103:
        case 104:
        case 105:
        case 106:
            return 1;
            break;
        default:
            DebugStop();
            break;
    }
    return 0;
}

template <typename MixedMaterial>
void TPZHDivErrorEstimateMaterial<MixedMaterial>:: ErrorsBC(TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors, TPZBndCondT<STATE> &bc){
    
    if(bc.Type()== 4){
        int H1functionposition = 0;
        H1functionposition = FirstNonNullApproxSpaceIndex(data);
        TPZVec<STATE> u_exact(1);
        TPZFMatrix<STATE> du_exact(3,1,0.);
        if(MixedMaterial::HasExactSol())
        {
            MixedMaterial::ExactSol()(data[H1functionposition].x,u_exact,du_exact);
        }
        else
        {
            return;
        }


    errors.Resize(NEvalErrors());
    errors.Fill(0.0);
    
    
    TPZFNMatrix<3,REAL> fluxreconstructed(3,1), fluxreconstructed2(3,1);
    TPZManVector<STATE,3> fluxfem(3);
        

    REAL normalsigmafem = 0.,normalsigmarec = 0.,urec=0.;;
    normalsigmafem = data[2].sol[0][0];// sigma.n
    urec = data[H1functionposition].sol[0][0];


    
    REAL u_D = 0.,g = 0.;
    REAL normflux = 0.;
        
        TPZManVector<STATE,3> fluxrec(this->Dimension());
        this->Solution(data,VariableIndex("FluxReconstructed"), fluxrec);
        
  //      std::cout<<"flux_rec "<<fluxrec[0]<<" , "<<fluxrec[1]<<"\n";
        
    
    TPZFNMatrix<9,REAL> PermTensor, InvPermTensor;
    TPZManVector<STATE> res(3);
    TPZFNMatrix<9, STATE> gradu(this->Dimension(), 1);
        u_D = u_exact[0];
        std::cout << "MUDEI O FLUXO DO PROGRAMA OLHE COM CUIDADO\n";
        
        
        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                
                normflux += data[2].normal[i]*PermTensor(i,j)*gradu(j,0);
            
            }
        }
        g = (-1)*normflux;
        
 //       std::cout<<"n_0 "<<data[2].normal[0]<<" n_1 "<<data[2].normal[1]<<"\n";


        
         
        
       
    REAL Km = bc.Val1()(0, 0);
    REAL InvKm = 1./Km;
 //   std::cout<<"Km "<<Km<<" InvKm "<<InvKm<<"\n";
    REAL errorEstimated =0.,errorReal = 0.;
        
        normalsigmarec = Km*(urec-u_D)+g;
        
//std::cout<<"normalsigmarec "<<normalsigmarec<<" normalsigmafem "<<normalsigmafem<<"\n";

    errorEstimated = InvKm * (normalsigmarec - normalsigmafem)* (normalsigmarec - normalsigmafem);
//std::cout<<"normflux "<<normflux<< " normalsigmafem "<<normalsigmafem<<"\n";
//        std::cout<<"----------"<<"\n";
    errorReal = InvKm * (normflux + normalsigmafem)* (normflux + normalsigmafem);
    errors[2] = errorReal;
    errors[3] = errorEstimated;
        
        
    }
}

template class TPZHDivErrorEstimateMaterial<TPZMixedDarcyFlow>;
template class TPZHDivErrorEstimateMaterial<TPZMixedElasticityND>;