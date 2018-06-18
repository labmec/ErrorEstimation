//
//  TPZMixedHDivErrorEstimate.cpp
//  ErrorEstimate
//
//  Created by Philippe Devloo on 20/04/18.
//

#include "TPZMixedHdivErrorEstimate.h"
#include "mixedpoisson.h"
#include "pzaxestools.h"


template<class MixedMat>
TPZMixedHDivErrorEstimate<MixedMat>::TPZMixedHDivErrorEstimate() : MixedMat()
{
    
}

template<class MixedMat>
TPZMixedHDivErrorEstimate<MixedMat>::TPZMixedHDivErrorEstimate(int matid, int dim) : MixedMat(matid,dim)
{
    
}

template<class MixedMat>
TPZMixedHDivErrorEstimate<MixedMat>::~TPZMixedHDivErrorEstimate()
{
    
}

template<class MixedMat>
TPZMixedHDivErrorEstimate<MixedMat>::TPZMixedHDivErrorEstimate(const TPZMixedHDivErrorEstimate &cp) : MixedMat(cp), fSignConvention(cp.fSignConvention)
{
    
}

template<class MixedMat>
TPZMixedHDivErrorEstimate<MixedMat> &TPZMixedHDivErrorEstimate<MixedMat>::operator=(const TPZMixedHDivErrorEstimate &copy)
{
    MixedMat::operator=(copy);
    fSignConvention = copy.fSignConvention;
    return *this;
}

template<class MixedMat>
void TPZMixedHDivErrorEstimate<MixedMat>::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
    MixedMat::FillDataRequirements(datavec);
    {
        int i = 1;
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
    }
    
}

template<class MixedMat>
int TPZMixedHDivErrorEstimate<MixedMat>::VariableIndex(const std::string &name)
{
    if(name == "FluxFem") return 40;
    if(name == "FluxReconstructed") return 41;
    if(name == "FluxExact") return 42;
    if(name == "PressureFem") return 43;
    if(name == "PressureReconstructed") return 44;
    if(name == "PressureExact") return 45;
    if(name == "PressureErrorExact") return 100;
    if(name == "PressureErrorEstimate") return 101;
    if(name == "EnergyErrorExact") return 102;
    if(name == "EnergyErrorEstimate") return 103;
    if(name == "PressureEffectivityIndex") return 104;
    if(name == "EnergyEffectivityIndex") return 105;
    return -1;
}

template<class MixedMat>
int TPZMixedHDivErrorEstimate<MixedMat>::NSolutionVariables(int var)
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
        case 100:
        case 101:
        case 102:
        case 103:
        case 104:
        case 105:
            return 1;
            break;
        default:
            DebugStop();
            break;
    }
    return 0;
}

/**
 * @brief It return a solution to multiphysics simulation.
 * @param datavec [in] Data material vector
 * @param var [in] number of solution variables. See  NSolutionVariables() method
 * @param Solout [out] is the solution vector
 */
template<class MixedMat>
void TPZMixedHDivErrorEstimate<MixedMat>::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    TPZFNMatrix<9,REAL> PermTensor = MixedMat::fTensorK;
    TPZFNMatrix<9,REAL> InvPermTensor = MixedMat::fInvK;
    //int rtens = 2*fDim;
    if(MixedMat::fPermeabilityFunction){
        PermTensor.Redim(3,3);
        InvPermTensor.Redim(3,3);
        TPZFNMatrix<3,STATE> resultMat;
        TPZManVector<STATE> res;
        MixedMat::fPermeabilityFunction->Execute(datavec[1].x,res,resultMat);
        for(int id=0; id<3; id++){
            for(int jd=0; jd<3; jd++){
                PermTensor(id,jd) = resultMat(id,jd);
                InvPermTensor(id,jd) = resultMat(id+3,jd);
            }
        }
    }

    STATE pressureexact = 0.;
    TPZManVector<STATE,2> pressvec(1,0.);
    TPZFNMatrix<9,STATE> gradu(3,1,0.), fluxinv(3,1);
    if(MixedMat::fForcingFunctionExact)
    {
        MixedMat::fForcingFunctionExact->Execute(datavec[0].x, pressvec,gradu);
        gradu.Resize(3, 1);
        gradu(2,0) = 0.;
    }
    PermTensor.Multiply(gradu, fluxinv);
    pressureexact = pressvec[0];
    switch (var)
    {
        case 40:
            for(int i=0; i<3; i++) Solout[i] = datavec[0].sol[0][i+3];
            break;
        case 41:
            for (int i=0; i<3; i++) Solout[i] = datavec[0].sol[0][i];
            break;
        case 42:
            for(int i=0; i<3; i++) Solout[i] = -fluxinv(i);
            break;
        case 43:
            Solout[0] = datavec[1].sol[0][1];
            break;
        case 44:
            Solout[0] = datavec[1].sol[0][0];
            break;
        case 45:
            Solout[0] = pressureexact;
            break;
        default:
            DebugStop();
    }
}


/// make a contribution to the error computation
template<class MixedMat>
void TPZMixedHDivErrorEstimate<MixedMat>::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);

    TPZManVector<STATE,3> fluxfem(3), fluxreconstructed(3), pressurefem(1), pressurereconstructed(1);
    
    for (int i=0; i<3; i++) {
        fluxreconstructed[i] = data[0].sol[0][i];
        fluxfem[i] = data[0].sol[0][i+3];
    }
    pressurereconstructed[0] = data[1].sol[0][0];
    pressurefem[0] = data[1].sol[0][1];
    
    TPZFNMatrix<9,REAL> PermTensor = MixedMat::fTensorK;
    TPZFNMatrix<9,REAL> InvPermTensor = MixedMat::fInvK;
    //int rtens = 2*fDim;
    if(MixedMat::fPermeabilityFunction){
        PermTensor.Redim(3,3);
        InvPermTensor.Redim(3,3);
        TPZFNMatrix<3,STATE> resultMat;
        TPZManVector<STATE> res;
        MixedMat::fPermeabilityFunction->Execute(data[1].x,res,resultMat);
        for(int id=0; id<3; id++){
            for(int jd=0; jd<3; jd++){
                PermTensor(id,jd) = resultMat(id,jd);
                InvPermTensor(id,jd) = resultMat(id+3,jd);
            }
        }
    }
    
    
    TPZFNMatrix<3,REAL> fluxexactneg;
    
    {
        TPZFNMatrix<9,REAL> gradpressure(3,1);
        for (int i=0; i<3; i++) {
            gradpressure(i,0) = du_exact[i];
        }
        PermTensor.Multiply(gradpressure,fluxexactneg);
    }
    REAL innerexact = 0.;
    REAL innerestimate = 0.;
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            innerexact += (fluxfem[i]+fluxexactneg(i,0))*InvPermTensor(i,j)*(fluxfem[j]+fluxexactneg(j,0));
            innerestimate += (fluxfem[i]-fluxreconstructed[i])*InvPermTensor(i,j)*(fluxfem[j]-fluxreconstructed[j]);
        }
    }
    errors[0] = (pressurefem[0]-u_exact[0])*(pressurefem[0]-u_exact[0]);
    errors[1] = (pressurefem[0]-pressurereconstructed[0])*(pressurefem[0]-pressurereconstructed[0]);
    errors[2] += innerexact;
    errors[3] += innerestimate;
}



template class TPZMixedHDivErrorEstimate<TPZMixedPoisson>;
