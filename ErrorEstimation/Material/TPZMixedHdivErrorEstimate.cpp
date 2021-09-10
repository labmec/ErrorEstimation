//
//  TPZMixedHDivErrorEstimate.cpp
//  ErrorEstimate
//
//  Created by Philippe Devloo on 20/04/18.
//

#include "TPZMixedHdivErrorEstimate.h"
#include "pzaxestools.h"

TPZMixedHDivErrorEstimate::TPZMixedHDivErrorEstimate() : TPZMixedDarcyFlow()
{
    
}

TPZMixedHDivErrorEstimate::TPZMixedHDivErrorEstimate(int matid, int dim) : TPZMixedDarcyFlow(matid,dim)
{
    
}

TPZMixedHDivErrorEstimate::TPZMixedHDivErrorEstimate(const TPZMixedDarcyFlow &copy) : TPZMixedDarcyFlow(copy)
{
    
}

TPZMixedHDivErrorEstimate::~TPZMixedHDivErrorEstimate() = default;

TPZMixedHDivErrorEstimate::TPZMixedHDivErrorEstimate(const TPZMixedHDivErrorEstimate &cp) : TPZMixedDarcyFlow(cp)
{
    
}

TPZMixedHDivErrorEstimate &TPZMixedHDivErrorEstimate::operator=(const TPZMixedHDivErrorEstimate &copy)
{
    TPZMixedDarcyFlow::operator=(copy);
    return *this;
}

void TPZMixedHDivErrorEstimate::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    TPZMixedDarcyFlow::FillDataRequirements(datavec);
    int i = 1;
    datavec[i].SetAllRequirements(false);
    datavec[i].fNeedsSol = true;
}

int TPZMixedHDivErrorEstimate::VariableIndex(const std::string &name) const
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
    if(name == "ResidualError") return 104;
    if(name == "PressureEffectivityIndex") return 105;
    if(name == "EnergyEffectivityIndex") return 106;
    if(name == "POrder") return 46;
   
    return -1;
}

int TPZMixedHDivErrorEstimate::NSolutionVariables(int var) const {
    switch (var) {
        case 40:
        case 41:
        case 42:
            return 3;
        case 43:
        case 44:
        case 45:
        case 46:
        case 100:
        case 101:
        case 102:
        case 103:
        case 104:
        case 105:
        case 106:
            return 1;
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
void
TPZMixedHDivErrorEstimate::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) {
    /**
     datavec[0]= Hdiv Resconstructed
     datavec[1]= Pressure Resconstructed
     datavec[2]= Hdiv FEM
     datavec[3]= Pressure FEM
    
     **/

    const auto perm = GetPermeability(datavec[1].x);

    TPZManVector<STATE, 2> pressexact(1, 0.);
    TPZFNMatrix<9, STATE> gradu(3, 1, 0.), fluxinv(3, 1);

    if (fExactSol) {
        this->fExactSol(datavec[0].x, pressexact, gradu);
    }

    for (int i = 0; i < 3; i++) {
        fluxinv(i, 0) = perm * gradu(i, 0);
    }

    int dim = this->fDim;

    switch (var)
    {
        case 40://FluxFem
            for(int i=0; i<dim; i++) Solout[i] = datavec[2].sol[0][i];
            break;
        case 41://FluxReconstructed
            for (int i=0; i<dim; i++) Solout[i] = datavec[0].sol[0][i];
            break;
        case 42:
            for(int i=0; i<dim; i++) Solout[i] = -fluxinv(i);
            break;
        case 43://PressureFem
            Solout[0] = datavec[3].sol[0][0];
            break;
        case 44://PressureReconstructed
            Solout[0] = datavec[1].sol[0][0];
            break;
        case 45:
            Solout[0] = pressexact[0];
            break;
        case 46:
            Solout[0] = datavec[1].p;
            break;
        default:
            DebugStop();
    }
}


/// make a contribution to the error computation
void TPZMixedHDivErrorEstimate::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {

    /**
     datavec[0]= Hdiv Reconstructed
     datavec[1]= Pressure Reconstructed
     datavec[2]= Hdiv FEM
     datavec[3]= Pressure FEM
     **/

    errors.Resize(NEvalErrors());
    errors.Fill(0.0);

    const int dim = this->fDim;

    TPZManVector<STATE,3> fluxfem(3), fluxreconstructed(3), pressurefem(1), pressurereconstructed(1);


    fluxreconstructed = data[0].sol[0];
    fluxfem = data[2].sol[0];
    STATE divsigmafem=data[2].divsol[0][0];


    TPZVec<STATE> divsigma(1,0.);

    TPZManVector<STATE, 1> u_exact(1);
    TPZFNMatrix<9, STATE> du_exact(3, 3);
    if(this->fExactSol){
        this->fExactSol(data[0].x,u_exact,du_exact);
        this->fForcingFunction(data[0].x,divsigma);
    }

    REAL residual = 0.;

    residual = (divsigma[0] - divsigmafem)*(divsigma[0] - divsigmafem);

    pressurereconstructed[0] = data[1].sol[0][0];
    pressurefem[0] = data[3].sol[0][0];

    const auto perm = GetPermeability(data[1].x);

    //sigma=-K grad(u)
    TPZFNMatrix<3, REAL> fluxexactneg(3, 1);

    {
        TPZFNMatrix<9, REAL> gradpressure(3, 1);
        for (int i = 0; i < 3; i++) {
            fluxexactneg(i, 0) = perm * du_exact[i];
        }
    }

    REAL innerexact = 0.;
    REAL innerestimate = 0.;
    for (int i = 0; i < dim; i++) {
        innerexact += (fluxfem[i] + fluxexactneg(i, 0)) * (1 / perm) * (fluxfem[i] + fluxexactneg(i, 0));
        innerestimate += (fluxfem[i] - fluxreconstructed[i]) * (1 / perm) * (fluxfem[i] - fluxreconstructed[i]);
    }
    errors[0] = (pressurefem[0]-u_exact[0])*(pressurefem[0]-u_exact[0]);//exact error pressure
    errors[1] = (pressurefem[0]-pressurereconstructed[0])*(pressurefem[0]-pressurereconstructed[0]);//error pressure reconstructed
    errors[2] = innerexact;//error flux exact
    errors[3] = innerestimate;//error flux reconstructed
    errors[4] = residual;
    for(int i = 0; i<5; i++)
    {
        if(std::isnan(errors[i])) DebugStop();
    }
}
