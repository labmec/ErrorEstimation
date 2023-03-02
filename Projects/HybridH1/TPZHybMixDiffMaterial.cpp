//
// Created by victor on 29/03/2021.
//

#include "TPZHybMixDiffMaterial.h"
#include "pzaxestools.h"

TPZHybMixDiffMaterial::TPZHybMixDiffMaterial(int matid, int dim) : TPZMixedPoisson(matid, dim)
{
}

TPZHybMixDiffMaterial::TPZHybMixDiffMaterial() : TPZMixedPoisson()
{
    
}

TPZHybMixDiffMaterial::TPZHybMixDiffMaterial(const TPZHybMixDiffMaterial &copy) : TPZMixedPoisson(copy)
{
}

TPZHybMixDiffMaterial::TPZHybMixDiffMaterial(const TPZMixedPoisson &copy) : TPZMixedPoisson(copy)
{
    
}

TPZHybMixDiffMaterial::TPZHybMixDiffMaterial(TPZMatLaplacianHybrid &matlaplacian)
{
    this->SetId(matlaplacian.Id());
    this->SetDimension(matlaplacian.Dimension());
    
    STATE Ks;
    Ks = matlaplacian.GetPermeability({0,0,0});
    this->SetConstantPermeability(Ks);
    
    if (matlaplacian.HasForcingFunction()) {
        int porder = matlaplacian.ForcingFunctionPOrder();
        this->SetForcingFunction(matlaplacian.ForcingFunction(),porder);
    }
    if(matlaplacian.HasExactSol())
    {
        int porder = matlaplacian.PolynomialOrderExact();
        this->SetExactSol(matlaplacian.ExactSol(), porder);
    }
}

TPZHybMixDiffMaterial::~TPZHybMixDiffMaterial()
{
    
}

TPZHybMixDiffMaterial &TPZHybMixDiffMaterial::operator=(const TPZHybMixDiffMaterial &copy)
{
    TPZMixedPoisson::operator=(copy);
    return *this;
}

void TPZHybMixDiffMaterial::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const {
    
    //fem solution for flux and potential
    for(int i =0 ; i< 3; i++){
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNormal = true;
    }
}

void TPZHybMixDiffMaterial::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const {
    for(int i =0 ; i< 3; i++){
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNormal = true;
    }
}

void TPZHybMixDiffMaterial::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors){
    /**
     datavec[0] HybH1 potential
     datavec[1] Mixed potential
     datavec[2] Mixed Flux
     
     error[0] - HybH1 potential exact error
     error[1] - Mixed potential exact error
     error[2] - HybH1 - Mixed error
     error[3] - HybH1 flux exact error (-k\nabla u)
     error[4] - Mixed flux exact error (\sigma)
     error[5] - HybH1 - Mixed flux
     
     **/
    
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);
    
    STATE hybP, mixP;
    
    hybP = data[0].sol[0][0];
    mixP = data[1].sol[0][0];
    
    TPZVec<STATE> u_exact(1,0);
    TPZFMatrix<STATE> du_exact(3,1,0);
    if(this->HasExactSol()){
        this->ExactSol()(data[0].x,u_exact,du_exact);
    }
    
    errors[0] = (hybP-u_exact[0])*(hybP-u_exact[0]);
    errors[1] = (mixP-u_exact[0])*(mixP-u_exact[0]);
    errors[2] = (hybP-mixP)*(hybP-mixP);
    
    TPZFNMatrix<3,REAL> gradHyb(3,1), mixF(3,1), hybF(3,1);
    
    gradHyb = data[0].dsol[0];
    for(int id=0 ; id<3; id++) {
        mixF(id,0) = data[2].sol[0][id];
    }
    
    gradHyb.Resize(3,1);
    
    TPZFNMatrix<9,REAL> PermTensor(3,3);
    TPZFNMatrix<9,REAL> InvPermTensor(3,3);
    STATE K = GetPermeability(data[0].x);
//    GetPermeabilities(data[0].x, PermTensor, InvPermTensor);
    PermTensor.Diagonal(K);
    InvPermTensor.Diagonal(1./K);
    
    PermTensor.Multiply(gradHyb,hybF);
    for(int ip = 0 ; ip < 3 ; ip++){
        hybF(ip,0) = -hybF(ip,0);
    }
    
    TPZFNMatrix<3,REAL> flux;
    
    {
        TPZFNMatrix<9,REAL> minusGradP(3,1);
        for (int i=0; i<3; i++) {
            minusGradP(i,0) = (-1.)*du_exact[i];
        }
        PermTensor.Multiply(minusGradP,flux);
    }
    
    STATE hybExactF = 0., mixExactF = 0., DiffF = 0.;
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            hybExactF += (hybF[i]-flux(i,0))*InvPermTensor(i,j)*(hybF[j]-flux(j,0));
            mixExactF += (mixF[i]-flux(i,0))*InvPermTensor(i,j)*(mixF[j]-flux(j,0));
            DiffF += (hybF[i]-mixF[i])*InvPermTensor(i,j)*(hybF[i]-mixF[i]);
        }
    }
    
    errors[3] = hybExactF;
    errors[4] = mixExactF;
    errors[5] = DiffF;
}


int TPZHybMixDiffMaterial::VariableIndex(const std::string &name) const
{
    if(name == "hybP") return 40;
    if(name == "mixP") return 41;
    if(name == "exactP") return 42;
    if(name == "DiffP") return 43;
    
    if(name == "hybF") return 50;
    if(name == "mixF") return 51;
    if(name == "exactF") return 52;
    if(name == "DiffF") return 53;
    
    return -1;
}

int TPZHybMixDiffMaterial::NSolutionVariables(int var) const
{
    switch (var) {
        case 40:
        case 41:
        case 42:
        case 43:
            return 3;
            break;
        case 50:
        case 51:
        case 52:
        case 53:
            return 1;
            break;
        default:
            DebugStop();
            break;
    }
    return 0;
}

void TPZHybMixDiffMaterial::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout)
{
    
    /**
     datavec[0] HybH1 potential
     datavec[1] Mixed potential
     datavec[2] Mixed Flux
     **/
    
    TPZFNMatrix<9,STATE> PermTensor(3,3);
    TPZFNMatrix<9,STATE> InvPermTensor(3,3);
    
    STATE K = GetPermeability(datavec[0].x);
//    GetPermeabilities(datavec[0].x, PermTensor, InvPermTensor);
    PermTensor.Diagonal(K);
    InvPermTensor.Diagonal(1./K);
    TPZManVector<STATE,2> pressexact(1,0.);
    TPZFNMatrix<9,STATE> gradu(3,1,0.), fluxinv(3,1);
    
    if(HasExactSol())
    {
        this->ExactSol()(datavec[0].x, pressexact,gradu);
    }
    
    PermTensor.Multiply(gradu, fluxinv);
    
    TPZFMatrix<REAL> &dsolaxes = datavec[0].dsol[0];
    TPZFNMatrix<9, REAL> dsol(3, 1);
    TPZFNMatrix<9, REAL> KGradsol(3, 1);
    TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, dsol, datavec[2].axes);
    
    PermTensor.Multiply(dsol, KGradsol);
    
    int dim=this->fDim;
    switch (var)
    {
        case 50://hybF
        {
            for (int i = 0; i < 3; i++) Solout[i]  = -KGradsol(i,0);
        }
            break;
            
        case 51:{// mixF
            for(int id=0 ; id<3; id++) Solout[id] = datavec[2].sol[0][id];
        }
            break;
            
        case 52://exactF
            for(int i=0; i<dim; i++) Solout[i] = -fluxinv(i);
            break;
            
        case 53://DiffF
            for(int i=0; i<dim; i++) Solout[i] = -KGradsol(i,0) - datavec[2].sol[0][i];
            break;
            
        case 40://hybP
            Solout[0] = datavec[0].sol[0][0];
            break;
            
        case 41://mixP
            Solout[0] = datavec[1].sol[0][0];
            break;
            
        case 42://exactP
            Solout[0] = pressexact[0];
            break;
            
        case 43://DiffP
            Solout[0] = datavec[0].sol[0][0] - datavec[1].sol[0][0];
            break;
        default:
            DebugStop();
    }
}
