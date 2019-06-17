//
//  TPZHDivErrorEstimateMaterial.cpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 11/06/19.
//

#include "TPZHDivErrorEstimateMaterial.h"
#include "pzaxestools.h"



TPZHDivErrorEstimateMaterial::TPZHDivErrorEstimateMaterial(int matid, int dim) : TPZMixedPoisson(matid,dim)
{
    
}

TPZHDivErrorEstimateMaterial::TPZHDivErrorEstimateMaterial() : TPZMixedPoisson()
{
    
}

TPZHDivErrorEstimateMaterial::TPZHDivErrorEstimateMaterial(const TPZMixedPoisson &copy) : TPZMixedPoisson(copy)
{
    
}

TPZHDivErrorEstimateMaterial::TPZHDivErrorEstimateMaterial(const TPZHDivErrorEstimateMaterial &copy) : TPZMixedPoisson(copy)
{
    
}

TPZHDivErrorEstimateMaterial::~TPZHDivErrorEstimateMaterial()
{
    
}

TPZHDivErrorEstimateMaterial &TPZHDivErrorEstimateMaterial::operator=(const TPZHDivErrorEstimateMaterial &copy)
{
    TPZMixedPoisson::operator=(copy);
    return *this;
}


void TPZHDivErrorEstimateMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    /**
     datavec[0] H1 mesh, local uk/grad v
     datavec[1] L2 mesh, restriction of local uk
     datavec[2] Hdiv mesh, sigmakE
     datavec[3] L2 mesh, ukE
     
     Implement the matrix
     |Sk Ck^T |  = |bk|
     |Ck  0   |    |uk|
     Sk = int_K graduk.gradv dx = int_K gradphi_i.gradphi_j dx
     CK = int_K uk dx = int_K phi_i dx
     bk = int_K sigmafem.gradv dx = int_K sigmafem.gradphi_i dx
     ck = int_K ukfem dx
   
     **/
    
    int dim = datavec[0].axes.Rows();
    //defining test functions
    // Setting the phis
    TPZFMatrix<REAL> &phiuk = datavec[0].phi;
    TPZFMatrix<REAL> &phirest = datavec[2].phi;// function of restriction term
    TPZFMatrix<REAL> &dphiukaxes = datavec[0].dphix;
//    TPZFMatrix<REAL> &dphiv = datavec[0].dphix;
    TPZFNMatrix<9,REAL> dphiuk(3,dphiukaxes.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiukaxes, dphiuk, datavec[0].axes);
    
   //  TPZFMatrix<REAL> &axes = datavec[0].axes;
    
    
    int nphiuk = phiuk.Rows();
    

    
    TPZFMatrix<STATE> solsigmafem(3,nphiuk),solukfem(1,1);
    solsigmafem.Zero();
    solukfem.Zero();
    
    //potetial fem
    solukfem(0,0) = datavec[3].sol[0][0];
    //flux fem
    for (int ip = 0; ip<3; ip++)
    {
        solsigmafem(ip,0) = datavec[2].sol[0][ip];
    }
    
   
    TPZFNMatrix<9,REAL> PermTensor = fTensorK;
    TPZFNMatrix<9,REAL> InvPermTensor = fInvK;
    
    

    TPZFMatrix<STATE> kgraduk(3,nphiuk,0.);
    
        
    for(int irow=0 ; irow<nphiuk; irow++){
        
        //K graduk
        for(int i=0; i< dim; i++){
            
            for(int jd=0; jd< dim;jd++){
                
                kgraduk(i,irow) += PermTensor(i,jd)*dphiuk(jd,irow);
                
            }
            //bk=int_k sigmaukfem.grad phi_i,here dphiuk is multiplied by axes
            
            ef(irow,0)+=weight*dphiuk(i,irow)*solsigmafem(i,0);
        }
        
        //matrix Sk= int_{K} K graduk.gradv
        for(int jcol=0; jcol<nphiuk;jcol++){
            
            for(int jd=0;  jd< dim; jd++)
            {
                ek(irow,jcol) +=weight*kgraduk(jd,irow)*dphiuk(jd,jcol);
            }
            
        }
        //Ck=int_{K} phi_i e Ck^t
        if(fNeumannLocalProblem){
            
            ek(irow,nphiuk)+= weight*phiuk(irow,0);
            ek(nphiuk,irow)+= weight*phiuk(irow,0);
            
        }
        
    }
    if(!fNeumannLocalProblem)
    {
        ek(nphiuk,nphiuk) += weight;
    }
    
    //muk = int_k ukfem
    ef(nphiuk,0)+= weight*solukfem(0,0);
    


    
}


void TPZHDivErrorEstimateMaterial::FillDataRequirements(TPZVec<TPZMaterialData > &datavec) {
    
    //fem solution for flux and potential
    datavec[2].SetAllRequirements(false);
    datavec[2].fNeedsSol = true;
    
    datavec[3].SetAllRequirements(false);
    datavec[3].fNeedsSol = true;
}

void TPZHDivErrorEstimateMaterial::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{
    /**
     datavec[0] H1 mesh, uh_reconstructed
     datavec[1] L2 mesh,
     datavec[2] Hdiv fem mesh, sigma_h
     datavec[3] L2 mesh fem, u_h
     
      error[0] - error computed with exact pressure
      error[1] - error computed with reconstructed pressure
      error[2] - energy error computed with exact solution
      error[3] - energy error computed with reconstructed solution
     **/
    
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);
    
    
    TPZManVector<STATE,3> fluxfem(3),fluxreconstructed(3), pressurefem(1), pressurereconstructed(1);
    
    for (int i=0; i<3; i++) {
        fluxfem[i] = data[2].sol[0][i];
       
    }
    
    TPZFMatrix<REAL> &dsolaxes = data[0].dsol[0];
    TPZFNMatrix<9,REAL> fluxrec(3,0);
    TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, fluxrec, data[0].axes);
    
    for(int id=0 ; id<3; id++) {
        fluxreconstructed[id] = fluxrec(id,0);
    }
    
    
    if(this->fForcingFunctionExact){
        
        this->fForcingFunctionExact->Execute(data[0].x,u_exact,du_exact);
    }
    
    pressurereconstructed[0] = data[0].sol[0][0];
  
    pressurefem[0] = data[3].sol[0][0];
    
    TPZFNMatrix<9,REAL> PermTensor = fTensorK;
    TPZFNMatrix<9,REAL> InvPermTensor = fInvK;
    
    
    TPZFNMatrix<3,REAL> fluxexactneg;
    
    //sigma=-K grad(u)
    
    {
        TPZFNMatrix<9,REAL> gradpressure(3,1);
        for (int i=0; i<3; i++) {
            gradpressure(i,0) = du_exact[i];
        }
        PermTensor.Multiply(gradpressure,fluxexactneg);
    }
    
    REAL innerexact = 0.;
    REAL innerestimate = 0.;
    
#ifdef PZDEBUG
    std::cout<<"flux fem "<<fluxfem<<std::endl;
    std::cout<<"flux reconst "<<fluxreconstructed<<std::endl;
    std::cout<<"-------"<<std::endl;
#endif
    
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            innerexact += (fluxfem[i]+fluxexactneg(i,0))*InvPermTensor(i,j)*(fluxfem[j]+fluxexactneg(j,0));//Pq esta somando: o fluxo fem esta + e o exato -
            innerestimate += (fluxfem[i]-fluxreconstructed[i])*InvPermTensor(i,j)*(fluxfem[j]-fluxreconstructed[j]);
        }
    }
    
#ifdef PZDEBUG
    std::cout<<"potential fem "<<pressurefem<<std::endl;
    std::cout<<"potential reconst "<<pressurereconstructed<<std::endl;
    std::cout<<"-------"<<std::endl;
#endif
    errors[0] = (pressurefem[0]-u_exact[0])*(pressurefem[0]-u_exact[0]);//exact error pressure
    errors[1] = (pressurefem[0]-pressurereconstructed[0])*(pressurefem[0]-pressurereconstructed[0]);//error pressure reconstructed
    errors[2] = innerexact;//error flux exact
    errors[3] = innerestimate;//error flux reconstructed

    
    
    
}


int TPZHDivErrorEstimateMaterial::VariableIndex(const std::string &name)
{
    if(name == "FluxFem") return 40;
    if(name == "FluxExact") return 42;
    if(name == "PressureFem") return 43;
    if(name == "PressureReconstructed") return 44;
    if(name == "FluxReconstructed") return 41;
    if(name == "PressureExact") return 45;
    if(name == "PressureErrorExact") return 100;
    if(name == "PressureErrorEstimate") return 101;
    if(name == "EnergyErrorExact") return 102;
    if(name == "EnergyErrorEstimate") return 103;
    if(name == "PressureEffectivityIndex") return 104;
    if(name == "EnergyEffectivityIndex") return 105;
    if(name == "POrder") return 46;
    return -1;
}


int TPZHDivErrorEstimateMaterial::NSolutionVariables(int var)
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

void TPZHDivErrorEstimateMaterial::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    /**
     datavec[0] H1 mesh, uh_reconstructed
     datavec[1] L2 mesh,
     datavec[2] Hdiv fem mesh, sigma_h
     datavec[3] L2 mesh fem, u_h
     **/
    
    TPZFNMatrix<9,REAL> PermTensor = fTensorK;
    TPZFNMatrix<9,REAL> InvPermTensor = fInvK;
    
    
    TPZManVector<STATE,2> pressexact(1,0.);
    TPZFNMatrix<9,STATE> gradu(3,1,0.), fluxinv(3,1);
    
    if(fForcingFunctionExact)
    {
        this->fForcingFunctionExact->Execute(datavec[0].x, pressexact,gradu);
       
    }
    
    PermTensor.Multiply(gradu, fluxinv);
    
    int dim=this->fDim;
    switch (var)
    {
        case 40://FluxFem
            for(int i=0; i<dim; i++) Solout[i] = datavec[2].sol[0][i];
            break;
        case 41:{//FluxReconstructed(???)
            
            TPZFMatrix<REAL> &dsolaxes = datavec[0].dsol[0];
            TPZFNMatrix<9,REAL> dsol(3,0);
            TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, dsol, datavec[0].axes);
            
            for(int id=0 ; id<fDim; id++) {
                Solout[id] = dsol(id,0);//derivate
            }
            //for (int i=0; i<dim; i++) Solout[i] = datavec[0].dsol[0][i];
        }
           
            
            break;
        case 42://flux exact
            for(int i=0; i<dim; i++) Solout[i] = -fluxinv(i);
            break;
        case 43://PressureFem
            Solout[0] = datavec[3].sol[0][0];
            break;
        case 44://PressureReconstructed
            Solout[0] = datavec[0].sol[0][0];
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

