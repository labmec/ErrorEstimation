//
//  TPZHDivErrorEstimateMaterial.cpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 11/06/19.
//

#include "TPZHDivErrorEstimateMaterial.h"
#include "pzaxestools.h"
#include "TPZAnalyticSolution.h"
#include "pzbndcond.h"


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

int TPZHDivErrorEstimateMaterial::IsH1Position(TPZVec<TPZMaterialData> &datavec){

    int nvec = datavec.NElements();
    int firstNoNullposition = -1;

    for(int ivec = 0; ivec < nvec ; ivec++){
        TPZMaterialData::MShapeFunctionType shapetype = datavec[ivec].fShapeType;
        if(shapetype != TPZMaterialData::EEmpty){
            firstNoNullposition = ivec;
            return firstNoNullposition;
        }

    }
    return firstNoNullposition;
}

void TPZHDivErrorEstimateMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    /**
   
     datavec[0] H1 mesh, local uk/grad v for Mark reconstruction and Empty for H1 reconstruction
     datavec[1] L2 mesh, restriction of local uk
     datavec[2] Hdiv mesh, sigmakE
     datavec[3] L2 mesh, ukE
     
     Implement the matrix
     |Sk Ck^T |  = |bk|
     |Ck  0   |    |uk|
     Sk = int_K graduk.gradv dx = int_K gradphi_i.gradphi_j dx
     CK = int_K uk dx = int_K phi_i dx
     bk = int_K sigmafem.gradv dx = int_K sigmafem.gradphi_i dx
     uk = int_K ukfem dx
   
     **/
    
    int H1functionposition = 0;
    H1functionposition = IsH1Position(datavec);
    
    int dim = datavec[H1functionposition].axes.Rows();
    //defining test functions
    // Setting the phis
    TPZFMatrix<REAL> &phiuk = datavec[H1functionposition].phi;
    TPZFMatrix<REAL> &dphiukaxes = datavec[H1functionposition].dphix;
    TPZFNMatrix<9,REAL> dphiuk(3,dphiukaxes.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiukaxes, dphiuk, datavec[H1functionposition].axes);
    
    
    int nphiuk = phiuk.Rows();
    
    TPZFMatrix<STATE> solsigmafem(3,nphiuk),solukfem(1,1);
    solsigmafem.Zero();
    solukfem.Zero();
    
    


        //potetial fem
        solukfem(0,0) = datavec[3].sol[0][0];
        //flux fem
        for (int ip = 0; ip<3; ip++){

            solsigmafem(ip,0) = datavec[2].sol[0][ip];
        }

    
    
    TPZFNMatrix<9,REAL> PermTensor;
    TPZFNMatrix<9,REAL> InvPermTensor;
    
    GetPermeabilities(datavec[1].x, PermTensor, InvPermTensor);
    
    
    TPZFMatrix<STATE> kgraduk(3,nphiuk,0.);
    
        
    for(int irow=0 ; irow<nphiuk; irow++){
        
        //K graduk
        for(int id=0; id< dim; id++){
            
            for(int jd=0; jd< dim;jd++){
                
                kgraduk(id,irow) += PermTensor(id,jd)*dphiuk(jd,irow);
                
            }
            //bk = (-1)*int_k sigmaukfem.grad phi_i,here dphiuk is multiplied by axes
            //the minus sign is necessary because we are working with sigma_h = - K grad u, Mark works with sigma_h = K grad u
            
            ef(irow,0)+=(-1.)*weight*dphiuk(id,irow)*solsigmafem(id,0);
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
    if(H1functionposition ==0){
    
        if(!fNeumannLocalProblem )
        {
            ek(nphiuk,nphiuk) += weight;
            ef(nphiuk,0)+= weight;
        }
        
        //muk = int_k ukfem
        
        else{
            
            ef(nphiuk,0)+= weight*solukfem(0,0);
        
        }
    }

    


    
}


void TPZHDivErrorEstimateMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
 /*
  Add Robin boundary condition for local problem
  ek+= <w,Km s_i>
  ef+= <w,Km*u_d + g + sigma_i.n>
  */
    int H1functionposition = 0;
    H1functionposition = IsH1Position(datavec);
    
    int dim = datavec[H1functionposition].axes.Rows();
    
    TPZFMatrix<REAL>  &phi_i = datavec[H1functionposition].phi;
    int nphi_i = phi_i.Rows();
    TPZFMatrix<STATE> solsigmafem(3,1);
    solsigmafem.Zero();
    
   
    TPZManVector<REAL,3> normal = datavec[2].normal;
    std::cout<<normal[0]<<", "<<normal[1]<<","<<normal[2]<<std::endl;

    REAL normalsigma = 0.;
    
    int nsol = datavec[2].sol[0].size();//ver pq nao é vetorial
    
    for (int ip = 0; ip< nsol; ip++){

        solsigmafem(ip,0) = datavec[2].sol[0][ip];
        normalsigma += datavec[2].sol[0][ip]*normal[ip];
    }
    
  
    
//tem que resover o UpdateBcValues para receber o g correto
    REAL g = bc.Val1()(1,0);//g
    REAL Km = bc.Val1()(0,0);//Km
    REAL u_D = 0.;
    REAL robinterm = 0.;
    if(bc.HasForcingFunction())
    {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9,STATE> gradu(dim,1);
        bc.ForcingFunction()->Execute(datavec[H1functionposition].x,res,gradu);
        u_D = res[0];
    

    }
    else{
        u_D = bc.Val2()(0,0);//usualmente updatebc coloca o valor exato no val2
    }

    
    if (bc.Type()==4) {
        robinterm = Km*u_D +g + normalsigma ;
            
            for(int iq = 0; iq < nphi_i; iq++) {
                //<w,Km*u_D+g+sigma_i*n>
                ef(iq,0) += robinterm*phi_i(iq,0)*weight;
                for (int jq = 0; jq < nphi_i; jq++) {
                    //<w,Km*s_i>
                    ek(iq,jq) += weight*Km*phi_i(iq,0)*phi_i(jq,0);
                }
            }
        
    }
    
    else{
        std::cout<<" Not implemented yet"<<std::endl;
        
    }
    
    
}


void TPZHDivErrorEstimateMaterial::FillDataRequirements(TPZVec<TPZMaterialData > &datavec) {

        //fem solution for flux and potential
        datavec[2].SetAllRequirements(false);
        datavec[2].fNeedsSol = true;
        datavec[2].fNeedsNormal = true;
        
        datavec[3].SetAllRequirements(false);
        datavec[3].fNeedsSol = true;
  

}

void TPZHDivErrorEstimateMaterial::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec){

    datavec[2].SetAllRequirements(false);
    datavec[2].fNeedsSol = true;
    datavec[2].fNeedsNormal = true;
        
    
}

void TPZHDivErrorEstimateMaterial::UpdateBCValues(TPZVec<TPZMaterialData> & datavec){
    DebugStop();
    
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
      error[4] - oscilatory data error
     **/
    
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);
    
    
    TPZManVector<STATE,3> fluxfem(3);
    STATE divsigmafem, pressurefem, pressurereconstructed;
    
    TPZFNMatrix<3,REAL> fluxreconstructed(3,1), fluxreconstructed2(3,1);
    

 
    fluxfem=data[2].sol[0];
    divsigmafem= data[2].divsol[0][0];
    
 

    int H1functionposition = 0;
    H1functionposition = IsH1Position(data);

   
    TPZFMatrix<REAL> &dsolaxes = data[H1functionposition].dsol[0];
    TPZFNMatrix<9,REAL> fluxrec(3,0);
    TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, fluxrec, data[H1functionposition].axes);
    
    for(int id=0 ; id<3; id++) {
        fluxreconstructed2(id,0) = (-1.)*fluxrec(id,0);
    }
    
    
    TPZVec<STATE> divsigma(1);
    
    if(this->fForcingFunctionExact){
        
        this->fForcingFunctionExact->Execute(data[H1functionposition].x,u_exact,du_exact);
    
        this->fForcingFunction->Execute(data[H1functionposition].x,divsigma);
    }
    
    REAL residual = 0.;
    
   // std::cout<<" divsigma[0] "<<divsigma[0]<<" divsigmafem "<<divsigmafem<<"\n";
    
    residual = (divsigma[0] - divsigmafem)*(divsigma[0] - divsigmafem);
   
    
    pressurereconstructed = data[H1functionposition].sol[0][0];
  
 
        pressurefem = data[3].sol[0][0];
    
    TPZFNMatrix<9,REAL> PermTensor;
    TPZFNMatrix<9,REAL> InvPermTensor;
    
    GetPermeabilities(data[2].x, PermTensor, InvPermTensor);
    
    TPZFNMatrix<3,REAL> fluxexactneg;
    
    //sigmarec = -K grad(urec)
    //  sigmak = -K graduk
    
    {
        TPZFNMatrix<9,REAL> gradpressure(3,1);
        for (int i=0; i<3; i++) {
            gradpressure(i,0) = du_exact[i];
        }
        PermTensor.Multiply(gradpressure,fluxexactneg);
    }
    
     PermTensor.Multiply(fluxreconstructed2,fluxreconstructed);
    
    
    REAL innerexact = 0.;
    REAL innerestimate = 0.;
    
    
    
#ifdef PZDEBUG2
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
    
#ifdef PZDEBUG2
    std::cout<<"potential fem "<<pressurefem<<std::endl;
    std::cout<<"potential reconst "<<pressurereconstructed<<std::endl;
    std::cout<<"-------"<<std::endl;
#endif
    errors[0] = (pressurefem-u_exact[0])*(pressurefem-u_exact[0]);//exact error pressure
    errors[1] = (pressurefem-pressurereconstructed)*(pressurefem-pressurereconstructed);//error pressure reconstructed
    errors[2] = innerexact;//error flux exact
    errors[3] = innerestimate;//error flux reconstructed
    errors[4] = residual; //||f - Proj_divsigma||

    
    
    
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
    if(name == "ResidualError") return 104;
    if(name == "PressureEffectivityIndex") return 105;
    if(name == "EnergyEffectivityIndex") return 106;
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

/**
 * @brief It return a solution to multiphysics simulation.
 * @param datavec [in] Data material vector
 * @param var [in] number of solution variables. See  NSolutionVariables() method
 * @param Solout [out] is the solution vector
 */

void TPZHDivErrorEstimateMaterial::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{

    /**
     datavec[0] H1 mesh, uh_reconstructed for Mark reconstruction and Empty for H1 reconstruction
     datavec[1] L2 mesh,
     datavec[2] Hdiv fem mesh, sigma_h
     datavec[3] L2 mesh fem, u_h
     **/
    
    int H1functionposition = 0;
    H1functionposition = IsH1Position(datavec);
    
    TPZFNMatrix<9,REAL> PermTensor;
    TPZFNMatrix<9,REAL> InvPermTensor;
    
    GetPermeabilities(datavec[2].x, PermTensor, InvPermTensor);
    
    
    TPZManVector<STATE,2> pressexact(1,0.);
    TPZFNMatrix<9,STATE> gradu(3,1,0.), fluxinv(3,1);
    
    if(fForcingFunctionExact)
    {
        this->fForcingFunctionExact->Execute(datavec[H1functionposition].x, pressexact,gradu);
       
    }
    
    PermTensor.Multiply(gradu, fluxinv);
    
    int dim=this->fDim;
    switch (var)
    {
        case 40://FluxFem
  
            for(int i=0; i<dim; i++) Solout[i] = datavec[2].sol[0][i];

            break;
        case 41:{//FluxReconstructed is grad U
            TPZFMatrix<REAL> &dsolaxes = datavec[H1functionposition].dsol[0];
            TPZFNMatrix<9,REAL> dsol(3,0);
             TPZFNMatrix<9,REAL> KGradsol(3,0);
            TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, dsol, datavec[H1functionposition].axes);
            
            PermTensor.Multiply(dsol,KGradsol);
            
            for(int id=0 ; id<fDim; id++) {
                Solout[id] = KGradsol(id,0);//dsol(id,0);//derivate
            }
            
            
 
        }
           
            
            break;
        case 42://flux exact
            for(int i=0; i<dim; i++) Solout[i] = -fluxinv(i);
            break;
        case 43://Pressure fem
            Solout[0] = datavec[3].sol[0][0];
       
            break;
        case 44://PressureReconstructed
            Solout[0] = datavec[H1functionposition].sol[0][0];

            break;
        case 45://pressureexact
            Solout[0] = pressexact[0];
            break;
        case 46://order p
            Solout[0] = datavec[1].p;
            break;

        default:
            DebugStop();
    }
}



