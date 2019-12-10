//
//  TPZMHMHDivErrorEstimateMaterial.cpp
//
//  Created by Denise Siqueira on 22/07/19.
//

#include "TPZMHMHDivErrorEstimatorMaterial.h"
#include "pzaxestools.h"
#include "TPZAnalyticSolution.h"
#include "pzbndcond.h"

TPZMHMHDivErrorEstimateMaterial::TPZMHMHDivErrorEstimateMaterial(int matid, int dim) : TPZMixedPoisson(matid,dim)
{
    
}

TPZMHMHDivErrorEstimateMaterial::TPZMHMHDivErrorEstimateMaterial() : TPZMixedPoisson()
{
    
}

TPZMHMHDivErrorEstimateMaterial::TPZMHMHDivErrorEstimateMaterial(const TPZMixedPoisson &copy) : TPZMixedPoisson(copy)
{
    
}

TPZMHMHDivErrorEstimateMaterial::TPZMHMHDivErrorEstimateMaterial(const TPZMHMHDivErrorEstimateMaterial &copy) : TPZMixedPoisson(copy)
{
    
}

TPZMHMHDivErrorEstimateMaterial::~TPZMHMHDivErrorEstimateMaterial()
{
    
}

TPZMHMHDivErrorEstimateMaterial &TPZMHMHDivErrorEstimateMaterial::operator=(const TPZMHMHDivErrorEstimateMaterial &copy)
{
    TPZMixedPoisson::operator=(copy);
    return *this;
}


void TPZMHMHDivErrorEstimateMaterial::FillDataRequirements(TPZVec<TPZMaterialData > &datavec) {
    
    //fem solution for flux and potential
    datavec[2].SetAllRequirements(false);
    datavec[2].fNeedsSol = true;
    
    datavec[3].SetAllRequirements(false);
    datavec[3].fNeedsSol = true;
}

int TPZMHMHDivErrorEstimateMaterial::IsH1Position(TPZVec<TPZMaterialData> &datavec){
    
    int nvec = datavec.NElements();
    int firstNoNullposition =0;
    
    for(int ivec = 0; ivec < nvec ; ivec++){
        TPZMaterialData::MShapeFunctionType shapetype = datavec[ivec].fShapeType;
        if(shapetype != TPZMaterialData::EEmpty){
            firstNoNullposition = ivec;
            return firstNoNullposition;
        }
        
    }
    
}

void TPZMHMHDivErrorEstimateMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    /**
     
     datavec[0] empty data
     datavec[1] H1 mesh, restriction of local uk
     datavec[2] Hdiv mesh, sigma_h
     datavec[3] L2 mesh, u_h
     
     Implement the matrix
     |Ak  |  = |bk|
     
     Ak = int_K K grads_i.gradv dx = int_K gradphi_i.gradphi_j dx
     bk = int_K K sigma_h.gradv dx = int_K - sigma_i. gradu_h dx
     
     **/
    
    int H1functionposition = IsH1Position(datavec);
    if(H1functionposition == 0){
       // TPZMixedPoisson:: Contribute(datavec, weight, ek, ef);
        ContributeHdiv(datavec, weight, ek, ef);
        return;
    }
    
    
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
    
    //flux fem
    for (int ip = 0; ip<3; ip++){
        
        solsigmafem(ip,0) = datavec[2].sol[0][ip];
    }
    
    
    //potetial fem
    solukfem(0,0) = datavec[3].sol[0][0];
    TPZFMatrix<REAL> dsolukfem(3,1,0),Kgradukfem(3,1,0);
    
    for(int i=0; i<3 ;i++){
        dsolukfem(i,0) = datavec[H1functionposition].dsol[0][i];
    }
    
   // dsolukfem = datavec[H1functionposition].dsol[0];
    
    TPZFNMatrix<9,REAL> PermTensor = fTensorK;
    TPZFNMatrix<9,REAL> InvPermTensor = fInvK;
    
    TPZFMatrix<STATE> kgraduk(3,nphiuk,0.);
    
//    for(int i=0; i< dim; i++){
//
//        for(int jd=0; jd< dim;jd++){
//            Kgradukfem(i,0) += PermTensor(i,jd)*dsolukfem(jd,0);
//
//
//        }
//    }
    
    
    
    
    for(int irow=0 ; irow<nphiuk; irow++){
        
        //K graduk
        for(int id=0; id< dim; id++){
            
            for(int jd=0; jd< dim;jd++){
                kgraduk(id,irow) += PermTensor(id,jd)*dphiuk(jd,irow);
                
                
            }
            ef(irow,0) += (-1.)*weight*solsigmafem(id,0)*dphiuk(id,irow);
        }
        
        
        
        //matrix Sk= int_{K} K grads_i.gradv
        for(int jcol=0; jcol<nphiuk;jcol++){
            
            for(int jd=0;  jd< dim; jd++)
            {
                ek(irow,jcol) += weight*kgraduk(jd,irow)*dphiuk(jd,jcol);
                //ef(irow,0) += weight*Kgradukfem(jd,irow)*dphiuk(jd,jcol);
            }
            
        }
        
    }
    
}

void TPZMHMHDivErrorEstimateMaterial::ContributeHdiv(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {

    
    STATE force = ff;
    if(fForcingFunction) {
        TPZManVector<STATE> res(1);
        fForcingFunction->Execute(datavec[1].x,res);
        force = res[0];
    }
    
    TPZFNMatrix<9,REAL> PermTensor;
    TPZFNMatrix<9,REAL> InvPermTensor;
    
    GetPermeability(datavec[1].x, PermTensor, InvPermTensor);
    
    // Setting the phis
    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &phip = datavec[1].phi;
    TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> &dphiP = datavec[1].dphix;
    TPZFNMatrix<9,REAL> dphiPXY(3,dphiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPXY, datavec[1].axes);
    
    
    int phrq, phrp;
    phrp = phip.Rows();
    phrq = datavec[0].fVecShapeIndex.NElements();
    
    int nactive = 0;
    for (int i=0; i<datavec.size(); i++) {
        if (datavec[i].fActiveApproxSpace) {
            nactive++;
        }
    }
#ifdef PZDEBUG
    if(nactive == 4)
    {
        int phrgb = datavec[2].phi.Rows();
        int phrub = datavec[3].phi.Rows();
        if(phrp+phrq+phrgb+phrub != ek.Rows())
        {
            DebugStop();
        }
    }else
    {
        if(phrp+phrq != ek.Rows())
        {
            DebugStop();
        }
    }
#endif
    //Calculate the matrix contribution for flux. Matrix A
    for(int iq=0; iq<phrq; iq++)
    {
        //ef(iq, 0) += 0.;
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3,REAL> ivec(3,1,0.);
        for(int id=0; id<3; id++){
            ivec(id,0) = datavec[0].fNormalVec(id,ivecind);
        }
    
        
        TPZFNMatrix<3,REAL> ivecZ(3,1,0.);
        TPZFNMatrix<3,REAL> jvecZ(3,1,0.);
        for (int jq=0; jq<phrq; jq++)
        {
            TPZFNMatrix<3,REAL> jvec(3,1,0.);
            int jvecind = datavec[0].fVecShapeIndex[jq].first;
            int jshapeind = datavec[0].fVecShapeIndex[jq].second;
            
            for(int id=0; id<3; id++){
                jvec(id,0) = datavec[0].fNormalVec(id,jvecind);
            }
            
            //dot product between Kinv[u]v
            jvecZ.Zero();
            for(int id=0; id<fDim; id++){
                for(int jd=0; jd<fDim; jd++){
                    jvecZ(id,0) += InvPermTensor(id,jd)*jvec(jd,0);
                }
            }
            //jvecZ.Print("mat1 = ");
            REAL prod1 = ivec(0,0)*jvecZ(0,0) + ivec(1,0)*jvecZ(1,0) + ivec(2,0)*jvecZ(2,0);
            ek(iq,jq) += fvisc*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod1;
            
        }
    }
    
    
    // Coupling terms between flux and pressure. Matrix B
    for(int iq=0; iq<phrq; iq++)
    {
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        
        TPZFNMatrix<3,REAL> ivec(3,1,0.);
        for(int id=0; id<3; id++){
            ivec(id,0) = datavec[0].fNormalVec(id,ivecind);
        }
        TPZFNMatrix<3,REAL> axesvec(3,1,0.);
        datavec[0].axes.Multiply(ivec,axesvec);
        
        REAL divwq = 0.;
        for(int iloc=0; iloc<fDim; iloc++)
        {
            divwq += axesvec(iloc,0)*dphiQ(iloc,ishapeind);
        }
        for (int jp=0; jp<phrp; jp++) {
            
            REAL fact = (-1.)*weight*phip(jp,0)*divwq;
            // Matrix B
            ek(iq, phrq+jp) += fact;
            
            // Matrix B^T
            ek(phrq+jp,iq) += fact;
        }
    }
    
    //termo fonte referente a equacao da pressao
    for(int ip=0; ip<phrp; ip++){
        ef(phrq+ip,0) += (-1.)*weight*force*phip(ip,0);
    }
    //
//    #ifdef LOG4CXX
//        if(logdata->isDebugEnabled())
//        {
//            std::stringstream sout;
//            sout<<"\n\n Matriz ek e vetor fk \n ";
//            ek.Print("ekmph = ",sout,EMathematicaInput);
//            ef.Print("efmph = ",sout,EMathematicaInput);
//            LOGPZ_DEBUG(logdata,sout.str());
//        }
//    #endif
    
}


void TPZMHMHDivErrorEstimateMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
     int H1functionposition = IsH1Position(datavec);
    if(H1functionposition == 0) return;

    TPZFMatrix<REAL>  &phi_u = datavec[1].phi;
    int phr_primal = phi_u.Rows();

    short in,jn;
    STATE v2[1];
    if(bc.HasForcingFunction()) {
        TPZManVector<STATE> res(1);
        TPZFNMatrix<3,STATE> dres(3,1);
        bc.ForcingFunction()->Execute(datavec[H1functionposition].x,res,dres);
        v2[0] = res[0];
    }
    else{
           v2[0] = bc.Val2()(0,0);
    }
//
    switch (bc.Type()) {
        case 0 :            // Dirichlet condition
            for(in = 0 ; in < phr_primal; in++) {
                ef(in,0) += (STATE)(gBigNumber* phi_u(in,0) * weight) * v2[0];
                for (jn = 0 ; jn < phr_primal; jn++) {
                    ek(in,jn) += gBigNumber * phi_u(in,0) * phi_u(jn,0) * weight;
                }
            }
            break;
        case 1 :            // Neumann condition
            std::cout<<"Not Implemented yet"<<std::endl;
            DebugStop();
            break;
        case 2 :        // mixed condition

            std::cout<<"Not Implemented yet"<<std::endl;
            DebugStop();

            break;
    }

}





void TPZMHMHDivErrorEstimateMaterial::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{
    /**
     datavec[0] Flux
     datavec[1] L2 mesh,
     datavec[2] flux fem
     datavec[3] L2 mesh fem
     
      error[0] - error computed with exact pressure
      error[1] - error computed with reconstructed pressure
      error[2] - energy error computed with exact solution
      error[3] - energy error computed with reconstructed solution
      error[4] - residual data error
     
     If the resconstructionis done using H1, the errors are: ||gradufem-gradurec||, ||gradufem-gradex|| and ||f- Proj divfem||
     **/
    
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);
    
    
    TPZManVector<STATE,3> fluxfem(3), fluxreconstructed(3), pressurefem(1), pressurereconstructed(1);


    int H1functionposition = 0;
    
    H1functionposition = IsH1Position(data);
    
    if(H1functionposition==0){
    
        fluxreconstructed = data[0].sol[0];
    }
    else{

        for(int i=0; i< 3;i++){
        fluxreconstructed[i] = data[1].dsol[0][i];
        }
        
        
    }
    
    fluxfem = data[2].sol[0];
    STATE divsigmafem=data[2].divsol[0][0];
    
    
    TPZVec<STATE> divsigma(1);
    
    if(this->fForcingFunctionExact){
        
        this->fForcingFunctionExact->Execute(data[H1functionposition].x,u_exact,du_exact);
        this->fForcingFunction->Execute(data[H1functionposition].x,divsigma);
    }
    
    
    
    REAL residual = 0.;
    residual = (divsigma[0] - divsigmafem)*(divsigma[0] - divsigmafem);
    

    
    pressurereconstructed[0] = data[1].sol[0][0];
    pressurefem[0] = data[3].sol[0][0];

    
    TPZFNMatrix<9,REAL> PermTensor = fTensorK;
    TPZFNMatrix<9,REAL> InvPermTensor = fInvK;
    
    
    
    TPZFNMatrix<3,REAL> fluxexactneg;
    
    //sigma=-K grad(u)
    
    {
        TPZFNMatrix<9,REAL> gradpressure(3,1);
       
        
        for (int i=0; i<fDim; i++) {
            gradpressure(i,0) = du_exact[i];
        }
        PermTensor.Multiply(gradpressure,fluxexactneg);
    }
    
    
    REAL innerexact = 0.;
    REAL innerestimate = 0.;
    for (int i=0; i<fDim; i++) {
        for (int j=0; j<fDim; j++) {
            innerexact += (fluxfem[i]+fluxexactneg(i,0))*InvPermTensor(i,j)*(fluxfem[j]+fluxexactneg(j,0));//Pq esta somando: o fluxo fem esta + e o exato -
            innerestimate += (fluxfem[i]-fluxreconstructed[i])*InvPermTensor(i,j)*(fluxfem[j]-fluxreconstructed[j]);
        }
    }
    errors[0] = (pressurefem[0]-u_exact[0])*(pressurefem[0]-u_exact[0]);//exact error pressure
    errors[1] = (pressurefem[0]-pressurereconstructed[0])*(pressurefem[0]-pressurereconstructed[0]);//error pressure reconstructed
    errors[2] = innerexact;//error flux exact
    errors[3] = innerestimate;//error flux reconstructed
    errors[4] = residual; //||f - Proj_divsigma||
    

    
    
    
}


int TPZMHMHDivErrorEstimateMaterial::VariableIndex(const std::string &name)
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


int TPZMHMHDivErrorEstimateMaterial::NSolutionVariables(int var)
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

void TPZMHMHDivErrorEstimateMaterial::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    /**
     datavec[0] Flux mesh
     datavec[1] L2 mesh,
     datavec[2] Hdiv fem mesh
     datavec[3] L2 mesh fem
     **/
    
    TPZFNMatrix<9,REAL> PermTensor = fTensorK;
    TPZFNMatrix<9,REAL> InvPermTensor = fInvK;
    
    
    TPZManVector<STATE,2> pressexact(1,0.);
    TPZFNMatrix<9,STATE> gradu(3,1,0.), fluxinv(3,1);
    
    int IsH1position = IsH1Position(datavec);
    
    
    if(fForcingFunctionExact)
    {
        this->fForcingFunctionExact->Execute(datavec[IsH1position].x, pressexact,gradu);
       
    }
    
    PermTensor.Multiply(gradu, fluxinv);
    
    switch (var)
    {
        case 40://FluxFem
            for(int i=0; i<fDim; i++) Solout[i] = datavec[2].sol[0][i];
            break;
        case 41://FluxReconstructed
            
            if(IsH1position == 0){
                for (int i=0; i<fDim; i++) Solout[i] = datavec[0].sol[0][i];
            }
            else{
                
                //FluxReconstructed is K gradU
                TPZFMatrix<REAL> &dsolaxes = datavec[IsH1position].dsol[0];
                TPZFNMatrix<9,REAL> dsol(3,0);
                TPZFNMatrix<9,REAL> KGradsol(3,0);
                TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, dsol, datavec[IsH1position].axes);
                
                PermTensor.Multiply(dsol,KGradsol);
                
                for(int id=0 ; id<fDim; id++) {
                    Solout[id] = KGradsol(id,0);
                }
                
            }
            break;
        case 42:
            for(int i=0; i<fDim; i++) Solout[i] = -fluxinv(i);
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
