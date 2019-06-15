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
    
}

