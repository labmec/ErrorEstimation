//
//  TPZHDivErrorEstimateMaterial.cpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 11/06/19.
//

#include "TPZHDivErrorEstimateMaterial.h"

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
    
    if(!fNeumannLocalProblem){
        
        TPZMatPoisson3d::Contribute(datavec[0], weight, ek,ef);
    }
    
    
    //defining test functions
    // Setting the phis
    TPZFMatrix<REAL> &phiuk = datavec[0].phi;
    TPZFMatrix<REAL> &phirest = datavec[2].phi;// function of restriction term
    TPZFMatrix<REAL> &dphiuk = datavec[0].dphix;
    
   //  TPZFMatrix<REAL> &axes = datavec[0].axes;
    
    
    int nphiuk = phiuk.Rows();
    

    
    TPZFMatrix<STATE> solsigmafem(3,nphiuk),solukfem(1,1);
    //potetial fem
    solukfem(0,0) = datavec[3].sol[0][0];
    //flux fem
    for (int ip = 0; ip<3; ip++)
    {
        solsigmafem(ip,0) = datavec[2].sol[0][ip];
    }
    
   
    TPZFNMatrix<3,REAL> PermTensor = fTensorK;
    TPZFNMatrix<3,REAL> InvPermTensor = fInvK;
    
    

    TPZFMatrix<STATE> kgraduk(3,nphiuk);
    
        
        for(int irow=0 ; irow<nphiuk; irow++){
            
             //K graduk
            for(int i=0; i< dphiuk.Rows(); i++){

                for(int j=0; j< dphiuk.Rows();j++){

                    kgraduk(i,irow)+=InvPermTensor(i,j)*dphiuk(j,irow);

                }
                //bk=int_k sigmaukfem.grad phi_i
                
                ef(irow,0)+=weight*dphiuk(i,irow)*solsigmafem(i,irow);
            }
            
                //matrix Sk= int_{K} K graduk.gradv
                for(int jcol=0; jcol<nphiuk;jcol++){
                    for(int jd=0;  jd< dphiuk.Rows(); jd++){
                    ek(irow,jcol) +=weight*kgraduk(jd,irow)*dphiuk(jd,jcol);

                    }
                    
                }
                //Ck=int_{K} phi_i e Ck^t
                        ek(irow,nphiuk+1)+= weight*phirest(irow,0);
                        ek(nphiuk+1,irow)+= weight*phirest(irow,0);
            
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
