//
//  SteklovProblem
//
//  Created by Denise De Siqueira on 27/11/20.
//

#include "TPZSteklovMaterial.h"
#include "pzaxestools.h"
#include "pzbndcond.h"


#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("TPZSteklovMaterial"));
#endif

TPZSteklovMaterial::TPZSteklovMaterial(int matid, int dim) : TPZMixedPoisson(matid,dim)
{
    
}

TPZSteklovMaterial::TPZSteklovMaterial() : TPZMixedPoisson()
{
    
}

TPZSteklovMaterial::TPZSteklovMaterial(const TPZMixedPoisson &copy) : TPZMixedPoisson(copy)
{
    
}

TPZSteklovMaterial::TPZSteklovMaterial(const TPZSteklovMaterial &copy) : TPZMixedPoisson(copy)
{
    
}

TPZSteklovMaterial::~TPZSteklovMaterial()
{
    
}

TPZSteklovMaterial &TPZSteklovMaterial::operator=(const TPZSteklovMaterial &copy)
{
    TPZMixedPoisson::operator=(copy);
    return *this;
}

void TPZSteklovMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    
    ef.Resize(ek.Rows(), ek.Cols());
    ef.Zero();
    
    STATE force = 0.;

    
    TPZFNMatrix<9,STATE> PermTensor;
    TPZFNMatrix<9,STATE> InvPermTensor;
    
    GetPermeabilities(datavec[1].x, PermTensor, InvPermTensor);
    
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
            ivec(id,0) = datavec[0].fDeformedDirections(id,ivecind);
        }
        
        
        TPZFNMatrix<3,REAL> ivecZ(3,1,0.);
        TPZFNMatrix<3,REAL> jvecZ(3,1,0.);
        for (int jq=0; jq<phrq; jq++)
        {
            TPZFNMatrix<3,REAL> jvec(3,1,0.);
            int jvecind = datavec[0].fVecShapeIndex[jq].first;
            int jshapeind = datavec[0].fVecShapeIndex[jq].second;
            
            for(int id=0; id<3; id++){
                jvec(id,0) = datavec[0].fDeformedDirections(id,jvecind);
            }
            
            //dot product between Kinv[u]v
            jvecZ.Zero();
            for(int id=0; id<3; id++){
                for(int jd=0; jd<3; jd++){
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
            ivec(id,0) = datavec[0].fDeformedDirections(id,ivecind);
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
//    for(int ip=0; ip<phrp; ip++){
//        ef(phrq+ip,0) += (-1.)*weight*force*phip(ip,0);
//    }
    

    #ifdef LOG4CXX
        if(logdata->isDebugEnabled())
        {
            std::stringstream sout;
            sout<<"\n\n Matriz ek e vetor fk \n ";
            ek.Print("ekmph = ",sout,EMathematicaInput);
            ef.Print("efmph = ",sout,EMathematicaInput);
            LOGPZ_DEBUG(logdata,sout.str());
        }
    #endif
    
}

void TPZSteklovMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    ef.Resize(ek.Rows(), ek.Cols());
    ef.Zero();
    
    
    int dim = Dimension();
    
    TPZFMatrix<REAL>  &phiQ = datavec[0].phi;
    int phrq = phiQ.Rows();

    REAL v2 = bc.Val2()(0,0);
    REAL v1 = bc.Val1()(0,0);
    REAL u_D = 0;
    REAL normflux = 0.;
    
    if(bc.HasForcingFunction())
    {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9,STATE> gradu(dim,1);
        bc.ForcingFunction()->Execute(datavec[0].x,res,gradu);
        TPZFNMatrix<9,STATE> PermTensor, InvPermTensor;
        GetPermeabilities(datavec[0].x, PermTensor, InvPermTensor);
        
        
        for(int i=0; i<3; i++)
        {
            for(int j=0; j<dim; j++)
            {
                normflux += datavec[0].normal[i]*PermTensor(i,j)*gradu(j,0);
            }
        }
        
        
        if(bc.Type() == 0)
        {
            v2 = res[0];
            u_D = res[0];
            normflux *= (-1.);
        }
        else if(bc.Type() == 1)
        {
            v2 = -normflux;

        }
        else
        {
            DebugStop();
        }
    }else
    {
        v2 = bc.Val2()(0,0);
    }

    switch (bc.Type()) {
        case 0 :        // Dirichlet condition
            //primeira equacao
            for(int iq=0; iq<phrq; iq++)
            {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq,0) += (-1.)*v2*phiQ(iq,0)*weight;
            }
            break;
            
        case 1 :            // Neumann condition
            //primeira equacao
            for(int iq=0; iq<phrq; iq++)
            {
                
                for (int jq=0; jq<phrq; jq++) {
                    ef(iq,jq)+= gBigNumber*phiQ(iq,0)*phiQ(jq,0)*weight;
                    
                    ek(iq,jq)+= gBigNumber*phiQ(iq,0)*phiQ(jq,0)*weight;
                }
            }
            break;
        
    }
    
}
