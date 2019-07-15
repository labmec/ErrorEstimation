//
//  TPZMatLaplacianHybrid.cpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 14/07/19.
//

#include "TPZMatLaplacianHybrid.h"
#include "pzbndcond.h"

TPZMatLaplacianHybrid::TPZMatLaplacianHybrid(int matid, int dim)
: TPZRegisterClassId(&TPZMatLaplacianHybrid::ClassId), TPZMatLaplacian(matid,dim)
{
    
}

TPZMatLaplacianHybrid::TPZMatLaplacianHybrid() :
TPZRegisterClassId(&TPZMatLaplacianHybrid::ClassId), TPZMatLaplacian()
{
    
}

TPZMatLaplacianHybrid::TPZMatLaplacianHybrid(const TPZMatLaplacian &copy) :
TPZRegisterClassId(&TPZMatLaplacianHybrid::ClassId), TPZMatLaplacian(copy)
{
    
}

TPZMatLaplacianHybrid::~TPZMatLaplacianHybrid()
{
    
}

TPZMatLaplacianHybrid &TPZMatLaplacianHybrid::operator=(const TPZMatLaplacianHybrid &copy)
{
    TPZMatLaplacian::operator=(copy);
    return *this;
}

TPZMaterial *TPZMatLaplacianHybrid::NewMaterial()
{
    return new TPZMatLaplacianHybrid(*this);
}

int TPZMatLaplacianHybrid::ClassId() const
{
    return Hash("TPZMatLaplacianHybrid") ^ TPZMatLaplacian::ClassId() << 1;
}


void TPZMatLaplacianHybrid::Write(TPZStream &buf, int withclassid) const
{
    TPZMatLaplacian::Write(buf,withclassid);
}

void TPZMatLaplacianHybrid::Read(TPZStream &buf, void *context)
{
    TPZMatLaplacian::Read(buf,context);
}


void TPZMatLaplacianHybrid::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    TPZFMatrix<REAL>  &phi = datavec[0].phi;
    TPZFMatrix<REAL> &dphi = datavec[0].dphix;
    TPZVec<REAL>  &x = datavec[0].x;
    //    TPZFMatrix<REAL> &axes = data.axes;
    //    TPZFMatrix<REAL> &jacinv = data.jacinv;
    int phr = phi.Rows();
    
    STATE fXfLoc = fXf;
    
    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        //TPZFMatrix<STATE> dres(Dimension(),1);
        //fForcingFunction->Execute(x,res,dres);       // dphi(i,j) = dphi_j/dxi
        fForcingFunction->Execute(x,res);
        fXfLoc = res[0];
    }
    
    STATE KPerm = fK;
    if (fPermeabilityFunction) {
        TPZFNMatrix<9,STATE> perm, invperm;
        TPZManVector<STATE,3> func;
        TPZFNMatrix<18,STATE> dfunc(6,3,0.);
        fPermeabilityFunction->Execute(x, func, dfunc);
        KPerm = dfunc(0,0);
    }
    
    //Equacao de Poisson
    for( int in = 0; in < phr; in++ ) {
        int kd;
        ef(in, 0) +=  (STATE)weight * fXfLoc * (STATE)phi(in,0);
        for( int jn = 0; jn < phr; jn++ ) {
            //ek(in,jn) += (STATE)weight*((STATE)(phi(in,0)*phi(jn,0)));
            for(kd=0; kd<fDim; kd++) {
                ek(in,jn) += (STATE)weight*(KPerm*(STATE)(dphi(kd,in)*dphi(kd,jn)));
            }
        }
    }
    for (int in =0; in < phr; in++) {
        ek(phr,in) += weight*phi(in,0);
        ek(in,phr) += weight*phi(in,0);
    }
    ek(phr,phr+1) -= weight;
    ek(phr+1,phr) -= weight;
    
    if (this->IsSymetric()){
        if ( !ek.VerifySymmetry(1.e-10) ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
}

void TPZMatLaplacianHybrid::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    TPZFMatrix<REAL>  &phi = datavec[0].phi;
    TPZFMatrix<REAL> &dphi = datavec[0].dphix;
    TPZVec<REAL>  &x = datavec[0].x;
    //    TPZFMatrix<REAL> &axes = data.axes;
    //    TPZFMatrix<REAL> &jacinv = data.jacinv;
    int phr = phi.Rows();
    
    STATE fXfLoc = fXf;
    
    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        //TPZFMatrix<STATE> dres(Dimension(),1);
        //fForcingFunction->Execute(x,res,dres);       // dphi(i,j) = dphi_j/dxi
        fForcingFunction->Execute(x,res);
        fXfLoc = res[0];
    }
    
    STATE KPerm = fK;
    if (fPermeabilityFunction) {
        TPZFNMatrix<9,STATE> perm, invperm;
        TPZManVector<STATE,3> func;
        TPZFNMatrix<18,STATE> dfunc(6,3,0.);
        fPermeabilityFunction->Execute(x, func, dfunc);
        KPerm = dfunc(0,0);
    }
    
    //Equacao de Poisson
    for( int in = 0; in < phr; in++ ) {
        int kd;
        ef(in, 0) +=  (STATE)weight * fXfLoc * (STATE)phi(in,0);
        for(kd=0; kd<fDim; kd++) {
            ef(in,0) -= (STATE)weight*(fK*(STATE)(dphi(kd,in)*datavec[0].dsol[0](kd,0)));
        }
    }
    ef(phr,0) += weight*(-datavec[0].sol[0][0]+ datavec[3].sol[0][0]);
    ef(phr+1,0) += weight*datavec[2].sol[0][0];

}

void TPZMatLaplacianHybrid::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    TPZFMatrix<REAL>  &phi_u = datavec[0].phi;
    TPZFMatrix<REAL>  &phi_flux = datavec[1].phi;
    //    TPZFMatrix<REAL> &axes = data.axes;
    int phr_primal = phi_u.Rows();
    int phr_hybrid = phi_flux.Rows();
    
    if(phr_hybrid)
    {
        std::cout << "Please implement me\n";
        DebugStop();
    }
    short in,jn;
    STATE v2[1];
    v2[0] = bc.Val2()(0,0);
    
    if(bc.HasForcingFunction()) {            // phi(in, 0) = phi_in                          // JORGE 2013 01 26
        TPZManVector<STATE> res(1);
        TPZFNMatrix<3,STATE> dres(3,1);
        bc.ForcingFunction()->Execute(datavec[0].x,res,dres);       // dphi(i,j) = dphi_j/dxi
        v2[0] = res[0];
    }
    
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
            for(in = 0 ; in < phr_primal; in++) {
                ef(in,0) += v2[0] * (STATE)(phi_u(in,0) * weight);
            }
            break;
        case 2 :        // mixed condition
            for(in = 0 ; in < phr_primal; in++) {
                ef(in, 0) += v2[0] * (STATE)(phi_u(in, 0) * weight);
                for (jn = 0 ; jn < phi_u.Rows(); jn++) {
                    ek(in,jn) += bc.Val1()(0,0) * (STATE)(phi_u(in,0) * phi_u(jn,0) * weight);     // peso de contorno => integral de contorno
                }
            }
            break;
    }

}

