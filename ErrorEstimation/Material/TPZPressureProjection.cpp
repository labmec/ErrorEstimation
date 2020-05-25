//
//  TPZPressureProjection.cpp
//  ErrorEstimationLib
//
//  Created by Denise De Siqueira on 21/05/20.
//

#include "TPZPressureProjection.h"
#include "TPZHDivErrorEstimateMaterial.h"
#include "TPZAnalyticSolution.h"
#include "pzbndcond.h"

TPZPressureProjection::TPZPressureProjection(int matid, int dim) : TPZHDivErrorEstimateMaterial(matid,dim)
{
    
}

TPZPressureProjection::TPZPressureProjection() : TPZHDivErrorEstimateMaterial()
{
    
}

TPZPressureProjection::TPZPressureProjection(const TPZHDivErrorEstimateMaterial &copy) : TPZHDivErrorEstimateMaterial(copy)
{
    
}


TPZPressureProjection::~TPZPressureProjection()
{
    
}

TPZPressureProjection &TPZPressureProjection::operator=(const TPZPressureProjection &copy)
{
    TPZHDivErrorEstimateMaterial::operator=(copy);
    return *this;
}


void TPZPressureProjection::ContributeBC(
    TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
    TPZFMatrix<STATE> &ef, TPZBndCond &bc) {

    /*
     Compute the L2 projection of potential into Robin boundary part
     ef= <InvKm(sigma.n+g)+u_d,phi_i>
     ek= <phi_i,phi_j>
     */
    int H1functionposition = FirstNonNullApproxSpaceIndex(datavec);
    int dim = datavec[H1functionposition].axes.Rows();

    TPZFMatrix<REAL> &phi_i = datavec[H1functionposition].phi;
    int nphi_i = phi_i.Rows();

    TPZManVector<REAL, 3> normal = datavec[2].normal;
    
    REAL normalsigma = datavec[2].sol[0][0];//sigma_h.n
 
    REAL u_D = 0.;
    REAL g = 0.;
    REAL normflux = 0.;
    
    
    if (bc.HasForcingFunction()) {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9, STATE> gradu(dim, 1);
        bc.ForcingFunction()->Execute(datavec[H1functionposition].x, res, gradu);
        TPZFNMatrix<9,REAL> PermTensor, InvPermTensor;
        GetPermeabilities(datavec[0].x, PermTensor, InvPermTensor);
        
        
        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                normflux += datavec[0].normal[i]*PermTensor(i,j)*gradu(j,0);
            }
        }
        
        g = (-1.)*normflux;
        u_D = res[0];
        
        

    }
    else {
        // usualmente updatebc coloca o valor exato no val2
        u_D = bc.Val2()(0, 0);
    }

    if (bc.Type() == 4 && !IsZero(bc.Val1()(0, 0))) {
        
        REAL Km = bc.Val1()(0, 0); // Km
        REAL InvKm = 1./Km;
        REAL robinterm = InvKm*(normalsigma+g)+u_D;

        for (int iq = 0; iq < nphi_i; iq++) {
            ef(iq, 0) += robinterm * phi_i(iq, 0) * weight;
            for (int jq = 0; jq < nphi_i; jq++) {
                ek(iq, jq) += weight *  phi_i(iq, 0) * phi_i(jq, 0);
            }
        }
    }
    else if (bc.Type() == 0 ||IsZero(bc.Val1()(0,0)) ) {
        
        for (int iq = 0; iq < nphi_i; iq++) {
            ef(iq, 0) += u_D * phi_i(iq, 0) * weight;
            for (int jq = 0; jq < nphi_i; jq++) {
                ek(iq, jq) += weight *  phi_i(iq, 0) * phi_i(jq, 0);
            }
        }
        
        
    }
    else{
        std::cout << "Projection is just for Dirichlet and Robin boundary conditions" << std::endl;
        
    }
}
