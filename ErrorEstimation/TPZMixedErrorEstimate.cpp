//
//  TPZMixedErrorEstimate.cpp
//  ErrorEstimate
//
//  Created by Philippe Devloo on 20/04/18.
//

#include "TPZMixedErrorEstimate.h"
#include "mixedpoisson.h"
#include "pzaxestools.h"


template<class MixedMat>
TPZMixedErrorEstimate<MixedMat>::TPZMixedErrorEstimate() : MixedMat()
{
    
}

template<class MixedMat>
TPZMixedErrorEstimate<MixedMat>::TPZMixedErrorEstimate(int matid, int dim) : MixedMat(matid,dim)
{
    
}

template<class MixedMat>
TPZMixedErrorEstimate<MixedMat>::~TPZMixedErrorEstimate()
{
    
}

template<class MixedMat>
TPZMixedErrorEstimate<MixedMat>::TPZMixedErrorEstimate(const TPZMixedErrorEstimate &cp) : MixedMat(cp), fSignConvention(cp.fSignConvention)
{
    
}

template<class MixedMat>
TPZMixedErrorEstimate<MixedMat> &TPZMixedErrorEstimate<MixedMat>::operator=(const TPZMixedErrorEstimate &copy)
{
    MixedMat::operator=(copy);
    fSignConvention = copy.fSignConvention;
    return *this;
}

template<class MixedMat>
void TPZMixedErrorEstimate<MixedMat>::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
    MixedMat::FillDataRequirements(datavec);
    {
        int i = 1;
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
    }
    
}

/**
 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
 * @param datavec [in] stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 */
template<class MixedMat>
void TPZMixedErrorEstimate<MixedMat>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    if (datavec[1].sol[0].size() != 3)
    {
        DebugStop();
    }
    TPZFNMatrix<100,STATE> efkeep(ef);
    MixedMat::Contribute(datavec,weight,ek,ef);
    ef = efkeep;
    TPZFMatrix<REAL> &phip = datavec[1].phi;
    int64_t phrp = phip.Rows();
    int64_t phrq = datavec[0].fVecShapeIndex.NElements();
    STATE force = MixedMat::ff;
    if(MixedMat::fForcingFunction) {
        TPZManVector<STATE> res(1);
        MixedMat::fForcingFunction->Execute(datavec[1].x,res);
        force = fSignConvention*res[0];
    }
    TPZFNMatrix<9,REAL> PermTensor = MixedMat::fTensorK;
    TPZFNMatrix<9,REAL> InvPermTensor = MixedMat::fInvK;
    //int rtens = 2*fDim;
    if(MixedMat::fPermeabilityFunction){
        PermTensor.Redim(3,3);
        InvPermTensor.Redim(3,3);
        TPZFNMatrix<3,STATE> resultMat;
        TPZManVector<STATE> res;
        MixedMat::fPermeabilityFunction->Execute(datavec[1].x,res,resultMat);
        for(int id=0; id<3; id++){
            for(int jd=0; jd<3; jd++){
                PermTensor(id,jd) = resultMat(id,jd);
                InvPermTensor(id,jd) = resultMat(id+3,jd);
            }
        }
    }
    

    TPZFNMatrix<3,REAL> fluxprimal;

    {
        REAL psival = datavec[1].sol[0][1];
        TPZFNMatrix<9,STATE> dsolprimal(3,3);
        TPZFNMatrix<9,REAL> gradpsi(3,1),gradpressure(3,1);
        TPZAxesTools<STATE>::Axes2XYZ(datavec[1].dsol[0], dsolprimal, datavec[1].axes);
        for (int i=0; i<3; i++) {
            gradpsi(i,0) = dsolprimal(i,1);
            gradpressure(i,0) = dsolprimal(i,2);
        }
        PermTensor.Multiply(gradpressure,fluxprimal);
        
        for (int64_t jq=0; jq<phrq; jq++)
        {
            // vector value
            TPZFNMatrix<3,REAL> jvec(3,1,0.);
            int jvecind = datavec[0].fVecShapeIndex[jq].first;
            int jshapeind = datavec[0].fVecShapeIndex[jq].second;
            REAL phival = datavec[0].phi(jshapeind);
            for(int id=0; id<3; id++){
                jvec(id,0) = datavec[0].fNormalVec(id,jvecind)*phival;
            }
            STATE inner = 0.;
            for(int i=0; i<3; i++) inner += jvec(i,0)*gradpressure(i,0);
            ef(jq) -= weight*inner*psival;
        }
        STATE inner2 = 0.;
        for (int i=0; i<3; i++) {
            inner2 += fluxprimal(i,0)*gradpsi(i,0);
        }
        for (int ip=0; ip<phrp; ip++) {
            ef(phrq+ip,0) += (-1.)*weight*force*phip(ip,0)*psival+weight*phip(ip,0)*inner2;
            //            ef(phrq+ip,0) += (-1.)*weight*phip(ip,0);
            
        }
    }

}

/// make a contribution to the error computation
template<class MixedMat>
void TPZMixedErrorEstimate<MixedMat>::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{
    errors.Resize(MixedMat::NEvalErrors());
    errors.Fill(0.0);
    TPZManVector<STATE,3> flux(3,0.), pressure(1,0.);
    this->Solution(data,MixedMat::VariableIndex("Flux"), flux);
    this->Solution(data,MixedMat::VariableIndex("Pressure"), pressure);

    TPZFNMatrix<9,REAL> PermTensor = MixedMat::fTensorK;
    TPZFNMatrix<9,REAL> InvPermTensor = MixedMat::fInvK;
    //int rtens = 2*fDim;
    if(MixedMat::fPermeabilityFunction){
        PermTensor.Redim(3,3);
        InvPermTensor.Redim(3,3);
        TPZFNMatrix<3,STATE> resultMat;
        TPZManVector<STATE> res;
        MixedMat::fPermeabilityFunction->Execute(data[1].x,res,resultMat);
        for(int id=0; id<3; id++){
            for(int jd=0; jd<3; jd++){
                PermTensor(id,jd) = resultMat(id,jd);
                InvPermTensor(id,jd) = resultMat(id+3,jd);
            }
        }
    }
    
    
    TPZFNMatrix<3,REAL> fluxprimalneg;
    
    {
        TPZFNMatrix<9,STATE> dsolprimal(3,3);
        TPZFNMatrix<9,REAL> gradpressure(3,1);
        TPZAxesTools<STATE>::Axes2XYZ(data[1].dsol[0], dsolprimal, data[1].axes);
        for (int i=0; i<3; i++) {
            gradpressure(i,0) = dsolprimal(i,2);
        }
        PermTensor.Multiply(gradpressure,fluxprimalneg);
    }
    REAL inner = 0.;
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            inner += (flux[i]+fluxprimalneg(i,0))*InvPermTensor(i,j)*(flux[j]+fluxprimalneg(j,0));
        }
    }
    errors[2] += inner;
}



template class TPZMixedErrorEstimate<TPZMixedPoisson>;
