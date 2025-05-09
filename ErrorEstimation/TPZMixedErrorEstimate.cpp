//
//  TPZMixedErrorEstimate.cpp
//  ErrorEstimate
//
//  Created by Philippe Devloo on 20/04/18.
//

#include "TPZMixedErrorEstimate.h"
#include "Material/DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZMaterialDataT.h"

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
void TPZMixedErrorEstimate<MixedMat>::FillDataRequirements( TPZVec<TPZMaterialDataT<STATE> > &datavec) const
{
    MixedMat::FillDataRequirements(datavec);
    
    for( int i = 0; i<5; i++){
        datavec[i].SetAllRequirements(true);
        datavec[i].fNeedsSol = true;
    }
        
    datavec[0].fNeedsHSize=true;

}

template<class MixedMat>
void TPZMixedErrorEstimate<MixedMat>::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const
{
    MixedMat::FillBoundaryConditionDataRequirements(type,datavec);
    datavec[2].fNeedsSol = true;
}

/**
 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
 * @param datavec [in] stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 */
template<class MixedMat>
void TPZMixedErrorEstimate<MixedMat>::Contribute(const TPZVec<TPZMaterialDataT<STATE> > &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    if (datavec.size() != 5  )
    {
        DebugStop();
    }
    
    TPZFNMatrix<100,STATE> efkeep(ef);
    
    if(datavec[Epatch].fNeedsSol == false || datavec[Eorigin].fNeedsSol == false){
        DebugStop();
    }
    
#ifdef PZ_LOG
    
    
#endif
    
    //REAL solpatch = datavec[2].sol[0][0];
    //TPZFMatrix<STATE> &gradH1 = datavec[3].dsol[0];
    TPZFMatrix<REAL> &phip = datavec[Epressure].phi;
    int64_t phrp = phip.Rows();
    int64_t phrq = datavec[Eflux].fVecShapeIndex.NElements();
    TPZFMatrix<STATE> eks(phrp+phrq,phrp+phrq,0.),efs(phrp+phrq,1,0.);

    MixedMat::Contribute(datavec,weight,eks,efs);
//    {
//        std::stringstream sout;
//        sout<<"\n\n Matriz ek e vetor fk \n ";
//        ek.Print("ekmph = ",sout,EMathematicaInput);
//    }
    ef = efkeep;
    ek.AddSub(0, 0, eks);
    STATE force = 0.;
    if(MixedMat::fForcingFunction) {
        TPZManVector<STATE> res(1);
        MixedMat::fForcingFunction(datavec[Epressure].x,res);
        force = fSignConvention*res[0];
    }
    STATE perm = MixedMat::GetPermeability(datavec[Epressure].x);
    
    TPZFNMatrix<3,REAL> fluxprimal;

    STATE soloriginal = datavec[Eorigin].sol[0][0];
    REAL psival = datavec[Epatch].sol[0][0];
    {
        //REAL psival = datavec[Epatch].sol[0][0]; //datavec[1].sol[0][1]??
        //TPZFNMatrix<9,STATE> dsolprimal(3,1);
        TPZFNMatrix<9,REAL> gradpsi(3,1),gradpressure(3,1);
        TPZAxesTools<STATE>::Axes2XYZ(datavec[Epatch].dsol[0], gradpsi, datavec[Epatch].axes); //datavec[1].axes
        TPZAxesTools<STATE>::Axes2XYZ(datavec[Eorigin].dsol[0], gradpressure, datavec[Eorigin].axes); //datavec[1].axes

        //STATE soloriginal = datavec[Eorigin].sol[0][0];
        
        fluxprimal = perm*gradpressure;
        
        for (int64_t jq=0; jq<phrq; jq++)
        {
            // vector value
            TPZFNMatrix<3,REAL> jvec(3,1,0.);

            for(int id=0; id<3; id++){
                jvec(id,0) = datavec[Eflux].fDeformedDirections(id,jq);
            }
            STATE inner = 0.;
            for(int i=0; i<3; i++) {
                //inner += jvec(i,0)*perm*gradpressure(i,0);
                inner += jvec(i,0)*gradpressure(i,0);
            }
            ef(jq) -= weight*psival*inner;
        }
        
        STATE inner2 = 0.;
        for (int i=0; i<3; i++) {
            inner2 += perm*gradpressure(i,0)*gradpsi(i,0);//fluxprimal(i,0)*gradpsi(i,0);

        }
        for (int ip=0; ip<phrp; ip++) {//because the signal, get the opposite of equilibrated flux
            ef(phrq+ip,0) += (-1.)*weight*psival*force*phip(ip,0) + weight*inner2*phip(ip,0);
            //            ef(phrq+ip,0) += (-1.)*weight*phip(ip,0);
        }
    }
        
    int nactive = 0;
    for (int i=0; i<datavec.size(); i++) {
        if (datavec[i].fActiveApproxSpace) {
            nactive++;
        }
    }

    if(nactive!=3) DebugStop();
    
    if(1){//for internal patches we insert additional blocks
        for(int j=0; j<phrp; j++){
            ek(phrq+j,phrq+phrp) += phip(j,0)*weight;
            ek(phrq+phrp,phrq+j) += phip(j,0)*weight;
        }
        ef(phrq+phrp,0) += 0;
        //ef(phrq+phrp,0) += weight*psival*soloriginal;
    }
}

/// make a contribution to the error computation
template<class MixedMat>
void TPZMixedErrorEstimate<MixedMat>::Errors(const TPZVec<TPZMaterialDataT<STATE> > &data, TPZVec<REAL> &errors)
{
    /*
     data[0] Flux
     data[1] pressure
     data[3] H1
     */
    errors.Resize(MixedMat::NEvalErrors());
    errors.Fill(0.0);
    
    TPZManVector<STATE,3> flux(3,0.), pressure(1,0.);
    //flux = equilibrate flux recovered
    this->Solution(data,MixedMat::VariableIndex("Flux"), flux);
    //this->Solution(data,MixedMat::VariableIndex("Pressure"), pressure);

    REAL hsize = data[Eflux].HSize;
    TPZManVector<REAL> errorsaux;
    errorsaux.Fill(0.0);
    MixedMat::Errors(data,errorsaux);
    errors[1] += errorsaux[2];
    errors[3] = 2*hsize*hsize*errorsaux[2]/(M_PI*M_PI);//for residual error
    
    STATE perm = MixedMat::GetPermeability(data[1].x);
    
    TPZFNMatrix<3,REAL> fluxprimalneg;
    {
        TPZFNMatrix<3,STATE> dsolprimal(3,1);
        TPZFNMatrix<3,REAL> gradpressureH1(3,1);
        TPZAxesTools<STATE>::Axes2XYZ(data[Eorigin].dsol[0], dsolprimal, data[Eorigin].axes);
        
        for (int i=0; i<3; i++) {
            gradpressureH1(i,0) = dsolprimal(i,0);
            fluxprimalneg = perm*gradpressureH1;
        }
    }
    
    REAL inner = 0.;
    for (int i=0; i<3; i++) {
        inner += (flux[i]+fluxprimalneg(i,0))*(flux[i]+fluxprimalneg(i,0))/perm;
    }
    errors[2] += inner; // Flux estimator

    //make a contribution to the partial error on each color uniquely
    if(0){
        TPZManVector<STATE,1> u_exact(1, 0);
        TPZFMatrix<STATE> du_exact(3, 1, 0);
                
        if (this->fExactSol) {
            this->fExactSol(data[Eflux].x, u_exact, du_exact);
        }
        
//        std::cout << "u_exact: " << u_exact[0] << std::endl;
//        std::cout << "du_exact: " << du_exact[0] << "," << du_exact[1] << std::endl;

        REAL psival = data[Epatch].sol[0][0]; //hat function
        
        REAL inner2 = 0.;
        for (int i=0; i<3; i++) {
            inner2 += (psival*du_exact(i,0)+flux[i])*(psival*du_exact(i,0)+flux[i]);
        }
        
        errors[4] += inner2; // for L2 norm
    }
}

template<class MixedMat>
void TPZMixedErrorEstimate<MixedMat>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {

    int dim = MixedMat::Dimension();

    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    int phrq = phiQ.Rows();

    REAL v2 = bc.Val2()[0];
    REAL v1 = bc.Val1()(0, 0);
    REAL u_D = 0;
    REAL normflux = 0.;

    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9, STATE> gradu(3, 1);
        bc.ForcingFunctionBC()(datavec[0].x, res, gradu);

        const STATE perm = MixedMat::GetPermeability(datavec[0].x);

        for (int i = 0; i < 3; i++) {
            normflux += datavec[0].normal[i] * perm * gradu(i, 0);
        }

        if (bc.Type() == 0 || bc.Type() == 4) {
            v2 = res[0];
            u_D = res[0];
            normflux *= (-1.);
        } else if (bc.Type() == 1 || bc.Type() == 2) {
            v2 = -normflux;
            if (bc.Type() == 2) {
                v2 = -res[0] + v2 / v1;
            }
        } else if (bc.Type() == 5) {
            v2 = res[0];
        } else {
            DebugStop();
        }
    } else {
        v2 = bc.Val2()[0];
    }
    
    STATE hatval = datavec[2].sol[0][0];

    switch (bc.Type()) {
        case 0 :        // Dirichlet condition
            for (int iq = 0; iq < phrq; iq++) {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq, 0) += (-1.) * v2 *hatval * phiQ(iq, 0) * weight;
            }
            break;
            
        case 1 :            // Neumann condition
            for (int iq = 0; iq < phrq; iq++) {
                ef(iq, 0) += TPZMaterial::fBigNumber * v2 * hatval * phiQ(iq, 0) * weight;
                for (int jq = 0; jq < phrq; jq++) {
                    
                    ek(iq, jq) += TPZMaterial::fBigNumber * phiQ(iq, 0) * phiQ(jq, 0) * weight;
                }
            }
            break;
        default:
            DebugStop();
    }
}

template class TPZMixedErrorEstimate<TPZMixedDarcyFlow>;
