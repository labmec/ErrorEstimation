//
//  TPZH1ErrorHybridH1EstimateMaterial.cpp
//  HybridH1vsMixed
//
//  Created by Philippe Devloo on 01/05/25.
//

#include "TPZH1ErrorHybridH1EstimateMaterial.h"

/**
 * @brief Returns an integer associated with a post-processing variable name
 * @param [in] name string containing the name of the post-processing variable. Ex: "Pressure".
 */
int TPZH1ErrorHybridH1EstimateMaterial::VariableIndex(const std::string &name) const {
    return TPZHybridDarcyFlow::VariableIndex(name);
}

/**
 * @brief Returns an integer with the dimension of a post-processing variable
 * @param [in] var index of the post-processing variable, according to TPZDarcyFlow::VariableIndex method.
 */
int TPZH1ErrorHybridH1EstimateMaterial::NSolutionVariables(int var) const {
    return TPZHybridDarcyFlow::TPZDarcyFlow::NSolutionVariables(var);
}
#include "pzaxestools.h"

//! @name Error
/** @{*/
/*!
  \brief Calculates the error at a given point x.
  \param[in] datavec input data
  \param[out] errors The calculated errors.
 */
void TPZH1ErrorHybridH1EstimateMaterial::Errors(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                                                TPZVec<REAL> &errors) {
    STATE KPerm = GetPermeability(datavec[0].x);
    STATE fXfLoc = 0;

    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        fForcingFunction(datavec[0].x, res);
        fXfLoc = res[0];
    }
    TPZManVector<STATE,3> exsol(1,0.);
    TPZFNMatrix<3,STATE> dsol(2, 1,0.);
    if(HasExactSol()) {
        this->ExactSol()(datavec[3].x,exsol,dsol);
    }
    errors.Fill(0.);
    errors[1] += (fXfLoc-datavec[0].sol[0][0])*(fXfLoc-datavec[0].sol[0][0]);
    errors[3] = errors[1]*datavec[0].HSize*datavec[0].HSize*2/(M_PI*M_PI);
    
    TPZFNMatrix<3,STATE> GraduH1xy(3,1,0.),GraduHybridxy(3,1,0.);
    TPZFMatrix<STATE> &dpressorigin = datavec[3].dsol[0];
    TPZFMatrix<STATE> &dH1Hybrid = datavec[1].dsol[0];

    TPZAxesTools<STATE>::Axes2XYZ(dpressorigin, GraduH1xy, datavec[3].axes);
    TPZAxesTools<STATE>::Axes2XYZ(dH1Hybrid, GraduHybridxy, datavec[1].axes);
    
    for (int d=0; d<2; d++) {
        errors[0] += (GraduH1xy(d,0)-dsol(d,0))*(GraduH1xy(d,0)-dsol(d,0));
        errors[4] += (GraduHybridxy(d,0)-dsol(d,0))*(GraduHybridxy(d,0)-dsol(d,0));
        errors[2] += (GraduH1xy(d,0)-GraduHybridxy(d,0))*(GraduH1xy(d,0)-GraduHybridxy(d,0));
    }
    errors[2] *= KPerm;
}

/**
 * @brief Returns an unique class identifier
 */
int TPZH1ErrorHybridH1EstimateMaterial::ClassId() const {
    return Hash("TPZH1ErrorHybridH1EstimateMaterial") ^ TPZHybridDarcyFlow::ClassId() << 1;

}


void TPZH1ErrorHybridH1EstimateMaterial::FillDataRequirements( TPZVec<TPZMaterialDataT<STATE> > &datavec) const
{
    TPZMatCombinedSpacesT<STATE>::FillDataRequirements(datavec);
    
    for( int i = 0; i<5; i++){
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
    }        
    datavec[0].fNeedsHSize=true;

}


/** @name Contribute */
/** @{ */
/**
 * @brief It computes a contribution to the stiffness matrix
 * and load vector at one integration point.
 * @param[in] datavec stores all input data
 * @param[in] weight is the weight of the integration rule
 * @param[out] ek is the element matrix
 * @param[out] ef is the rhs vector
 */


void TPZH1ErrorHybridH1EstimateMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                        REAL weight,TPZFMatrix<STATE> &ek,
                                                    TPZFMatrix<STATE> &ef) {
    TPZFMatrix<REAL> &phires = datavec[0].phi;
    TPZFMatrix<REAL> &dphi = datavec[1].dphix;
    TPZFMatrix<REAL> &phi = datavec[1].phi;
    TPZFNMatrix<40> dphix;
    TPZAxesTools<STATE>::Axes2XYZ(dphi, dphix, datavec[1].axes);
    TPZVec<REAL>  &x = datavec[1].x;
    STATE hatval = datavec[2].sol[0][0];
    TPZFMatrix<STATE> &dhatval = datavec[2].dsol[0];
    TPZFNMatrix<3> dhatvalxy(3,1);
    TPZAxesTools<STATE>::Axes2XYZ(dhatval, dhatvalxy, datavec[2].axes);
    STATE pressorigin = datavec[3].sol[0][0];
    TPZFMatrix<STATE> &dpressorigin = datavec[3].dsol[0];
    TPZFNMatrix<3> dpressoriginxy(3, 1);
    TPZAxesTools<STATE>::Axes2XYZ(dpressorigin, dpressoriginxy, datavec[3].axes);

    int nphiRes = phires.Rows();
    int nphiH1 = phi.Rows();
    
    if(ek.Rows() != nphiRes+nphiH1+1) DebugStop();
    
    STATE fXfLoc = 0;

    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        fForcingFunction(x, res);
        fXfLoc = res[0];
    }
    // contribution of the residual error estimator
    ek.AddContribution(0,0, phires, 0, phires, 1, weight);
    TPZFMatrix<REAL> onebyone(1,1,1.);
    ef.AddContribution(0, 0, phires, 0, onebyone, 0, fXfLoc*weight*hatval);
    
    // contribution of the hybrid H1 contribution
    STATE KPerm = GetPermeability(datavec[0].x);
    ek.AddContribution(nphiRes, nphiRes, dphix, 1, dphix, 0, weight*KPerm);
    STATE gradhatgradorig = 0.;
    for(int i=0; i<3; i++) gradhatgradorig += KPerm*dhatvalxy(i,0)*dpressoriginxy(i,0);

    TPZFNMatrix<1> forcemat(1, 1, fXfLoc*hatval-gradhatgradorig);
    ef.AddContribution(nphiRes,0, phi, 0, forcemat, 0, weight);
    
    // computing KPerm \nabla psi_a u
    TPZFNMatrix<3,STATE> uhAgradPsi(3, 1,0.);
    for(int d=0; d<3; d++) uhAgradPsi(d,0) = KPerm*pressorigin*dhatvalxy(d,0);
    ef.AddContribution(nphiRes, 0, dphix, 1, uhAgradPsi, 0, weight);
    
    // contribution of the lagrange multiplier
    TPZFMatrix<STATE> &phiconst = datavec[Eaveragepressure].phi;
    if(phiconst.Rows() != 1) DebugStop();
    ek.AddContribution(nphiRes+nphiH1, nphiRes, phiconst, 0, phi, 1, weight);
    ek.AddContribution(nphiRes, nphiH1+nphiRes, phi, 0, phiconst, 0, weight);
    ef(nphiRes+nphiH1,0) += weight*hatval*pressorigin;
}
/**@}*/

void TPZH1ErrorHybridH1EstimateMaterial::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const
{
    //TPZHybridDarcyFlow::FillBoundaryConditionDataRequirements(type,datavec);
    datavec[2].fNeedsSol = true;
}

/** @name ContributeBC
    @ingroup Contribute*/
/**@{*/
/**
 * @brief It computes a contribution to the stiffness matrix
 * and load vector at one BC integration point.
 * @param[in] datavec stores all input data
 * @param[in] weight is the weight of the integration rule
 * @param[out] ek is the element matrix
 * @param[out] ef is the rhs vector
 * @param[in] bc is the boundary condition material
 */
void TPZH1ErrorHybridH1EstimateMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                          REAL weight, TPZFMatrix<STATE> &ek,
                          TPZFMatrix<STATE> &ef,
                                                      TPZBndCondT<STATE> &bc) {
    int dim = Dimension();

    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    int64_t phrq = phiQ.Rows();
    STATE hatval = datavec[2].sol[0][0];

    REAL v2 = bc.Val2()[0];
    REAL v1 = bc.Val1()(0, 0);
    REAL u_D = 0;
    REAL normflux = 0.;

    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9, STATE> gradu(3, 1);
        bc.ForcingFunctionBC()(datavec[0].x, res, gradu);

        const STATE perm = GetPermeability(datavec[0].x);

        for (int i = 0; i < 3; i++) {
            normflux += datavec[0].normal[i] * perm * gradu(i, 0);
        }

        if (bc.Type() == 0) {
            v2 = res[0];
            u_D = res[0];
            normflux *= (-1.);
        } else if (bc.Type() == 1) {
            v2 = -normflux;
        } else {
            DebugStop();
        }
    } else {
        v2 = bc.Val2()[0];
    }

    switch (bc.Type()) {
        case 0 :        // Dirichlet condition
            for (int iq = 0; iq < phrq; iq++) {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq, 0) += (-1.) * v2 * phiQ(iq, 0) * hatval * weight;
            }
            break;
        default:
            DebugStop();
    }
}
