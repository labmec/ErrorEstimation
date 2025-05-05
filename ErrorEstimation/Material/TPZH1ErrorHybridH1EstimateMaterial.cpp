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
    
}

/**
 * @brief Returns an integer with the dimension of a post-processing variable
 * @param [in] var index of the post-processing variable, according to TPZDarcyFlow::VariableIndex method.
 */
int TPZH1ErrorHybridH1EstimateMaterial::NSolutionVariables(int var) const {
    
}

//! @name Error
/** @{*/
/*!
  \brief Calculates the error at a given point x.
  \param[in] datavec input data
  \param[out] errors The calculated errors.
 */
void TPZH1ErrorHybridH1EstimateMaterial::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data,
                                                TPZVec<REAL> &errors) {
    
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

#include "pzaxestools.h"

void TPZH1ErrorHybridH1EstimateMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                        REAL weight,TPZFMatrix<STATE> &ek,
                                                    TPZFMatrix<STATE> &ef) {
    std::cout << "Favor me implementar\n";
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

    int nphi = phi.Rows();
    
    STATE fXfLoc = 0;

    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        fForcingFunction(x, res);
        fXfLoc = res[0];
    }
    STATE KPerm = GetPermeability(datavec[0].x);
    ek.AddContribution(0, 0, dphix, 1, dphix, 0, weight*KPerm*hatval);
    STATE gradhatgradorig = 0.;
    for(int i=0; i<3; i++) gradhatgradorig += dhatvalxy(i,0)*dpressoriginxy(i,0);

    TPZFNMatrix<1> forcemat(1, 1, fXfLoc*hatval-gradhatgradorig);
    ef.AddContribution(0,0, phi, 0, forcemat, 0, weight);
    
    TPZFMatrix<STATE> &phiconst = datavec[Eaveragepressure].phi;
    if(phiconst.Rows() != 1) DebugStop();
    ek.AddContribution(nphi, 0, phiconst, 0, phi, 1, weight);
    ek.AddContribution(0, nphi, phi, 0, phiconst, 0, weight);
}
/**@}*/

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
    std::cout << "Favor me implementar BC\n";
    int dim = Dimension();

    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    int64_t phrq = phiQ.Rows();

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

    switch (bc.Type()) {
        case 0 :        // Dirichlet condition
            for (int iq = 0; iq < phrq; iq++) {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq, 0) += (-1.) * v2 * phiQ(iq, 0) * weight;
            }
            break;
        default:
            DebugStop();
    }
}
