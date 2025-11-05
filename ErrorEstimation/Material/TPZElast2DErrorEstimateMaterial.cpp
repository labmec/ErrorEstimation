//
//  TPZElast2DErrorEstimateMaterial.cpp
//  HybridH1vsMixed
//
//  Created by Philippe Devloo on 01/05/25.
//

#include "TPZElast2DErrorEstimateMaterial.h"

/**
 * @brief Returns an integer associated with a post-processing variable name
 * @param [in] name string containing the name of the post-processing variable. Ex: "Pressure".
 */
int TPZElast2DErrorEstimateMaterial::VariableIndex(const std::string &name) const {
    if(name == "Partition") return 31;
    if(name == "DistFlux") return 32;
    if(name == "SolH1Hat") return 33;
    if(name == "H1_Error") return 100;
    if(name == "Hybrid_Error") return 101;
    if(name == "Estimated_Error") return 102;
    if(name == "EffectivityIndex") return 103;
    return TPZHybridElasticity2D::VariableIndex(name);
}

/**
 * @brief Returns an integer with the dimension of a post-processing variable
 * @param [in] var index of the post-processing variable, according to TPZDarcyFlow::VariableIndex method.
 */
int TPZElast2DErrorEstimateMaterial::NSolutionVariables(int var) const {
    if(var == 31) return 1;
    if(var == 32) return 1;
    if(var == 33) return 2;
    if(var == 100) return 1;
    if(var == 101) return 1;
    if(var == 102) return 1;
    if(var == 103) return 1;

    return TPZHybridElasticity2D::NSolutionVariables(var);
}

/** @brief Returns the solution associated with a given index
    based on the finite element approximation.
    @param[in] datavec Stores all the input data.
    @param[in] var Index of the queried solution
    @param[out] sol FEM Solution at the integration point
*/
void TPZElast2DErrorEstimateMaterial::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                      int var, TPZVec<STATE> &sol) {
    if(var < 15) {
        TPZElasticity2D::Solution(datavec[1],var,sol);
    }
    // hat function
    if(var == 31) sol = datavec[2].sol[0][0];
    if(var == 32) sol = datavec[4].sol[0][0];
    TPZManVector<STATE,3> solhat(datavec[Eorigin].sol[0]);
    for(int i=0; i<2; i++) solhat[i] *= datavec[2].sol[0][0];
    if(var == 33) sol =  solhat;
}


#include "pzaxestools.h"

//! @name Error
/** @{*/
/*!
  \brief Calculates the error at a given point x.
  \param[in] datavec input data
  \param[out] errors The calculated errors.
 */
void TPZElast2DErrorEstimateMaterial::Errors(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                                                TPZVec<REAL> &errors) {
    TPZManVector<STATE,2> fXfLoc(2,0.);

    if(fForcingFunction) {
        fForcingFunction(datavec[0].x, fXfLoc);
    }
    TPZManVector<STATE,3> solexact(2,0.);
    TPZFNMatrix<3,STATE> dsolexact(2, 2,0.);
    if(HasExactSol()) {
        this->ExactSol()(datavec[0].x,solexact,dsolexact);
    }
    errors.Fill(0.);
    TPZFNMatrix<9,STATE> D(3, 3);
    ComputeDMatrix(fE_def, fnu_def, D);
    
    TPZFNMatrix<3,STATE> GraduH1xy(3,1,0.),GraduHybridxy(3,1,0.);
    TPZFMatrix<STATE> &dH1 = datavec[Eorigin].dsol[0];
    TPZFMatrix<STATE> &dH1Hybrid = datavec[Epressure].dsol[0];

    TPZAxesTools<STATE>::Axes2XYZ(dH1, GraduH1xy, datavec[Eorigin].axes);
    TPZAxesTools<STATE>::Axes2XYZ(dH1Hybrid, GraduHybridxy, datavec[Epressure].axes);
    TPZFNMatrix<3,STATE> dsolexactV(3, 1,0.), GraduH1xyV(3,1,0.),GraduHybridxyV(3,1,0.);
    
    for(int i=0; i<2; i++) {
        dsolexactV(i) = dsolexact(i,i);
        GraduH1xyV(i) = GraduH1xy(i,i);
        GraduHybridxyV(i) = GraduHybridxy(i,i);
        dsolexactV(2) += dsolexact(i,1-i);
        GraduH1xyV(2) += GraduH1xy(i,1-i);
        GraduHybridxyV(2) += GraduHybridxy(i,1-i);
    }
    TPZFNMatrix<3,STATE> sigmaexactV(3,1),sigmaH1(3,1),sigmaHybrid(3,1);
    D.Multiply(dsolexactV, sigmaexactV);
    D.Multiply(GraduH1xyV, sigmaH1);
    D.Multiply(GraduHybridxyV, sigmaHybrid);
    for (int d=0; d<3; d++) {
        errors[0] += (GraduH1xyV(d,0)-dsolexactV(d,0))*(sigmaH1(d,0)-sigmaexactV(d,0));
        errors[1] += (GraduHybridxyV(d,0)-dsolexactV(d,0))*(sigmaHybrid(d,0)-sigmaexactV(d,0));
        errors[2] += (GraduH1xyV(d,0)-GraduHybridxyV(d,0))*(sigmaH1(d,0)-sigmaHybrid(d,0));
    }
}

/**
 * @brief Returns an unique class identifier
 */
int TPZElast2DErrorEstimateMaterial::ClassId() const {
    return Hash("TPZElast2DErrorEstimateMaterial") ^ TPZHybridElasticity2D::ClassId() << 1;

}


void TPZElast2DErrorEstimateMaterial::FillDataRequirements( TPZVec<TPZMaterialDataT<STATE> > &datavec) const
{
    TPZMatCombinedSpacesT<STATE>::FillDataRequirements(datavec);
    
    for( int i = 0; i<5; i++){
        datavec[i].SetAllRequirements(false);
//        datavec[i].fNeedsSol = true;
    }
    datavec[Epatch].fNeedsSol = true;
    datavec[Eorigin].fNeedsSol = true;
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


void TPZElast2DErrorEstimateMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                        REAL weight,TPZFMatrix<STATE> &ek,
                                                    TPZFMatrix<STATE> &ef) {
    TPZFMatrix<REAL> &phires = datavec[Eflux].phi;
    TPZFMatrix<REAL> &dphi = datavec[Epressure].dphix;
    TPZFMatrix<REAL> &phi = datavec[Epressure].phi;
    TPZFNMatrix<40> dphix;
    TPZAxesTools<STATE>::Axes2XYZ(dphi, dphix, datavec[1].axes);
    TPZVec<REAL>  &x = datavec[Epressure].x;
    STATE hatval = datavec[Epatch].sol[0][0];
    TPZFMatrix<STATE> &dhatval = datavec[Epatch].dsol[0];
    TPZFNMatrix<3> dhatvalxy(3,1);
    TPZAxesTools<STATE>::Axes2XYZ(dhatval, dhatvalxy, datavec[2].axes);
    
    auto &displacement = datavec[Eorigin].sol[0];
    TPZFMatrix<STATE> &graddisp = datavec[Eorigin].dsol[0];
    TPZFNMatrix<3,STATE> graddispxy(3, 1);
    TPZAxesTools<STATE>::Axes2XYZ(graddisp, graddispxy, datavec[Eorigin].axes);
    TPZFMatrix<STATE> &dphixyz = dphix;
    
    TPZFNMatrix<3,STATE> exactderiv(3,2);
    TPZManVector<STATE,1> exactsol(2,0);
    if(fExactSol) {
        fExactSol(x,exactsol,exactderiv);
    }
    TPZFNMatrix<9,REAL> D(3, 3), solvoight(3,1),sigmavoight;
    ComputeDMatrix(fE_def, fnu_def, D);
    solvoight(0,0) = graddisp(0,0);
    solvoight(1,0) = graddisp(1,1);
    solvoight(2,0) = graddisp(1,0)+graddisp(0,1);
    D.Multiply(solvoight, sigmavoight);

    
    int64_t nphiRes = phires.Rows();
    int64_t nphiH1 = phi.Rows();
    
    if(ek.Rows() != 2*nphiRes+2*nphiH1+3) DebugStop();
    
    TPZFMatrix<STATE> dphixyV(3,2*nphiH1,0),DB(3,2*nphiH1,0);
    // contribution of the hybrid H1 contribution
    for(int i=0; i<nphiH1; i++) {
        dphixyV(0,2*i) = dphix(0,i);
        dphixyV(1,2*i+1) = dphix(1,i);
        dphixyV(2,2*i) = dphix(1,i);
        dphixyV(2,2*i+1) = dphix(0,i);
    }
    D.Multiply(dphixyV, DB);
    ek.AddContribution(nphiRes*2, nphiRes*2, dphix, 1, DB, 0, weight);
    
    // contribution of the residual error estimator
    TPZFNMatrix<60,REAL> phires2(2*nphiRes,2,0.);
    for(int i=0; i<nphiRes; i++) {
        phires2(2*i,0) = phires(i,0);
        phires2(2*i+1,1) = phires(i,0);
    }
    ek.AddContribution(0,0, phires2, 0, phires2, 1, weight);

    // contribution of the lagrange multiplier
    TPZManVector<REAL,3> xcenter = datavec[Eaveragepressure].XCenter;
    {
        int64_t phr = nphiH1;
        int64_t fi = 2*nphiRes;
        for (int in =0; in < nphiH1; in++) {
            ek(fi+2*phr,fi+2*in) += weight*phi(in,0);//lambda*phi
            ek(fi+2*phr+1,fi+2*in+1) = weight*phi(in,0);
            ek(fi+2*phr+2,fi+2*in) += -weight*(x[1]-xcenter[1])*phi(in,0);
            ek(fi+2*phr+2,fi+2*in+1) += weight*(x[0]-xcenter[0])*phi(in,0);
            ek(fi+2*in,fi+2*phr) += weight*phi(in,0);//lambda*phi
            ek(fi+2*in+1,fi+2*phr+1) = weight*phi(in,0);
            ek(fi+2*in,fi+2*phr+2) += -weight*(x[1]-xcenter[1])*phi(in,0);
            ek(fi+2*in+1,fi+2*phr+2) += weight*(x[0]-xcenter[0])*phi(in,0);
        }
//        TPZFNMatrix<9,REAL> rigcouple(3, 3,0.);
//        rigcouple(0,0) = weight;
//        rigcouple(1,1) = weight;
//        rigcouple(0,2) = -(x[1]-xcenter[1])*weight;
//        rigcouple(1,2) = (x[0]-xcenter[0])*weight;
//        rigcouple(2,0) = rigcouple(0,2);
//        rigcouple(2,1) = rigcouple(1,2);
//        rigcouple(2,2) = (x[0]-xcenter[0])*(x[0]-xcenter[0])+(x[1]-xcenter[1])*(x[1]-xcenter[1])*weight;
//        for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
//            ek(fi+2*phr+i,fi+2*phr+3+j) += -rigcouple(i,j);
//            ek(fi+2*phr+3+i,fi+2*phr+j) += -rigcouple(j,i);
//        }
    }

    TPZFMatrix<STATE> &phiconst = datavec[Eaveragepressure].phi;
    if(phiconst.Rows() != 1) DebugStop();
    ek.AddContribution(nphiRes+nphiH1, nphiRes, phiconst, 0, phi, 1, weight);
    ek.AddContribution(nphiRes, nphiH1+nphiRes, phi, 0, phiconst, 0, weight);

    TPZManVector<STATE,3> XfLoc(2,0.);

    if(fForcingFunction) {          ;
        fForcingFunction(x, XfLoc);
    }
    TPZFMatrix<REAL> twobyone(2,1,1.);
    for(int i=0; i<2; i++) twobyone(i) = XfLoc[i];
    ef.AddContribution(0, 0, phires2, 0, twobyone, 0, weight*hatval);
    
    // compute sigma grad psia *
    TPZManVector<STATE,3> sigmagradhat(2,0.);
    sigmagradhat[0] = sigmavoight[0]*dhatvalxy[0]+sigmavoight[2]*dhatvalxy[1];
    sigmagradhat[1] = sigmavoight[2]*dhatvalxy[0]+sigmavoight[1]*dhatvalxy[1];

//    std::cout << "x " << x << " gradhatgradorig " << gradhatgradorig << std::endl;
    extern bool Print;
    if(Print) {
        std::cout << "x " << x << " dhatvalxy ";
        for(int i=0; i<2; i++) std::cout << dhatvalxy(i,0) << " ";
//        std::cout << "\nx " << x << " dhatval ";
//        for(int i=0; i<2; i++) std::cout << dhatval(i,0) << " ";
        std::cout << std::endl;
    }

    TPZFNMatrix<2, STATE> forcemat(2, 1);
    TPZFNMatrix<60,REAL> phi2(2*nphiH1,2,0.);
    for(int i=0; i<nphiH1; i++) {
        phi2(2*i,0) = phi(i);
        phi2(2*i+1,1) = phi(i);
    }
    for(int i=0; i<2; i++) forcemat(i,0) = (XfLoc[i]*hatval-sigmagradhat[i]);
    ef.AddContribution(2*nphiRes ,0, phi2, 0, forcemat, 0, weight);
    
    extern std::complex<STATE> integrateF;
//    integrateF += weight;
//    integrateF += gradhatgradorig*weight;
    integrateF += XfLoc[0]*hatval*weight;
//    integrateF += exactderiv(1,0)*weight;
    int dim = dphix.Rows();
    // computing KPerm \nabla psi_a u
    TPZFNMatrix<3,STATE> BPsia(3, 1,0.);
    BPsia(0,0) = dhatvalxy[0]*displacement[0];
    BPsia(1,0) = dhatvalxy[1]*displacement[1];
    BPsia(2,0) = dhatvalxy[1]*displacement[0]+dhatvalxy[0]*displacement[1];
    TPZFNMatrix<3,STATE> sigpsia(3,1);
    D.Multiply(BPsia, sigpsia);
//    for(int d=0; d<3; d++) uhAgradPsi(d,0) = KPerm*exactsol[0]*dhatvalxy(d,0);
//    for(int d=0; d<3; d++) uhAgradPsi(d,0) = KPerm*pressorigin*dhatvalxy(d,0);
    ef.AddContribution(2*nphiRes, 0, dphixyV, 1, sigpsia, 0, weight);
//    std::cout << "Norm dphixyz " << Norm(dphixyz) << " Norm uhAgradPsi " << Norm(uhAgradPsi) << std::endl;
    int64_t first = 2*nphiRes+2*nphiH1;
    ef(first) += weight*hatval*displacement[0];
    ef(first+1) += weight*hatval*displacement[1];
    ef(first+2) += weight*hatval*(-(x[1]-xcenter[1])*displacement[0]+(x[0]-xcenter[0])*displacement[1]);
}

/**
 * @brief It computes a contribution to the stiffness matrix
 * and load vector at one integration point.
 * @param[in] datavec stores all input data
 * @param[in] weight is the weight of the integration rule
 * @param[out] ek is the element matrix
 * @param[out] ef is the rhs vector
 */
void TPZElast2DErrorEstimateMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                        const TPZFMatrix<CSTATE> &phi, const TPZFMatrix<CSTATE> &dphix,
                        REAL weight,
                                                    TPZFMatrix<CSTATE> &ef, TPZFMatrix<CSTATE> &RBM) {
    int64_t nphiH1 = phi.Rows();
    if(nphiH1 == 0) return;
#ifdef PZDEBUG
    if(phi.Cols() != 2) DebugStop();
    if(dphix.Rows() != 4) DebugStop();
    if(dphix.Cols() != nphiH1) DebugStop();
#endif
    int nstate = NStateVariables();
    if(nstate != 2) DebugStop();
    TPZVec<REAL>  &x = datavec[Epatch].x;
    // contribution of the hybrid H1 contribution
    STATE hatval = datavec[Epatch].sol[0][0];
    TPZFMatrix<STATE> &dhatval = datavec[Epatch].dsol[0];
    TPZFNMatrix<3> dhatvalxy(3,1);
    TPZAxesTools<STATE>::Axes2XYZ(dhatval, dhatvalxy, datavec[2].axes);
    auto &displacement = datavec[Eorigin].sol[0];
    TPZFMatrix<STATE> &graddisp = datavec[Eorigin].dsol[0];
    TPZFNMatrix<9,STATE> graddispxy(3, 2);
    TPZAxesTools<STATE>::Axes2XYZ(graddisp, graddispxy, datavec[Eorigin].axes);
    TPZFMatrix<CSTATE> dphixyz(4,dphix.Cols());
//    datavec[Eorigin].axes.Print("Axes ",std::cout);
    // rotate the derivative for the x and y displacement separately
    {
        for(int istate=0; istate<nstate; istate++) {
            TPZFNMatrix<40,CSTATE> dphist(2, nphiH1), dphistxy(3,nphiH1);
            for(int iphi=0; iphi<nphiH1; iphi++) {
                dphist(0,iphi) = dphix(2*istate,iphi);
                dphist(1,iphi) = dphix(2*istate+1,iphi);
            }
            TPZAxesTools<CSTATE>::Axes2XYZ(dphist, dphistxy, datavec[Eorigin].axes);
            for(int iphi=0; iphi<nphiH1; iphi++) {
                dphixyz(2*istate,iphi) = dphistxy(0,iphi);
                dphixyz(2*istate+1,iphi) = dphistxy(1,iphi);
            }
        }
    }
        
    
    TPZFNMatrix<3,STATE> exactderiv(3,2);
    TPZManVector<STATE,1> exactsol(2,0);
    if(fExactSol) {
        fExactSol(x,exactsol,exactderiv);
    }
    TPZFNMatrix<9,REAL> D(3, 3), solvoight(3,1),sigmavoight;
    ComputeDMatrix(fE_def, fnu_def, D);
    solvoight(0,0) = graddispxy(0,0);
    solvoight(1,0) = graddispxy(1,1);
    solvoight(2,0) = graddispxy(1,0)+graddispxy(0,1);
    D.Multiply(solvoight, sigmavoight);

    TPZFNMatrix<9,REAL> solexactvoight(3,1),sigmaexactvoight;
    solexactvoight(0,0) = exactderiv(0,0);
    solexactvoight(1,0) = exactderiv(1,1);
    solexactvoight(2,0) = exactderiv(1,0)+exactderiv(0,1);
    D.Multiply(solexactvoight, sigmaexactvoight);

    bool use_exact = false;
    // compute sigma grad psia *
    TPZManVector<STATE,3> sigmagradhat(2,0.);
    if(!use_exact) {
        sigmagradhat[0] = sigmavoight[0]*dhatvalxy[0]+sigmavoight[2]*dhatvalxy[1];
        sigmagradhat[1] = sigmavoight[2]*dhatvalxy[0]+sigmavoight[1]*dhatvalxy[1];
    } else {
        sigmagradhat[0] = sigmaexactvoight[0]*dhatvalxy[0]+sigmaexactvoight[2]*dhatvalxy[1];
        sigmagradhat[1] = sigmaexactvoight[2]*dhatvalxy[0]+sigmaexactvoight[1]*dhatvalxy[1];

    }
    TPZManVector<STATE,3> XfLoc(2,0.);

    if(fForcingFunction) {
        fForcingFunction(x, XfLoc);
//        std::cout << "x " << x << " XfLoc " << XfLoc << " pressorigin " << pressorigin << std::endl;
    }
//    std::cout << "x " << x << " gradhatgradorig " << gradhatgradorig << std::endl;
    extern bool Print;

    TPZFNMatrix<1, CSTATE> forcemat(2, 1);
    for(int i=0; i<2; i++) forcemat(i,0) = (XfLoc[i]*hatval-sigmagradhat[i]);
    ef.AddContribution(0 ,0, phi, 0, forcemat, 0, weight);
    if(Print) {
        std::cout << "x " << x << " dhatvalxy ";
        for(int i=0; i<2; i++) std::cout << dhatvalxy(i,0) << " ";
        std::cout << std::endl;
        std::cout << "forcemat " << forcemat << std::endl;
//        std::cout << "\nx " << x << " dhatval ";
//        for(int i=0; i<2; i++) std::cout << dhatval(i,0) << " ";
    }

    extern std::complex<STATE> integrateF;
//    integrateF += weight;
//    integrateF += gradhatgradorig*weight;
    integrateF += XfLoc[0]*hatval*weight;
//    integrateF += exactderiv(1,0)*weight;
    int dim = dphix.Rows();
    // computing KPerm \nabla psi_a u
    TPZFNMatrix<3,STATE> BPsia(3, 1,0.);
    if(!use_exact) {
        BPsia(0,0) = dhatvalxy[0]*displacement[0];
        BPsia(1,0) = dhatvalxy[1]*displacement[1];
        BPsia(2,0) = dhatvalxy[1]*displacement[0]+dhatvalxy[0]*displacement[1];
    } else {
        BPsia(0,0) = dhatvalxy[0]*exactsol[0];
        BPsia(1,0) = dhatvalxy[1]*exactsol[1];
        BPsia(2,0) = dhatvalxy[1]*exactsol[0]+dhatvalxy[0]*exactsol[1];
    }
    TPZFNMatrix<3,STATE> sigpsia(3,1);
    D.Multiply(BPsia, sigpsia);
    TPZFNMatrix<3,CSTATE> Csigpsia(3,1, 0.);
    for(int i=0; i<3; i++) Csigpsia(i,0) = sigpsia(i,0);
    TPZFNMatrix<60,CSTATE> dphixyzVoight(3, nphiH1);
    for(int i=0; i<nphiH1; i++) {
        dphixyzVoight(0,i) = dphixyz(0,i);
        dphixyzVoight(1,i) = dphixyz(3,i);
        dphixyzVoight(2,i) = dphixyz(1,i) + dphixyz(2,i);
    }
    if(Print) {
        std::cout << "sigpsia ";
        for(int i=0; i<3; i++) std::cout << sigpsia(i) << " ";
        std::cout << std::endl;
    }
//    for(int d=0; d<3; d++) uhAgradPsi(d,0) = KPerm*exactsol[0]*dhatvalxy(d,0);
//    for(int d=0; d<3; d++) uhAgradPsi(d,0) = KPerm*pressorigin*dhatvalxy(d,0);
    ef.AddContribution(0, 0, dphixyzVoight, 1, Csigpsia, 0, weight);
//    std::cout << "Norm dphixyz " << Norm(dphixyz) << " Norm uhAgradPsi " << Norm(uhAgradPsi) << std::endl;
    int numRBM = RBM.Cols();
    if(numRBM) {
        TPZVec<REAL> &xcenter = datavec[Eaveragepressure].XCenter;
        for(int i=0; i<nphiH1; i++) {
            RBM(i,0) += phi(i,0)*weight;
            RBM(i,1) += phi(i,1)*weight;
            RBM(i,2) += (-(x[1]-xcenter[1])*phi(i,0)+(x[0]-xcenter[0])*phi(i,1))*weight;
        }
    }
}
/**@}*/

void TPZElast2DErrorEstimateMaterial::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const
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
void TPZElast2DErrorEstimateMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                          REAL weight, TPZFMatrix<STATE> &ek,
                          TPZFMatrix<STATE> &ef,
                                                      TPZBndCondT<STATE> &bc) {
    int dim = Dimension();

    TPZFMatrix<REAL> &phiQ = datavec[1].phi;
    auto &x = datavec[1].x;
    int64_t phrq = phiQ.Rows();
    if(2*phrq != ef.Rows()) DebugStop();
    STATE hatval = datavec[2].sol[0][0];

    TPZManVector<REAL,3> v2 = bc.Val2();
    REAL u_D = 0;

    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9, STATE> gradu(3, 2);
        bc.ForcingFunctionBC()(x,res,gradu);
        if (bc.Type() == 0) {
            v2 = res;
        } else if (bc.Type() == 1) {
            TPZFNMatrix<9, STATE> graduV(3,1), D(3,3), sigmaV(3,1);
            graduV(0,0)=gradu(0,0); graduV(1,0) = gradu(1,1);
            graduV(2,0) = gradu(0,1)+gradu(1,0);
            D.Multiply(graduV, sigmaV);
            TPZManVector<REAL,3> normflux(3,0.);
            auto &normal = datavec[0].normal;
            normflux[0] = sigmaV[0]*normal[0]+sigmaV[2]*normal[1];
            normflux[1] = sigmaV[2]*normal[0]+sigmaV[1]*normal[1];
            v2 = normflux;
        } else {
            DebugStop();
        }
    }

    switch (bc.Type()) {
        case 0 :        // Dirichlet condition
            for (int iq = 0; iq < phrq; iq++) {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(2*iq, 0) += (1.) * v2[0] * phiQ(iq, 0) * hatval * weight;
                ef(2*iq+1, 0) += (1.) * v2[1] * phiQ(iq, 0) * hatval * weight;
            }
            break;
        case 1:
            for (int iq = 0; iq <phrq; iq++) {
                for(int jq = 0; jq <phrq; jq++) {
                    ek(2*iq,2*jq) += fBigNumber* phiQ(iq, 0) * phiQ(jq, 0) * weight;
                    ek(2*iq+1,2*jq+1) += fBigNumber* phiQ(iq, 0) * phiQ(jq, 0) * weight;
                }
                ef(2*iq) += fBigNumber* phiQ(iq, 0) * v2[0] * weight;
                ef(2*iq+1) += fBigNumber* phiQ(iq, 0) * v2[1] * weight;
            }
            break;
        default:
            DebugStop();
    }
}

