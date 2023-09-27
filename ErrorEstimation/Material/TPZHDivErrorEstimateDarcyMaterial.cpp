/* 
 * File:   TPZHDivErrorEstimateDarcyMaterial.cpp
 * Author: quinelato
 * 
 * Created on September 6, 2023, 10:50 PM
 */

#include "TPZHDivErrorEstimateDarcyMaterial.h"
#include "Util/pzaxestools.h"

void TPZHDivErrorEstimateDarcyMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    /**
   
     datavec[0] H1 mesh, local uk/grad v for Mark's reconstruction and Empty for H1 reconstruction
     datavec[1] L2 mesh, restriction of local uk
     datavec[2] Hdiv mesh, sigmakE
     datavec[3] L2 mesh, ukE
     
     Implements the matrix
     |Sk Ck^T |  = |bk|
     |Ck  0   |    |uk|
     Sk = int_K graduk.gradv dx = int_K gradphi_i.gradphi_j dx
     CK = int_K uk dx = int_K phi_i dx
     bk = int_K sigmafem.gradv dx = int_K sigmafem.gradphi_i dx
     uk = int_K ukfem dx
   
     **/

    int H1functionIdx = FirstNonNullApproxSpaceIndex(datavec);

    int dim = datavec[H1functionIdx].axes.Rows();
    //defining test functions
    // Setting the phis
    TPZFMatrix<REAL> &phiuk = datavec[H1functionIdx].phi;
    TPZFMatrix<REAL> &dphiukaxes = datavec[H1functionIdx].dphix;
    TPZFNMatrix<9, REAL> dphiuk(3, dphiukaxes.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiukaxes, dphiuk, datavec[H1functionIdx].axes);

    int nphiuk = phiuk.Rows();

    TPZFMatrix<STATE> solsigmafem(3, nphiuk), solukfem(1, 1);
    solsigmafem.Zero();
    solukfem.Zero();

    //potential fem
    solukfem(0, 0) = datavec[3].sol[0][0];
    //flux fem
    for (int ip = 0; ip < 3; ip++) {
        solsigmafem(ip, 0) = datavec[2].sol[0][ip];
    }

    TPZFNMatrix<9, REAL> PermTensor(3, 3);
    TPZFNMatrix<9, REAL> InvPermTensor(3, 3);
    STATE perm = this->GetPermeability(datavec[1].x);
    PermTensor.Diagonal(perm);
    InvPermTensor.Diagonal(1. / perm);

    //K graduk
    TPZFMatrix<STATE> kgraduk(3, nphiuk, 0.);
    PermTensor.Multiply(dphiuk, kgraduk);
    for (int irow = 0; irow < nphiuk; irow++) {
        for (int id = 0; id < dim; id++) {
            // bk = (-1)*int_k sigmaukfem.grad phi_i, here dphiuk is multiplied by axes
            // the minus sign is necessary because we are working with sigma_h = - K grad u, Mark works with sigma_h = K grad u
            ef(irow, 0) += (-1.) * weight * dphiuk(id, irow) * solsigmafem(id, 0);
        }

        //matrix Sk= int_{K} K graduk.gradv
        for (int jcol = 0; jcol < nphiuk; jcol++) {
            for (int jd = 0; jd < dim; jd++) {
                ek(irow, jcol) += weight * kgraduk(jd, irow) * dphiuk(jd, jcol);
            }
        }
        //Ck=int_{K} phi_i e Ck^t
        if (fNeumannLocalProblem) {
            ek(irow, nphiuk) += weight * phiuk(irow, 0);
            ek(nphiuk, irow) += weight * phiuk(irow, 0);
        }
    }
    if (H1functionIdx == 0) {
        if (!fNeumannLocalProblem) {
            ek(nphiuk, nphiuk) += weight;
            ef(nphiuk, 0) += weight;
        }
            //muk = int_k ukfem
        else {
            ef(nphiuk, 0) += weight * solukfem(0, 0);
        }
    }
}

void TPZHDivErrorEstimateDarcyMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                                TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {
    /*
     Add Robin boundary condition for local problem
     ek+= <w,Km s_i>
     ef+= <w,Km*u_d - g + sigma_i.n>
    */
    int H1functionIdx = FirstNonNullApproxSpaceIndex(datavec);
    int dim = datavec[H1functionIdx].axes.Rows();

    TPZFMatrix<REAL> &phi_i = datavec[H1functionIdx].phi;
    int nphi_i = phi_i.Rows();

    TPZFMatrix<STATE> solsigmafem(3, 1);
    solsigmafem.Zero();

    TPZManVector<REAL, 3> normal = datavec[2].normal;

    REAL normalsigma = datavec[2].sol[0][0];

    REAL u_D;
    REAL g = 0.;
    REAL normflux = 0.;
    TPZFNMatrix<9,REAL> PermTensor(3,3);
    TPZFNMatrix<9,REAL> InvPermTensor(3,3);
    STATE perm = this->GetPermeability(datavec[1].x);
    PermTensor.Diagonal(perm);
    InvPermTensor.Diagonal(1./perm);
    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9, STATE> gradu(dim, 1);
        bc.ForcingFunctionBC()(datavec[H1functionIdx].x, res, gradu);
        u_D = res[0];

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                normflux += datavec[2].normal[i] * PermTensor(i, j) * gradu(j, 0);
            }
        }
        g = (-1.) * normflux;
    } else {
        // u_D is usually stored in val2(0, 0)
        u_D = bc.Val2()[0];
    }

    switch (bc.Type()) {
        case (4): {
            REAL Km = bc.Val1()(0, 0);
            REAL robinterm = (Km * u_D - g + normalsigma);

            for (int iq = 0; iq < nphi_i; iq++) {
                //<w,Km*u_D-g+sigma_i*n>
                ef(iq, 0) += robinterm * phi_i(iq, 0) * weight;
                for (int jq = 0; jq < nphi_i; jq++) {
                    //<w,Km*s_i>
                    ek(iq, jq) += weight * Km * phi_i(iq, 0) * phi_i(jq, 0);
                }
            }
            break;
        }

        case (0): {
            for (int iq = 0; iq < nphi_i; iq++) {
                ef(iq, 0) += this->fBigNumber * u_D * phi_i(iq, 0) * weight;
                for (int jq = 0; jq < nphi_i; jq++) {
                    ek(iq, jq) += this->fBigNumber * weight * phi_i(iq, 0) * phi_i(jq, 0);
                }
            }
            break;
        }

        default: {
            // std::cout << " This material does not implement BC Type " << bc.Type()<< std::endl;
            break;
        }
    }
}

void TPZHDivErrorEstimateDarcyMaterial::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {    
    /**
     datavec[0] H1 mesh, uh_reconstructed
     datavec[1] L2 mesh,
     datavec[2] Hdiv fem mesh, sigma_h
     datavec[3] L2 mesh fem, u_h
     
      error[0] - error computed with exact pressure
      error[1] - error computed with reconstructed pressure
      error[2] - energy error computed with exact solution
      error[3] - energy error computed with reconstructed solution
      error[4] - oscilatory data error
     **/
    
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);
    
    TPZManVector<STATE,3> fluxfem(3);
    STATE divsigmafem, pressurefem, pressurereconstructed;
    divsigmafem = 0.;
    pressurefem = 0.;
    pressurereconstructed = 0.;
    
    TPZFNMatrix<3,REAL> fluxreconstructed(3,1), fluxreconstructed2(3,1);
    fluxfem = data[2].sol[0];
    divsigmafem= data[2].divsol[0][0];
    
    STATE divtest=0.;
    for (int j=0; j < this->Dimension(); j++) {
        divtest += data[2].dsol[0](j,j);
    }

    int H1functionposition = 0;
    H1functionposition = FirstNonNullApproxSpaceIndex(data);
    
    TPZVec<STATE> divsigma(1);
    divsigma[0]=0.;

    TPZManVector<STATE, 1> u_exact(1);
    TPZFNMatrix<9, STATE> du_exact(3, 3);
    if(this->fExactSol){
        
        this->fExactSol(data[H1functionposition].x,u_exact,du_exact);
    
    }
    if(this->HasForcingFunction())
    {
        this->ForcingFunction()(data[H1functionposition].x,divsigma);

    }
    
    REAL residual = 0.;

    residual = (divsigma[0] - divsigmafem)*(divsigma[0] - divsigmafem);
   
    
    pressurereconstructed = data[H1functionposition].sol[0][0];
  
 
        pressurefem = data[3].sol[0][0];
    
    
    TPZFNMatrix<9,REAL> PermTensor(3,3);
    TPZFNMatrix<9,REAL> InvPermTensor(3,3);
    STATE perm = TPZMixedDarcyFlow::GetPermeability(data[1].x);
    PermTensor.Diagonal(perm);
    InvPermTensor.Diagonal(1./perm);
    
    TPZFNMatrix<3,REAL> fluxexactneg;
    
    //sigmarec = -K grad(urec)
    //  sigmak = -K graduk
    
    {
        TPZFNMatrix<9,REAL> gradpressure(3,1);
        for (int i=0; i<3; i++) {
            gradpressure(i,0) = du_exact[i];
        }
        PermTensor.Multiply(gradpressure,fluxexactneg);
    }
    
    
    
    TPZFMatrix<REAL> &dsolaxes = data[H1functionposition].dsol[0];
    TPZFNMatrix<9,REAL> fluxrec(3,0);
    TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, fluxrec, data[H1functionposition].axes);
    
    for(int id=0 ; id<3; id++) {
        fluxreconstructed2(id,0) = (-1.)*fluxrec(id,0);
    }
    
     PermTensor.Multiply(fluxreconstructed2,fluxreconstructed);
    
    
    REAL innerexact = 0.;
    REAL innerestimate = 0.;
    
    
    
#ifdef ERRORESTIMATION_DEBUG2
    std::cout<<"flux fem "<<fluxfem<<std::endl;
    std::cout<<"flux reconst "<<fluxreconstructed<<std::endl;
    std::cout<<"-------"<<std::endl;
#endif
    
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            innerexact += (fluxfem[i]+fluxexactneg(i,0))*InvPermTensor(i,j)*(fluxfem[j]+fluxexactneg(j,0));//Pq esta somando: o fluxo fem esta + e o exato -
            innerestimate += (fluxfem[i]-fluxreconstructed[i])*InvPermTensor(i,j)*(fluxfem[j]-fluxreconstructed[j]);
        }
    }
    
#ifdef ERRORESTIMATION_DEBUG2
    std::cout<<"potential fem "<<pressurefem<<std::endl;
    std::cout<<"potential reconst "<<pressurereconstructed<<std::endl;
    std::cout<<"-------"<<std::endl;
#endif
    errors[0] = (pressurefem-u_exact[0])*(pressurefem-u_exact[0]);//exact error pressure
    errors[1] = (pressurefem-pressurereconstructed)*(pressurefem-pressurereconstructed);//error pressure reconstructed
    errors[2] = innerexact;//error flux exact
    errors[3] = innerestimate;//error flux reconstructed
    errors[4] = residual; //||f - Proj_divsigma||

#ifdef LOG4CXX
    if(logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Coord: " << data[H1functionIdx].x[0] << ", " << data[H1functionIdx].x[1] << ", "
             << data[H1functionIdx].x[2] << '\n';
        sout << "PressureReconstructed = " << pressurereconstructed << "\n";
        sout << "FluxReconstructed = " << fluxreconstructed[0] << ", " << fluxreconstructed[1] << ", "
             << fluxreconstructed[2] << "\n";
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

/**
 * @brief It return a solution to multiphysics simulation.
 * @param datavec [in] Data material vector
 * @param var [in] number of solution variables. See  NSolutionVariables() method
 * @param Solout [out] is the solution vector
 */
void TPZHDivErrorEstimateDarcyMaterial::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) {

    /**
     datavec[0] H1 mesh, uh_reconstructed for Mark reconstruction and Empty for H1 reconstruction
     datavec[1] L2 mesh,
     datavec[2] Hdiv fem mesh, sigma_h
     datavec[3] L2 mesh fem, u_h
     **/

    int H1functionIdx = FirstNonNullApproxSpaceIndex(datavec);

    TPZFNMatrix<9, REAL> PermTensor(3, 3);
    TPZFNMatrix<9, REAL> InvPermTensor(3, 3);
    STATE perm = TPZMixedDarcyFlow::GetPermeability(datavec[1].x);
    PermTensor.Diagonal(perm);
    InvPermTensor.Diagonal(1. / perm);



    TPZManVector<STATE, 2> pressexact(1, 0.);
    TPZFNMatrix<9, STATE> gradu(3, 1, 0.), fluxinv(3, 1);

    if (TPZMixedDarcyFlow::HasExactSol()) {
        this->ExactSol()(datavec[H1functionIdx].x, pressexact, gradu);
    }

    PermTensor.Multiply(gradu, fluxinv);
    int dim = this->fDim;
    switch (var) {
        case 40://FluxFem
            for (int i = 0; i < 3; i++) {
                Solout[i] = datavec[2].sol[0][i];
            }
            break;
        case 41:
        {//FluxReconstructed is grad U
            TPZFMatrix<REAL> &dsolaxes = datavec[H1functionIdx].dsol[0];
            TPZFNMatrix<9, REAL> dsol(3, 0);
            TPZFNMatrix<9, REAL> KGradsol(3, 0);
            TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, dsol, datavec[H1functionIdx].axes);

            PermTensor.Multiply(dsol, KGradsol);

            for (int id = 0; id<this->Dimension(); id++) {
                Solout[id] = -KGradsol(id, 0); //dsol(id,0);//derivate
            }
        }
            break;
        case 42://flux exact
            for (int i = 0; i < dim; i++) Solout[i] = -fluxinv(i);
            break;
        case 43://Pressure fem
            Solout[0] = datavec[3].sol[0][0];
            break;
        case 44://PressureReconstructed
            Solout[0] = datavec[H1functionIdx].sol[0][0];
            break;
        case 45://pressureexact
            Solout[0] = pressexact[0];
            break;
        case 46://order p
            Solout[0] = datavec[1].p;
            break;
        default:
            DebugStop();
    }
}