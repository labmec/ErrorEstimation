/* 
 * File:   TPZHDivErrorEstimateElasticityMaterial.cpp
 * Author: quinelato
 * 
 * Created on September 6, 2023, 10:50 PM
 */

#include "TPZHDivErrorEstimateElasticityMaterial.h"
#include "Util/pzaxestools.h"

TPZHDivErrorEstimateElasticityMaterial::TPZHDivErrorEstimateElasticityMaterial() {
}

TPZHDivErrorEstimateElasticityMaterial::~TPZHDivErrorEstimateElasticityMaterial() {
}

void TPZHDivErrorEstimateElasticityMaterial::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {
    /**
     datavec[0] H1 mesh, uh_reconstructed
     datavec[1] L2 mesh,
     datavec[2] Hdiv fem mesh, sigma_h
     datavec[3] L2 mesh fem, u_h
     
     error[0] - error computed with exact displacement (|| u_fem-u_exact ||)
     error[1] - error computed with reconstructed displacement  (|| u_exact-u_rec ||)
     error[2] = || u_rec - u_fem ||
     error[3] - energy error computed with exact solution  (|| sigma - sigma_fem ||_{C})
     error[4] -  energy error computed with reconstructed displacement  (|| sigma_fem - A epsilon(u_rec)||_{C})
     error[5] - oscilatory data error (|| f - Proj_divsigma ||)
     **/
    
    std:: cout<<"Computing Aposteriori Error Estimation for Linear Elasticity"<<std::endl;
    

    
    int dim= this-> fDimension;
    int nerrors = 6;
    errors.Resize(6);
    errors.Fill(0.0);
    
    TPZFNMatrix<9, STATE> stressfem(dim, dim, 0.);
    TPZManVector<double,3> displacementreconstructed(3,0);
    TPZManVector<double,3> displacementfem(3,0);
    
    TPZManVector<STATE> divsigma(dim,0.);
    TPZManVector<STATE> divsigmafem(dim,0.);
    
    
    for (int i = 0; i < dim; i++) {
        divsigmafem[i] = data[2].divsol[0][i];
    }

    
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            stressfem(i, j) = data[2].sol[0][j + i * 3];
            
        }
    }
    
    


    
    
    STATE divtest=0.;
   
    for (int j=0; j < dim; j++) {
        divtest += data[2].dsol[0](j,j);
    }

    int H1functionposition = 0;
    H1functionposition = FirstNonNullApproxSpaceIndex(data);
    
  
    


    TPZManVector<STATE,3> u_exact(3,0.);
    TPZFNMatrix<9, STATE> du_exact(3, 3);
    if(this->fExactSol){
        
        this->fExactSol(data[H1functionposition].x,u_exact,du_exact);
    
    }
    if(this->HasForcingFunction())
    {
        this->ForcingFunction()(data[H1functionposition].x,divsigma);

    }
    
    REAL residual = 0.;

    for (int idf = 0; idf < dim; idf++) {
        residual += (divsigma[idf] - divsigmafem[idf])*(divsigma[idf] - divsigmafem[idf]);
    }
    
    displacementreconstructed = data[H1functionposition].sol[0];
    displacementfem = data[3].sol[0];
    
    std::cout<<" stressfem ----\n"<<stressfem;
    stressfem.Print(std::cout);
    
    std::cout<<" displacement fem ----\n"<<displacementfem;
    displacementfem.Print(std::cout);
    
    std::cout<<" displacement reconstructed ----\n"<<displacementreconstructed;
    displacementreconstructed.Print(std::cout);

    
    /// calculo do erro de sigma na norma energia || sigma_fem-sigma_ex||_C
    int nstate = fDimension;
    int matdim = nstate*nstate;
    TPZManVector<STATE, 9> Sigma_fem(matdim, 0.), sigma_exactV(matdim, 0.), eps_exactV(matdim, 0.), EPSZV(matdim, 0.);

    TPZFNMatrix<9, STATE> sigma(nstate, nstate, 0.), eps(nstate, nstate, 0.), grad(nstate, nstate, 0.);
    TPZFNMatrix<9, STATE> eps_exact(nstate, nstate, 0.);
    TPZFNMatrix<9, STATE> eps_reconstructed(nstate, nstate, 0.);


    ToVoigt(stressfem, Sigma_fem);
 
//eps(exact displacement)
    eps_reconstructed(0, 0) = du_exact(0, 0);
    eps_reconstructed(1, 0) = eps_reconstructed(0, 1) = 0.5 * (du_exact(0, 1) + du_exact(1, 0));
    eps_reconstructed(1, 1) = du_exact(1, 1);
    
    //eps(reconstructed displacement)
   
    const auto &dudxreconstructed = data[H1functionposition].dsol[0];
    const auto &axes = data[2].axes;
   
    TPZFNMatrix<6,STATE> du(3,3);
    TPZAxesTools<STATE>::Axes2XYZ(dudxreconstructed,du,axes);
    eps_reconstructed(0, 0) =dudxreconstructed(0, 0);
    eps_reconstructed(1, 0) = eps_reconstructed(0, 1) = 0.5 * (dudxreconstructed(0, 1) + dudxreconstructed(1, 0));
    eps_reconstructed(1, 1) = dudxreconstructed(1, 1);

    
    ToVoigt(eps_exact, eps_exactV);
    
    TPZManVector<REAL, 3> x = data[2].x;
    TElasticityAtPoint elast(fE_const,fnu_const);
    

    if(TPZMixedElasticityND::fElasticity)
    {
        //TPZManVector<REAL,3> result(2);
        TPZManVector<STATE, 3> result(2);
        TPZFNMatrix<4,STATE> Dres(0,0);
        fElasticity(x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        TElasticityAtPoint modify(E,nu);
        elast = modify;
    }
    //compute C(sigma) ?
    ComputeStressVector(eps_exactV, sigma_exactV, elast);

    errors[3] = 0.;
    for (int i = 0; i < matdim; i++) {
        //L2_Error: for the stress tensor (sigma)
       // errors[0] += (Sigma_fem[i] - sigma_exactV[i])*(Sigma_fem[i] - sigma_exactV[i]);
        
        //Energy_Error: for the stress tensor (sigma)
        errors[3] += (Sigma_fem[i] - sigma_exactV[i])*(EPSZV[i] - eps_exactV[i]);
    }
    
    /// Como calcular do erro estimado na norma energia?  || sigma_fem - Aeps(u_rec)||_C
    
//
    
#ifdef ERRORESTIMATION_DEBUG2
    std::cout<<"displacement fem "<<displacementfem<<std::endl;
    std::cout<<"displacement reconst "<<displacementreconstructed<<std::endl;
    std::cout<<"-------"<<std::endl;
#endif
    
    //exact error displacement
    errors[0] = 0.;
    for (int idf = 0; idf < dim; idf++) {
        errors[0] += (displacementfem[idf] - u_exact[idf])*(displacementfem[idf] - u_exact[idf]);
    }
    
    //exact error displacement reconstructed
    errors[1]=0;
    for (int idf = 0; idf < dim; idf++) {
        errors[1] += (displacementreconstructed[idf] - u_exact[idf])*(displacementreconstructed[idf] - u_exact[idf]);
    }
    // error displacement reconstructed and displacement fem
    errors[2]=0;
    for (int idf = 0; idf < dim; idf++) {
        errors[2] += (displacementreconstructed[idf] - displacementfem[idf])*(displacementreconstructed[idf] - displacementfem[idf]);
    }

    //||f - Proj_divsigma||
    errors[5] = residual;

#ifdef LOG4CXX
    if(logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Coord: " << data[H1functionIdx].x[0] << ", " << data[H1functionIdx].x[1] << ", "
             << data[H1functionIdx].x[2] << '\n';
        sout << "DisplacementReconstructed = " << displacementreconstructed << "\n";
        sout << "DisplacementFem = " << displacementfem[0] << ", " << displacementfem[1] << ", "
             << displacementfem[2] << "\n";
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
void TPZHDivErrorEstimateElasticityMaterial::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) {

    /**
     datavec[0] H1 mesh, uh_reconstructed
     datavec[1] L2 mesh,
     datavec[2] Hdiv fem mesh, sigma_h
     datavec[3] L2 mesh fem, u_h
     **/

    int H1functionIdx = FirstNonNullApproxSpaceIndex(datavec);

    TPZFNMatrix<9, REAL> PermTensor(3, 3);
    TPZFNMatrix<9, REAL> InvPermTensor(3, 3);
    STATE perm = 0; // this->GetPermeability(data[1].x);
    PermTensor.Diagonal(perm);
    InvPermTensor.Diagonal(1. / perm);
    int dim = this->Dimension();


    TPZManVector<STATE, 2> pressexact(dim, 0.);
    TPZFNMatrix<9, STATE> gradu(3, dim, 0.), fluxinv(3, 1);

    if (TPZMixedElasticityND::HasExactSol()) {
        this->ExactSol()(datavec[H1functionIdx].x, pressexact, gradu);
    }

    //PermTensor.Multiply(gradu, fluxinv);
    TPZFNMatrix<9, STATE> StressExac(3, 3, 0.), StressFem(3, 3, 0.), eps(3, 3, 0.);
    
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < 3; j++) {
            StressFem(i, j) = datavec[2].sol[0][j + i * 3];
        }
    }
    
    switch (var) {
        case 40://Stress FEM
        {
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < 3; j++) {
                    StressFem(i, j) = datavec[2].sol[0][j + i * 3];
                }
            }
            
            Solout[0] = StressFem(0, 0);//StressFemX
            Solout[1] = StressFem(1, 1);//StressFemY
        }
            break;

            
        case 41://StressExact
        {
                TPZVec<STATE> u_exact(fDimension,0.);
                TPZFMatrix<STATE> du_exact(fDimension,fDimension,0.);
                if (this->fExactSol) {
                    this->fExactSol(datavec[2].x, u_exact, du_exact);
                }
                // std::cout << "duexact = " << du_exact << std::endl;
                // For 2D only
                Solout[0] = du_exact(0,0);//Sigma x
                Solout[1] = du_exact(1,1);//Sigma y
                return;
            }
            break;
            
        case 43://displacement Exact
        { TPZVec<STATE> u_exact(fDimension,0.);
            TPZFMatrix<STATE> du_exact(fDimension,fDimension,0.);
            if (this->fExactSol) {
                this->fExactSol(datavec[2].x, u_exact, du_exact);
            }
            for (int idf = 0; idf < dim; idf++) {
                Solout[idf] = u_exact[idf];
            }
            return;
        }
            break;
        case 42://displacement fem
            for (int i = 0; i < dim; i++) {
                    Solout[i] = datavec[3].sol[0][i];
            }
            break;
            
        case 44://displacement Reconstructed
            Solout[0] = datavec[H1functionIdx].sol[0][0];
            Solout[1] = datavec[H1functionIdx].sol[0][1];
            
            break;
        case 45://order p
            Solout[0] = datavec[H1functionIdx].p;
            break;
        default:
            DebugStop();
    }
}

void TPZHDivErrorEstimateElasticityMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    /**
   
     datavec[0] H1 mesh, local uk/grad v for Mark reconstruction and Empty for H1 reconstruction
     datavec[1] L2 mesh, restriction of local uk
     datavec[2] Hdiv mesh, sigmakE
     datavec[3] L2 mesh, ukE
     
     Implement the matrix
     |Sk Ck^T |  = |bk|
     |Ck  0   |    |uk|
     Sk = int_K graduk.gradv dx = int_K gradphi_i.gradphi_j dx
     CK = int_K uk dx = int_K phi_i dx
     bk = int_K sigmafem.gradv dx = int_K sigmafem.gradphi_i dx
     uk = int_K ukfem dx
   
     **/
    
    int H1functionposition = 0;
    H1functionposition = FirstNonNullApproxSpaceIndex(datavec);
    
    int dim = datavec[H1functionposition].axes.Rows();
    //defining test functions
    // Setting the phis
    TPZFMatrix<REAL> &phiuk = datavec[H1functionposition].phi;
    TPZFMatrix<REAL> &dphiukaxes = datavec[H1functionposition].dphix;
    TPZFNMatrix<9,REAL> dphiuk(3,dphiukaxes.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiukaxes, dphiuk, datavec[H1functionposition].axes);
    
    
    int nphiuk = phiuk.Rows();
    
    TPZFMatrix<STATE> solsigmafem(3,nphiuk),solukfem(1,1);
    solsigmafem.Zero();
    solukfem.Zero();
    
    
        //potetial fem
        solukfem(0,0) = datavec[3].sol[0][0];
        //flux fem
        for (int ip = 0; ip<3; ip++){

            solsigmafem(ip,0) = datavec[2].sol[0][ip];
        }

    
    
    TPZFNMatrix<9,REAL> PermTensor(3,3);
    TPZFNMatrix<9,REAL> InvPermTensor(3,3);
    STATE perm = 0; //this->GetPermeability(datavec[1].x);
    PermTensor.Diagonal(perm);
    InvPermTensor.Diagonal(1./perm);
        
    
    TPZFMatrix<STATE> kgraduk(3,nphiuk,0.);
    
        
    for(int irow=0 ; irow<nphiuk; irow++){
        
        //K graduk
        for(int id=0; id< dim; id++){
            
            for(int jd=0; jd< dim;jd++){
                
                kgraduk(id,irow) += PermTensor(id,jd)*dphiuk(jd,irow);
                
            }
            //bk = (-1)*int_k sigmaukfem.grad phi_i,here dphiuk is multiplied by axes
            //the minus sign is necessary because we are working with sigma_h = - K grad u, Mark works with sigma_h = K grad u
            
            ef(irow,0)+=(-1.)*weight*dphiuk(id,irow)*solsigmafem(id,0);
        }
        
        //matrix Sk= int_{K} K graduk.gradv
        for(int jcol=0; jcol<nphiuk;jcol++){
            
            for(int jd=0;  jd< dim; jd++)
            {
                ek(irow,jcol) +=weight*kgraduk(jd,irow)*dphiuk(jd,jcol);
            }
            
        }
        //Ck=int_{K} phi_i e Ck^t
        if(fNeumannLocalProblem){
            
            ek(irow,nphiuk)+= weight*phiuk(irow,0);
            ek(nphiuk,irow)+= weight*phiuk(irow,0);
            
        }
        
    }
    if(H1functionposition ==0){
    
        if(!fNeumannLocalProblem )
        {
            ek(nphiuk,nphiuk) += weight;
            ef(nphiuk,0)+= weight;
        }
        
        //muk = int_k ukfem
        
        else{
            
            ef(nphiuk,0)+= weight*solukfem(0,0);
        
        }
    }
    
}

void TPZHDivErrorEstimateElasticityMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                                TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {

    /*
     Add Robin boundary condition for local problem
     ek+= <w,Km s_i>
     ef+= <w,Km*u_d - g + sigma_i.n>
    */
    int H1functionposition = FirstNonNullApproxSpaceIndex(datavec);
    int dim = datavec[H1functionposition].axes.Rows();

    TPZFMatrix<REAL> &phi_i = datavec[H1functionposition].phi;
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
    STATE perm = 0; // GetPermeability(datavec[1].x);
    PermTensor.Diagonal(perm);
    InvPermTensor.Diagonal(1./perm);

    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9, STATE> gradu(dim, 1);
        bc.ForcingFunctionBC()(datavec[H1functionposition].x, res, gradu);
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
                ef(iq, 0) += fBigNumber * u_D * phi_i(iq, 0) * weight;
                for (int jq = 0; jq < nphi_i; jq++) {
                    ek(iq, jq) += fBigNumber * weight * phi_i(iq, 0) * phi_i(jq, 0);
                }
            }
            break;
        }

        default: {
            // std::cout << " This material not implement BC Type " << bc.Type()<< std::endl;
            break;
        }
    }
}



int TPZHDivErrorEstimateElasticityMaterial::VariableIndex(const std::string &name)const
{
    if(name == "StressFem") return 40;
    
    if(name == "StressExact") return 41;
    
    if(name == "DisplacementFem") return 42;
    if(name == "DisplacementExact") return 43;
    if(name == "DisplacementReconstructed") return 44;
    if(name == "POrder") return 45;
    if(name == "DisplacementErrorExact") return 100;
    if(name == "DisplacementErrorEstimate") return 101;
    if(name == "EnergyErrorExact") return 102;
    if(name == "EnergyErrorEstimate") return 103;
    if(name == "ResidualError") return 104;
    if(name == "DisplacementEffectivityIndex") return 105;
    if(name == "EnergyEffectivityIndex") return 106;
    
     
    return -1;
}
