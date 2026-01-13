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

int TPZHDivErrorEstimateElasticityMaterial::NEvalErrors() const {
    /*
     * errors[0] - error computed with exact displacement (|| u_fem-u_exact ||)
     * errors[1] - error computed with reconstructed displacement  (|| u_exact-u_rec ||)
     * errors[2] = || u_rec - u_fem ||
     * errors[3] - energy error computed with exact solution  (|| sigma - sigma_fem ||_{C})
     * errors[4] -  energy error computed with reconstructed displacement  (|| sigma_fem - A epsilon(u_rec)||_{C})
     * errors[5] - oscilatory data error (|| f - Proj_divsigma ||)
     * TODO: errors[6] - antisymmetric part (||0.5*(sigma_fem-sigma_fem^T)||_{C})
     * TODO: error[7] - ||u_femH1-u_ex||
     * error[8] = | sigma_fem - A epsilon(u_femH1)||_{C}
     * error[9] = | sigma - A epsilon(u_femH1) ||_{C})
     **/
    return 10;
}

void TPZHDivErrorEstimateElasticityMaterial::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {
    /**
     data[0] H1 mesh, uh_reconstructed
     data[1] L2 mesh,
     data[2] Hdiv fem mesh, sigma_h
     data[3] L2 mesh fem, u_h
     
     error[0] - error computed with exact displacement (|| u_fem-u_exact ||)
     error[1] - error computed with reconstructed displacement  (|| u_exact-u_rec ||)
     error[2] - energy error computed with exact solution  (|| sigma_ex - sigma_fem ||_{C})
     error[3] - energy error computed with reconstructed displacement  (|| sigma_fem - A epsilon(u_rec)||_{C})
     error[4] = || u_rec - u_fem ||
     error[5] - oscilatory data error (|| f - Proj_divsigma ||)
     TODO: 
     errors[6] - antisymmetric part (||0.5*(sigma_fem-sigma_fem^T)||_{C})
     error[7] = ||u_femH1-u_ex||
     error[8] = | sigma_fem - A epsilon(u_femH1)||_{C}
     error[9] = | sigma_ex - A epsilon(u_femH1) ||_{C})
     **/

    //std::cout << "Computing Aposteriori Error Estimation for Linear Elasticity" << std::endl;

    int dim = this->Dimension();
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);
//stressfem = sigmafem
    TPZFNMatrix<9, STATE> stressfem(dim, dim, 0.);
    
    for (unsigned int i = 0; i < dim; i++) {
        for (unsigned int j = 0; j < dim; j++) {
            stressfem(i, j) = data[2].sol[0][j + i * 3];
        }
    }
    
//Symmetric part of stressfem
    
        TPZFNMatrix<9, STATE> stressSym(dim, dim, 0.0);

        for (unsigned int i = 0; i < dim; i++) {
            for (unsigned int j = 0; j < dim; j++) {
                stressSym(i, j) = 0.5 * (stressfem(i, j) + stressfem(j, i));
            }
        }

  //  std::cout<<"stress "<<stressfem<<std::endl;
    
    
   // std::cout<<"stressSym "<<stressSym<<std::endl;
    
    
    //Antisymmetric part of stressfem (stressfemAS)
    TPZFNMatrix<9, STATE> stressfemAS(dim, dim, 0.);
    for (unsigned int i = 0; i < dim; i++) {
        for (unsigned int j = 0; j < dim; j++) {
            if (i==j) continue;
            stressfemAS(i,j) = 0.5*(stressfem(i, j)-stressfem(j, i));
        }
    }

   // std::cout<<"stressAS "<<stressfemAS<<std::endl;

    //
    TPZFNMatrix<9, STATE> rotfem(dim, dim, 0.);
    rotfem(0,1) = data[4].sol[0][0];
    rotfem(1,0) = (-1.)*data[4].sol[0][0];
    
  //  std::cout<<"rot= "<<rotfem<<std::endl;
    
    TPZManVector<STATE> divstressfem(dim, 0.);
    divstressfem.Fill(0);
    for (int i = 0; i < dim; i++) {
        divstressfem[i] = data[2].divsol[0][i];
    }

    int H1functionposition = 0;
    H1functionposition = FirstNonNullApproxSpaceIndex(data);


    TPZManVector<STATE> divstress(dim, 0.);
    
    divstress.Fill(0);

    TPZManVector<STATE, 3> u_exact(dim, 0.);
    TPZFNMatrix<9, STATE> du_exact(dim, dim, 0.);
    if (this->fExactSol) {
        this->fExactSol(data[H1functionposition].x, u_exact, du_exact);
    }
    if (this->HasForcingFunction()) {
        this->ForcingFunction()(data[H1functionposition].x, divstress);
        for (unsigned int i = 0; i < dim; ++i) {
            divstress[i] *=-1;
        }
    }

    REAL residual = 0.;
    for (unsigned int idf = 0; idf < dim; idf++) {
        residual = (divstress[idf] - divstressfem[idf])*(divstress[idf] - divstressfem[idf]);
    }
    
    //||f - Proj_divsigma||
    errors[5] = residual;
    
  //  std::cout<<"residual error "<<residual<<std::endl;

    //TODO:add the new error for displacement in H1
    TPZManVector<STATE, 3> displacementfemH1(dim, 0.);
    displacementfemH1 = data[5].sol[0];
    
    TPZManVector<STATE, 3> displacementreconstructed(dim, 0.);
    displacementreconstructed = data[H1functionposition].sol[0];
    

    

    TPZManVector<STATE, 3> displacementfem(dim, 0.);
    displacementfem = data[3].sol[0];
    

    /// calculo do erro de sigma na norma energia || sigma_fem-sigma_ex||_C
    int nstate = dim;
    int matdim = nstate*nstate;
    TPZManVector<STATE, 9> stress_femV(matdim, 0.), sigma_exactV(matdim, 0.), eps_exactV(matdim, 0.), EPSZV(matdim, 0.),stressfemAS_V(matdim,0.),stressSym_V(matdim,0.),rotfem_V(matdim,0.);
    TPZFNMatrix<9, STATE> sigma(nstate, nstate, 0.), eps(nstate, nstate, 0.), grad(nstate, nstate, 0.);
    TPZFNMatrix<9, STATE> eps_exact(nstate, nstate, 0.);
    TPZFNMatrix<9, STATE> eps_reconstructed(nstate, nstate, 0.);
    TPZFNMatrix<9, STATE> assym_reconstructed(nstate, nstate, 0.);

    ToVoigt(stressfem, stress_femV);
    ToVoigt(stressSym, stressSym_V);
    ToVoigt(stressfemAS, stressfemAS_V);
    ToVoigt(rotfem, rotfem_V);
    

    //eps(exact displacement)
    eps_exact(0, 0) = du_exact(0, 0);
    eps_exact(1, 0) = eps_exact(0, 1) = 0.5 * (du_exact(0, 1) + du_exact(1, 0));
    eps_exact(1, 1) = du_exact(1, 1);

    //eps(reconstructed displacement)
    const auto &dudxreconstructed = data[H1functionposition].dsol[0];
    const auto &axes = data[H1functionposition].axes;

    TPZFNMatrix<6, STATE> du(3, 3),gradS(3, 3);
    TPZAxesTools<STATE>::Axes2XYZ(dudxreconstructed, du, axes);
    eps_reconstructed(0, 0) = du(0, 0);
    eps_reconstructed(1, 0) = eps_reconstructed(0, 1) = 0.5 * (du(0, 1) + du(1, 0));
    eps_reconstructed(1, 1) = du(1, 1);
    
    //std::cout<<"eps_reconstructed "<<eps_reconstructed<<std::endl;
    
    gradS(0,0)=du(0,0);
    gradS(0,1)=du(0,1);
    gradS(1,0)=du(1,0);
    gradS(1,1)=du(1,1);
    TPZManVector<STATE, 9> gradS_V(matdim, 0.);
    ToVoigt(gradS, gradS_V);
    
    TPZManVector<STATE, 9> AgradS_V(matdim, 0.);
    
    
   
   
    
    
    //R(reconstructed displacement)=0.5(gradS-gradS^T)
    assym_reconstructed(0, 0) = 0.;
    assym_reconstructed(1, 0) =  0.5 * (du(1, 0) - du(0, 1));
    assym_reconstructed(0, 1) =(-1.)*assym_reconstructed(1, 0);
    assym_reconstructed(1, 1) = 0.;
    
    
    //operacoes com H1
    
    const auto &dudxh1 = data[5].dsol[0];
    const auto &axesh1 = data[5].axes;
    TPZFNMatrix<6, STATE> duh1(3, 3);
    TPZFNMatrix<9, STATE> eps_h1(nstate, nstate, 0.);
    
    TPZAxesTools<STATE>::Axes2XYZ(dudxh1, duh1, axesh1);
    eps_h1(0, 0) = duh1(0, 0);
    eps_h1(1, 0) = eps_h1(0, 1) = 0.5 * (duh1(0, 1) + duh1(1, 0));
    eps_h1(1, 1) = duh1(1, 1);
    
    TPZManVector<STATE, 9> eps_h1V(matdim, 0.);
    TPZManVector<STATE, 9> sigma_h1V(matdim, 0.);
    ToVoigt(eps_h1, eps_h1V);
    
    
    
    
    
    
    
#ifdef ERRORESTIMATION_DEBUG2
    std::cout<<"eps exact = "<<eps_exact(0,0)<<" , "<<eps_exact(1,0)<<","<<eps_exact(1,1)<<"\n";
    std::cout<<"eps rec = "<<eps_reconstructed(0,0)<<" , "<<eps_reconstructed(1,0)<<","<<eps_reconstructed(1,1)<<"\n";
    std::cout<<"----\n";
    
#endif

    ToVoigt(eps_exact, eps_exactV);

    TPZManVector<REAL, 3> x = data[2].x;

    TElasticityAtPoint elast(fE_const, fnu_const);
   // std::cout<<"E = "<<fE_const<<std::endl;
   // std::cout<<"nu = "<<fnu_const<<std::endl;

    if (TPZMixedElasticityND::fElasticity) {
        //TPZManVector<REAL,3> result(2);
        TPZManVector<STATE, 3> result(2);
        TPZFNMatrix<4, STATE> Dres(0, 0);
        fElasticity(x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        TElasticityAtPoint modify(E, nu);
        elast = modify;
        
        //std::cout<<"E = "<<E<<std::endl;
       // std::cout<<"nu = "<<nu<<std::endl;
    }
    
    
    
    
    /// || sigma - sigma_fem||_C^2 = (C(sigma - sigma_fem), sigma - sigma_fem) = (eps - Csigma_fem, sigma - sigma_fem)
    /// compute C(sigma)
    ComputeStressVector(eps_exactV, sigma_exactV, elast);
    TPZManVector<STATE, 9> Csigma_femV(matdim, 0.),CsigmaSym_femV(matdim,0.);
    ComputeDeformationVector(stress_femV, Csigma_femV, elast);
    ComputeDeformationVector(stressSym_V, CsigmaSym_femV, elast);
    
    TPZManVector<STATE, 9> part1(matdim, 0.);
     TPZManVector<STATE, 9> part2(matdim, 0.);
     for (unsigned int i = 0; i < matdim; ++i) {
         part1[i] = eps_exactV[i]-Csigma_femV[i];
         part2[i] = sigma_exactV[i]-stress_femV[i];
     }
     errors[2] =TPZMixedElasticityND::InnerVec(part1, part2);
 
    /// || sigma_fem^S - Aeps(u_rec)||_C^2 = (Csigma_fem^S - eps(u_rec), sigma_fem^S - Aeps(u_rec))
  
    ComputeStressVector(gradS_V, AgradS_V, elast);
    
    
    TPZManVector<STATE, 9> Sigma_reconstructed(matdim, 0.);
    TPZManVector<STATE, 9> sigma_reconstructedV(matdim, 0.), eps_reconstructedV(matdim, 0.) ;
    ToVoigt(eps_reconstructed, eps_reconstructedV);
   
    ComputeStressVector(eps_reconstructedV, sigma_reconstructedV, elast);
    
    part1.Fill(0.);
    part2.Fill(0.);
    

    for (unsigned int i = 0; i < matdim; ++i) {
       // part1[i] = CsigmaSym_femV[i] - eps_reconstructedV[i];
      //  part2[i] = stressSym_V[i] - sigma_reconstructedV[i];
        
//        part1[i] = Csigma_femV[i] - gradS_V[i];
//        part2[i] = stress_femV[i] - AgradS_V[i];
        part1[i] = Csigma_femV[i] - eps_reconstructedV[i];
        part2[i] = stress_femV [i] - sigma_reconstructedV[i];

    }
    errors[3] = TPZMixedElasticityND::InnerVec(part1, part2);
    
    
    part1.Fill(0.);
    part2.Fill(0.);
    

    //||sigma_femAS-AR(urec)||_C=(Csigma_femAS-R(urec), sigma_femAS-AR(urec))
    TPZManVector<STATE, 9> assym_reconstructedV(matdim, 0.),ARot_reconstructedV(matdim,0.),Csigma_femAS_V(matdim,0.),CRot_femV(matdim,0.) ;
    ToVoigt(assym_reconstructed, assym_reconstructedV);
    
   
    
    ComputeStressVector(assym_reconstructedV, ARot_reconstructedV, elast);
    ComputeDeformationVector(rotfem_V, CRot_femV, elast);
    
    
    
    
    ComputeDeformationVector(stressfemAS_V, Csigma_femAS_V, elast);
    
//    std::cout<<"------- "<<std::endl;
//    std::cout<<"Csigma_femV - "<<Csigma_femV<<std::endl;
//    std::cout<<"Csigma_femAS_V - "<<Csigma_femAS_V<<std::endl;
//    std::cout<<"CsigmaSym_femV - "<<CsigmaSym_femV<<std::endl;
//    std::cout<<"assym_reconstructedV - "<<assym_reconstructedV<<std::endl;
//    std::cout<<"ARot_reconstructedV - "<<ARot_reconstructedV<<std::endl;
//    std::cout<<"eps_reconstructedV - "<<eps_reconstructedV<<std::endl;
//    std::cout<<"sigma_reconstructedV - "<<sigma_reconstructedV<<std::endl;
//    std::cout<<"rotfem_V - "<<rotfem_V<<std::endl;
//    std::cout<<"CRot_femV - "<<CRot_femV<<std::endl;
    
    
    for (unsigned int i = 0; i < matdim; ++i) {
        
       // std::cout<<"Csigma_femAS_V[i]= "<<Csigma_femAS_V[i]<< " assym_reconstructedV[i]= "<<assym_reconstructedV[i]<<std::endl;
        part1[i] = (Csigma_femAS_V[i] - assym_reconstructedV[i]);
        //std::cout<<"------"<<std::endl;
      //  std::cout<<"stressfemAS_V[i]= "<<stressfemAS_V[i]<< " ARot_reconstructedV[i]= "<<ARot_reconstructedV[i]<<std::endl;
        part2[i] = (stressfemAS_V[i] -ARot_reconstructedV[i]);
        

    }

    errors[6] = TPZMixedElasticityND::InnerVec(part1, part2);
    
//calculando somente a norma C de sigma_femAS
 //   errors[6] = TPZMixedElasticityND::InnerVec(Csigma_femAS_V, stressfemAS_V);

    
    
    /// || sigma_fem - Aeps(u_h1)||_C^2 = (Csigma_fem - eps(u_h1), sigma_fem - Aeps(u_h1))
                            /// = (Csigma_femV - eps_h1V, stress_femV - sigma_h1V)
    

    
    ComputeStressVector(eps_h1V, sigma_h1V, elast);
    
    
    //std:: cout<<"---Csigma_femV "<<Csigma_femV<<"eps_h1 "<<eps_h1V<<"\n";
    //std:: cout<<"---stress_femV "<<stress_femV<<"sigma_h1V "<<sigma_h1V<<"\n";
    
    part1.Fill(0.);
    part2.Fill(0.);
    
    
    for (unsigned int i = 0; i < matdim; ++i) {
        part1[i] = Csigma_femV[i] - eps_h1V[i];
        part2[i] = stress_femV[i] - sigma_h1V[i];

    }
    errors[8] = TPZMixedElasticityND::InnerVec(part1, part2);
    
    //
    part1.Fill(0.);
    part2.Fill(0.);
    
    for (unsigned int i = 0; i < matdim; ++i) {
        part1[i] = eps_exactV[i] - eps_h1V[i];
        part2[i] = sigma_exactV[i] - sigma_h1V[i];

    }
    errors[9] = TPZMixedElasticityND::InnerVec(part1, part2);
    
    

#ifdef ERRORESTIMATION_DEBUG2
    std::cout << "stress fem " << stress_femV << std::endl;
    std::cout << "stress reconst " << sigma_reconstructedV << std::endl;
    std::cout << "stress exact " << sigma_exactV << std::endl;
    std::cout << "-------" << std::endl;
#endif

    //exact error displacement
    errors[0] = 0.;
    errors[7] = 0.;
//    std::cout << " displacement fem = " << displacementfem<<std::endl;
//    std::cout << "---displacement femH1= " << displacementfemH1<<std::endl;
//    std::cout << "----displacement reconstructed= " << displacementreconstructed<<std::endl;
//    std::cout << "displacement exact= " << u_exact<<std::endl;
    for (int idf = 0; idf < dim; idf++) {
        errors[0] += (displacementfem[idf] - u_exact[idf])*(displacementfem[idf] - u_exact[idf]);
        errors[7]+= (displacementfemH1[idf] - u_exact[idf])*(displacementfemH1[idf] - u_exact[idf]);
    }
  
    
    //exact error displacement reconstructed
    errors[1] = 0;
    for (int idf = 0; idf < dim; idf++) {
        errors[1] += (displacementreconstructed[idf] - u_exact[idf])*(displacementreconstructed[idf] - u_exact[idf]);
    }
    // error displacement reconstructed and displacement fem
    errors[4] = 0;
    for (int idf = 0; idf < dim; idf++) {
        errors[4] += (displacementreconstructed[idf] - displacementfem[idf])*(displacementreconstructed[idf] - displacementfem[idf]);
    }


    
   // std::cout<<"residual "<<residual<<std::endl;
    
  

#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Coord: " << data[H1functionposition].x[0] << ", " << data[H1functionposition].x[1] << ", "
                << data[H1functionposition].x[2] << '\n';
        sout << "DisplacementReconstructed = " << displacementreconstructed << "\n";
        sout << "DisplacementFem = " << displacementfem[0] << ", " << displacementfem[1] << ", "
                << displacementfem[2] << "\n";
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}


/**
 * @brief Returns the solution of the multiphysics simulation.
 * @param datavec [in] Data material vector
 * @param var [in] number of solution variables. See NSolutionVariables() method
 * @param Solout [out] The solution vector
 */
void TPZHDivErrorEstimateElasticityMaterial::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) {

    /**
     datavec[0] H1 mesh, uh_reconstructed for Mark reconstruction and Empty for H1 reconstruction
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


    TPZManVector<STATE, 2> u_exact(dim, 0.);
    u_exact.Fill(0);
    TPZFNMatrix<9, STATE> gradu(3, dim, 0.), fluxinv(3, 1);

    auto x = datavec[H1functionIdx].x;
    

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
                for (int j = 0; j < dim; j++) {
                    StressFem(i, j) = datavec[2].sol[0][j + i * 3];
                }
            }

            Solout.Resize(3);
            Solout[0] = StressFem(0, 0); //StressFemXX
            Solout[1] = StressFem(1, 1); //StressFemYY
            Solout[2] = StressFem(0, 1); //StressFemXY
        }
            break;


        case 41://StressExact
        {

            unsigned int nstate = fDimension*fDimension;
            TPZVec<STATE> u_exact(fDimension, 0.);
            TPZFMatrix<STATE> du_exact(fDimension, fDimension, 0.);
            if (this->fExactSol) {
                this->fExactSol(x, u_exact, du_exact);
            }
            TPZFMatrix<STATE> eps_exact(fDimension, fDimension, 0.);
            eps_exact(0, 0) = du_exact(0, 0);
            eps_exact(1, 0) = eps_exact(0, 1) = 0.5 * (du_exact(0, 1) + du_exact(1, 0));
            eps_exact(1, 1) = du_exact(1, 1);
            
            TPZVec<STATE> eps_exactV(nstate, 0.);
            ToVoigt(eps_exact, eps_exactV);

            TElasticityAtPoint elast(fE_const, fnu_const);
            if (TPZMixedElasticityND::fElasticity) {
                TPZManVector<STATE, 3> result(2);
                TPZFNMatrix<4, STATE> Dres(0, 0);
                fElasticity(x, result, Dres);
                REAL E = result[0];
                REAL nu = result[1];
                TElasticityAtPoint modify(E, nu);
                elast = modify;
            }

            TPZVec<STATE> sigma_exactV(nstate, 0.);
            ComputeStressVector(eps_exactV, sigma_exactV, elast);
            // std::cout << "duexact = " << du_exact << std::endl;
            // For 2D only
            Solout.Resize(3);
            Solout[0] = sigma_exactV[Exx]; //Sigma XX
            Solout[1] = sigma_exactV[Eyy]; //Sigma YY
            Solout[2] = sigma_exactV[Exy]; //Sigma XY
            return;
        }
            break;

        case 43://displacement Exact
        {

            TPZVec<STATE> u_exact(fDimension, 0.);
            TPZFMatrix<STATE> du_exact(fDimension, fDimension, 0.);
            if (this->fExactSol) {
                this->fExactSol(datavec[2].x, u_exact, du_exact);
            }
            //std::cout << "uexact = " << u_exact << std::endl;

            for (int idf = 0; idf < dim; idf++) {
                Solout[idf] =u_exact[idf];
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
        case 46://Stress reconstructed
            //eps(reconstructed displacement)
            {
                TElasticityAtPoint elast(fE_const, fnu_const);
                if (TPZMixedElasticityND::fElasticity) {
                    TPZManVector<STATE, 3> result(2);
                    TPZFNMatrix<4, STATE> Dres(0, 0);
                    fElasticity(x, result, Dres);
                    REAL E = result[0];
                    REAL nu = result[1];
                    TElasticityAtPoint modify(E, nu);
                    elast = modify;
                }
                TPZFNMatrix<9, STATE> eps_reconstructed(fDimension, fDimension, 0.);
                const auto &dudxreconstructed = datavec[H1functionIdx].dsol[0];
                const auto &axes = datavec[H1functionIdx].axes;
                TPZFNMatrix<6, STATE> du(3, 3);
                TPZAxesTools<STATE>::Axes2XYZ(dudxreconstructed, du, axes);
                eps_reconstructed(0, 0) = du(0, 0);
                eps_reconstructed(1, 0) = eps_reconstructed(0, 1) = 0.5 * (du(0, 1) + du(1, 0));
                eps_reconstructed(1, 1) = du(1, 1);
                TPZManVector<STATE, 9> sigma_reconstructedV(fDimension*fDimension, 0.), eps_reconstructedV(fDimension*fDimension, 0.);
                ToVoigt(eps_reconstructed, eps_reconstructedV);
                ComputeStressVector(eps_reconstructedV, sigma_reconstructedV, elast);
                Solout[0] = sigma_reconstructedV[0];
                Solout[1] = sigma_reconstructedV[1];
                Solout[2] = sigma_reconstructedV[2];
            }
            break;
            case 47://eps(u_exact)
            {
                
                unsigned int nstate = fDimension*fDimension;
                TPZVec<STATE> u_exact(fDimension, 0.);
                TPZFMatrix<STATE> du_exact(fDimension, fDimension, 0.);
                if (this->fExactSol) {
                    this->fExactSol(x, u_exact, du_exact);
                }
                TPZFMatrix<STATE> eps_exact(fDimension, fDimension, 0.);
                TPZFNMatrix<6, STATE> du(3, 3);
                const auto &axes = datavec[H1functionIdx].axes;
                TPZAxesTools<STATE>::Axes2XYZ(du_exact, du, axes);
                
                eps_exact(0, 0) = du(0, 0);
                eps_exact(1, 0) = eps_exact(0, 1) = 0.5 * (du(0, 1) + du(1, 0));
                eps_exact(1, 1) = du(1, 1);
                
                TPZVec<STATE> eps_exactV(nstate, 0.);
                Solout.Resize(3, 0.);
                
                Solout[Exx] = eps_exact(0, 0);
                Solout[Exy] = eps_exact(1, 0);
                Solout[Eyy] = eps_exact(1, 1);
                
                
            }
                break;
                
            case 48://eps(u_rec)
            {
                
                TElasticityAtPoint elast(fE_const, fnu_const);
                if (TPZMixedElasticityND::fElasticity) {
                    TPZManVector<STATE, 3> result(2);
                    TPZFNMatrix<4, STATE> Dres(0, 0);
                    fElasticity(x, result, Dres);
                    REAL E = result[0];
                    REAL nu = result[1];
                    TElasticityAtPoint modify(E, nu);
                    elast = modify;
                }
                TPZFNMatrix<9, STATE> eps_reconstructed(fDimension, fDimension, 0.);
                const auto &dudxreconstructed = datavec[H1functionIdx].dsol[0];
                const auto &axes = datavec[H1functionIdx].axes;
                TPZFNMatrix<6, STATE> du(3, 3);
                TPZAxesTools<STATE>::Axes2XYZ(dudxreconstructed, du, axes);
                eps_reconstructed(0, 0) = du(0, 0);
                eps_reconstructed(1, 0) = eps_reconstructed(0, 1) = 0.5 * (du(0, 1) + du(1, 0));
                eps_reconstructed(1, 1) = du(1, 1);
                
                Solout.Resize(3, 0.);
                Solout[Exx] = eps_reconstructed(0, 0);
                Solout[Exy] = eps_reconstructed(1, 0);
                Solout[Eyy] = eps_reconstructed(1, 1);
                
                
            }
            break;
        case 49:
        {
            for (int i = 0; i < dim; i++) {
                Solout[i] = datavec[5].sol[0][i];
            }
        }
            break;
        case 50:
        {
            std::cout<<"to be implemented"<<std::endl;
        }
        default:
            DebugStop();
    }
}

void TPZHDivErrorEstimateElasticityMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    /**
   
     datavec[0] H1 mesh, local uk/grad v for Mark reconstruction and Empty for H1 reconstruction
     datavec[1] L2 mesh, restriction of local uk
     datavec[2] Hdiv mesh, sigmakE
     datavec[3] L2 mesh, ukE
     
     Implement the matrix
     (eps (u_rec),eps(v)) = (C sigma_fem,eps(v))
     
     int eps(u_rec): eps(v) dx = int C sigma_fem : eps(v) dx
   
     **/

    int dim = this->fDimension;
    int matdim = dim*dim;
    TPZManVector<REAL, 3> x = datavec[2].x;

    TElasticityAtPoint elast(fE_const, fnu_const);
    if (TPZMixedElasticityND::fElasticity) {
        TPZManVector<STATE, 3> result(2);
        TPZFNMatrix<4, STATE> Dres(0, 0);
        fElasticity(x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        TElasticityAtPoint modify(E, nu);
        elast = modify;
    }
    TPZFNMatrix<9, STATE> stressfem(dim, dim, 0.);
    for (unsigned int i = 0; i < dim; i++) {
        for (unsigned int j = 0; j < dim; j++) {
            stressfem(i, j) = datavec[2].sol[0][j + i * dim];
        }
    }
 
    //std::cout<<"stressfem= "<<stressfem<<std::endl;
    
    //compute C(sigma)
    TPZManVector<STATE, 9> stress_femV(matdim, 0.);
    ToVoigt(stressfem, stress_femV);

    TPZManVector<STATE, 9> Csigma_femV(matdim, 0.);
    ComputeDeformationVector(stress_femV, Csigma_femV, elast);
    
    TPZFNMatrix<9, STATE> Csigma_fem(dim, dim, 0.);
    FromVoigt(Csigma_femV, Csigma_fem);
    
    
    
    //
    TPZFNMatrix<9, STATE> rotfem(dim, dim, 0.);
    rotfem(0,1) = datavec[4].sol[0][0];
    rotfem(1,0) = (-1.)*datavec[4].sol[0][0];
    
  //  std::cout<<"rot= "<<rotfem<<std::endl;
    
    
    TPZFNMatrix<9, STATE> rhs_term(dim,dim,0.);
    for (unsigned int i = 0; i < dim; i++) {
        for (unsigned int j = 0; j < dim; j++) {
            rhs_term(i, j) = stressfem(i, j)+rotfem(i,j);
        }
    }
   // std::cout<<"rhs_term= "<<rhs_term<<std::endl;

    
//    //defining test functions
//    // Setting the phis
    int H1functionposition = 0;
    H1functionposition = FirstNonNullApproxSpaceIndex(datavec);
    TPZFMatrix<REAL> &phiuk = datavec[H1functionposition].phi;
    TPZFMatrix<REAL> &dphiukaxes = datavec[H1functionposition].dphix;
    TPZFNMatrix<9, REAL> dphiuk(3, dphiukaxes.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiukaxes, dphiuk, datavec[H1functionposition].axes);
    
    const auto nphiuk = phiuk.Rows();
    const auto &axes = datavec[2].axes;
    TPZVec<TPZManVector<STATE,9>> eps_phiukV(nphiuk, TPZManVector<STATE, 9>(matdim, 0.));
    for(int in = 0; in < nphiuk; in++ ) {
	TPZFNMatrix<4,STATE> du(2,2);
        du(0,0) = dphiuk(0,in)*axes(0,0)+dphiuk(1,in)*axes(1,0);//dvx
        du(1,0) = dphiuk(0,in)*axes(0,1)+dphiuk(1,in)*axes(1,1);//dvy
        


//        TPZFNMatrix<9, STATE> eps_phiuk(dim, dim, 0.);
//        eps_phiuk(0, 0) = du(0, 0);
//        eps_phiuk(1, 0) = eps_phiuk(0, 1) = 0.5 * (du(0, 1) + du(1, 0));
//        eps_phiuk(1, 1) = du(1, 1);
//
//        ToVoigt(eps_phiuk, eps_phiukV[in]);

        ef(2*in, 0) +=weight * (Csigma_fem[Exx]*du(0,0) + 0.5*(Csigma_fem[Exy]+Csigma_fem[Eyx])*du(1,0));
        ef(2*in+1, 0) += weight * (Csigma_fem[Eyy]*du(1,0) + 0.5*(Csigma_fem[Exy]+Csigma_fem[Eyx])*du(0,0));

        //ef(2*in, 0) += weight * (Csigma_fem[Exx]*du(0,0) + (Csigma_fem[Exy] + rotfem(0, 1))*du(1,0));
        //ef(2*in+1, 0) += weight * (Csigma_fem[Eyy]*du(1,0) + (Csigma_fem[Eyx] + rotfem(1, 0))*du(0,0));

        for(int jn = 0; jn < nphiuk; jn++ ) {
            du(0,1) = dphiuk(0,jn)*axes(0,0)+dphiuk(1,jn)*axes(1,0);//dux
            du(1,1) = dphiuk(0,jn)*axes(0,1)+dphiuk(1,jn)*axes(1,1);//duy
        
            ek(2*in,2*jn) += weight * (du(0,1)*du(0,0) + 0.5*du(1,1)*du(1,0));
            ek(2*in,2*jn+1) += weight * 0.5*du(0,1)*du(1,0);
            ek(2*in+1,2*jn) += weight * 0.5*du(1,1)*du(0,0);
            ek(2*in+1,2*jn+1) += weight * (du(1,1)*du(1,0) + 0.5*du(0,1)*du(0,0));

//            ek(2*in,2*jn) += weight * (du(0,1)*du(0,0) + du(1,1)*du(1,0));
//            ek(2*in+1,2*jn+1) += weight * (du(0,1)*du(0,0) + du(1,1)*du(1,0));


        }

    }
    
}

// void TPZHDivErrorEstimateElasticityMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
//     /**
   
//      datavec[0] H1 mesh, local uk/grad v for Mark reconstruction and Empty for H1 reconstruction
//      datavec[1] L2 mesh, restriction of local uk
//      datavec[2] Hdiv mesh, sigmakE
//      datavec[3] L2 mesh, ukE
     
//      Implement the matrix
//      (eps (u_rec),eps(v)) = (C sigma_fem,eps(v))
     
//      int eps(u_rec): eps(v) dx = int C sigma_fem : eps(v) dx
   
//      **/

//     int dim = this->fDimension;
//     int matdim = dim*dim;
//     TPZManVector<REAL, 3> x = datavec[2].x;

//     TElasticityAtPoint elast(fE_const, fnu_const);
//     if (TPZMixedElasticityND::fElasticity) {
//         TPZManVector<STATE, 3> result(2);
//         TPZFNMatrix<4, STATE> Dres(0, 0);
//         fElasticity(x, result, Dres);
//         REAL E = result[0];
//         REAL nu = result[1];
//         TElasticityAtPoint modify(E, nu);
//         elast = modify;
//     }
//     TPZFNMatrix<9, STATE> stressfem(dim, dim, 0.);
//     for (unsigned int i = 0; i < dim; i++) {
//         for (unsigned int j = 0; j < dim; j++) {
//             stressfem(i, j) = datavec[2].sol[0][j + i * dim];
//         }
//     }
//     //compute C(sigma)
//     TPZManVector<STATE, 9> stress_femV(matdim, 0.);
//     ToVoigt(stressfem, stress_femV);

//     TPZManVector<STATE, 9> Csigma_femV(matdim, 0.);
//     ComputeDeformationVector(stress_femV, Csigma_femV, elast);
    
//     TPZFNMatrix<9, STATE> Csigma_fem(dim, dim, 0.);
//     FromVoigt(Csigma_femV, Csigma_fem);

    
// //    //defining test functions
// //    // Setting the phis
//     int H1functionposition = 0;
//     H1functionposition = FirstNonNullApproxSpaceIndex(datavec);
//     TPZFMatrix<REAL> &phiuk = datavec[H1functionposition].phi;
//     TPZFMatrix<REAL> &dphiukaxes = datavec[H1functionposition].dphix;
//     TPZFNMatrix<9, REAL> dphiuk(3, dphiukaxes.Cols());
//     TPZAxesTools<REAL>::Axes2XYZ(dphiukaxes, dphiuk, datavec[H1functionposition].axes);
    
//     const auto nphiuk = phiuk.Rows();
//     const auto &axes = datavec[2].axes;
//     TPZVec<TPZManVector<STATE,9>> eps_phiukV(nphiuk, TPZManVector<STATE, 9>(matdim, 0.));
//     for(int in = 0; in < nphiuk; in++ ) {
// 	TPZFNMatrix<4,STATE> du(2,2);
//         du(0,0) = dphiuk(0,in)*axes(0,0)+dphiuk(1,in)*axes(1,0);//dvx
//         du(1,0) = dphiuk(0,in)*axes(0,1)+dphiuk(1,in)*axes(1,1);//dvy
        
// //        TPZFNMatrix<9, STATE> eps_phiuk(dim, dim, 0.);
// //        eps_phiuk(0, 0) = du(0, 0);
// //        eps_phiuk(1, 0) = eps_phiuk(0, 1) = 0.5 * (du(0, 1) + du(1, 0));
// //        eps_phiuk(1, 1) = du(1, 1);
// //
// //        ToVoigt(eps_phiuk, eps_phiukV[in]);
//         ef(2*in, 0) += weight * (Csigma_fem[Exx]*du(0,0) + 0.5*(Csigma_fem[Exy]+Csigma_fem[Eyx])*du(1,0));
//         ef(2*in+1, 0) += weight * (Csigma_fem[Eyy]*du(1,0) + 0.5*(Csigma_fem[Exy]+Csigma_fem[Eyx])*du(0,0));
//         for(int jn = 0; jn < nphiuk; jn++ ) {
//             du(0,1) = dphiuk(0,jn)*axes(0,0)+dphiuk(1,jn)*axes(1,0);//dux
//             du(1,1) = dphiuk(0,jn)*axes(0,1)+dphiuk(1,jn)*axes(1,1);//duy
            
//             ek(2*in,2*jn) += weight * (du(0,1)*du(0,0) + 0.5*du(1,1)*du(1,0));
//             ek(2*in,2*jn+1) += weight * 0.5*du(0,1)*du(1,0);
//             ek(2*in+1,2*jn) += weight * 0.5*du(1,1)*du(0,0);
//             ek(2*in+1,2*jn+1) += weight * (du(1,1)*du(1,0) + 0.5*du(0,1)*du(0,0));
//         }        
//     }
    
// }

void TPZHDivErrorEstimateElasticityMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
        TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {

    /*
     Add Dirichlet boundary condition for local problem
     ek+= <w,Km s_i>
     ef+= <w,Km*u_d - g + sigma_i.n>
     */
    int H1functionposition = FirstNonNullApproxSpaceIndex(datavec);
    int dim = datavec[H1functionposition].axes.Rows();

    TPZFMatrix<REAL> &phi_i = datavec[H1functionposition].phi;
    int nphi_i = phi_i.Rows();

   // std::cout<<" aplicando BC "<<std::endl;
    
    const TPZVec<REAL> u_D = bc.Val2();
    TPZFNMatrix<9, STATE> g = bc.Val1();

    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9, STATE> gradu(dim, 1);
        bc.ForcingFunctionBC()(datavec[H1functionposition].x, res, gradu);

        u_D[0]=res[0];
        u_D[1]=res[1];


    } else {
        // u_D is usually stored in val2(0, 0)
        u_D[0] = bc.Val2()[0];
        u_D[1] = bc.Val2()[1];
    }

    int nstate = 2;

    
    switch (bc.Type()) {


        case (0):
        {
            for(int in = 0 ; in < nphi_i; in++) {
   
                ef(2*in,0)   += fBigNumber * u_D[0] * phi_i(in,0) * weight;        // forced v2 displacement
                ef(2*in+1,0) += fBigNumber * u_D[1] * phi_i(in,0) * weight;        // forced v2 displacement
                for (int jn = 0 ; jn < nphi_i; jn++)
                {
                    ek(2*in,2*jn)     += fBigNumber * phi_i(in,0) *phi_i(jn,0) * weight;
                    ek(2*in+1,2*jn+1) += fBigNumber * phi_i(in,0) *phi_i(jn,0) * weight;
                }
            }
        }

//        default:
//        {
//             std::cout << " This material not implement BC Type " << bc.Type()<< std::endl;
//            break;
//        }
    }
}

int TPZHDivErrorEstimateElasticityMaterial::VariableIndex(const std::string &name)const {
    if (name == "StressFem") return 40;

    if (name == "StressExact") return 41;

    if (name == "DisplacementFem") return 42;
    if (name == "DisplacementExact") return 43;
    if (name == "DisplacementReconstructed") return 44;
    if (name == "POrder") return 45;
    if (name == "DisplacementErrorExact") return 100;
    
    if (name == "DisplacementErrorEstimate") return 101;
    if (name == "EnergyErrorExact") return 102;
    if (name == "EnergyErrorEstimate") return 103;
    if (name == "ResidualError") return 104;
    if (name == "EnergyH1Error") return 109;
    
    if (name == "DisplacementEffectivityIndex") return 110;
    if (name == "EnergyEffectivityIndex") return 111;
    if (name == "LocalErrorIndicator") return 112;
    
    
    
    if (name == "StressReconstructed") return 46;
    if (name == "EpsExact") return 47;
    if (name == "EpsRec") return 48;
    if (name == "DisplacementFemH1") return 49;
    if (name == "ErrorDispFemH1Rec") return 50;
 


    return -1;
}

int TPZHDivErrorEstimateElasticityMaterial::NSolutionVariables(int var) const
{
    switch (var) {
        case 40:
        case 41:
        case 42:
        case 43:
        case 44:
        case 46:
        case 47:
        case 48:
        case 49:
            return 3;
            break;
            
        case 100:
        case 101:
        case 102:
        case 103:
        case 104:
        case 107:
        case 108:
        case 109:
        case 110:
        case 111:
        case 112:
        case 113:
            
        case 45:
        case 50:
            return 1;
            break;
        default:
            DebugStop();
            break;
    }
    return 0;
}

void TPZHDivErrorEstimateElasticityMaterial::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const {

    
        //fem solution for flux and potential
        datavec[0].SetAllRequirements(false);
        datavec[0].fNeedsSol = true;
        datavec[0].fNeedsNormal = true;

        datavec[1].SetAllRequirements(false);
        datavec[1].fNeedsSol = true;
        datavec[1].fNeedsNormal = true;

        datavec[2].SetAllRequirements(false);
        datavec[2].fNeedsSol = true;
        datavec[2].fNeedsNormal = true;

        datavec[4].SetAllRequirements(false);
        datavec[4].fNeedsSol = true;
    

}
