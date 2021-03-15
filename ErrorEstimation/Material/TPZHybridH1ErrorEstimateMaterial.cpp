//
// Created by victor on 07/10/2020.
//

#include "TPZHybridH1ErrorEstimateMaterial.h"
#include "pzaxestools.h"
#include "TPZAnalyticSolution.h"
#include "pzbndcond.h"


#ifdef LOG4CXX
    static LoggerPtr loggerF(Logger::getLogger("DebuggingF"));
#endif


TPZHybridH1ErrorEstimateMaterial::TPZHybridH1ErrorEstimateMaterial(int matid, int dim) : TPZMixedPoisson(matid, dim)
{
}

TPZHybridH1ErrorEstimateMaterial::TPZHybridH1ErrorEstimateMaterial() : TPZMixedPoisson()
{

}

TPZHybridH1ErrorEstimateMaterial::TPZHybridH1ErrorEstimateMaterial(const TPZHybridH1ErrorEstimateMaterial &copy) : TPZMixedPoisson(copy)
{
}

TPZHybridH1ErrorEstimateMaterial::TPZHybridH1ErrorEstimateMaterial(const TPZMixedPoisson &copy) : TPZMixedPoisson(copy)
{

}

TPZHybridH1ErrorEstimateMaterial::TPZHybridH1ErrorEstimateMaterial(TPZMatLaplacianHybrid &matlaplacian)
{
    this->SetId(matlaplacian.Id());
    this->SetDimension(matlaplacian.Dimension());

    TPZFNMatrix<9,STATE> K,invK;
    matlaplacian.GetPermeability(K);
    matlaplacian.GetInvPermeability(invK);
    this->SetPermeabilityTensor(K,invK);

    if (matlaplacian.HasForcingFunction()) {
        this->SetForcingFunctionExact(matlaplacian.ForcingFunctionExact());
        this->SetForcingFunction(matlaplacian.ForcingFunction());
    }
}

TPZHybridH1ErrorEstimateMaterial::~TPZHybridH1ErrorEstimateMaterial()
{

}

TPZHybridH1ErrorEstimateMaterial &TPZHybridH1ErrorEstimateMaterial::operator=(const TPZHybridH1ErrorEstimateMaterial &copy)
{
    TPZMixedPoisson::operator=(copy);
    return *this;
}

void TPZHybridH1ErrorEstimateMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    /**
     find (sigma_h,p_h,g_h) satisfying
        (invK.sigma_h,v_h) - (p_h,div v_h) = 0;
        -(div sigma_h,w_h) + (g_h,w_h) = (f,v);
        (p_h,t_h) = (u_avg,t_h),
    for all (v_h,w_h,t_h)
     **/

    if(datavec[0].fActiveApproxSpace){
        TPZFNMatrix<9,REAL> PermTensor;
        TPZFNMatrix<9,REAL> InvPermTensor;
        int dim = datavec[0].axes.Rows();
        TPZFMatrix<STATE> gradSol(3, 1,0);

        GetPermeabilities(datavec[0].x, PermTensor, InvPermTensor);

        // Setting the phis
        TPZFMatrix<REAL> &phiQ = datavec[0].phi;
        TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
        TPZFMatrix<REAL> &phip = datavec[1].phi;

        STATE force = ff;
        if(fForcingFunction) {
            TPZManVector<STATE> res(1);
            fForcingFunction->Execute(datavec[1].x,res);
            force = res[0];
        }


        TPZFMatrix<REAL> &phiuk = datavec[0].phi;
        TPZFMatrix<REAL> &dphiukaxes = datavec[0].dphix;
        TPZFNMatrix<9,REAL> dphiuk(3,dphiukaxes.Cols());
        TPZAxesTools<REAL>::Axes2XYZ(dphiukaxes, dphiuk, datavec[0].axes);

        int phrq,phrp;
        phrq = datavec[0].fVecShapeIndex.NElements();
        phrp = phip.Rows();

        TPZFNMatrix<15, STATE> &dsolpaxes = datavec[1].dsol[0];
        TPZFNMatrix<9, REAL> dsolp(2, dphiukaxes.Cols());
        TPZAxesTools<REAL>::Axes2XYZ(dsolpaxes, dsolp, datavec[1].axes);

        for (int ip = 0; ip < dim; ip++) {
            gradSol(ip, 0) = dsolp.Get(ip, 0);
        }

        for(int iq=0; iq<phrq; iq++)
        {
            //ef(iq, 0) += 0.;
            int ivecind = datavec[0].fVecShapeIndex[iq].first;
            int ishapeind = datavec[0].fVecShapeIndex[iq].second;
            TPZFNMatrix<3,REAL> ivec(3,1,0.);
            for(int id=0; id<3; id++){
                ivec(id,0) = datavec[0].fDeformedDirections(id,ivecind);
            }

            TPZFNMatrix<3,REAL> ivecZ(3,1,0.);
            TPZFNMatrix<3,REAL> jvecZ(3,1,0.);
            for (int jq=0; jq<phrq; jq++)
            {

                TPZFNMatrix<3,REAL> jvec(3,1,0.);
                int jvecind = datavec[0].fVecShapeIndex[jq].first;
                int jshapeind = datavec[0].fVecShapeIndex[jq].second;

                for(int id=0; id<3; id++){
                    jvec(id,0) = datavec[0].fDeformedDirections(id,jvecind);
                }

                //dot product between Kinv[u]v
                jvecZ.Zero();
                for(int id=0; id<3; id++){
                    for(int jd=0; jd<3; jd++){
                        jvecZ(id,0) += InvPermTensor(id,jd)*jvec(jd,0);
                    }
                }
                /// (invK.sigma_h,v)_K
                REAL prod1 = ivec(0,0)*jvecZ(0,0) + ivec(1,0)*jvecZ(1,0) + ivec(2,0)*jvecZ(2,0);
                ek(iq,jq) += fvisc*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod1;
            }
        }

        // Coupling terms between flux and pressure. Matrix B
        for(int iq=0; iq<phrq; iq++)
        {
            int ivecind = datavec[0].fVecShapeIndex[iq].first;
            int ishapeind = datavec[0].fVecShapeIndex[iq].second;

            TPZFNMatrix<3,REAL> ivec(3,1,0.);
            for(int id=0; id<3; id++){
                ivec(id,0) = datavec[0].fDeformedDirections(id,ivecind);
                //ivec(1,0) = datavec[0].fDeformedDirections(1,ivecind);
                //ivec(2,0) = datavec[0].fDeformedDirections(2,ivecind);
            }
            TPZFNMatrix<3,REAL> axesvec(3,1,0.);
            datavec[0].axes.Multiply(ivec,axesvec);

            REAL divwq = 0.;
            for(int iloc=0; iloc<fDim; iloc++)
            {
                divwq += axesvec(iloc,0)*dphiQ(iloc,ishapeind);
            }
            for (int jp=0; jp<phrp; jp++) {
                /// - (p_h,div v_h)
                REAL fact = (-1.)*weight*phip(jp,0)*divwq;
                // Matrix B
                ek(iq, phrq+jp) += fact;

                // Matrix B^T
                ek(phrq+jp,iq) += fact;

            }
        }
        //termo fonte referente a equacao da pressao
        for(int ip=0; ip<phrp; ip++){
            /// (f,v)
            ef(phrq+ip,0) += (-1.)*weight*force*phip(ip,0);
        }

        for(int ip=0; ip<phrp; ip++){
            /// (g_h,w_h)
            ek(phrq+ip,phrq+phrp) += phip(ip,0)*weight;
            ek(phrq+phrp,phrq+ip) += phip(ip,0)*weight;
        }
        ef(phrq+phrp,0) += weight*datavec[3].sol[0][0];
    }

    /**
     Implement

     1) (K grad s_h, grad v)_K = (f,v):
            int_K phi_i.phi_j dx = int_K f phi_i  dx;
     **/
    else {

        int H1functionposition = 1;

        int dim = datavec[H1functionposition].axes.Rows();
        //defining test functions
        // Setting the phis
        TPZFMatrix<REAL> &phiuk = datavec[H1functionposition].phi;
        TPZFMatrix<REAL> &dphiukaxes = datavec[H1functionposition].dphix; //(2xnphiuk)
        TPZFNMatrix<9, REAL> dphiuk(2, dphiukaxes.Cols());
        TPZAxesTools<REAL>::Axes2XYZ(dphiukaxes, dphiuk, datavec[H1functionposition].axes); //(3xnphiuk)
        TPZFNMatrix<15, STATE> &dsolaxes = datavec[H1functionposition].dsol[0];
        TPZFNMatrix<9, REAL> dsol(2, dphiukaxes.Cols());
        TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, dsol, datavec[H1functionposition].axes);

        int nphiuk = phiuk.Rows();

        TPZFMatrix<STATE> gradSol(dim, 1,0), kGradSol(dim, 1,0), solukfem(1, 1,0);

        TPZVec<STATE> divsigma(1);
        divsigma[0] = 0;

        if (!this->fForcingFunction) DebugStop();
        this->fForcingFunction->Execute(datavec[H1functionposition].x, divsigma);

        TPZFNMatrix<9, REAL> PermTensor;
        TPZFNMatrix<9, REAL> InvPermTensor;

        GetPermeabilities(datavec[1].x, PermTensor, InvPermTensor);

        //potetial fem
        /*solukfem(0, 0) = datavec[3].sol[0][0];
        for (int ip = 0; ip < dim; ip++) {
           gradSol(ip, 0) = dsol.Get(ip, 0);
        }
        for (int id = 0; id < dim; id++) {
            for (int jd = 0; jd < dim; jd++) {
                kGradSol(id, 0) += PermTensor(id, jd) * gradSol(jd, 0);
            }
        }*/

        TPZFMatrix<STATE> kgraduk(dim, nphiuk, 0.);



        for (int irow = 0; irow < nphiuk; irow++) {

            for (int id = 0; id < dim; id++) {
                for (int jd = 0; jd < dim; jd++) {
                    kgraduk(id, irow) += PermTensor(id, jd) * dphiuk(jd, irow);
                }
                /// ... = (grad u_h, grad v_h)
                //ef(irow,0) +=weight*dphiuk(id,irow)*kGradSol(id,0);
            }
            ///... = (f , v_h)
            ef(irow, 0) += weight * phiuk(irow, 0) * divsigma[0];

            //matrix Sk= int_{K} K graduk.gradv
            for (int jcol = 0; jcol < nphiuk; jcol++) {

                for (int jd = 0; jd < dim; jd++) {
                    ek(irow, jcol) += weight * kgraduk(jd, irow) * dphiuk(jd, jcol);
                }

            }
        }
#ifdef LOG4CXX
        if (loggerF->isDebugEnabled()) {
            std::stringstream ss;
            ss << "X = [" << datavec[1].x[0] << "," << datavec[1].x[1] << "," << datavec[1].x[2] << "]\n";
            ss << "f = " << divsigma[0] << "\n";
            ef.Print("EF", ss, EMathematicaInput);
            LOGPZ_DEBUG(loggerF, ss.str())
        }
#endif
    }
}

void TPZHybridH1ErrorEstimateMaterial::ContributeBC(
        TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
        TPZFMatrix<STATE> &ef, TPZBndCond &bc) {

    if(!datavec[0].fActiveApproxSpace) {
        /*
         Add Robin boundary condition for local problem
         ek+= <w,Km s_i>
         ef+= <w,Km*u_d - g + sigma_i.n>
         */
        int H1functionposition = 1;
        int dim = datavec[H1functionposition].axes.Rows();

        TPZFMatrix<REAL> &phi_i = datavec[H1functionposition].phi;
        int nphi_i = phi_i.Rows();

        TPZFMatrix<STATE> solsigmafem(3, 1);
        solsigmafem.Zero();

        TPZManVector<REAL, 3> normal = datavec[2].normal;

        REAL normalsigma = 0.;
        normalsigma = datavec[2].sol[0][0];

        REAL u_D = 0.;
        REAL g = 0.;
        //REAL normflux = 0.;



        if (bc.HasForcingFunction()) {
            TPZManVector<STATE> res(3);
            TPZFNMatrix<9, STATE> gradu(dim, 1);
            bc.ForcingFunction()->Execute(datavec[H1functionposition].x, res, gradu);
            u_D = res[0];
            g = normalsigma;


            /*TPZFNMatrix<9,REAL> PermTensor, InvPermTensor;
            GetPermeabilities(datavec[0].x, PermTensor, InvPermTensor);


            for(int i=0; i<3; i++)
            {
                for(int j=0; j<3; j++)
                {

                    normflux += datavec[2].normal[i]*PermTensor(i,j)*gradu(j,0);
                }
            }*/





        } else {
            // usualmente updatebc coloca o valor exato no val2
            u_D = bc.Val2()(0, 0);
        }


        switch (bc.Type()) {
            case (4): {

                REAL Km = bc.Val1()(0, 0); // Km
                REAL InvKm = 1. / Km;
                // std::cout<< " g "<< g<< " normalsigma "<<normalsigma<<"\n";
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
                    ef(iq, 0) += gBigNumber * u_D * phi_i(iq, 0) * weight;
                    for (int jq = 0; jq < nphi_i; jq++) {
                        ek(iq, jq) += gBigNumber * weight * phi_i(iq, 0) * phi_i(jq, 0);
                    }
                }

                break;
            }

            default: {

                std::cout << " This material not implement BC Type " << bc.Type() << std::endl;
                break;
            }
        }
    }
}


void TPZHybridH1ErrorEstimateMaterial::FillDataRequirements(TPZVec<TPZMaterialData > &datavec) {

    if(!datavec[0].fActiveApproxSpace) {
        //fem solution for flux and potential
        datavec[0].SetAllRequirements(false);
        datavec[0].fNeedsSol = true;
        datavec[0].fNeedsNormal = true;

        datavec[2].SetAllRequirements(false);
        datavec[2].fNeedsSol = true;
        datavec[2].fNeedsNormal = true;

        datavec[3].SetAllRequirements(false);
        datavec[3].fNeedsSol = true;

        datavec[4].SetAllRequirements(false);
        datavec[4].fNeedsSol = true;
    }
    else{
        datavec[3].SetAllRequirements(false);
        datavec[3].fNeedsSol = true;
        datavec[3].fNeedsNormal = true;
    }

}

void TPZHybridH1ErrorEstimateMaterial::FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData > &datavec){
    if(!datavec[0].fActiveApproxSpace) {
        datavec[0].SetAllRequirements(false);
        datavec[0].fNeedsSol = true;
        datavec[0].fNeedsNormal = true;

        datavec[2].SetAllRequirements(false);
        datavec[2].fNeedsSol = true;
        datavec[2].fNeedsNormal = true;

        datavec[4].SetAllRequirements(false);
    }
    else{
        datavec[3].SetAllRequirements(false);
        datavec[3].fNeedsSol = true;
        datavec[3].fNeedsNormal = true;
    }
}

void TPZHybridH1ErrorEstimateMaterial::UpdateBCValues(TPZVec<TPZMaterialData> & datavec){
    DebugStop();

}

void TPZHybridH1ErrorEstimateMaterial::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{
    /**
     datavec[0] H1 mesh, uh_reconstructed
     datavec[1] L2 mesh,
     datavec[2] Hdiv fem mesh, sigma_h
     datavec[3] L2 mesh fem, u_h

      error[0] - error computed with exact pressure
      error[1] - error computed with reconstructed pressure
      error[2] - energy error computed with exact solution
      error[3] - energy error computed with reconstructed flux
      error[4] - energy error computed with reconstructed potential
      error[5] - oscilatory data error

     **/

    errors.Resize(NEvalErrors());
    errors.Fill(0.0);


    STATE divsigmarec, pressurefem, pressurereconstructed, forceProj;

    TPZFNMatrix<3,REAL> fluxreconstructed(3,1), fluxreconstructed2(3,1), gradreconstructed(3,1);
    TPZFMatrix<REAL> gradfem,fluxfem(3,1);

    gradfem=data[3].dsol[0];

    gradfem.Resize(3,1);


    int H1functionposition = 1;

    TPZVec<STATE> divsigma(1);

    if(this->fForcingFunctionExact){

        this->fForcingFunctionExact->Execute(data[H1functionposition].x,u_exact,du_exact);

        this->fForcingFunction->Execute(data[H1functionposition].x,divsigma);
    }

    REAL residual = 0.,altResidual = 0.;
    divsigmarec= data[0].divsol[0][0];
    residual = (divsigma[0] - divsigmarec)*(divsigma[0] - divsigmarec);

    forceProj = data[4].sol[0][0];
    altResidual = forceProj - divsigma[0];

    pressurereconstructed = data[H1functionposition].sol[0][0];


    pressurefem = data[3].sol[0][0];

    TPZFNMatrix<9,REAL> PermTensor;
    TPZFNMatrix<9,REAL> InvPermTensor;

    GetPermeabilities(data[1].x, PermTensor, InvPermTensor);

    PermTensor.Resize(3,3);
    InvPermTensor.Resize(3,3);


    TPZFNMatrix<3,REAL> fluxexactneg;

    {
        TPZFNMatrix<9,REAL> gradpressure(3,1);
        for (int i=0; i<3; i++) {
            gradpressure(i,0) = (-1.)*du_exact[i];
            gradfem(i,0) = (-1.)* gradfem(i,0);
        }
        PermTensor.Multiply(gradpressure,fluxexactneg);
    }
    PermTensor.Multiply(gradfem,fluxfem);

    TPZFMatrix<REAL> &dsolaxes = data[H1functionposition].dsol[0];
    TPZFNMatrix<9,REAL> fluxrec(dsolaxes.Rows(),0);


    for(int ip = 0 ; ip < 3 ; ip++){
        fluxreconstructed(ip,0) = data[0].sol[0][ip];
    }

    TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, fluxrec, data[H1functionposition].axes);
    for(int id=0 ; id<3; id++) {
        fluxreconstructed2(id,0) = (-1.)*fluxrec(id,0);
    }
    PermTensor.Multiply(fluxreconstructed2,gradreconstructed);

    //data[H1functionposition].axes.Print(std::cout);
    //dsolaxes.Print(std::cout);
    //fluxrec.Print(std::cout);


    REAL innerexact = 0.;
    REAL innerestimate = 0.;
    REAL gradinnerestimate = 0.;
    REAL npz =0.;



#ifdef PZDEBUG2
    std::cout<<"flux fem "<<fluxfem<<std::endl;
    std::cout<<"flux reconst "<<fluxreconstructed<<std::endl;
    std::cout<<"-------"<<std::endl;
#endif



    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            innerexact += (fluxfem[i]-fluxexactneg(i,0))*InvPermTensor(i,j)*(fluxfem[j]-fluxexactneg(j,0));//Pq esta somando: o fluxo fem esta + e o exato -
            innerestimate += (fluxfem[i]-fluxreconstructed[i])*InvPermTensor(i,j)*(fluxfem[j]-fluxreconstructed[j]);
            gradinnerestimate += (fluxfem[i]-gradreconstructed[i])*InvPermTensor(i,j)*(fluxfem[j]-gradreconstructed[j]);
            npz += (gradreconstructed[i]-fluxreconstructed[i])*InvPermTensor(i,j)*(gradreconstructed[i]-fluxreconstructed[i]);
        }
    }

#ifdef PZDEBUG2
    std::cout<<"potential fem "<<pressurefem<<std::endl;
    std::cout<<"potential reconst "<<pressurereconstructed<<std::endl;
    std::cout<<"-------"<<std::endl;
#endif
    errors[0] = (pressurefem-u_exact[0])*(pressurefem-u_exact[0]);//exact error pressure
    errors[1] = (pressurefem-pressurereconstructed)*(pressurefem-pressurereconstructed);//error pressure reconstructed
    errors[2] = innerexact;//error flux exact
    errors[3] = gradinnerestimate; // NFC: ||grad(u_h-s_h)||
    errors[4] = residual; //||f - Proj_divsigma||
    errors[5] = innerestimate;//NF: ||grad(u_h)+sigma_h)||
    errors[6] = altResidual*altResidual;
}


int TPZHybridH1ErrorEstimateMaterial::VariableIndex(const std::string &name)
{
    if(name == "FluxFem") return 40;
    if(name == "FluxExact") return 42;
    if(name == "PressureFem") return 43;
    if(name == "PressureReconstructed") return 44;
    if(name == "FluxSigmaReconstructed") return 39;
    if(name == "FluxReconstructed") return 41;
    if(name == "PressureExact") return 45;
    if(name == "PressureErrorExact") return 100;
    if(name == "PressureErrorEstimate") return 101;
    if(name == "EnergyErrorExact") return 102;
    if(name == "EnergyErrorEstimate") return 103;
    if(name == "ResidualError") return 104;
    if(name == "PressureEffectivityIndex") return 107;
    if(name == "EnergyEffectivityIndex") return 108;
    if(name == "POrder") return 46;

    return -1;
}


int TPZHybridH1ErrorEstimateMaterial::NSolutionVariables(int var)
{
    switch (var) {
        case 39:
        case 40:
        case 41:
        case 42:
            return 3;
            break;
        case 43:
        case 44:
        case 45:
        case 46:
        case 47:
        case 100:
        case 101:
        case 102:
        case 103:
        case 104:
        case 105:
        case 106:
        case 107:
        case 108:
            return 1;
            break;
        default:
            DebugStop();
            break;
    }
    return 0;
}

/**
 * @brief It return a solution to multiphysics simulation.
 * @param datavec [in] Data material vector
 * @param var [in] number of solution variables. See  NSolutionVariables() method
 * @param Solout [out] is the solution vector
 */

void TPZHybridH1ErrorEstimateMaterial::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{

    /**
     datavec[0] H1 mesh, uh_reconstructed for Mark reconstruction and Empty for H1 reconstruction
     datavec[1] L2 mesh,
     datavec[2] Hdiv fem mesh, sigma_h
     datavec[3] L2 mesh fem, u_h
     **/

    int H1functionposition = 1;

    TPZFNMatrix<9,REAL> PermTensor;
    TPZFNMatrix<9,REAL> InvPermTensor;

    GetPermeabilities(datavec[0].x, PermTensor, InvPermTensor);


    TPZManVector<STATE,2> pressexact(1,0.);
    TPZFNMatrix<9,STATE> gradu(3,1,0.), fluxinv(3,1);

    if(fForcingFunctionExact)
    {
        this->fForcingFunctionExact->Execute(datavec[H1functionposition].x, pressexact,gradu);

    }

    PermTensor.Multiply(gradu, fluxinv);

    int dim=this->fDim;
    switch (var)
    {
        case 40://FluxFem
        {
            TPZFMatrix<REAL> &dsolaxes = datavec[3].dsol[0];
            TPZFNMatrix<9, REAL> dsol(3, 0);
            TPZFNMatrix<9, REAL> KGradsol(3, 0);
            TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, dsol, datavec[3].axes);

            PermTensor.Multiply(dsol, KGradsol);


            for (int i = 0; i < 3; i++) Solout[i]  = -KGradsol(i,0);
        }
            break;
        //Flux reconstrucion
        case 39: // sigma_h
        case 41:{// grad s_h
            TPZFMatrix<REAL> &dsolaxes = datavec[H1functionposition].dsol[0];
            TPZFNMatrix<9,REAL> dsol(3,0);
            TPZFNMatrix<9,REAL> KGradsol(3,0);
            TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, dsol, datavec[H1functionposition].axes);

            PermTensor.Multiply(dsol,KGradsol);

            if(var == 39){
                for(int id=0 ; id<3; id++) {
                    Solout[id] = datavec[0].sol[0][id];
                }
            }else{
                for(int id=0 ; id<3; id++) {
                    Solout[id] = -KGradsol(id,0);
                }
            }
        }


            break;
        case 42://flux exact
            for(int i=0; i<dim; i++) Solout[i] = -fluxinv(i);
            break;
        case 43://Pressure fem
            Solout[0] = datavec[3].sol[0][0];

            break;
        case 44://PressureReconstructed
            Solout[0] = datavec[H1functionposition].sol[0][0];

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

void TPZHybridH1ErrorEstimateMaterial:: ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors, TPZBndCond &bc){

    if(bc.Type()== 4){


        errors.Resize(NEvalErrors());
        errors.Fill(0.0);


        TPZFNMatrix<3,REAL> fluxreconstructed(3,1), fluxreconstructed2(3,1);
        TPZManVector<STATE,3> fluxfem(3);

        int H1functionposition = 1;

        REAL normalsigmafem = 0.,normalsigmarec = 0.,urec=0.;;
        normalsigmafem = data[2].sol[0][0];// sigma.n
        urec = data[H1functionposition].sol[0][0];



        REAL u_D = 0.,g = 0.;
        REAL normflux = 0.;

        TPZManVector<STATE,3> fluxrec(fDim);
        this->Solution(data,VariableIndex("FluxReconstructed"), fluxrec);

        std::cout<<"flux_rec "<<fluxrec[0]<<" , "<<fluxrec[1]<<"\n";


        TPZFNMatrix<9,REAL> PermTensor, InvPermTensor;
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9, STATE> gradu(this->Dimension(), 1);

        if (bc.HasForcingFunction()) {
            bc.ForcingFunction()->Execute(data[H1functionposition].x, res, gradu);
            GetPermeabilities(data[0].x, PermTensor, InvPermTensor);
            u_D = res[0];


        } else {
            // usualmente updatebc coloca o valor exato no val2
            u_D = bc.Val2()(0, 0);
        }


        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {

                normflux += data[2].normal[i]*PermTensor(i,j)*gradu(j,0);

            }
        }
        g = (-1)*normflux;

        std::cout<<"n_0 "<<data[2].normal[0]<<" n_1 "<<data[2].normal[1]<<"\n";






        REAL Km = bc.Val1()(0, 0);
        REAL InvKm = 1./Km;
        std::cout<<"Km "<<Km<<" InvKm "<<InvKm<<"\n";
        REAL errorEstimated =0.,errorReal = 0.;

        normalsigmarec = Km*(urec-u_D)+g;

//    std::cout<<"normalsigmarec "<<normalsigmarec<<"\n";
//    std::cout<<"normalsigmafem "<<normalsigmafem<<"\n";
//    std::cout<<"----------"<<"\n";
        errorEstimated = InvKm * (normalsigmarec - normalsigmafem)* (normalsigmarec - normalsigmafem);
        errorReal = InvKm * (normflux - normalsigmafem)* (normflux - normalsigmafem);
        errors[2] = errorReal;
        errors[3] = errorEstimated;
    }


}