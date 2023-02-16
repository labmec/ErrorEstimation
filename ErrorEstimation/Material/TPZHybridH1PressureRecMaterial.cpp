//
// Created by victor on 14/02/23.
//

#include "TPZHybridH1PressureRecMaterial.h"

//
// Created by victor on 07/10/2020.
//

#include "TPZHybridH1ErrorEstimateMaterial.h"
#include "pzaxestools.h"
#include "TPZAnalyticSolution.h"
#include "TPZMaterialDataT.h"


#ifdef LOG4CXX
static LoggerPtr loggerF(Logger::getLogger("DebuggingF"));
#endif

TPZHybridH1PressureRecMaterial::TPZHybridH1PressureRecMaterial(int matid, int dim) : TPZMixedPoisson(matid, dim)
{
}

TPZHybridH1PressureRecMaterial::TPZHybridH1PressureRecMaterial() : TPZMixedPoisson()
{

}

TPZHybridH1PressureRecMaterial::TPZHybridH1PressureRecMaterial(const TPZHybridH1PressureRecMaterial &copy) : TPZMixedPoisson(copy)
{
}

TPZHybridH1PressureRecMaterial::TPZHybridH1PressureRecMaterial(const TPZMixedPoisson &copy) : TPZMixedPoisson(copy)
{

}

TPZHybridH1PressureRecMaterial::TPZHybridH1PressureRecMaterial(TPZMatLaplacianHybrid &matlaplacian)
{
    this->SetId(matlaplacian.Id());
    this->SetDimension(matlaplacian.Dimension());

    STATE perm = matlaplacian.GetPermeability({0.,0.,0.});
    SetConstantPermeability(perm);

    if (matlaplacian.HasForcingFunction()) {
        TPZMatErrorCombinedSpaces<STATE> *me = this;
        TPZMatErrorCombinedSpaces<STATE> *lapl = &matlaplacian;
        auto funcpt = lapl->ExactSol();
        me->SetExactSol(funcpt,lapl->PolynomialOrderExact());
        this->SetForcingFunction(matlaplacian.ForcingFunction(),matlaplacian.ForcingFunctionPOrder());
    }
}

TPZHybridH1PressureRecMaterial::~TPZHybridH1PressureRecMaterial()
{

}

TPZHybridH1PressureRecMaterial &TPZHybridH1PressureRecMaterial::operator=(const TPZHybridH1PressureRecMaterial &copy)
{
    TPZMixedPoisson::operator=(copy);
    return *this;
}

void TPZHybridH1PressureRecMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    /**
     Implement

     1) (K grad s_h, grad v)_K = (f,v):
            int_K phi_i.phi_j dx = int_K f phi_i  dx;
     **/
        int H1functionposition = 1;

        int dim = datavec[H1functionposition].axes.Rows();
        //defining test functions
        // Setting the phis
        TPZFMatrix<REAL> &phiuk = datavec[H1functionposition].phi;
        TPZFMatrix<REAL> &dphiukaxes = datavec[H1functionposition].dphix; //(2xnphiuk)
        TPZFNMatrix<9, REAL> dphiuk(2, dphiukaxes.Cols());
        TPZAxesTools<REAL>::Axes2XYZ(dphiukaxes, dphiuk, datavec[H1functionposition].axes); //(3xnphiuk)
        TPZFMatrix<STATE> &dsolaxes = datavec[3].dsol[0];
        TPZFNMatrix<9, REAL> dsol(2, dphiukaxes.Cols());
        TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, dsol, datavec[H1functionposition].axes);

        int nphiuk = phiuk.Rows();

        TPZFMatrix<STATE> gradSol(dim, 1,0), kGradSol(dim, 1,0), solukfem(1, 1,0);

        TPZVec<STATE> divsigma(1);
        divsigma[0] = 0;

        if (!this->fForcingFunction) DebugStop();
        this->fForcingFunction(datavec[H1functionposition].x, divsigma);

        TPZFNMatrix<9, REAL> PermTensor(3,3);
        TPZFNMatrix<9, REAL> InvPermTensor(3,3);
        auto perm = GetPermeability(datavec[H1functionposition].x);
        PermTensor.Diagonal(perm);
        InvPermTensor.Diagonal(1./perm);


        //potetial fem
        solukfem(0, 0) = datavec[3].sol[0][0];
        for (int ip = 0; ip < dim; ip++) {
            gradSol(ip, 0) = dsol.Get(ip, 0);
        }
        for (int id = 0; id < dim; id++) {
            for (int jd = 0; jd < dim; jd++) {
                kGradSol(id, 0) += PermTensor(id, jd) * gradSol(jd, 0);
            }
        }

        TPZFMatrix<STATE> kgraduk(dim, nphiuk, 0.);



        for (int irow = 0; irow < nphiuk; irow++) {

            for (int id = 0; id < dim; id++) {
                for (int jd = 0; jd < dim; jd++) {
                    kgraduk(id, irow) += PermTensor(id, jd) * dphiuk(jd, irow);
                }
                /// ... = (grad u_h, grad v_h)
                ef(irow,0) +=weight*dphiuk(id,irow)*kGradSol(id,0);
            }
            ///... = (f , v_h)
            //ef(irow, 0) += weight * phiuk(irow, 0) * divsigma[0];

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

void TPZHybridH1PressureRecMaterial::ContributeBC(
    const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
    TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {

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



        if (bc.HasForcingFunctionBC()) {
            TPZManVector<STATE> res(3);
            TPZFNMatrix<9, STATE> gradu(dim, 1);
            bc.ForcingFunctionBC()(datavec[H1functionposition].x, res, gradu);
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
            u_D = bc.Val2()[0];
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
                ef(iq, 0) += fBigNumber * u_D * phi_i(iq, 0) * weight;
                for (int jq = 0; jq < nphi_i; jq++) {
                    ek(iq, jq) += fBigNumber * weight * phi_i(iq, 0) * phi_i(jq, 0);
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


void TPZHybridH1PressureRecMaterial::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const {

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

void TPZHybridH1PressureRecMaterial::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const{

        datavec[0].SetAllRequirements(false);
        datavec[0].fNeedsSol = true;
        datavec[0].fNeedsNormal = true;

        datavec[2].SetAllRequirements(false);
        datavec[2].fNeedsSol = true;
        datavec[2].fNeedsNormal = true;

        datavec[4].SetAllRequirements(false);

}

void TPZHybridH1PressureRecMaterial::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors)
{
    /**
     datavec[0] H1 mesh, uh_reconstructed
     datavec[1] L2 mesh,
     datavec[2] Hdiv fem mesh, sigma_h
     datavec[3] L2 mesh fem, u_h
     datavec[4] L2 projection

      error[0] - error computed with exact pressure
      error[1] - error computed with reconstructed pressure
      error[2] - energy error computed with exact solution
      error[3] - energy error computed with reconstructed flux
      error[4] - energy error computed with reconstructed potential
      error[5] - oscilatory data error

     **/
    if(!ExactSol()) DebugStop();
    TPZVec<STATE> u_exact(1);
    TPZFMatrix<STATE> du_exact(3,1,0.);
    ExactSol()(data[1].x,u_exact,du_exact);

    errors.Resize(NEvalErrors());
    errors.Fill(0.0);


    STATE divsigmarec, pressurefem, pressurereconstructed, forceProj;

    TPZFNMatrix<3,REAL> fluxreconstructed(3,1), fluxreconstructed2(3,1), gradreconstructed(3,1);
    TPZFMatrix<REAL> gradfem,fluxfem(3,1);

    gradfem=data[3].dsol[0];

    gradfem.Resize(3,1);


    int H1functionposition = 1;

    TPZVec<STATE> divsigma(1);

    if(this->fForcingFunction){

        this->fForcingFunction(data[H1functionposition].x,divsigma);
    }

    REAL residual = 0.,altResidual = 0.;
    divsigmarec= data[0].divsol[0][0];
    residual = (divsigma[0] - divsigmarec)*(divsigma[0] - divsigmarec);

    forceProj = data[4].sol[0][0];
    altResidual = forceProj - divsigma[0];

    pressurereconstructed = data[H1functionposition].sol[0][0];


    pressurefem = data[3].sol[0][0];

    TPZFNMatrix<9,REAL> PermTensor(3,3);
    TPZFNMatrix<9,REAL> InvPermTensor(3,3);
    auto perm = GetPermeability(data[1].x);

    PermTensor.Diagonal(perm);
    InvPermTensor.Diagonal(1./perm);


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



#ifdef ERRORESTIMATION_DEBUG2
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

#ifdef ERRORESTIMATION_DEBUG2
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
    //errors[5] = npz;
    errors[6] = altResidual*altResidual;
}



int TPZHybridH1PressureRecMaterial::VariableIndex(const std::string &name) const
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
    if(name == "NCIndex") return 103;
    if(name == "NRIndex") return 104;
    if(name == "NFIndex") return 105;
    if(name == "PressureEffectivityIndex") return 107;
    if(name == "EnergyEffectivityIndex") return 108;
    if(name == "EnergyErrorEstimated") return 109;
    if(name == "POrder") return 46;

    return -1;
}


int TPZHybridH1PressureRecMaterial::NSolutionVariables(int var) const
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
    case 100:
    case 101:
    case 102:
    case 103:
    case 104:
    case 105:
    case 106:
    case 107:
    case 108:
    case 109:
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

void TPZHybridH1PressureRecMaterial::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout)
{

    /**
     datavec[0] H1 mesh, uh_reconstructed for Mark reconstruction and Empty for H1 reconstruction
     datavec[1] L2 mesh,
     datavec[2] Hdiv fem mesh, sigma_h
     datavec[3] L2 mesh fem, u_h
     **/

    int H1functionposition = 1;

    TPZFNMatrix<9,REAL> PermTensor(3,3);
    TPZFNMatrix<9,REAL> InvPermTensor(3,3);

    auto perm = GetPermeability(datavec[0].x);
    PermTensor.Diagonal(perm);
    InvPermTensor.Diagonal(1./perm);

    TPZManVector<STATE,2> pressexact(1,0.);
    TPZFNMatrix<9,STATE> gradu(3,1,0.), fluxinv(3,1);

    if(fExactSol)
    {
        this->fExactSol(datavec[H1functionposition].x, pressexact,gradu);

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
        //TPZMixedPoisson::Solution(datavec,var, Solout);
        DebugStop();
    }
}

