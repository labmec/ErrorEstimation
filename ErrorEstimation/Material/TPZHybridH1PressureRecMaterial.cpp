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
     * datavec[0]: H1-conform reconstructed potential
     * datavec[1]: Broken-H1 simulation result
     *
     Implement
     1) (K grad s_h, grad v)_K = (f,v):
            int_K phi_i.phi_j dx = int_K f phi_i  dx;
     **/

        int dim = datavec[fH1ReconstructionPosition].axes.Rows();
        //defining test functions
        // Setting the phis
        TPZFMatrix<REAL> &phiuk = datavec[fH1ReconstructionPosition].phi;
        TPZFMatrix<REAL> &dphiukaxes = datavec[fH1ReconstructionPosition].dphix; //(2xnphiuk)
        TPZFNMatrix<9, REAL> dphiuk(2, dphiukaxes.Cols());
        TPZAxesTools<REAL>::Axes2XYZ(dphiukaxes, dphiuk, datavec[fH1ReconstructionPosition].axes); //(3xnphiuk)
        TPZFMatrix<STATE> &dsolaxes = datavec[fFEMbrokenH1Position].dsol[0];
        TPZFNMatrix<9, REAL> dsol(2, dphiukaxes.Cols());
        TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, dsol, datavec[fH1ReconstructionPosition].axes);

        int nphiuk = phiuk.Rows();

        TPZFMatrix<STATE> gradSol(dim, 1,0), kGradSol(dim, 1,0), solukfem(1, 1,0);

        TPZFNMatrix<9, REAL> PermTensor(3,3);
        TPZFNMatrix<9, REAL> InvPermTensor(3,3);
        auto perm = GetPermeability(datavec[fH1ReconstructionPosition].x);
        PermTensor.Diagonal(perm);
        InvPermTensor.Diagonal(1./perm);


        //potetial fem
        solukfem(0, 0) = datavec[fFEMbrokenH1Position].sol[0][0];
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

void TPZHybridH1PressureRecMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
                                                  TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {

    /*
     Add Robin boundary condition for local problem
     ek+= <w,Km s_i>
     ef+= <w,Km*u_d - g + sigma_i.n>
     */

    int dim = datavec[fH1ReconstructionPosition].axes.Rows();

    TPZFMatrix<REAL> &phi_i = datavec[fH1ReconstructionPosition].phi;
    int nphi_i = phi_i.Rows();

    TPZFMatrix<STATE> solsigmafem(3, 1);
    solsigmafem.Zero();

    TPZManVector<REAL, 3> normal = datavec[fLagrangeCoeffPosition].normal;

    REAL normalsigma = 0.;
    normalsigma = datavec[fLagrangeCoeffPosition].sol[0][0];

    REAL u_D = 0.;
    REAL g = 0.;
    //REAL normflux = 0.;



    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9, STATE> gradu(dim, 1);
        bc.ForcingFunctionBC()(datavec[fH1ReconstructionPosition].x, res, gradu);
        u_D = res[0];
        g = normalsigma;

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

        datavec[fFEMbrokenH1Position].SetAllRequirements(false);
        datavec[fFEMbrokenH1Position].fNeedsSol = true;
}

void TPZHybridH1PressureRecMaterial::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const{

        datavec[fLagrangeCoeffPosition].SetAllRequirements(false);
        datavec[fLagrangeCoeffPosition].fNeedsSol = true;
        datavec[fLagrangeCoeffPosition].fNeedsNormal = true;
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
    ExactSol()(data[fH1ReconstructionPosition].x,u_exact,du_exact);

    errors.Resize(NEvalErrors());
    errors.Fill(0.0);

    STATE divsigmarec, u_h, s_h, forceProj;

    TPZFNMatrix<3,REAL> fluxreconstructed(3,1), KGradSh(3,1), gradreconstructed(3,1);
    TPZFMatrix<REAL> gradfemaxes(3,1),gradfem(3,1),KGradUh(3,1);

    gradfemaxes=data[fFEMbrokenH1Position].dsol[0];
    TPZAxesTools<REAL>::Axes2XYZ(gradfemaxes,gradfem,data[fFEMbrokenH1Position].axes);

    s_h = data[fH1ReconstructionPosition].sol[0][0];
    u_h = data[fFEMbrokenH1Position].sol[0][0];

    auto perm = GetPermeability(data[fH1ReconstructionPosition].x);
    auto invperm = 1./perm;

    TPZFNMatrix<3,REAL> KGradU(3,1);
    {
        for (int i=0; i<3; i++) {
            KGradU(i,0) = perm*du_exact[i];
            KGradUh(i,0) = perm*gradfem(i,0);
        }
    }

    TPZFMatrix<REAL> &gradshaxes = data[fH1ReconstructionPosition].dsol[0];
    TPZFNMatrix<9,REAL> gradsh(gradshaxes.Rows(),0);
    TPZAxesTools<REAL>::Axes2XYZ(gradshaxes, gradsh, data[fH1ReconstructionPosition].axes);
    for(int id=0 ; id<3; id++) {
        KGradSh(id,0) = perm*gradsh(id,0);
    }

    REAL FEMexactError = 0.;
    REAL FEMreconstructionError = 0.;
    REAL reconstructionExactError =0.;

#ifdef ERRORESTIMATION_DEBUG2
    std::cout<<"flux fem "<<fluxfem<<std::endl;
    std::cout<<"flux reconst "<<fluxreconstructed<<std::endl;
    std::cout<<"-------"<<std::endl;
#endif

    for (int i=0; i<3; i++) {
        FEMexactError += (KGradUh[i]-KGradU(i,0))*invperm*(KGradUh[i]-KGradU(i,0));
        reconstructionExactError += (KGradSh[i]-KGradU[i])*invperm*(KGradSh[i]-KGradU[i]);
        FEMreconstructionError += (KGradUh[i]-KGradSh[i])*invperm*(KGradUh[i]-KGradSh[i]);
    }

#ifdef ERRORESTIMATION_DEBUG2
    std::cout<<"potential fem "<<pressurefem<<std::endl;
    std::cout<<"potential reconst "<<pressurereconstructed<<std::endl;
    std::cout<<"-------"<<std::endl;
#endif
    
    errors[0] = FEMexactError;             // ||grad(u_h-u)||
    errors[1] = reconstructionExactError;  // ||grad(s_h-u)||
    errors[2] = FEMreconstructionError;    // ||grad(u_h-s_h)||
}



int TPZHybridH1PressureRecMaterial::VariableIndex(const std::string &name) const
{
    if(name == "FluxFem") return 40;
    if(name == "FluxExact") return 42;
    if(name == "PressureFem") return 43;
    if(name == "PressureReconstructed") return 44;
    if(name == "GradReconstructed") return 41;
    if(name == "PressureExact") return 45;
    if(name == "GradFEMerror") return 100;
    if(name == "GradReconstructionH1Error") return 101;
    if(name == "GradFEMreconstructionsH1Error") return 102;
    if(name == "POrder") return 46;

    return -1;
}


int TPZHybridH1PressureRecMaterial::NSolutionVariables(int var) const
{
    switch (var) {
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

    auto perm = GetPermeability(datavec[0].x);

    TPZManVector<STATE,2> pressexact(1,0.);
    TPZFNMatrix<9,STATE> gradu(3,1,0.), KGradU(3,1);

    if(fExactSol)
    {
        this->fExactSol(datavec[fH1ReconstructionPosition].x, pressexact,gradu);

    }

    for(int i=0; i<3;i++)
        KGradU = perm*gradu;

    int dim=this->fDim;
    switch (var)
    {
    case 40://FluxFem
    {
        TPZFMatrix<REAL> &dsolaxes = datavec[fFEMbrokenH1Position].dsol[0];
        TPZFNMatrix<9, REAL> dsol(3, 0);
        TPZFNMatrix<9, REAL> KGradsol(3, 0);
        TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, dsol, datavec[fFEMbrokenH1Position].axes);

        for(int i=0;i<3;i++)
            KGradsol = perm*dsol;

        for (int i = 0; i < 3; i++) Solout[i]  = -KGradsol(i,0);
    }
    break;
    //Flux reconstrucion
    case 41:// grad s_h
    {
        TPZFMatrix<REAL> &dsolaxes = datavec[fH1ReconstructionPosition].dsol[0];
        TPZFNMatrix<9, REAL> dsol(3, 0);
        TPZFNMatrix<9, REAL> KGradsol(3, 0);
        TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, dsol, datavec[fH1ReconstructionPosition].axes);

        for (int i = 0; i < 3; i++)
            KGradsol = dsol * perm;

        for (int id = 0; id < 3; id++) {
            Solout[id] = -KGradsol(id, 0);
        }
    }
    break;
    case 42://flux exact
        for(int i=0; i<dim; i++) Solout[i] = -KGradU(i);
        break;
    case 43://Pressure fem
        Solout[0] = datavec[fFEMbrokenH1Position].sol[0][0];

        break;
    case 44://PressureReconstructed
        Solout[0] = datavec[fH1ReconstructionPosition].sol[0][0];

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

int TPZHybridH1PressureSingleSpace::NSolutionVariables(int var) const{
     if (var == 0) {
        return 1;
    } else if (var == 1) {
        return 3;
    } else {
        DebugStop();
        return 0;
    }
}

void TPZHybridH1PressureSingleSpace::Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<STATE> &sol) {
    if (var == 0) {
        sol = data.sol[0];
    } else if (var == 1) {
        for (int i = 0; i < 3; i++) {
            sol[i] += data.dsol[0].GetVal(i, 0);
        }
    } else {
        DebugStop();
    }
}