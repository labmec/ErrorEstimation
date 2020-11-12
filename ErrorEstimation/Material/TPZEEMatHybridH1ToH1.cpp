//
// Created by victor on 07/10/2020.
//

#include "TPZEEMatHybridH1ToH1.h"
#include "pzaxestools.h"
#include "TPZAnalyticSolution.h"
#include "pzbndcond.h"


TPZEEMatHybridH1ToH1::TPZEEMatHybridH1ToH1(int matid, int dim) : TPZMatLaplacianHybrid(matid,dim)
{

}

TPZEEMatHybridH1ToH1::TPZEEMatHybridH1ToH1() : TPZMatLaplacianHybrid()
{

}

TPZEEMatHybridH1ToH1::TPZEEMatHybridH1ToH1(const TPZMatLaplacianHybrid &copy) : TPZMatLaplacianHybrid(copy)
{

}

TPZEEMatHybridH1ToH1::TPZEEMatHybridH1ToH1(const TPZEEMatHybridH1ToH1 &copy) : TPZMatLaplacianHybrid(copy)
{

}

TPZEEMatHybridH1ToH1::~TPZEEMatHybridH1ToH1()
{

}

TPZEEMatHybridH1ToH1 &TPZEEMatHybridH1ToH1::operator=(const TPZEEMatHybridH1ToH1 &copy)
{
    TPZMatLaplacianHybrid::operator=(copy);
    return *this;
}

int TPZEEMatHybridH1ToH1::IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const{
    int order = 0;
    if (fForcingFunction) {
        order = fForcingFunction->PolynomialOrder();
    }

    int pmax = 0;
    for (int ip = 0; ip < elPMaxOrder.size(); ip++) {
        if (elPMaxOrder[ip] > pmax) pmax = elPMaxOrder[ip];
    }
    pmax += 1; // HDiv simulations use an additional integration order. In order to test if HybridH1 and HDiv are equivalents, both integration orders shall be equivalent.

    int integrationorder = 2 * pmax;
    if (pmax < order) {
        integrationorder = order + pmax;
    }

    return integrationorder;
}

int TPZEEMatHybridH1ToH1::FirstNonNullApproxSpaceIndex(TPZVec<TPZMaterialData> &datavec){

    int nvec = datavec.NElements();
    int firstNoNullposition = -1;

    for(int ivec = 0; ivec < nvec ; ivec++){
        TPZMaterialData::MShapeFunctionType shapetype = datavec[ivec].fShapeType;
        if(shapetype != TPZMaterialData::EEmpty){
            firstNoNullposition = ivec;
            return firstNoNullposition;
        }

    }
    return firstNoNullposition;
}

void TPZEEMatHybridH1ToH1::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
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
    TPZFMatrix<REAL> &dphiukaxes = datavec[H1functionposition].dphix; //(2xnphiuk)
    TPZFNMatrix<9,REAL> dphiuk(2,dphiukaxes.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiukaxes, dphiuk, datavec[H1functionposition].axes); //(3xnphiuk)
    TPZFNMatrix<15, STATE> &dsolaxes = datavec[H1functionposition].dsol[0];
    TPZFNMatrix<9,REAL> dsol(2,dphiukaxes.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, dsol, datavec[H1functionposition].axes);

    int nphiuk = phiuk.Rows();

    TPZFMatrix<STATE> gradSol(dim,1),kGradSol(dim,1),solukfem(1,1);
    gradSol.Zero();
    solukfem.Zero();
    kGradSol.Zero();

    TPZVec<STATE> divsigma(1); divsigma[0] = 0;

    if(this->fForcingFunction){
        this->fForcingFunction->Execute(datavec[H1functionposition].x,divsigma);
    }

    //potetial fem
    solukfem(0,0) = datavec[3].sol[0][0];
    for(int ip = 0 ; ip <dim ; ip++){
        gradSol(ip,0) = dsol.Get(ip,0);
    }

    TPZFNMatrix<9,REAL> PermTensor;
    TPZFNMatrix<9,REAL> InvPermTensor;

    GetPermeabilities(datavec[1].x, PermTensor, InvPermTensor);

    for (int id = 0; id < dim; id++) {
        for (int jd = 0; jd < dim; jd++) {
            kGradSol(id, 0) += PermTensor(id, jd) * gradSol(jd, 0);
        }
    }

    TPZFMatrix<STATE> kgraduk(dim,nphiuk,0.);

    for(int irow=0 ; irow<nphiuk; irow++){
        //K graduk
        for(int id=0; id< dim; id++){
            for(int jd=0; jd< dim;jd++){

                kgraduk(id,irow) += PermTensor(id,jd)*dphiuk(jd,irow);

            }
            //bk = (-1)*int_k sigmaukfem.grad phi_i,here dphiuk is multiplied by axes

            //ef(irow,0)+=weight*dphiuk(id,irow)*kGradSol(id,0);
        }
        ef(irow,0)+= weight*phiuk(irow,0)*divsigma[0];

        //matrix Sk= int_{K} K graduk.gradv
        for(int jcol=0; jcol<nphiuk;jcol++){

            for(int jd=0;  jd< dim; jd++)
            {
                ek(irow,jcol) +=weight*kgraduk(jd,irow)*dphiuk(jd,jcol);
            }

        }
    }
}

void TPZEEMatHybridH1ToH1::ContributeBC(
        TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
        TPZFMatrix<STATE> &ef, TPZBndCond &bc) {

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
        case (4):{

            REAL Km = bc.Val1()(0, 0); // Km
            REAL InvKm = 1./Km;
            // std::cout<< " g "<< g<< " normalsigma "<<normalsigma<<"\n";
            REAL robinterm = (Km * u_D - g+ normalsigma);

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




        case( 0):{

        for (int iq = 0; iq < nphi_i; iq++) {
            ef(iq, 0) += gBigNumber*u_D * phi_i(iq, 0) * weight;
            for (int jq = 0; jq < nphi_i; jq++) {
                ek(iq, jq) += gBigNumber*weight *  phi_i(iq, 0) * phi_i(jq, 0);
            }
        }

            break;
        }

        default:{

            std::cout << " This material not implement BC Type " << bc.Type()<< std::endl;
            break;
        }
    }


}


void TPZEEMatHybridH1ToH1::FillDataRequirements(TPZVec<TPZMaterialData > &datavec) {

    //fem solution for flux and potential
    datavec[2].SetAllRequirements(false);
    datavec[2].fNeedsSol = true;
    datavec[2].fNeedsNormal = true;

    datavec[3].SetAllRequirements(false);
    datavec[3].fNeedsSol = true;


}

void TPZEEMatHybridH1ToH1::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec){

    datavec[2].SetAllRequirements(false);
    datavec[2].fNeedsSol = true;
    datavec[2].fNeedsNormal = true;


}

void TPZEEMatHybridH1ToH1::UpdateBCValues(TPZVec<TPZMaterialData> & datavec){
    DebugStop();

}

// Todo: Suport sigma errors
void TPZEEMatHybridH1ToH1::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{
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


    //TPZManVector<STATE,3> fluxfem(3);
    STATE divsigmafem, pressurefem, pressurereconstructed;

    TPZFNMatrix<3,REAL> fluxreconstructed(3,1), fluxreconstructed2(3,1);
    TPZFMatrix<REAL> fluxfem, gradfem;

    //fluxfem=data[2].sol[0];
    //divsigmafem= data[2].divsol[0][0];

    gradfem=data[3].dsol[0];

    gradfem.Resize(3,1);

    //divsigmafem= data[2].divsol[0][0];



    int H1functionposition = 1;
    H1functionposition = FirstNonNullApproxSpaceIndex(data);


    TPZVec<STATE> divsigma(1);

    if(this->fForcingFunctionExact){

        this->fForcingFunctionExact->Execute(data[H1functionposition].x,u_exact,du_exact);

        this->fForcingFunction->Execute(data[H1functionposition].x,divsigma);
    }

    REAL residual = 0.;

    //std::cout<<" divsigma[0] "<<divsigma[0]<<" divsigmafem "<<divsigmafem<<"\n";

    //residual = (divsigma[0] - divsigmafem)*(divsigma[0] - divsigmafem);


    pressurereconstructed = data[H1functionposition].sol[0][0];


    pressurefem = data[3].sol[0][0];

    TPZFNMatrix<9,REAL> PermTensor;
    TPZFNMatrix<9,REAL> InvPermTensor;

    GetPermeabilities(data[1].x, PermTensor, InvPermTensor);

    PermTensor.Resize(3,3);
    InvPermTensor.Resize(3,3);


    TPZFNMatrix<3,REAL> fluxexactneg;

    // sigmarec = -K grad(urec)
    //  sigmak = -K graduk

    {
        TPZFNMatrix<9,REAL> gradpressure(3,1);
        for (int i=0; i<3; i++) {
            gradpressure(i,0) = (-1.)*du_exact[i];
            gradfem(i,0) = (-1.)* gradfem(i,0);
        }
        PermTensor.Multiply(gradpressure,fluxexactneg);
        PermTensor.Multiply(gradfem,fluxfem);
    }



    TPZFMatrix<REAL> &dsolaxes = data[H1functionposition].dsol[0];
    TPZFNMatrix<9,REAL> fluxrec(dsolaxes.Rows(),0);
    TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, fluxrec, data[H1functionposition].axes);

    for(int id=0 ; id<3; id++) {
        fluxreconstructed2(id,0) = (-1.)*fluxrec(id,0);
    }

    PermTensor.Multiply(fluxreconstructed2,fluxreconstructed);


    REAL innerexact = 0.;
    REAL innerestimate = 0.;



#ifdef PZDEBUG2
    std::cout<<"flux fem "<<fluxfem<<std::endl;
    std::cout<<"flux reconst "<<fluxreconstructed<<std::endl;
    std::cout<<"-------"<<std::endl;
#endif



    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            innerexact += (fluxfem[i]-fluxexactneg(i,0))*InvPermTensor(i,j)*(fluxfem[j]-fluxexactneg(j,0));//Pq esta somando: o fluxo fem esta + e o exato -
            innerestimate += (fluxfem[i]-fluxreconstructed[i])*InvPermTensor(i,j)*(fluxfem[j]-fluxreconstructed[j]);
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
    errors[3] = innerestimate;//error flux reconstructed
    errors[4] = residual; //||f - Proj_divsigma||




}


int TPZEEMatHybridH1ToH1::VariableIndex(const std::string &name)
{
    if(name == "FluxFem") return 40;
    if(name == "FluxExact") return 42;
    if(name == "PressureFem") return 43;
    if(name == "PressureReconstructed") return 44;
    if(name == "FluxReconstructed") return 41;
    if(name == "PressureExact") return 45;
    if(name == "PressureErrorExact") return 100;
    if(name == "PressureErrorEstimate") return 101;
    if(name == "EnergyErrorExact") return 102;
    if(name == "EnergyErrorEstimate") return 103;
    if(name == "ResidualError") return 104;
    if(name == "PressureEffectivityIndex") return 105;
    if(name == "EnergyEffectivityIndex") return 106;
    if(name == "POrder") return 46;

    return -1;
}


int TPZEEMatHybridH1ToH1::NSolutionVariables(int var)
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
        case 47:
        case 100:
        case 101:
        case 102:
        case 103:
        case 104:
        case 105:
        case 106:
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

void TPZEEMatHybridH1ToH1::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{

    /**
     datavec[0] H1 mesh, uh_reconstructed for Mark reconstruction and Empty for H1 reconstruction
     datavec[1] L2 mesh,
     datavec[2] Hdiv fem mesh, sigma_h
     datavec[3] L2 mesh fem, u_h
     **/

    int H1functionposition = 0;
    H1functionposition = FirstNonNullApproxSpaceIndex(datavec);

    TPZFNMatrix<9,REAL> PermTensor;
    TPZFNMatrix<9,REAL> InvPermTensor;

    GetPermeabilities(datavec[2].x, PermTensor, InvPermTensor);


    TPZManVector<STATE,2> pressexact(1,0.);
    TPZFNMatrix<9,STATE> gradu(3,1,0.), fluxinv(3,1);

    if(fForcingFunctionExact)
    {
        this->fForcingFunctionExact->Execute(datavec[H1functionposition].x, pressexact,gradu);

    }

    PermTensor.Redim(3,3);

    PermTensor.Multiply(gradu, fluxinv);

    int dim=this->fDim;
    switch (var)
    {
        case 40://FluxFem

            for(int i=0; i<3; i++) Solout[i] = datavec[2].sol[0][i];

            break;
        case 41:{//FluxReconstructed is grad U
            TPZFMatrix<REAL> &dsolaxes = datavec[H1functionposition].dsol[0];
            TPZFNMatrix<9,REAL> dsol(3,0);
            TPZFNMatrix<9,REAL> KGradsol(3,0);
            TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, dsol, datavec[H1functionposition].axes);

            PermTensor.Multiply(dsol,KGradsol);

            for(int id=0 ; id<fDim; id++) {
                Solout[id] = KGradsol(id,0);//dsol(id,0);//derivate
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

void TPZEEMatHybridH1ToH1:: ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors,TPZBndCond &bc){

    if(bc.Type()== 4){


        errors.Resize(NEvalErrors());
        errors.Fill(0.0);


        TPZFNMatrix<3,REAL> fluxreconstructed(3,1), fluxreconstructed2(3,1);
        TPZManVector<STATE,3> fluxfem(3);

        int H1functionposition = 0;
        H1functionposition = FirstNonNullApproxSpaceIndex(data);

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