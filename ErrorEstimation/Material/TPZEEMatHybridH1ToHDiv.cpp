//
//  TPZHybridH1ErrorEstimateMaterial.cpp
//  ErrorEstimation
//

#include "TPZEEMatHybridH1ToHDiv.h"
#include "pzaxestools.h"
#include "TPZAnalyticSolution.h"
#include "pzbndcond.h"


TPZEEMatHybridH1ToHDiv::TPZEEMatHybridH1ToHDiv(int matid, int dim) : TPZMixedPoisson(matid,dim)
{

}

TPZEEMatHybridH1ToHDiv::TPZEEMatHybridH1ToHDiv() : TPZMixedPoisson()
{

}

TPZEEMatHybridH1ToHDiv::TPZEEMatHybridH1ToHDiv(TPZMixedPoisson &copy) :TPZMixedPoisson(copy)
{

}

TPZEEMatHybridH1ToHDiv::TPZEEMatHybridH1ToHDiv(const TPZEEMatHybridH1ToHDiv &copy) : TPZMixedPoisson(copy)
{

}

TPZEEMatHybridH1ToHDiv::~TPZEEMatHybridH1ToHDiv()
{

}

TPZEEMatHybridH1ToHDiv &TPZEEMatHybridH1ToHDiv::operator=(const TPZEEMatHybridH1ToHDiv &copy)
{
    TPZMixedPoisson::operator=(copy);
    return *this;
}

TPZEEMatHybridH1ToHDiv::TPZEEMatHybridH1ToHDiv(TPZMatLaplacianHybrid matlaplacian){
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

int TPZEEMatHybridH1ToHDiv::FirstNonNullApproxSpaceIndex(TPZVec<TPZMaterialData> &datavec){

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

void TPZEEMatHybridH1ToHDiv::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    /**
     Implement
     (invK.sigma_h,v)_K = -(u_h,\nabla.v), ou
     (invK.sigma_h,v)_K = (f,v)
     int_K phi_i.phi_j dx = int_K u_h phi_i  dx
     **/

    TPZFNMatrix<9,REAL> PermTensor;
    TPZFNMatrix<9,REAL> InvPermTensor;

    GetPermeabilities(datavec[0].x, PermTensor, InvPermTensor);

    // Setting the phis
    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;

    TPZFMatrix<REAL> &phiuk = datavec[0].phi;
    TPZFMatrix<REAL> &dphiukaxes = datavec[0].dphix;
    TPZFNMatrix<9,REAL> dphiuk(3,dphiukaxes.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiukaxes, dphiuk, datavec[0].axes);

    int phrq;
    phrq = datavec[0].fVecShapeIndex.NElements();

    phiuk.Print(std::cout);
    dphiuk.Print(std::cout);

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
            /// (invK.sigma_h,v)_K

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
            //jvecZ.Print("mat1 = ");

            REAL prod1 = ivec(0,0)*jvecZ(0,0) + ivec(1,0)*jvecZ(1,0) + ivec(2,0)*jvecZ(2,0);
            ek(iq,jq) += fvisc*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod1;
        }

        /// -(u_h,\nabla.v)
        TPZFNMatrix<3,REAL> axesvec(3,1,0.);
        datavec[0].axes.Multiply(ivec,axesvec);

        REAL divwq = 0.;
        for(int iloc=0; iloc<fDim; iloc++)
        {
            divwq += axesvec(iloc,0)*dphiQ(iloc,ishapeind);
        }
        ef(iq,0) += (-1.)*weight*divwq*datavec[1].sol[0][0];
    }
}

void TPZEEMatHybridH1ToHDiv::ContributeBC(
        TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
        TPZFMatrix<STATE> &ef, TPZBndCond &bc) {

}


void TPZEEMatHybridH1ToHDiv::FillDataRequirements(TPZVec<TPZMaterialData > &datavec) {

    //fem solution for flux and potential
    /*datavec[2].SetAllRequirements(false);
    datavec[2].fNeedsSol = true;
    datavec[2].fNeedsNormal = true;

    datavec[3].SetAllRequirements(false);
    datavec[3].fNeedsSol = true;
*/

}

void TPZEEMatHybridH1ToHDiv::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec){

    datavec[2].SetAllRequirements(false);
    datavec[2].fNeedsSol = true;
    datavec[2].fNeedsNormal = true;


}

void TPZEEMatHybridH1ToHDiv::UpdateBCValues(TPZVec<TPZMaterialData> & datavec){
    DebugStop();

}

void TPZEEMatHybridH1ToHDiv::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
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


    TPZManVector<STATE,3> fluxfem(3);
    STATE divsigmafem, pressurefem, pressurereconstructed;

    TPZFNMatrix<3,REAL> fluxreconstructed(3,1), fluxreconstructed2(3,1);



    fluxfem=data[2].sol[0];
    divsigmafem= data[2].divsol[0][0];



    int H1functionposition = 0;
    H1functionposition = FirstNonNullApproxSpaceIndex(data);


    TPZVec<STATE> divsigma(1);

    if(this->fForcingFunctionExact){

        this->fForcingFunctionExact->Execute(data[H1functionposition].x,u_exact,du_exact);

        this->fForcingFunction->Execute(data[H1functionposition].x,divsigma);
    }

    REAL residual = 0.;

    // std::cout<<" divsigma[0] "<<divsigma[0]<<" divsigmafem "<<divsigmafem<<"\n";

    residual = (divsigma[0] - divsigmafem)*(divsigma[0] - divsigmafem);


    pressurereconstructed = data[H1functionposition].sol[0][0];


    pressurefem = data[3].sol[0][0];

    TPZFNMatrix<9,REAL> PermTensor;
    TPZFNMatrix<9,REAL> InvPermTensor;

    GetPermeabilities(data[2].x, PermTensor, InvPermTensor);

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



#ifdef PZDEBUG2
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


int TPZEEMatHybridH1ToHDiv::VariableIndex(const std::string &name)
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


int TPZEEMatHybridH1ToHDiv::NSolutionVariables(int var)
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

void TPZEEMatHybridH1ToHDiv::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
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

void TPZEEMatHybridH1ToHDiv:: ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors,TPZBndCond &bc){

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


