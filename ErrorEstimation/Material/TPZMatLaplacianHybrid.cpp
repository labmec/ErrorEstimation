//
//  TPZMatLaplacianHybrid.cpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 14/07/19.
//

#include "TPZMatLaplacianHybrid.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "pzaxestools.h"

TPZMatLaplacianHybrid::TPZMatLaplacianHybrid(int matid, int dim)
: TPZRegisterClassId(&TPZMatLaplacianHybrid::ClassId), TPZMatLaplacian(matid,dim)
{
    
}

TPZMatLaplacianHybrid::TPZMatLaplacianHybrid() :
TPZRegisterClassId(&TPZMatLaplacianHybrid::ClassId), TPZMatLaplacian()
{
    
}

TPZMatLaplacianHybrid::TPZMatLaplacianHybrid(const TPZMatLaplacianHybrid &copy) : TPZDarcyFlow(copy),
TPZMatCombinedSpacesT<STATE>(copy), TPZMatErrorCombinedSpaces<STATE>(copy), TPZMatError<STATE>(copy) {
//    TPZMatError<STATE>::operator=(copy);
//    const TPZMatErrorCombinedSpaces *errcs = &copy;
}

TPZMatLaplacianHybrid::~TPZMatLaplacianHybrid()
{
    
}

TPZMatLaplacianHybrid &TPZMatLaplacianHybrid::operator=(const TPZMatLaplacianHybrid &copy)
{
    TPZMatLaplacian::operator=(copy);
    TPZMatCombinedSpacesT<STATE>::operator=(copy);
    TPZMatErrorCombinedSpaces<STATE>::operator=(copy);
//    fDim = copy.fDim;
    return *this;
}

TPZMaterial *TPZMatLaplacianHybrid::NewMaterial() const
{
    return new TPZMatLaplacianHybrid(*this);
}

int TPZMatLaplacianHybrid::ClassId() const
{
    return Hash("TPZMatLaplacianHybrid") ^ TPZMatLaplacian::ClassId() << 1;
}


void TPZMatLaplacianHybrid::Write(TPZStream &buf, int withclassid) const
{
    TPZMatLaplacian::Write(buf,withclassid);
}

void TPZMatLaplacianHybrid::Read(TPZStream &buf, void *context)
{
    TPZMatLaplacian::Read(buf,context);
}

int TPZMatLaplacianHybrid::VariableIndex(const std::string &name) const
{
    
    if(name == "Pressure") return 44;
    if(name == "PressureExact") return 45;
    
    if(name == "Flux") return 7;
    if(name == "ExactFlux") return 10;
    
    if(name == "ExactFluxShiftedOrigin") return 23;

    if(name == "Kperm") return 46;
    
    return -1;
}

int TPZMatLaplacianHybrid::NSolutionVariables(int var) const {
    if(var == 44 || var==45 || var==46) return 1;
    if(var == 7 || var==10 || var == 23) return fDim;
    
    else{
        
        return TPZMatLaplacian::NSolutionVariables(var);
        return 0;
        
    }
}

int TPZMatLaplacianHybrid::IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const{
    int order = 0;
    if (fForcingFunction) {
        order = fForcingFunctionPOrder;
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




void TPZMatLaplacianHybrid::Contribute(const TPZVec<TPZMaterialDataT<STATE> > &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    /**
     datavec[1] L2 mesh (phi's)
     datavec[0] Hdiv mesh,
     datavec[2] Interface Mesh
     datavec[3] Interface Mesh
     
     Implement the matrix
     |Sk Ck^T |  = |f1|
     |Ck  0   |    |f2|
     Sk = int_K K graduk.gradv dx = int_K K gradphi_i.gradphi_j dx
     CK = int_partialK lambda_k*uk dx = int_K phi_i dx
     f1 = int_K f*v dx = int_K f*phi_j dx
     ck = int_partialK phi_i*mu_j dx
     f2 = int_partialK g*mu_j dx
     
     **/
    
    TPZFMatrix<REAL>  &phi = datavec[1].phi;
    TPZFMatrix<REAL> &dphiaxes = datavec[1].dphix;
    TPZFNMatrix<9, REAL> dphi(3, dphiaxes.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiaxes, dphi, datavec[1].axes);
    TPZVec<REAL>  &x = datavec[1].x;
    
    int phr = phi.Rows();
    
    STATE fXfLoc = 0.;
    
    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        //TPZFMatrix<STATE> dres(Dimension(),1);
        //fForcingFunction->Execute(x,res,dres);       // dphi(i,j) = dphi_j/dxi
        fForcingFunction(x,res);
        fXfLoc = res[0];
    }
    
    STATE KPerm = GetPermeability(x);
    
    //Equacao de Poisson
    for( int in = 0; in < phr; in++ ) {
        int kd;
        ef(in, 0) +=  (STATE)weight * fXfLoc * (STATE)phi(in,0);
        
        //matrix Sk
        for( int jn = 0; jn < phr; jn++ ) {
            for(kd=0; kd<fDim; kd++) {
                ek(in,jn) += (STATE)weight*(KPerm*(STATE)(dphi(kd,in)*dphi(kd,jn)));
            }
        }
    }
    // u_average and g contributions
    for (int in =0; in < phr; in++) {
        ek(phr,in) += weight*phi(in,0);
        ek(in,phr) += weight*phi(in,0);
    }
    
    ek(phr,phr+1) -= weight;
    ek(phr+1,phr) -= weight;
    
#ifdef PZDEBUG
    if ( ek.VerifySymmetry(1.e-10) == SymProp::NonSym) std::cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << std::endl;
#endif
    
}

void TPZMatLaplacianHybrid::Contribute(const TPZVec<TPZMaterialDataT<STATE> > &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    TPZFMatrix<REAL>  &phi = datavec[1].phi;
    TPZFMatrix<REAL> &dphi = datavec[1].dphix;
    TPZVec<REAL>  &x = datavec[1].x;
    //    TPZFMatrix<REAL> &axes = data.axes;
    //    TPZFMatrix<REAL> &jacinv = data.jacinv;
    int phr = phi.Rows();
    
    STATE fXfLoc = 0.;
    
    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        //TPZFMatrix<STATE> dres(Dimension(),1);
        //fForcingFunction->Execute(x,res,dres);       // dphi(i,j) = dphi_j/dxi
        fForcingFunction(x,res);
        fXfLoc = res[0];
    } 
    
    STATE KPerm = GetPermeability(x);
    
    //Equacao de Poisson
    for( int in = 0; in < phr; in++ ) {
        int kd;
        ef(in, 0) +=  (STATE)weight * fXfLoc * (STATE)phi(in,0);
        for(kd=0; kd<fDim; kd++) {
            ef(in,0) -= (STATE)weight*(KPerm*(dphi(kd,in)*datavec[1].dsol[0](kd,0)));
        }
    }
    ef(phr,0) += weight*(-datavec[1].sol[0][0]+ datavec[3].sol[0][0]);
    ef(phr+1,0) += weight*datavec[2].sol[0][0];
}

void TPZMatLaplacianHybrid::ContributeBC(const TPZVec<TPZMaterialDataT<STATE> > &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc) {
    TPZFMatrix<REAL> &phi_u = datavec[1].phi;
    TPZFMatrix<REAL> &phi_flux = datavec[0].phi;
    //    TPZFMatrix<REAL> &axes = data.axes;
    int phr_primal = phi_u.Rows();
    int phr_hybrid = phi_flux.Rows();
    bool primal = true;
    TPZManVector<REAL, 3> x(3);
    if (phr_hybrid) {
        primal = false;
        x = datavec[0].x;
    } else {
        x = datavec[1].x;
    }
    short in, jn;
    STATE v2[1];
    v2[0] = bc.Val2()[0];
    REAL normflux = 0.;
    int dim = Dimension();
    int dimBC = Dimension()-1;

    if (bc.HasForcingFunctionBC()) {            // phi(in, 0) = phi_in                          // JORGE 2013 01 26
        TPZManVector<STATE> res(1);
        TPZFNMatrix<3, STATE> dres(3, 1);
        bc.ForcingFunctionBC()(x, res, dres);       // dphi(i,j) = dphi_j/dxi
        v2[0] = res[0];
        TPZManVector<REAL, 3> normal(3, 0), vec(3), axes1(3), axes2(3), zaxes(3, 0);
        zaxes[2] = 1;
        bool unitary = true;
        if (!primal && bc.Type() == 1) {
            TPZFNMatrix<9, STATE> PermTensor(3,3), InvPermTensor(3,3);
            STATE KPerm = GetPermeability(datavec[0].x);
            PermTensor.Diagonal(KPerm);
            for (int i = 0; i < 3; i++) {
                axes1[i] = datavec[0].axes(0, i);
            }
            switch (dimBC) {
                case 0: //normal = axes
                    std::cout << "Implement me\n";
                    DebugStop();
                    break;
                case 1: // axes x (0,0,1)
                    this->VectorialProd(axes1, zaxes, normal, unitary);
                    break;
                case 2:
                    for (int i = 0; i < 3; i++) {
                        axes2[i] = datavec[0].axes(1, i);
                    }
                    this->VectorialProd(axes1, axes2, normal, unitary);
                    break;
                default:
                    DebugStop();
            }

            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    normflux += normal[i] * PermTensor(i, j) * dres(j, 0);
                }
            }
            v2[0] = -normflux;
        }
    }

        if (primal) {
            switch (bc.Type()) {
                case 0 :            // Dirichlet condition
                    for (in = 0; in < phr_primal; in++) {
                        ef(in, 0) += (STATE) (fBigNumber * phi_u(in, 0) * weight) * v2[0];
                        for (jn = 0; jn < phr_primal; jn++) {
                            ek(in, jn) += fBigNumber * phi_u(in, 0) * phi_u(jn, 0) * weight;
                        }
                    }
                    break;
                case 1 :            // Neumann condition
                    for (in = 0; in < phr_primal; in++) {
                        ef(in, 0) += v2[0] * (STATE) (phi_u(in, 0) * weight);
                    }
                    break;
                case 2 :        // mixed condition
                    for (in = 0; in < phr_primal; in++) {
                        ef(in, 0) += v2[0] * (STATE) (phi_u(in, 0) * weight);
                        for (jn = 0; jn < phi_u.Rows(); jn++) {
                            ek(in, jn) += bc.Val1()(0, 0) * (STATE) (phi_u(in, 0) * phi_u(jn, 0) *
                                                                     weight);     // peso de contorno => integral de contorno
                        }
                    }
                    break;
                default:
                    DebugStop();
            }
        } else {
            switch (bc.Type()) {
                case 0 :            // Dirichlet condition
                    for (in = 0; in < phr_hybrid; in++) {
                        ef(in, 0) += v2[0] * (STATE) (phi_flux(in, 0) * weight);
                    }
                    break;
                case 1 :            // Neumann condition
                    for (in = 0; in < phr_hybrid; in++) {
                        ef(in, 0) += (STATE) (fBigNumber * phi_flux(in, 0) * weight) * v2[0];
                        for (jn = 0; jn < phr_hybrid; jn++) {
                            ek(in, jn) += fBigNumber * phi_flux(in, 0) * phi_flux(jn, 0) * weight;
                        }
                    }
                    break;
                case 2 :        // mixed condition
                    DebugStop();
                    for (in = 0; in < phr_hybrid; in++) {
                        ef(in, 0) += v2[0] * (STATE) (phi_flux(in, 0) * weight);
                        for (jn = 0; jn < phi_flux.Rows(); jn++) {
                            ek(in, jn) += 1. / bc.Val1()(0, 0) * (STATE) (phi_flux(in, 0) * phi_flux(jn, 0) *
                                                                          weight);     // peso de contorno => integral de contorno
                        }
                    }
                    break;
            }
        }
}

void TPZMatLaplacianHybrid::VectorialProd(TPZVec<REAL> & ivec, TPZVec<REAL> & jvec, TPZVec<REAL> & kvec, bool unitary)
{
    kvec.Resize(3);
    kvec[0] =  ivec[1]*jvec[2] - ivec[2]*jvec[1];
    kvec[1] = -ivec[0]*jvec[2] + ivec[2]*jvec[0];
    kvec[2] =  ivec[0]*jvec[1] - ivec[1]*jvec[0];

    if(unitary)
    {
        REAL size = 0.;
        int i;
        for(i = 0; i < 3; i++)size += kvec[i] * kvec[i];
        size = sqrt(size);
        //if(size <= 1.e-9)PZError << "\nTPZInterpolationSpace::VectorialProd - null result\n";
        for(i = 0; i < 3; i++)kvec[i] /= size;
    }
}

void TPZMatLaplacianHybrid::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout)
{
    /**
     datavec[1] L2 mesh
     datavec[0] Hdiv mesh,
     datavec[2] Interface Mesh
     datavec[3] Interface Mesh
     **/

    if(var == 0)
    {
        TPZMatLaplacian::Solution(datavec[1],var,Solout);
        return;
    }

    STATE perm = GetPermeability(datavec[1].x);

    if(var==46){
        Solout[0] = perm;
        return;
    }

    TPZManVector<STATE,2> pressexact(1,0.);
    TPZFNMatrix<9,STATE> grad(fDim,1,0.);

    if(TPZMatErrorCombinedSpaces<STATE>::HasExactSol())
    {
        TPZMatErrorCombinedSpaces<STATE>::ExactSol()(datavec[1].x, pressexact,grad);
    }

    switch (var)
    {
        case 44://PressureFem
            Solout[0] = datavec[1].sol[0][0];
            break;
        case 45://Pressure Exact
            Solout[0] = pressexact[0];
            break;
        case 7:
        case 10:
        case 23:
            TPZMatErrorSingleSpace<STATE>::SetExactSol(TPZMatErrorCombinedSpaces<STATE>::ExactSol(), TPZMatErrorCombinedSpaces<STATE>::PolynomialOrderExact());
            TPZMatLaplacian::Solution(datavec[1],var,Solout);
            if(var==10){
                for(int i=0;i<Solout.size();i++){
                    Solout[i]*=-perm;
                }
            }
            break;
        default:
            DebugStop();
    }
}

void TPZMatLaplacianHybrid::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors)
{
    if(!TPZMatErrorCombinedSpaces<STATE>::HasExactSol()) return;
    
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);
    
    TPZManVector<STATE> u_exact(1);
    TPZFNMatrix<9,STATE> du_exact,Kgradu;
    
    
    TPZMatErrorCombinedSpaces<STATE>::ExactSol()(data[1].x,u_exact,du_exact);

    REAL pressure = data[1].sol[0][0];

    // errors[0] norm L2 || u ||_l2
    
    errors[0] = (pressure-u_exact[0])*(pressure-u_exact[0]);//exact error pressure
    
    // errors[1] Semi norm H1 || grad u ||_l2
    
    TPZManVector<STATE,3> sol(1),dsol(3,0.);

    STATE perm = GetPermeability(data[1].x);
    STATE invperm = 1./perm;;
    TPZFMatrix<REAL> &dsolaxes = data[1].dsol[0];
    TPZFNMatrix<9,REAL> gradUh(3,0),KgradUh(3,0);
    TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, gradUh, data[1].axes);
    KgradUh= gradUh*perm;
    Kgradu = du_exact*perm;
    
    for(int id=0; id<fDim; id++) {
        REAL diff = fabs(du_exact(id,0) - gradUh(id,0));
        errors[1]  += diff*diff;
    }
    
    errors[2] = errors[0] +errors[1];
    
    
    
    REAL energy = 0.;
    for (int i=0; i<fDim; i++) {
            energy += (KgradUh(i,0) - Kgradu(i,0))*invperm*(KgradUh(i,0) - Kgradu(i,0));
    }
    
    errors[3] = energy;
}