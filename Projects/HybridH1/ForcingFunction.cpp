//
// Created by victor on 22/07/2021.
//

#include "ForcingFunction.h"

void LinearFunc( const TPZVec<REAL> &x, TPZVec<REAL> &g, TPZFMatrix<REAL> &dg)
{
    dg.Resize(3,1);
    dg(0,0) = 0;
    dg(1,0) = 0.5-x[0];
    dg(2,0) = 0;
}

void SingularityExact( const TPZVec<REAL> &x, TPZVec<REAL> &g, TPZFMatrix<REAL> &dg)
{
    double theta = atan2(x[1], x[0]);//theta=arctan(y/x)
    auto thetaval = shapeFAD::val(theta);
    if (thetaval < (0.)) theta += 2. * M_PI;

    double r = sqrt(x[0]*x[0]+x[1]*x[1]);

    // Verification to avoid numerical errors when x > 0 and y = 0
    if (x[0] > 0 && x[1] < 1e-15 && x[1] > -1.e-15) {
        g[0] = 0.;
    }

    else {
        double factor = pow(r,2./3.);//pow(r,TVar (2.)/TVar (3.))-pow(r,TVar (2.));//
        g[0] = factor * (sin( (2.) * theta / (3.)));
    }

    dg.Resize(3,1);
    //double denominator = 3*pow(x[0],1./3.);//3*pow(x[0]*x[0]+x[1]*x[1],2./3.);
    //if (denominator < 1e-15 ) {
    //   denominator = 1e-15;
    //}
    double denominator = 3.*pow(x[0]*x[0]+x[1]*x[1],2./3.);
    if (denominator < 1e-15 &&  denominator > -1e-15) {
        if (denominator >= 0.) denominator = 1e-15;
        else denominator = -1e-15;
    }
    dg(0,0) = 2.*(-x[1]*cos(2.*theta/3.)+x[0]*sin(2.*theta/3.))/denominator;
    dg(1,0) = 2.*(x[0]*cos(2.*theta/3.)+x[1]*sin(2.*theta/3.))/denominator;
    //dg(0,0) = 2.*sin(2./3.*theta)/denominator;
    //dg(1,0) = 2./3.*pow(x[0]*x[0]+x[1]*x[1],1./3.)*cos(2./3.*theta);
    dg(2,0) = 0;

    //dg(0,0) *=-1.;
    //dg(1,0) *=-1.;
}

void SingularityForcingFunction( const TPZVec<REAL> &x, TPZVec<REAL> &f){
    f[0] = 0.;
}

/*void SingularityExact( const TPZVec<REAL> &x, TPZVec<REAL> &g, TPZFMatrix<REAL> &dg)
{
    double x2 = x[0]*x[0];
    double y2 = x[1]*x[1];
    double r = sqrt(x2+y2);
    g.Resize(1);
    g[0] = pow(r,1./4.)*(1.-x2)*(1.-y2);

    dg.Resize(3,1);
    dg(0,0) = x[0]*(-1+y2)*(-1+201*x2+200*y2)/(100*pow(x2+y2,199/200));
    dg(1,0) = x[1]*(-1+x2)*(-1+200*x2+201*y2)/(100*pow(x2+y2,199/200));
    dg(2,0) = 0;
}

void SingularityForcingFunction( const TPZVec<REAL> &x, TPZVec<REAL> &f){
    f.Resize(1);
    double x2 = x[0]*x[0];
    double y2 = x[1]*x[1];
    double x2y2 = x2+y2;
    double x2y2pow = pow(x2y2,199./200.);

    f[0] = (-1-20000*x2*x2+40401*y2-20000*y2*y2+x2*(40401-40801*y2))/(10000*x2y2pow);
}*/

 /*
 * void SingularityExact( const TPZVec<REAL> &x, TPZVec<REAL> &g, TPZFMatrix<REAL> &dg)
{
    double x2 = x[0]*x[0];
    double y2 = x[1]*x[1];
    double r = sqrt(x2+y2);
    g.Resize(1);
    g[0] = pow(r,1./100.)*(1.-x2)*(1.-y2);

    dg.Resize(3,1);
    dg(0,0) = x[0]*(-1+y2)*(-1+201*x2+200*y2)/(100*pow(x2+y2,199/200));
    dg(1,0) = x[1]*(-1+x2)*(-1+200*x2+201*y2)/(10*pow(x2+y2,199/200));
    dg(2,0) = 0;
}

void SingularityForcingFunction( const TPZVec<REAL> &x, TPZVec<REAL> &f){
    f.Resize(1);
    double x2 = x[0]*x[0];
    double y2 = x[1]*x[1];
    double x2y2 = x2+y2;
    double x2y2pow = pow(x2y2,199./200.);

    f[0] = (-1-20000*x2*x2+40401*y2-20000*y2*y2+x2*(40401-40801*y2))/(10000*x2y2pow);
}
void SingularityExact( const TPZVec<REAL> &x, TPZVec<REAL> &g, TPZFMatrix<REAL> &dg)
{
    double x2 = x[0]*x[0];
    double y2 = x[1]*x[1];
    double r = sqrt(x2+y2);
    g.Resize(1);
    g[0] = pow(r,2./3.)*(1.-x2)*(1.-y2);

    dg.Resize(3,1);
    dg(0,0) = 2*x[0]*(-1+y2)*(-1+4*x2+3*y2)/(3*pow(x2+y2,2/3));
    dg(1,0) = 2*x[1]*(-1+x2)*(-1+3*x2+4*y2)/(3*pow(x2+y2,2/3));
    dg(2,0) = 0;
}

void SingularityForcingFunction( const TPZVec<REAL> &x, TPZVec<REAL> &f){
    f.Resize(1);
    double x2 = x[0]*x[0];
    double y2 = x[1]*x[1];
    double x2y2 = x2+y2;
    double x2y2pow = pow(x2y2,2./3.);

    f[0] = -2*(2+9*x2*x2-32*y2+9*y2*y2+4*x2*(-8+11*y2))/(9*x2y2pow);
}*/