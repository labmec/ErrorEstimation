/**
 * @file
 * @brief Contains the implementation of the TPZNacaProfile methods. 
 */
#include "tpznacaprofile.h"

#include "pzreal.h"
#include "pzvec.h"
#include "TPZGeoElement.h"
#include "pzshapequad.h"
#include "pzrefquad.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"

#include <iostream>
#include <fstream>

using namespace pzgeom;
using namespace pzshape;
using namespace pzrefine;

TPZNacaProfile::TPZNacaProfile()
{
}

TPZNacaProfile::TPZNacaProfile(REAL cord, int FourDigits, REAL angle, TPZVec<REAL> &x0) :
    TPZBlendNACA(cord, FourDigits, angle, x0)
{
    fX0[0] = x0[0];
    fX0[1] = x0[1];
    fX0[2] = x0[2];
    fP = P<REAL>();
    fM = M<REAL>();
    fTT = TT<REAL>();
}

TPZNacaProfile::~TPZNacaProfile()
{
}

template<class Type>
Type TPZNacaProfile::P()
{
    int aux = fFourDigits/100;
    aux -= ((int)(aux/10))*10;
    return (Type)(aux/10.);
}
template REAL TPZNacaProfile::P();

template<class Type>
Type TPZNacaProfile::M()
{
    int aux = fFourDigits/1000;
    return (Type)(aux/100.)*fCord;
}
template REAL TPZNacaProfile::M();


template<class toto>
toto TPZNacaProfile::TT()
{
    int aux = fFourDigits - ((int)(fFourDigits/100))*100;
    toto result;
    result = (toto)(aux/100.)*fCord;
    // std::cout << "TT = " << result << std::endl;
    return result;
}
template REAL TPZNacaProfile::TT();


template <class toto>
toto TPZNacaProfile::yc(toto t) const
{
    toto x = xp(t);
    return TPZBlendNACA::yc(x);
}

template REAL TPZNacaProfile::yc(REAL x) const;

//*******
template <class toto>
toto TPZNacaProfile::tgphi(toto t) const
{
    toto x = xp(t);
    return TPZBlendNACA::tgphi(x);
}
template REAL TPZNacaProfile::tgphi(REAL t) const;

//*******
template <class toto>
toto TPZNacaProfile::dtgphi(toto t) const
{
    toto x = xp(t);
    toto result = TPZBlendNACA::dtgphi(x)*dxp(t);
    return result;
}
template REAL TPZNacaProfile::dtgphi(REAL t) const;

//********
template <class toto>
toto TPZNacaProfile::yt(toto t) const
{
    toto x = t*t/fCord;
    toto t2 = t/fCord;
    toto aux = t2*t2;
    const REAL a0 =  1.4845,
	a1 = -0.6300,
	a2 = -1.7580,
	a3 =  1.4215,
//	a4 = -0.5075;
	a4 = -0.518;
	
    toto val = fTT * (a0*t2 + a1*aux + a2*aux*aux + a3*aux*aux*aux + a4*aux*aux*aux*aux);
    // std::cout << " t " << t << " x " << x << " yt = " << val << std::endl;
    // TPZBlendNACA::yt(x);
    return val;
}

//********
template <class toto>
toto TPZNacaProfile::dyt(toto t) const
{
    toto t2 = t/fCord;
    toto aux = t2*t2;
    const REAL a0 =  1.4845,
	a1 = -0.6300,
	a2 = -1.7580,
	a3 =  1.4215,
//	a4 = -0.5075;
	a4 = -0.518;
	
    return fTT/fCord * (a0 + (a1 + 2.*a2*aux + 3.*a3*aux*aux + 4.*a4*aux*aux*aux)*2.*t2);
}

template REAL TPZNacaProfile::dyt(REAL x) const;

//********
template <class toto>
toto TPZNacaProfile::xu(toto t) const
{
    toto x = xp(t);
    toto val = x-yt(t)*sin(atan(tgphi(t)));
    return val;
    // toto x = xp(t);
    // std::cout << "t = " << t << " x = " << x << std::endl;
    // toto val = TPZBlendNACA::xu(x);
    // std::cout << "t = " << t << " x = " << x << " xu = " << val << std::endl;
}

template REAL TPZNacaProfile::xu(REAL t) const;
//*********
template <class toto>
toto TPZNacaProfile::xl(toto t) const
{
    toto x = xp(t);
    return x+yt(t)*sin(atan(tgphi(t)));
}
template REAL TPZNacaProfile::xl(REAL t) const;


//********
template <class toto>
toto TPZNacaProfile::dxu(toto t) const
{
    toto x = xp(t);
    toto dx = dxp(t);
    toto tg = tgphi(t);
    toto val = dx-(dyt(t)*sin(atan(tg))+
        yt(t)/(1+tg*tg)/sqrt(1.+tg*tg)*dtgphi(t));
    return val;

    // toto x = xp(t);
    // toto dxuv = 0.;
    // if(x > 0.) {
    //     dxuv = TPZBlendNACA::dxu(x);
    // }
    // auto dxpv = dxp(t);
    // if(dxpv == 0.) dxuv = 0.;
    // auto resp = dxuv*dxpv;
    // // std::cout << " x = " << x << " dxu = " << dxuv << " dxp = " << dxpv << std::endl;
    // return resp;
}
template REAL TPZNacaProfile::dxu(REAL t) const;
//*********
template <class toto>
toto TPZNacaProfile::dxl(toto t) const
{
    toto x = xp(t);
    toto dx = dxp(t);
    toto tg = tgphi(t);
    toto val = dx+dyt(t)*sin(atan(tg))+
        yt(t)/(1+tg*tg)/sqrt(1.+tg*tg)*dtgphi(t);
    return val;
}
template REAL TPZNacaProfile::dxl(REAL t) const;



template <class toto>
toto TPZNacaProfile::yu(toto t) const
{
    return yc(t) + yt(t)*cos(atan(tgphi(t)));
}
template REAL TPZNacaProfile::yu(REAL x) const;
//********

template <class toto>
toto TPZNacaProfile::dyu(toto t) const
{
    auto loctgphi = tgphi(t);
    auto resp =  dyc(t)+dyt(t)*cos(atan(loctgphi))
        -yt(t)*loctgphi/(1+loctgphi*loctgphi)/sqrt(1.+loctgphi*loctgphi)*dtgphi(t);
    return resp;
}
template REAL TPZNacaProfile::dyu(REAL t) const;
//********

//********
template <class toto>
toto TPZNacaProfile::yl(toto t) const
{
    return yc(t) - yt(t)*cos(atan(tgphi(t)));
}
template REAL TPZNacaProfile::yl(REAL t) const;

//********
template <class toto>
toto TPZNacaProfile::dyl(toto t) const
{
    toto tg = tgphi(t);
    return dyc(t) - dyt(t)*cos(atan(tg))+
        yt(t)*tg/(1+tg*tg)/sqrt(1.+tg*tg)*dtgphi(t);
}
template REAL TPZNacaProfile::dyl(REAL t) const;


//********

template<class toto>
void TPZNacaProfile::NearestParameter(TPZVec<REAL> &pt, int &uplow, int maxPt,toto& Result)
{
    REAL distminlow = 20.*fCord;
    REAL distminup  = 20.*fCord;
    REAL distlow,distup;
    REAL ptl[2],ptu[2];
    int ip,maxp=maxPt;
    REAL par,parlow = (REAL)0.0;
	REAL parup = (REAL)0.0;
    for(ip=0; ip<=maxp; ip++)
    {
        par = ip*fCord/maxp;
        ptu[0] = xua(par);
        ptu[1] = yua(par);
        ptl[0] = xla(par);
        ptl[1] = yla(par);
        distlow = (ptl[0]-pt[0])*(ptl[0]-pt[0])+(ptl[1]-pt[1])*(ptl[1]-pt[1]);
        distup  = (ptu[0]-pt[0])*(ptu[0]-pt[0])+(ptu[1]-pt[1])*(ptu[1]-pt[1]);
        if(distlow < distminlow) parlow = par;
        if(distup < distminup) parup = par;
        distminlow = distminlow < distlow ? distminlow : distlow;
        distminup = distminup < distup ? distminup : distup;
    }
    REAL delpar = (0.1L)/maxp;
    if(distminlow < distminup)
    {
        uplow = -1;
        REAL distprev = distminlow;
        par = parlow;
        while(fabs(delpar) > 0.00001/maxp)
        {
            ptl[0] = xla(par+delpar);
            ptl[1] = yla(par+delpar);
            distlow = (ptl[0]-pt[0])*(ptl[0]-pt[0])+(ptl[1]-pt[1])*(ptl[1]-pt[1]);
            if(distlow < distprev)
            {
                par += delpar;
                distprev = distlow;
            }
            else if (delpar > 0.)
            {
                delpar *= -1.;
            }
            else
            {
                delpar *= -0.1;
            }
        }
    } 
    else 
    {
        uplow = 1;
        REAL distprev = distminup;
        par = parup;
        while(fabs(delpar) > 0.001/maxp)
        {
            ptu[0] = xua(par+delpar);
            ptu[1] = yua(par+delpar);
            distup = (ptu[0]-pt[0])*(ptu[0]-pt[0])+(ptu[1]-pt[1])*(ptu[1]-pt[1]);
            if(distup < distprev)
            {
                par += delpar;
                distprev = distup;
            }
            else if (delpar > 0.)
            {
                delpar *= -1.;
            }
            else
            {
                delpar *= -0.1;
            }
        }
    }
    Result = par;
}
template void TPZNacaProfile::NearestParameter(TPZVec<REAL> &pt, int &uplow, int maxPt,REAL& Result);

///************************************ public methods
template <class toto>
toto TPZNacaProfile::xua(toto t) const
{

    return fX0[0]+(xu(t)-fCord/2.)*cos(fAngle) + yu(t) * sin(fAngle) + fCord/2.;
}
template REAL TPZNacaProfile::xua(REAL t) const;

//***************

template <class toto>
toto TPZNacaProfile::yua(toto t) const
{
    toto val = fX0[1]+yu(t)*cos(fAngle) - (xu(t)-fCord/2.) * sin(fAngle);
    // std::cout << "yua = " << val << std::endl;
    return val;
}
template REAL TPZNacaProfile::yua(REAL t) const;

//************
template <class toto>
toto TPZNacaProfile::xla(toto t) const
{
    return fX0[0]+(xl(t)-fCord/2.)*cos(fAngle) + yl(t) * sin(fAngle) + fCord/2.;
}
template REAL TPZNacaProfile::xla(REAL t) const;

//************
template <class toto>
toto TPZNacaProfile::yla(toto t) const
{
    return fX0[1]+yl(t)*cos(fAngle) - (xl(t)-fCord/2.) * sin(fAngle);
}
template REAL TPZNacaProfile::yla(REAL) const;

template <class toto>
toto TPZNacaProfile::dxua(toto t) const
{
    toto dxuval = dxu(t);
    toto dyuval = dyu(t);
    // std::cout << " t " << t << " dxu " << dxuval << " dyu " << dyuval << std::endl;
    return fX0[0]+(dxuval*cos(fAngle) + dyuval * sin(fAngle));
}
template REAL TPZNacaProfile::dxua(REAL t) const;

//***************

template <class toto>
toto TPZNacaProfile::dyua(toto t) const
{
    return dyu(t)*cos(fAngle) - (dxu(t)) * sin(fAngle);
}
template REAL TPZNacaProfile::dyua(REAL t) const;

//************
template <class toto>
toto TPZNacaProfile::dxla(toto t) const
{
    return (dxl(t))*cos(fAngle) + dyl(t) * sin(fAngle);
}
template REAL TPZNacaProfile::dxla(REAL t) const;

//************
template <class toto>
toto TPZNacaProfile::dyla(toto t) const
{
    return dyl(t)*cos(fAngle) - (xl(t)) * sin(fAngle);
}
template REAL TPZNacaProfile::dyla(REAL) const;
//template REAL TPZNacaProfile::yla(REAL x);

///************************************
template <class toto>
void TPZNacaProfile::ProjectPoint(TPZVec<toto> &pt, int maxPt)
{
    int uplow;
    
    toto par;
    NearestParameter(pt,uplow, maxPt,par);
    
    if(uplow == 0)
    {
        pt[0] = xla(par);
        pt[1] = yla(par);
    }
    else
    {
        pt[0] = xua(par);
        pt[1] = yua(par);
    }
}
template<> void TPZNacaProfile::ProjectPoint(TPZVec<REAL> &pt, int maxPt);

int TPZNacaProfile::ClassId() const{
    return Hash("TPZNacaProfile") ^ pzgeom::TPZNodeRep<2,pztopology::TPZLine>::ClassId() << 1;
}

#include "tpzgeoelmapped.h"

/**
 * Creates a geometric element according to the toto of the father element
 */

// TPZGeoEl *TPZNacaProfile::CreateGeoElement(TPZGeoMesh &mesh, MElementtoto toto,
// 										 TPZVec<int64_t>& nodeindexes,
// 										 int matid,
// 										 int64_t& index)
// {
// 	return CreateGeoElementMapped(mesh,toto,nodeindexes,matid,index);
// }
//template<> class TPZGeoElRefPattern<pzgeom::TPZNacaProfile>;


template class
TPZRestoreClass< TPZGeoElRefPattern<TPZNacaProfile>>;


