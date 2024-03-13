/**
 * @file
 * @brief Contains the TPZNacaProfile class. It is a special map.
 */

#ifndef TPZNACAPROFILE_H
#define TPZNACAPROFILE_H

#include <math.h>
#include "pzreal.h"
#include "pzvec.h"
#include "pzeltype.h"
#include "pznoderep.h"
#include "tpzline.h"
#include "tpzblendnaca.h"
#include "pzgeoel.h"
#include "pzgnode.h"

class TPZGeoEl;
class TPZGeoMesh;

/**
 * @author Philippe Devloo
 * @ingroup geometry
 * @brief Special map to NACA. \ref geometry "Geometry"
 */

namespace pzgeom
{

class TPZNacaProfile  : public TPZBlendNACA
{
public:

    REAL firstpos = 0.; // first position of the parametrization
    REAL lastpos = -1.; // last position of the parametrization
    int uplow = 1; // upper or lower profile >0 upper <0 lower

	/** @brief Default constructor */
    TPZNacaProfile();
    /** @brief Constructor */
    TPZNacaProfile(REAL cord, int FourDigits, REAL angle, TPZVec<REAL> &x0);

    TPZNacaProfile(TPZVec<int64_t> &nodeindexes) : TPZRegisterClassId(&TPZNacaProfile::ClassId),
        TPZBlendNACA(nodeindexes) {

    }

    TPZNacaProfile(const TPZNacaProfile &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZNacaProfile::ClassId),
        TPZBlendNACA(cp){
            DebugStop();
        }
        /** @brief Copy constructor with map of nodes */
    TPZNacaProfile(const TPZNacaProfile &cp,std::map<int64_t,int64_t> & gl2lcNdMap) : 
        TPZRegisterClassId(&TPZNacaProfile::ClassId),
        TPZBlendNACA(cp,gl2lcNdMap){
            DebugStop();
        
        }
	/** @brief Default destructor */
    ~TPZNacaProfile();
	
    /// with attack angle
    /// superior profile
    template <class Type>
    Type xua(Type x) const;
    template <class Type>
    Type yua(Type x) const;
	
    template <class Type>
    Type dxua(Type x) const;
    template <class Type>
    Type dyua(Type x) const;
	
    /// inferior profile
    template <class Type>
    Type xla(Type x) const;
    template <class Type>
    Type yla(Type x) const;
    
    /// inferior profile
    template <class Type>
    Type dxla(Type x) const;
    template <class Type>
    Type dyla(Type x) const;
    
    template <class Type>
    void ProjectPoint(TPZVec<Type> &pt, int maxPt = 1000);
	
public:
    
    // declare the methods needed to implement the interface of TPZGeoElRefPattern<TPZNacaProfile>
    

    virtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord) override
    {
        //Dont have
    }
    
    int ClassId() const override;

	/** @brief Creates a geometric element according to the type of the father element */
	// static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
	// 								  TPZVec<int64_t>& nodeindexes,
	// 								  int matid,
	// 								  int64_t& index);
private:
	
    template<class Type>
    Type P();
    template<class Type>
    Type M();
    template<class Type>
    Type TT();
	
    template <class Type>
    Type xp(Type t) const{
        return t*t/fCord;
    }

    template <class Type>
    Type dxp(Type t) const {
        return 2.*t/fCord;
    }
    /** @brief Mean line for the wing */
    template <class Type>
    Type yc(Type x) const;

    /** @brief Tg of the angle of the mean line */
    template <class Type>
    Type tgphi(Type x) const;
	
    /** @brief Derivative of the Tg of the angle of the mean line */
    template <class Type>
    Type dtgphi(Type x) const;
	
    /** @brief Thickness */
    template <class Type>
    Type yt(Type x) const;
	
    /** @brief Derivative of the thickness */
    template <class toto>
    toto dyt(toto x) const;

    /** @brief Superior profile */
    template <class Type>
    Type xu(Type x) const;
    template <class Type>
    // Derivative of the superior profile
    Type dxu(Type x) const;

    template <class Type>
    Type yu(Type x) const;
    template <class Type>
    Type dyu(Type x) const;
	
    /** @brief Inferior profile */
    template <class Type>
    Type xl(Type x) const;
    template <class Type>
    Type yl(Type x) const;
    /** @brief Inferior profile */
    template <class Type>
    Type dxl(Type x) const;
    template <class Type>
    Type dyl(Type x) const;
    
public:
	template <class Type>
    void NearestParameter(TPZVec<REAL> &pt, int &uplow, int maxPt,Type &);

public:

    void Read(TPZStream& buf, void* context) override {
        DebugStop();
    }
    void Write(TPZStream &buf, int withclassid) const override
    {
        DebugStop();
    }

    void Initialize(TPZGeoEl *refel)
    {
        auto firstnode = refel->NodePtr(0);
        TPZManVector<REAL,3> firstcoord(3,0.), lastcoord(3,0.);
        firstnode->GetCoordinates(firstcoord);
        int uplowfirst, uplowlast;
        NearestParameter(firstcoord,uplowfirst,1000,firstpos);
        auto lastnode = refel->NodePtr(1);
        lastnode->GetCoordinates(lastcoord);
        NearestParameter(lastcoord,uplowlast,1000,lastpos);
        REAL valfirst = (uplowfirst-firstpos)*(lastpos-uplowfirst);
        REAL vallast = (uplowlast-firstpos)*(lastpos-uplowlast);
        if(fabs(valfirst) > fabs(vallast))
        {
            uplow = uplowfirst;
        }
        else
        {
            uplow = uplowlast;
        }
        if(uplow == 1) {
            REAL xv = xu(firstpos);
            REAL yv = yu(firstpos);
            if(fabs(xv-firstcoord[0]) > 1.e-6 || fabs(yv-firstcoord[1]) > 1.e-6)
            {
                DebugStop();
            }
        } else {
            REAL xv = xl(firstpos);
            REAL yv = yl(firstpos);
            if(fabs(xv-firstcoord[0]) > 1.e-6 || fabs(yv-firstcoord[1]) > 1.e-6)
            {
                DebugStop();
            }
        }
        if(uplow == 1) {
            REAL xv = xu(lastpos);
            REAL yv = yu(lastpos);
            if(fabs(xv-lastcoord[0]) > 1.e-6 || fabs(yv-lastcoord[1]) > 1.e-6)
            {
                DebugStop();
            }
        } else {
            REAL xv = xl(lastpos);
            REAL yv = yl(lastpos);
             if(fabs(xv-lastcoord[0]) > 1.e-6 || fabs(yv-lastcoord[1]) > 1.e-6)
            {
                DebugStop();
            }
        }

    }
    
    template<class T>
    void X(TPZFMatrix<REAL> &coord,TPZVec<T> &loc,TPZVec<T> &result) const
    {
        T par = firstpos + (lastpos-firstpos)*(loc[0]+1.)/2.;
        result[2] = 0.;
        // std::cout << " firstpos " << firstpos << " lastpos " << lastpos << " loc[0] " << loc[0] << " par " << par << " uplow " << uplow << "\n";
        if(uplow == 1)
        {
            result[0] = xua(par);
            result[1] = yua(par);
        }
        else
        {
            result[0] = xla(par);
            result[1] = yla(par);
        }
        // std::cout << "result " << result << std::endl;
    }   

    template<class T>
    void GradX(TPZFMatrix<REAL> &coord, TPZVec<T> &par, TPZFMatrix<T> &gradx) const
    {
        T locpar = firstpos + (lastpos-firstpos)*(par[0]+1.)/2.;
        gradx.Resize(3,1);
        gradx(2,0) = 0.;
        REAL mult = 2./(lastpos-firstpos);
        if(uplow == 1)
        {
            gradx(0,0) = dxua(locpar)*mult;
            gradx(1,0) = dyua(locpar)*mult;
        }
        else
        {
            gradx(0,0) = dxla(locpar)*mult;
            gradx(1,0) = dyla(locpar)*mult;
        }
    }
    static std::string TypeName() { return "NACA";}
    
    static bool IsLinearMapping(int side)
    {
        return false;
    }

};
};

#endif
