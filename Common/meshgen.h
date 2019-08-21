
#ifndef MESHGEN
#define MESHGEN

#include "pzvec.h"
#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzfunction.h"
#include "pzcompel.h"
#include "TPZAnalyticSolution.h"
#include <string>

class TPZGeoMesh;

struct TRunConfig;

TPZGeoMesh *MalhaGeomFredQuadrada(int nelx, int nely, TPZVec<REAL> &x0, TPZVec<REAL> &x1, TPZVec<int64_t> &coarseindices, int ndiv);

/// Solve the problem composed of a multiphysics mesh composed of compmeshes - applies to MHM and MHM-H(div)
void SolveProblem(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<TPZAutoPointer<TPZCompMesh> > compmeshes, TPZAnalyticSolution &analytic, std::string prefix, TRunConfig config);

/// Solve the problem composed of a multiphysics mesh composed of compmeshes - applies to MHM and MHM-H(div)
void SolveParabolic(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<TPZAutoPointer<TPZCompMesh> > compmeshes, std::string prefix, TRunConfig config);

struct TRunConfig
{
    int nelxcoarse = -1;
    int nelycoarse = -1;
    int numHDivisions = 0;
    int pOrderInternal = 1;
    int hdivmaismais = 1;
    int numDivSkeleton = 0;
    int pOrderSkeleton = 1;
    int Hybridize = 0;
    int Condensed = 1;
    int LagrangeMult = 0;
    int newline = 0;
    int n_threads = 0;
    bool MHM_HDiv_Elast = false;
    
    /// directory where the files will be stored
    std::string dir_name = ".";
    
    /// number of equations when not condensing anything
    int64_t fGlobalSystemSize = -1;
    /// number of equations considering local condensation
    int64_t fGlobalSystemWithLocalCondensationSize = -1;
    /// number of equations of the global system
    int64_t fNumeq = -1;

    REAL fDeltaT = 1.;
    
    /// number of timesteps
    int64_t nTimeSteps = 10;
    
    std::ostream &InlinePrint(std::ostream &out)
    {
//        out << "nelxCoarse " << nelxcoarse << " nelyCoarse " << nelycoarse << " numHDiv " << numHDivisions << " porderInternal " << pOrderInternal << " numDivSkeleton " << numDivSkeleton
//        << " porderSkeleton " << pOrderSkeleton << " Hybridize " << Hybridize << " Condensed " << Condensed << " LagrangeMult " << LagrangeMult
//        << " sysnocondense " << fGlobalSystemSize << " syslocalcondense " << fGlobalSystemWithLocalCondensationSize << " neq " << fNumeq;
        // <inputs>: n_dx n_dy k_skeleton m_div k_subelement
        out << "n_dx " << nelxcoarse << " n_dy " << nelycoarse  << " k_skeleton " << pOrderSkeleton << " m_div " << numHDivisions << " k_subelement " << pOrderInternal << " upscaling_dof " << fNumeq << " total_dof " << fGlobalSystemSize;
        return out;
    }
    std::ostream &MathematicaInlinePrint(std::ostream &out)
    {
        out << "nelxCoarse, " << nelxcoarse << ", nelyCoarse, " << nelycoarse << " ,numHDiv, " << numHDivisions << " ,porderInternal, " << pOrderInternal << " ,numDivSkeleton, " << numDivSkeleton
        << " ,porderSkeleton, " << pOrderSkeleton << " ,Hybridize, " << Hybridize << " ,Condensed, " << Condensed << " ,LagrangeMult, " << LagrangeMult
        << " ,sysnocondense, " << fGlobalSystemSize << " ,syslocalcondense, " << fGlobalSystemWithLocalCondensationSize << " ,neq, " << fNumeq;
        return out;
    }
    
    std::ostream &ConfigPrint(std::ostream &out)
    {
        out << nelxcoarse << "x" << nelycoarse << "_HSkel" << numDivSkeleton << "_pSkel" << pOrderSkeleton << "_HDiv" << numHDivisions << "_pInt" << pOrderInternal;
        return out;
    }
};

typedef void (ExactFunc)(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu);

#include <iostream>
//#include <boost/numeric/odeint.hpp>       // odeint function definitions

//using namespace boost::numeric::odeint;

// Defining a shorthand for the type of the mathematical state
typedef std::vector< double > state_type;

//static TPZManVector<REAL,3> g_qsi(2.0);
//static int64_t g_elindex;
//static TPZGeoMesh *g_gmesh = 0;

struct TStreamLinePoint
{
    TPZManVector<REAL,4> fx;
    REAL fTime;
    TStreamLinePoint(): fx(4,0.), fTime(0.) {}
    TStreamLinePoint(const TStreamLinePoint &cp) : fx(cp.fx), fTime(cp.fTime){}
    TStreamLinePoint &operator=(const TStreamLinePoint &cp)
    {
        fx = cp.fx;
        fTime = cp.fTime;
        return *this;
    }
};

struct TStreamLine
{
    TPZGeoMesh *fGMesh;
    TPZManVector<REAL,3> fQsi;
    int64_t fElIndex;
    int fVarIndex = 0;
    
    TStreamLine() : fGMesh(0), fQsi(0), fElIndex(0), fVarIndex(0) {}
    
    TStreamLine(TPZGeoMesh *gmesh) : fGMesh(gmesh), fElIndex(0)
    {
        fQsi.resize(gmesh->Dimension());
    }
    
    TStreamLine(const TStreamLine &cp) : fGMesh(cp.fGMesh), fQsi(cp.fQsi), fElIndex(cp.fElIndex),
    fVarIndex(cp.fVarIndex) {}
    
    TStreamLine &operator=(const TStreamLine &cp)
    {
        fGMesh = cp.fGMesh;
        fQsi = cp.fQsi;
        fElIndex = cp.fElIndex;
        fVarIndex = cp.fVarIndex;
        return *this;
    }
    
    void operator()(const state_type &x , state_type &dxdt , const double t )
    {
        
        //    std::cout << "x = " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl;
        TPZManVector<REAL,3> xel(3,0.);
        for(int i=0; i<3; i++) xel[i] = x[i];
        TPZGeoEl *gel = fGMesh->FindElementCaju(xel, fQsi, fElIndex, 2);
        //    std::cout << " elindex " << g_elindex << " qsi " << g_qsi << std::endl;
        if(gel)
        {
            TPZCompEl *cel = gel->Reference();
            TPZManVector<STATE,3> flux(3,0.);
            cel->Solution(fQsi, fVarIndex, flux);
            //    std::cout << " flux " << flux << std::endl;
            REAL fluxnorm = sqrt(flux[0]*flux[0]+flux[1]*flux[1]+flux[2]*flux[2]);
            if(fluxnorm > 1.e-8)
            {
                for(int i=0; i<3; i++) dxdt[i] = flux[i] / fluxnorm;
                dxdt[3] = 1./fluxnorm;
            }
            else
            {
                for(int i=0; i<3; i++) dxdt[i] = flux[i];
                dxdt[3] = 0.;
            }
        }
        else
        {
            for(int i=0; i<4; i++) dxdt[i] = 0.;
        }
    }

    
};

struct TStreamLineData
{
    std::shared_ptr<TPZStack<TStreamLinePoint>> fStreamLine;
    
    TStreamLineData() : fStreamLine(new TPZStack<TStreamLinePoint>){}
    
    TStreamLineData(const TStreamLineData &cp) : fStreamLine(cp.fStreamLine){}
    
    TStreamLineData &operator=(const TStreamLineData &cp){
        fStreamLine = cp.fStreamLine;
        return *this;
    }
    
    void operator()( const state_type &x, const double t )
    {
        TStreamLinePoint pt;
        for(int i=0; i<4; i++) pt.fx[i] = x[i];
        pt.fTime = t;
        std::cout << "Adding time " << t << " x " << pt.fx << " size " << fStreamLine->size() << std::endl;
        fStreamLine->Push(pt);
    }
    
    void Print(std::ostream &out)
    {
        for(int i=0; i<fStreamLine->size(); i++)
        {
            out << " time " << (*fStreamLine)[i].fTime << " x " << (*fStreamLine)[i].fx << std::endl;
        }
    }
    
    operator TPZFMatrix<REAL> (){
        TPZFMatrix<REAL> result(fStreamLine->size(),5);
        for(int i=0; i<fStreamLine->size(); i++)
        {
            result(i,0) = (*fStreamLine)[i].fTime;
            for(int ic = 0; ic<4; ic++) result(i,ic+1) = (*fStreamLine)[i].fx[ic];
        }
        return result;
    }
};

TPZGeoEl *FindEntry(TPZGeoMesh *gmesh);

TStreamLineData ComputeStreamLine(TPZCompMesh *fluxmesh, TPZVec<REAL> &xco);

#include "TPZAnalyticSolution.h"

#ifndef TPZANALYTICSOLUTION

/*
struct TAnalyticSolution
{
    
    
    virtual TPZAutoPointer<TPZFunction<STATE> > ForcingFunction() = 0;
    
    virtual TPZAutoPointer<TPZFunction<STATE> > ValueFunction() = 0;
    
    virtual TPZAutoPointer<TPZFunction<STATE> > ConstitutiveLawFunction() = 0;
    
    virtual ExactFunc *Exact() = 0;
    
    virtual ~TAnalyticSolution()
    {
        
    }
};


#ifdef _AUTODIFF

struct TElasticityExample1 : public TAnalyticSolution
{
    static void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force)
    {
        TPZManVector<REAL,3> locforce(2);
        DivSigma(x, locforce);
        force[0] = -locforce[0];
        force[1] = -locforce[1];
    }
    
    virtual TPZAutoPointer<TPZFunction<STATE> > ForcingFunction();
    
    virtual TPZAutoPointer<TPZFunction<STATE> > ValueFunction();
    
    virtual TPZAutoPointer<TPZFunction<STATE> > ConstitutiveLawFunction();
    
    static void GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu);
    
    virtual ExactFunc *Exact()
    {
        return GradU;
    }
    
    virtual ~TElasticityExample1()
    {
        
    }
    
    template<class TVar>
    static void Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma);
    
    template<class TVar>
    static void DivSigma(const TPZVec<TVar> &x, TPZVec<TVar> &divsigma);
    
    template<class TVar>
    static void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp);
    
    static void Dirichlet(const TPZVec<REAL> &x, TPZVec<STATE> &disp)
    {
        TPZManVector<REAL,3> disploc(2,0.);
        uxy(x,disploc);
        for(int i=0; i<2; i++) disp[i] = disploc[i];
    }
    
    template<class TVar>
    static void graduxy(const TPZVec<TVar> &x, TPZFMatrix<TVar> &grad);

    template<class TVar>
    static void Elastic(const TPZVec<TVar> &x, TVar &Elast, TVar &nu);

    static void ElasticDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv);

};

struct TLaplaceExample1 : public TAnalyticSolution
{
    virtual TPZAutoPointer<TPZFunction<STATE> > ForcingFunction();
    
    virtual TPZAutoPointer<TPZFunction<STATE> > ValueFunction();
    
    virtual TPZAutoPointer<TPZFunction<STATE> > ConstitutiveLawFunction();
    
    virtual ~TLaplaceExample1()
    {
        
    }
    
    static void GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu);
    
    virtual ExactFunc *Exact()
    {
        return GradU;
    }

    template<class TVar>
    static void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp);
    
    template<class TVar>
    static void graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad);
    
    template<class TVar>
    static void Permeability(const TPZVec<TVar> &x, TVar &Elast);

    static void PermeabilityDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv);
    
    template<class TVar>
    static void Sigma(const TPZVec<TVar> &x, TPZVec<TVar> &sigma);
    
    template<class TVar>
    static void DivSigma(const TPZVec<TVar> &x, TVar &divsigma);
    
    static void Dirichlet(const TPZVec<REAL> &x, TPZVec<STATE> &disp)
    {
        TPZManVector<REAL,3> disploc(2,0.);
        uxy(x,disploc);
        for(int i=0; i<1; i++) disp[i] = disploc[i];
    }
    
    static void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force)
    {
        REAL locforce;
        DivSigma(x, locforce);
        force[0] = locforce;
    }

};

struct TLaplaceExampleSmooth : public TAnalyticSolution
{
    virtual TPZAutoPointer<TPZFunction<STATE> > ForcingFunction();
    
    virtual TPZAutoPointer<TPZFunction<STATE> > ValueFunction();
    
    virtual TPZAutoPointer<TPZFunction<STATE> > ConstitutiveLawFunction();
    
    virtual ~TLaplaceExampleSmooth()
    {
        
    }
    
    static void GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu);
    
    virtual ExactFunc *Exact()
    {
        return GradU;
    }
    
    template<class TVar>
    static void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp);
    
    template<class TVar>
    static void graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad);
    
    template<class TVar>
    static void Permeability(const TPZVec<TVar> &x, TVar &Elast);
    
    static void PermeabilityDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv);
    
    template<class TVar>
    static void Sigma(const TPZVec<TVar> &x, TPZVec<TVar> &sigma);
    
    template<class TVar>
    static void DivSigma(const TPZVec<TVar> &x, TVar &divsigma);
    
    static void Dirichlet(const TPZVec<REAL> &x, TPZVec<STATE> &disp)
    {
        TPZManVector<REAL,3> disploc(2,0.);
        uxy(x,disploc);
        for(int i=0; i<1; i++) disp[i] = disploc[i];
    }
    
    static void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force)
    {
        REAL locforce;
        DivSigma(x, locforce);
        force[0] = locforce;
    }
    
};


 #endif
 
 */

#endif

#endif
