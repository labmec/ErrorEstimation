//
//  New_main.cpp
//  ErrorEstimateHDiv
//
//  Created by Philippe Devloo on 14th March 2020.
//

#include "pzlog.h"
#include "pzgenericshape.h"
#include "pzshapepoint.h"

extern int64_t counter;

#include "boost/date_time/posix_time/posix_time.hpp"

template<class TSHAPE>
void TimeComparaison();
template<class TSHAPE>
void TimeComparaison(TPZVec<int> &orders);

int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
//    TPZManVector<int,27> orders(19,5);
//    TimeComparaison<pzshape::TPZShapeCube>(orders);
//
//    return 1;
    
    TimeComparaison<pzshape::TPZShapeLinear>();
    TimeComparaison<pzshape::TPZShapeTriang>();
    TimeComparaison<pzshape::TPZShapeQuad>();
//    TimeComparaison<pzshape::TPZShapePoint>();


    TimeComparaison<pzshape::TPZShapeTetra>();
    TimeComparaison<pzshape::TPZShapeCube>();
    TimeComparaison<pzshape::TPZShapePrism>();
    TimeComparaison<pzshape::TPZShapePiram>();
    

    return 0;
    
}

template<class TSHAPE>
void TimeComparaison()
{
    const int ncorner = TSHAPE::NCornerNodes;
    const int nsides = TSHAPE::NSides;
    TPZManVector<int,27> orders(nsides-ncorner,1);
    
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::cout << "orders = " << orders << std::endl;
    TimeComparaison<TSHAPE>(orders);
    orders.Fill(2);
    std::cout << "orders = " << orders << std::endl;
    TimeComparaison<TSHAPE>(orders);
    orders.Fill(1);
    if(nsides-ncorner > 1)
    {
        orders[nsides-ncorner-1] = 3;
        std::cout << "orders = " << orders << std::endl;
        TimeComparaison<TSHAPE>(orders);
    }
    if(TSHAPE::Dimension > 1)
    {
        int nsides1 = TSHAPE::NumSides(1);
        orders.Fill(1);
        for(int i=0; i<nsides1; i++) orders[i] = 5;
        std::cout << "orders = " << orders << std::endl;
        TimeComparaison<TSHAPE>(orders);
    }
    if(TSHAPE::Dimension > 2)
    {
        int nsides1 = TSHAPE::NumSides(1);
        int nsides2 = TSHAPE::NumSides(2);
        orders.Fill(1);
        for(int i=nsides1; i<nsides1+nsides2; i++) orders[i] = 5;
        std::cout << "orders = " << orders << std::endl;
        TimeComparaison<TSHAPE>(orders);
    }
    orders.Fill(5);
    std::cout << "orders = " << orders << std::endl;
    TimeComparaison<TSHAPE>(orders);
}

template<class TSHAPE>
void TimeComparaison(TPZVec<int> &orders)
{
    const int ncorner = TSHAPE::NCornerNodes;
    const int nsides = TSHAPE::NSides;

    TParDefs par;
    TPZManVector<TPZTransform<REAL>,27> transforms;
    TPZManVector<int64_t,27> ids = {8,7,6,5,4,3,2,1,0};
    ComputeTransforms<TSHAPE>(ids, par.transvec);
    TPZManVector<REAL,3> pt = {1./sqrt(5.),1./M_PI,1./(4.*M_PI)};
    pt.Resize(TSHAPE::Dimension);
    par.orders = orders;
    int numshape = TSHAPE::NShapeF(orders);
    TPZManVector<int,27> nshape(TSHAPE::NSides-TSHAPE::NCornerNodes);
        
    for(int is = ncorner; is< nsides; is++)
    {
        nshape[is-ncorner] = TSHAPE::NConnectShapeF(is, orders[is-ncorner]);
    }
    std::cout << "nshape = " << nshape << std::endl;
    par.nshape = nshape;
    //    par.phi.Resize(numshape, 1);
    //    par.dphi.Resize(3, numshape);
    TPZFNMatrix<100,REAL> phi2(numshape,1,0.), dphi2(3,numshape,0.);
    TPZFNMatrix<100,REAL> phi1(numshape,1,0.), dphi1(3,numshape,0.);
    boost::posix_time::ptime::time_system_type::time_duration_type timenew, timeold;
    boost::posix_time::ptime tsim0 = boost::posix_time::microsec_clock::local_time();
    TSHAPE::Shape(pt, ids, orders, phi2, dphi2);
    Shape<TSHAPE>(pt, par, phi1, dphi1);
    TSHAPE::ShapeCorner(pt, phi1, dphi1);
    boost::posix_time::ptime tsim01 = boost::posix_time::microsec_clock::local_time();
    for(int i=0; i<1; i++)
    {
        counter = 0;
        boost::posix_time::ptime tsim1 = boost::posix_time::microsec_clock::local_time();
        for(int j=0; j<100000; j++)
                TSHAPE::Shape(pt, ids, orders, phi2, dphi2);
        std::cout << "Chebyshev counter " << counter << std::endl;
        counter = 0;
        boost::posix_time::ptime tsim2 = boost::posix_time::microsec_clock::local_time();
        for(int j=0; j<100000; j++)
            Shape<TSHAPE>(pt, par, phi1, dphi1);
        std::cout << "Chebyshev counter " << counter << std::endl;
        boost::posix_time::ptime tsim3 = boost::posix_time::microsec_clock::local_time();
        timeold += tsim2-tsim1;
        timenew += tsim3-tsim2;
    }
//    std::cout << "Sum wall time = " << tsim01 - tsim0 << " s" << std::endl;
    std::cout << "Total wall time of Shape old = " << timeold << " s" << std::endl;
    std::cout << "Total wall time of Shape new = " << timenew << " s" << std::endl;
    //    std::cout << "Phi 1 " << phi1 << std::endl;
    phi2 -= phi1;
    dphi2 -= dphi1;
    std::cout << Norm(phi2) << " " << Norm(dphi2) << std::endl;

}
