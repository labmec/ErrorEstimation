//
//  New_main.cpp
//  ErrorEstimateHDiv
//
//  Created by Philippe Devloo on 14th March 2020.
//

#include "pzlog.h"
#include "pzgenericshape.h"
#include "pzshapepoint.h"
#include "TPZGeoMeshTools.h"
#include "TPZMixedDarcyFlow.h"
#include "pzelchdiv.h"
#include "tpzautopointer.h"

extern int64_t counter;

#include "boost/date_time/posix_time/posix_time.hpp"

TPZAutoPointer<TPZCompMesh> GenerateCompMesh(MMeshType type, int order);

template<class TSHAPE>
TPZCompElHDiv<TSHAPE> *FindElement(TPZAutoPointer<TPZCompMesh> &cmesh);

template<class TSHAPE>
void TimeComparaisonRequiredData();
template<class TSHAPE>
void TimeComparaisonRequiredData(int order);


MElementType ElType(MMeshType meshtype)
{
    switch (meshtype) {
        case MMeshType::EOneDimensional:
            return EOned;
            break;
        case MMeshType::ETriangular:
            return ETriangle;
            break;
        case MMeshType::EQuadrilateral:
            return EQuadrilateral;
        case MMeshType::ETetrahedral:
            return ETetraedro;
        case MMeshType::EHexahedral:
            return ECube;
        case MMeshType::EPrismatic:
            return EPrisma;
        case MMeshType::EPyramidal:
            return EPiramide;
        case MMeshType::EHexaPyrMixed:
        default:
            DebugStop();
    }
    return MElementType::ENoType;
}

MMeshType MeshType(MElementType eltype)
{
    switch (eltype) {
        case EOned:
            return MMeshType::EOneDimensional;
        case ETriangle:
            return MMeshType::ETriangular;
        case EQuadrilateral:
            return MMeshType::EQuadrilateral;
        case ECube:
            return MMeshType::EHexahedral;
        case ETetraedro:
            return MMeshType::ETetrahedral;
        case EPiramide:
            return MMeshType::EPyramidal;
        case EPrisma:
            return MMeshType::EPrismatic;
            break;
            
        default:
            DebugStop();
            break;
    }
    return MMeshType::ENoType;
}

int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    TimeComparaisonRequiredData<pzshape::TPZShapeLinear>();
    TimeComparaisonRequiredData<pzshape::TPZShapeTriang>();
    TimeComparaisonRequiredData<pzshape::TPZShapeQuad>();
    TimeComparaisonRequiredData<pzshape::TPZShapeTetra>();
    TimeComparaisonRequiredData<pzshape::TPZShapeCube>();
    TimeComparaisonRequiredData<pzshape::TPZShapePrism>();
    

    return 0;
    
}

template<class TSHAPE>
void TimeComparaisonRequiredData()
{
    int order = 1;
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::cout << "order = " << order << std::endl;
    TimeComparaisonRequiredData<TSHAPE>(order);
    order = 2;
    std::cout << "order = " << order << std::endl;
    TimeComparaisonRequiredData<TSHAPE>(order);
    order = 3;
    std::cout << "order = " << order << std::endl;
    TimeComparaisonRequiredData<TSHAPE>(order);
}

template<class TSHAPE>
void TimeComparaisonRequiredData(int order)
{
    TPZAutoPointer<TPZCompMesh> cmesh = GenerateCompMesh(MeshType(TSHAPE::Type()), order);
    TPZCompElHDiv<TSHAPE> *cel = FindElement<TSHAPE>(cmesh);
    if(!cel) DebugStop();
    int nc = cel->NConnects();
    int nsides = cel->Reference()->NSides();
    TPZMixedDarcyFlow *mixed = dynamic_cast<TPZMixedDarcyFlow*>(cmesh->FindMaterial(1));
    if(!mixed) DebugStop();
    TPZManVector<int64_t,27> first_old(nsides+1),first_new(nsides+1);
    cel->FirstShapeIndex(first_old);
    cel->FirstShapeIndex(order+1, first_new);
//    std::cout << "before " << cel->AllSideOrient() << std::endl;
    cel->AllSideOrient().Fill(1);
//    std::cout << "after " << cel->AllSideOrient() << std::endl;
    for(int i=0; i<=nc; i++) if(first_old[i] != first_new[i]) DebugStop();
    for(int i=0; i<nc; i++) if(cel->NConnectShapeF(i, order) != cel->NConnectShapeFNew(i, order)) DebugStop();
    TPZVec<TPZMaterialData> data_old(2),data_new(2);
    data_old[1].fShapeType = TPZMaterialData::EScalarShape;
    data_new[1].fShapeType = TPZMaterialData::EScalarShape;
    TPZManVector<REAL,3> qsi(TSHAPE::Dimension);
    TSHAPE::RandomPoint(qsi);
    boost::posix_time::ptime::time_system_type::time_duration_type timenew, timeold;
    boost::posix_time::ptime tsim0 = boost::posix_time::microsec_clock::local_time();
    cel->InitMaterialData(data_old[0]);
    cel->ComputeRequiredData(data_old[0],qsi);
    cel->InitMaterialDataNew(data_new[0]);
    cel->ComputeRequiredDataNew(data_new[0],qsi);
    int nflux = data_new[0].fVecShapeIndex.size();
    TPZFMatrix<REAL> ek_old(nflux,nflux,0.), ef(nflux,0.), ek_new(nflux,nflux,0.);
    REAL weight = 1.;
    boost::posix_time::ptime tsim01 = boost::posix_time::microsec_clock::local_time();
    mixed->m_fast = true;
    boost::posix_time::ptime tsim1 = boost::posix_time::microsec_clock::local_time();
    for(int i=0; i<10000; i++)
    {
//        cel->InitMaterialDataNew(data_new[0]);
        cel->ComputeRequiredDataNew(data_new[0],qsi);
//        mixed->Contribute(data_new, weight, ek_new, ef);
    }
    boost::posix_time::ptime tsim2 = boost::posix_time::microsec_clock::local_time();
    mixed->m_fast = false;
    for(int i=0; i<10000; i++)
    {
//        cel->InitMaterialData(data_old[0]);
        cel->ComputeRequiredData(data_old[0],qsi);
//        mixed->Contribute(data_old, weight, ek_old, ef);
    }
    boost::posix_time::ptime tsim3 = boost::posix_time::microsec_clock::local_time();
    if(data_old[0].fVecShapeIndex.size() != data_new[0].fVecShapeIndex.size()) DebugStop();
    int nvec = data_old[0].fVecShapeIndex.size();
    int err = 0;
    for (int ivec = 0; ivec<nvec ; ivec++) {
        bool compareabs = false;
        if(TSHAPE::Type() == EQuadrilateral && ivec > nvec-cel->Connect(4).NShape())
        {
            compareabs = true;
        }
        int ind_old = data_old[0].fVecShapeIndex[ivec].first;
        int ind_new = data_new[0].fVecShapeIndex[ivec].first;
        int shape_ind_old = data_old[0].fVecShapeIndex[ivec].second;
        int shape_ind_new = data_new[0].fVecShapeIndex[ivec].second;
        if(shape_ind_old != shape_ind_new) DebugStop();
        TPZManVector<int,3> vec_old(TSHAPE::Dimension),vec_new(TSHAPE::Dimension);
        for(int i=0; i<TSHAPE::Dimension; i++)
        {
            vec_old[i] = data_old[0].fDeformedDirections(i,ind_old);
            vec_new[i] = data_new[0].fDeformedDirections(i,ind_new);
        }
        for(int i=0; i<TSHAPE::Dimension; i++)
        {
            if(compareabs==false && !IsZero(vec_old[i]-vec_new[i]))
            {
                err = 1;
            }
            if(compareabs==true && !IsZero(abs(vec_old[i])-abs(vec_new[i]))) DebugStop();
        }
        if(err) DebugStop();
    }
//    ek_old -= ek_new;
//    REAL err_ek = Norm(ek_old);
//    std::cout << "Error ek " << err_ek << std::endl;
    
    timenew += tsim2-tsim1;
    timeold += tsim3-tsim2;
//    std::cout << "Sum wall time = " << tsim01 - tsim0 << " s" << std::endl;
    std::cout << "Total wall time of ComputeRequiredData old = " << timeold << " s" << std::endl;
    std::cout << "Total wall time of ComputeRequiredData new = " << timenew << " s" << std::endl;

}

TPZAutoPointer<TPZCompMesh> GenerateCompMesh(MMeshType eltype, int order)
{
    TPZGeoMesh *gmesh = 0;
    bool createBoundEls = false;
    TPZManVector<int> nDivs(3,1), matids(1,1);
    REAL meshsize = 1.;
    gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(meshsize, matids, nDivs, eltype, createBoundEls);
    if(!gmesh) DebugStop();
    TPZAutoPointer<TPZGeoMesh> gmeshauto(gmesh);
    TPZCompMesh *cmesh = new TPZCompMesh(gmeshauto);
    cmesh->SetDefaultOrder(order);
    int matid = matids[0];
    TPZMixedDarcyFlow *mixed = new TPZMixedDarcyFlow(matid,gmesh->Dimension());
    mixed->SetPermeability(1.);
    cmesh->InsertMaterialObject(mixed);
    int bcid = -1;
    int bctype = 0;
    TPZFNMatrix<1> val1(1,1,0.), val2(1,1,1.);
    TPZBndCond *bc = mixed->CreateBC(mixed, bcid, bctype, val1, val2);
    cmesh->InsertMaterialObject(bc);
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild();
    return cmesh;
}

template<class TSHAPE>
TPZCompElHDiv<TSHAPE> *FindElement(TPZAutoPointer<TPZCompMesh> &cmesh)
{
    int64_t nelem = cmesh->NElements();
    for (int64_t el=0; el<nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZCompElHDiv<TSHAPE> *intel = dynamic_cast<TPZCompElHDiv<TSHAPE> *>(cel);
        if(intel) return intel;
    }
    return NULL;
}
