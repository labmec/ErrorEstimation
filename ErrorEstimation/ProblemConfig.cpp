//
// Created by Gustavo on 15/07/18.
//

#include <iostream>
#include "ProblemConfig.h"
#include "pzintel.h"

void ProblemConfig::ApplyDivision()
{
    for(auto &itlist : fElIndexDivide) {
        for(auto eleindex : itlist) {
            TPZStack<TPZGeoEl *> subels;
            gmesh->Element(eleindex)->Divide(subels);
        }
        ApplyTwoonOneRestraint();
        DivideEncircledElements();
        DivideBoundaryElements();
    }
}

static void AllSons(TPZGeoElSide gelside, TPZStack<TPZGeoElSide> &sons)
{
    if(!gelside.Element()->HasSubElement()) return;
    int64_t nsons = sons.size();
    gelside.GetSubElements2(sons);
    int64_t nsons2 = nsons+gelside.NSubElements();
    for(int64_t is = nsons; is<nsons2; is++) {
        if(sons[is].HasSubElement()) AllSons(sons[is],sons);
    }
    if(0 && sons.size() == nsons2) {
        std::cout << "gelside level " << gelside.Element()->Level() << std::endl;
        int64_t nsons = sons.size();
        for(int64_t is = 0; is<nsons; is++) {
            std::cout << "is " << is << " elindex " << sons[is].Element()->Index() << " side " << sons[is].Side() << " level " << sons[is].Element()->Level() << std::endl;
        }
    }
}

static int MaxLevel(TPZGeoElSide gelside)
{
    TPZStack<TPZGeoElSide> sons;
    AllSons(gelside,sons);
    int64_t nsons = sons.size();
    int maxlevel = 0;
    for (int64_t s = 0; s<nsons; s++) {
        int sonlevel = sons[s].Element()->Level();
        maxlevel = maxlevel < sonlevel ? sonlevel : maxlevel;
    }
    return maxlevel;
}

void ProblemConfig::ApplyTwoonOneRestraint()
{
    int64_t nel = gmesh->NElements();
    int dim = gmesh->Dimension();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel || gel->Dimension() != dim || gel->HasSubElement()) continue;
        int mylev = gel->Level();
        for (int side = gel->NCornerNodes(); side < gel->NSides()-1; side++) {
            if(gel->HasSubElement()) break;
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while(neighbour != gelside) {
                int maxlev = MaxLevel(neighbour);
                if(maxlev > mylev+1) {
                    TPZStack<TPZGeoEl *> subs;
                    std::cout << "gel index " << gel->Index() << " level " << gel->Level() << " maxlev " << maxlev << " needs divide " << std::endl;
                    gel->Divide(subs);
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
        }
    }
}

void ProblemConfig::DivideEncircledElements()
{
    int64_t nel = gmesh->NElements();
    int dim = gmesh->Dimension();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel || gel->Dimension() != dim || gel->HasSubElement()) continue;
        int firstside = gel->FirstSide(dim-1);
        int nsides = gel->NSides(dim-1);
        int numneighdivided = 0;
        for (int side = firstside; side < gel->NSides()-1; side++) {
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while(neighbour != gelside)
            {
                auto neighel = neighbour.Element();
                if(neighel->Dimension() == dim && neighel->HasSubElement()) {
                    numneighdivided++;
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
        }
        if(numneighdivided >= nsides-1) {
            TPZStack<TPZGeoEl *> subels;
            std::cout << "Element " << el << " divided because surrounded by divided elements\n";
            gel->Divide(subels);
        }
    }
}

void ProblemConfig::DivideBoundaryElements()
{
    int64_t nel = gmesh->NElements();
    int dim = gmesh->Dimension();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel || gel->Dimension() == dim || gel->HasSubElement()) continue;
        TPZGeoElSide gelside(gel);
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside) {
            if(neighbour.Element()->HasSubElement()) {
                TPZStack<TPZGeoEl *> subs;
                gel->Divide(subs);
                break;
            }
            neighbour = neighbour.Neighbour();
        }
    }
}

static int CelOrder(TPZCompEl* cel){
    int ncon = cel->NConnects();
    int celorder = cel->Connect(ncon-1).Order();
    return celorder;
}

void ProblemConfig::AdjustH1PorderDistrib(){
    TPZCompMesh *cmeshH1 = gmesh->Reference();
    int64_t nels2 = cmeshH1->NElements();
    bool change = true;
    while (change) {
        change = false;
        for(int64_t el = 0; el < nels2; el++){
            TPZCompEl* cel = cmeshH1->Element(el);
            int celorder = CelOrder(cel);
            TPZGeoEl* gel = cel->Reference();
            
            for (int side = 0; side < gel->NSides()-1; side++) {
                TPZGeoElSide gelside(gel,side);
                TPZCompElSide celsidebig = gelside.LowerLevelCompElementList2(1);
                
                if(celsidebig){
                    int orderbig = CelOrder(celsidebig.Element());
                    if(orderbig > celorder + 1){
                        TPZInterpolatedElement* intel = dynamic_cast<TPZInterpolatedElement*>(cel);
                        intel->PRefine(orderbig-1);
                        celorder = orderbig - 1;
                        change = true;
                        //break;
                    }
                }
                
                TPZInterpolatedElement* intel = dynamic_cast<TPZInterpolatedElement*>(cel);
                
                TPZGeoElSide neighbour = gelside.Neighbour();
                while(neighbour != gelside){
                    TPZGeoEl* gelneigh = neighbour.Element();
                    TPZCompEl* celneigh = gelneigh->Reference();
                    
                    if(celneigh){
                        TPZInterpolatedElement* intelneigh = dynamic_cast<TPZInterpolatedElement*>(celneigh);
                        
                        int orderneigh = CelOrder(celneigh);
                        
                        if(celorder > orderneigh + 1) {
                            intelneigh->PRefine(celorder - 1);
                            change = true;
                            //break;
                        }
                        
                        if(orderneigh > celorder + 1) {
                            intel->PRefine(orderneigh - 1);
                            change = true;
                            celorder = orderneigh - 1;
                            
                        }
                    }
                    neighbour = neighbour.Neighbour();
                }
                TPZStack<TPZCompElSide> celsidestack;
                gelside.HigherLevelCompElementList2(celsidestack, 1, 1);
                for(int i = 0; i < celsidestack.size(); i++ ){
                    TPZCompElSide celside_i = celsidestack[i];
                    TPZGeoElSide gelside_i = celside_i.Reference();
                    if(gelside_i.Dimension() != cmeshH1->Dimension()-1){
                        continue;
                    }
                    int order_i = CelOrder(celside_i.Element());
                    if (order_i > celorder + 1){
                        intel->PRefine(order_i - 1);
                        change = true;
                        celorder = order_i - 1;
                    }
                    if (celorder > order_i + 1){
                        TPZInterpolatedElement* intel_i = dynamic_cast<TPZInterpolatedElement*>(celside_i.Element());
                        intel_i->PRefine(celorder - 1);
                        change = true;
                    }
                }
            }
        }
    }
}


// To make p refinement incrementing or decrement
void ProblemConfig::PorderIncrement() {
    TPZCompMesh *cmeshH1 = gmesh->Reference();
    if (!fElIndexPplus.size() || !cmeshH1) {
        return;
    }
    int64_t nels = fElIndexPplus.size();
    
    for(auto &itlist : fElIndexPplus) {
        for(auto eleindex : itlist) {
            if (eleindex.first < 0) continue;
            TPZGeoEl* gel = gmesh->ElementVec()[eleindex.first];
            TPZCompEl* cel = gel->Reference();
            
            if (!cel || cel->Dimension() != cmeshH1->Dimension()){
                continue;
            }
                        
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(cel);
            if(!intel) continue;
            //int order = cel->GetgOrder();
            intel->PRefine(eleindex.second);
            
        }
    }
    
    //Adjustments to smooth the distribution of polynomial orders
    //AdjustH1PorderDistrib();
    
    //cmeshH1->AdjustBoundaryElements();
    cmeshH1->CleanUpUnconnectedNodes();
    
}
