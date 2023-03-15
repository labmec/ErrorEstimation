//
// Created by Gustavo on 15/07/18.
//

#include <iostream>
#include "ProblemConfig.h"

void ProblemConfig::ApplyDivision()
{
    for(auto &itlist : fElIndexDivide) {
        for(auto eleindex : itlist) {
            TPZStack<TPZGeoEl *> subels;
            gmesh->Element(eleindex)->Divide(subels);
        }
        ApplyTwoonOneRestraint();
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

