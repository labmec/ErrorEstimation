/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZDarcyMHMHDivErrorEstimator.cpp
 * Author: quinelato
 * 
 * Created on September 27, 2023, 11:47 AM
 */

#include "TPZDarcyMHMHDivErrorEstimator.h"

TPZDarcyMHMHDivErrorEstimator::TPZDarcyMHMHDivErrorEstimator(TPZMultiphysicsCompMesh &originalMesh, TPZMHMixedMeshControl *mhm, bool postProcWithHDiv): TPZMHMHDivErrorEstimator<TPZMixedDarcyFlow>(originalMesh, mhm, postProcWithHDiv) {
}

TPZDarcyMHMHDivErrorEstimator::~TPZDarcyMHMHDivErrorEstimator() {
}

/// compute the pressure weights and material weights
// fills in the data structure fPressureweights and fMatid_weights
void TPZDarcyMHMHDivErrorEstimator::ComputePrimalWeights() {
    std::cout << "Computing pressure weights\n";
    TPZCompMesh *primalMesh = fPostProcMesh.MeshVector()[1];
    const int dim = primalMesh->Dimension();
    const int64_t nel = primalMesh->NElements();
    fPrimalWeights.Resize(nel, 0);
    fMatid_weights.clear();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = primalMesh->Element(el);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        TPZMaterial *mat = this->fOriginal->FindMaterial(matid);
        if (matid == fPrimalSkeletonMatId || matid == fHybridizer.fLagrangeInterface) {
            fPrimalWeights[el] = 0.;
            fMatid_weights[matid] = 0.;
            continue;
        }
        if (!mat) DebugStop();

        TPZBndCondT<STATE> *bcmat = dynamic_cast<TPZBndCondT<STATE> *>(mat);
        if (bcmat) {
            switch(bcmat->Type()) {
                case 0: // Dirichlet BC
                    fPrimalWeights[el] = 1.e12;
                    fMatid_weights[matid] = 1.e12;
                    break;
                case 1: // Neumann BC
                    fPrimalWeights[el] = 0;
                    fMatid_weights[matid] = 0;
                    break;
                case 4: // Robin BC, weight = Km
                    fPrimalWeights[el] = bcmat->Val1()(0, 0);
                    fMatid_weights[matid] = bcmat->Val1()(0, 0);
                    break;
                default:
                    DebugStop();
            }
        } else {
            TPZMixedDarcyFlow *mixpoisson = dynamic_cast<TPZMixedDarcyFlow *>(mat);
            if (!mixpoisson) DebugStop();

            REAL perm;
            TPZVec<REAL> xi(gel->Dimension(), 0.);
            gel->CenterPoint(gel->NSides() - 1, xi);
            TPZVec<REAL> x(3, 0.);
            gel->X(xi, x);
            perm = mixpoisson->GetPermeability(x);

            if (IsZero(perm)) DebugStop();
            this->fPrimalWeights[el] = perm;
            fMatid_weights[matid] = perm;
        }
    }
    std::cout << "Finished computing pressure weights\n";
}