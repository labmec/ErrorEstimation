/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZElasticityMHMHDivErrorEstimator.cc
 * Author: quinelato
 * 
 * Created on September 27, 2023, 12:00 PM
 */

#include "TPZElasticityMHMHDivErrorEstimator.h"
#include "Common/pzerror.h"

TPZElasticityMHMHDivErrorEstimator::TPZElasticityMHMHDivErrorEstimator(TPZMultiphysicsCompMesh &originalMesh, TPZMHMixedMeshControl *mhm, bool postProcWithHDiv) : TPZMHMHDivErrorEstimator<TPZMixedElasticityND>(originalMesh, mhm, postProcWithHDiv) {
}

TPZElasticityMHMHDivErrorEstimator::~TPZElasticityMHMHDivErrorEstimator() {
}

void TPZElasticityMHMHDivErrorEstimator::ComputePrimalWeights() {
    std::cout << "Computing displacement weights\n";
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
            TPZMixedElasticityND *mixedElasticityMaterial = dynamic_cast<TPZMixedElasticityND *>(mat);
            if (!mixedElasticityMaterial) DebugStop();

            STATE weight;
            TPZVec<REAL> xi(gel->Dimension(), 0.);
            gel->CenterPoint(gel->NSides() - 1, xi);
            TPZVec<REAL> x(3, 0.);
            gel->X(xi, x);
            // DebugStop();
            weight = std::norm(mixedElasticityMaterial->GetMaxComplianceEigenvalue(x));

            if (IsZero(weight)) DebugStop();
            this->fPrimalWeights[el] = weight;
            fMatid_weights[matid] = weight;
        }
    }
    std::cout << "Finished computing displacement weights\n";
}


void TPZElasticityMHMHDivErrorEstimator::PostProcessing(TPZAnalysis &an, std::string &out) {

    TPZMaterial *mat = fPostProcMesh.FindMaterial(1);
    int varindex = -1;
    if (mat) varindex = mat->VariableIndex("DisplacementFem");
    if (varindex != -1) {
        TPZStack<std::string> scalnames, vecnames;
        if (fExact) {
            vecnames.Push("DisplacementExact");
            scalnames.Push("DisplacementErrorExact");
            scalnames.Push("EnergyErrorExact");
            scalnames.Push("DisplacementEffectivityIndex");
            scalnames.Push("EnergyEffectivityIndex");
            vecnames.Push("StressExact");
        }
        vecnames.Push("DisplacementFem");
        vecnames.Push("DisplacementReconstructed");
        scalnames.Push("DisplacementErrorEstimate");
        scalnames.Push("EnergyErrorEstimate");
        vecnames.Push("StressFem");
        scalnames.Push("POrder");
        
        int dim = fPostProcMesh.Reference()->Dimension();

        an.DefineGraphMesh(dim, scalnames, vecnames, out);
        an.PostProcess(0, dim);
    }
    else {
        std::cout << __PRETTY_FUNCTION__ << "\nPost Processing variable not found!\n";
    }
}
