//
//  TPZSBFemElementGroupPostProcess.cpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 19/09/25.
//

#include "TPZSBFemElementGroupPostProcess.h"
#include "pzmultiphysicselement.h"
#include "TPZH1ErrorHybridH1EstimateMaterial.h"
#include "TPZElast2DErrorEstimateMaterial.h"

TPZSBFemElementGroupPostProcess::~TPZSBFemElementGroupPostProcess() {
    
}

void TPZSBFemElementGroupPostProcess::CalcStiff(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef) {
    // reorganize the connectindices to correspond to the sequence of the refered group
    // put the rigid body mode at the end
    ReorganizeConnectOrder();
    TPZCompEl::InitializeElementMatrix(ek,ef);
    // copy the stiffness from the referred element group
    fReferred->ContributeStiffness(ek.fMat);
    // if there is no boundary, add a line representing the integral of the shape functions
    int64_t nr = ek.fMat.Rows();
    
    int nstate = fReferred->NState();
    int dim = Mesh()->Dimension();
    int nc = NConnects();
    int lastconnectsize = Connect(nc-1).NDof();
    int numRBM = 1;
    if(fHasBoundary) {
        if(lastconnectsize != 0) DebugStop();
//        ek.fMat(nr-1,nr-1) = 1.;
        numRBM = 0;
    } else {
        if(nstate == 1 && lastconnectsize != 1) DebugStop();
        else if(nstate == 2 && dim == 2 && lastconnectsize != 3) DebugStop();
        if(nstate == 1) numRBM = 1;
        if(nstate == 2) numRBM = 3;
    }
    // compute the rhs by taking the contributions of the SBFemVolume elements
    // loop over the SBFemVolume elements
        // integrate the rhs corresponding to SBFem eigenvectors
        // integrate the rhs corresponding to the bubble functions
    // compute the eigenvalue rhs by multiplying by phiinv
    // compute the bubble rhs by multiplying by AMat
    int64_t numeig = fReferred->NumEigenValues();
    TPZFMatrix<CSTATE> rhssbfem(numeig,1,0.), rhsbubble(fReferred->NumEigenValuesBubble(),1,0.),
        RBM(numeig,numRBM,0.);
    ComputeRhsRBM(rhssbfem,rhsbubble, RBM);
    {
        auto eigval = fReferred->EigenValues();
        TPZFMatrix<CSTATE> temp;
        int transpose = 1;
        fReferred->PhiInverse().Multiply(rhssbfem, temp,transpose);
        for(int i=0; i<numeig; i++) ef.fMat(i,0) = temp(i,0).real();

//        if(numRBM) {
//            RBM.Print("RBM = ", std::cout , EMathematicaInput);
//        }
        fReferred->PhiInverse().Multiply(RBM, temp,transpose);
//        if(numRBM) {
//            temp.Print("PhiRBM = ", std::cout , EMathematicaInput);
//        }
        for(int i=0; i<numeig; i++) for(int c=0; c<numRBM; c++) {
            ek.fMat(i,numeig+c) = temp(i,c).real();
            ek.fMat(numeig+c,i) = temp(i,c).real();
        }
        fReferred->MatBubble().Multiply(rhsbubble, fRhsBubble,transpose);
        if(0) {
            std::ofstream out("HatDiag.txt");
            fReferred->Phi().Print("eigvec = ",out, EMathematicaInput);
            fReferred->PhiInverse().Print("eivecPZInv = ",out,EMathematicaInput);
            out << "Eigenvalues " << eigval << std::endl;
            rhssbfem.Print("rhssbfem ",out);
            rhsbubble.Print("rhsbubble ",out);
            fRhsBubble.Print("fRhsBubble", out);
            out << "eigvalBubble = "<< fReferred->EigenValuesBubble() << std::endl;
            fReferred->PhiBubble().Print("phibubble = ",out,EMathematicaInput);
            fReferred->MatBubble().Print("matbubble = ",out,EMathematicaInput);
        }
    }
//    fRhsBubble.Zero();
    
}

/// reorganize the connect indexes to correspond to the original sbfem group
void TPZSBFemElementGroupPostProcess::ReorganizeConnectOrder() {
    int GroupNCAt = fReferred->NConnects();
    int GroupNCMF = NConnects();
    std::map<int64_t,int64_t> OldtoNew;
    const TPZVec<TPZCompEl *> &elvecAtomic = fReferred->GetElGroup();
    const TPZVec<TPZCompEl *> &elvecMF = this->GetElGroup();
    if(elvecAtomic.size() != elvecMF.size()) DebugStop();
    std::set<int64_t> constant;
    for(int64_t ivol=0; ivol < elvecAtomic.size(); ivol++) {
        TPZCompEl *celAt = elvecAtomic[ivol];
        TPZCompEl *celMF = elvecMF[ivol];
        int ncAt = celAt->NConnects();
        int ncMF = celMF->NConnects();
        for(int ic = 0; ic<ncAt; ic++) {
            int64_t cindexAt = celAt->ConnectIndex(ic);
            int64_t cindexMF = celMF->ConnectIndex(ic);
            if(OldtoNew.find(cindexAt) == OldtoNew.end()) {
                OldtoNew[cindexAt] = cindexMF;
            } else {
                if(OldtoNew[cindexAt] != cindexMF) DebugStop();
            }
        }
        for(int ic = ncAt; ic < ncMF; ic++) {
            int64_t cindexMF = celMF->ConnectIndex(ic);
            constant.insert(cindexMF);
        }
    }
    TPZManVector<int64_t> connectindexes(GroupNCMF);
    for(int64_t i = 0; i<GroupNCAt; i++) {
        int64_t cindexAt = fReferred->ConnectIndex(i);
        int64_t cindexMF = OldtoNew[cindexAt];
        connectindexes[i] = cindexMF;
    }
    if(constant.size() != GroupNCMF-GroupNCAt) DebugStop();
    auto it = constant.begin();
    for(int64_t i = GroupNCAt; i< GroupNCMF; i++) {
        connectindexes[i] = *it;
        it++;
    }
    ReorderConnects(connectindexes);
    
}

#include "pztrnsform.h"

/// cpmpute the right hand side contribution of the hybrid h1 reconstruction for sbfem volume elements
void TPZSBFemElementGroupPostProcess::ComputeRhsRBM(TPZFMatrix<CSTATE> &rhssbfem, TPZFMatrix<CSTATE> &rhsbubble, TPZFMatrix<CSTATE> &RBM) {
    TPZElementMatrixT<CSTATE> ef(Mesh(),TPZElementMatrix::EF);
    int dim = Mesh()->Dimension();
    int nstate = fReferred->NState();
    InitializeElementMatrix(ef);
    int grouporder = 0;
    int nc = NConnects();
    for(int ic = 0; ic<nc; ic++) {
        int corder = Connect(ic).Order();
        grouporder = grouporder < corder ? corder : grouporder;
    }
    int64_t nelgrp = fElGroup.size();
    for(int el = 0; el<nelgrp; el++) {
        TPZCompEl *cel = fElGroup[el];
        TPZGeoEl *gel = cel->Reference();
        TPZMultiphysicsElement *mcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if(!mcel) DebugStop();
        int nelmp = mcel->ElementVec().size();
        TPZMaterial *mat = mcel->Material();
        TPZH1ErrorHybridH1EstimateMaterial *errmat1 = dynamic_cast<TPZH1ErrorHybridH1EstimateMaterial *>(mat);
        TPZElast2DErrorEstimateMaterial *errmat2 = dynamic_cast<TPZElast2DErrorEstimateMaterial *>(mat);
        if(!errmat1 && !errmat2) DebugStop();
        int nstate = mat->NStateVariables();
        TPZManVector<TPZMaterialDataT<STATE>> datavec(nelmp);
        mcel->InitMaterialData(datavec);
        extern std::complex<STATE> integrateF;
        std::complex<STATE> storeF = integrateF;
        integrateF = 0.;
        
        int sbfemindex = TPZH1ErrorHybridH1EstimateMaterial::Epressure;
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(mcel->Element(sbfemindex));
        if(!sbfem) DebugStop();
        
        TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(gel->NSides() - 1, 1);
        int maxIntOrder = intrule->GetMaxOrder();
//        TPZManVector<int, 3> maxorder(Dimension(), 2*grouporder+8);
        TPZManVector<int, 3> maxorder(Dimension(), maxIntOrder);
        intrule->SetOrder(maxorder);

        int64_t nfunc = sbfem->Phi().Cols();
        int64_t nfuncbubble = sbfem->PhiBubble().Cols();
        TPZFMatrix<CSTATE> phieig(nfunc,nstate),dphixeig(dim*nstate,nfunc);
        TPZFMatrix<CSTATE> phibubble(nfuncbubble,nstate),dphixbubble(dim*nstate,nfuncbubble);
        TPZManVector<TPZTransform<> > trvec(nelmp);
        for(auto &it : trvec) it = TPZTransform<>(dim);

        TPZManVector<REAL,3> qsi(dim,0.1);
        REAL weight = 1;
        int64_t npts = intrule->NPoints();
        for(int ip = 0; ip<npts; ip++) {
//            extern bool Print;
//            if(ip == 0) Print = true;
//            else Print = false;
            intrule->Point(ip, qsi, weight);
            sbfem->ComputeEigenShape(qsi, phieig, dphixeig);
            mcel->ComputeRequiredData(qsi, trvec, datavec);
            REAL detjac = datavec[2].detjac;
            weight *= detjac;
//            std::cout << " Norm dphixeig " << Norm(dphixeig) << std::endl;
            if(errmat1) {
                errmat1->Contribute(datavec, phieig, dphixeig, weight, rhssbfem, RBM);
            } else if(errmat2) {
                errmat2->Contribute(datavec, phieig, dphixeig, weight, rhssbfem, RBM);
            }
//            Print = false;
            TPZFMatrix<CSTATE> RBMFake;
            sbfem->ComputeBubbleShape(qsi, phibubble, dphixbubble);
            if(errmat1) {
                errmat1->Contribute(datavec, phibubble, dphixbubble, weight, rhsbubble, RBMFake);
            } else {
                errmat2->Contribute(datavec, phibubble, dphixbubble, weight, rhsbubble, RBMFake);
            }
//            std::cout << "norm rhsbubble " << Norm(rhsbubble) << std::endl;
        }
        integrateF /= 2.;
//        std::cout << "For sbfem " << sbfem->Index() << " integrate = " << integrateF << std::endl;
        integrateF += storeF;
        
    }
//    extern std::complex<STATE> integrateF;
//    std::cout << "Completed integral " << integrateF << std::endl;
}

void TPZSBFemElementGroupPostProcess::LoadSolution() {
    // assume the solution has been transferred to the atomic meshes
    // we need to gather the solution of the connects
    // set the coefficients of the referred element
    // compute the value of the bubble functions using the bubble rhs
    int nc = fReferred->NConnects();
    int64_t ncoef = fReferred->NumEigenValues();
    
    TPZFNMatrix<100, std::complex<double> > uh_local(ncoef, fMesh->Solution().Cols(),0.);
    
    TPZFMatrix<STATE> &meshSol = fMesh->Solution();
    
    int count = 0;
    for (int ic=0; ic<nc; ic++) {
        TPZConnect &c = Connect(ic);
        int nshape = c.NShape();
        int nstate = c.NState();
        int blsize = nshape*nstate;
        int64_t seqnum = c.SequenceNumber();
        int64_t pos = fMesh->Block().Position(seqnum);
        for (int seq=0; seq < blsize; seq++) {
            for (int j = 0; j < uh_local.Cols(); j++)
            {
                uh_local(count+seq,j) = meshSol(pos+seq,j);
            }
        }
        count += blsize;
    }
    // this will compute the coefficients that allow us to post process the error
    // and copy the internal coefficients to the connects
    fReferred->LoadSolution(uh_local,fRhsBubble);
}
