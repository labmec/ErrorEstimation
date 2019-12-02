//
//  TPZMFSolutionTransfer.cpp
//  LinearTracer
//
//  Created by Jose on 11/21/19.
//

#include "TPZMFSolutionTransfer.h"


void TPZMFSolutionTransfer::Match::TransferFromMultiphysics(TPZCompMesh * cmesh){
    
    TPZBlock<STATE> &blockMF = cmesh->Block();
    TPZBlock<STATE>* blockAto = fblockTarget.first;
    int seqMF = fblockTarget.second;
    int seqAto = fblocknumber;
    int blocksizeAto = blockAto->Size(seqAto);
    int blocksizeMF = blockMF.Size(seqMF);
    if(blocksizeMF != blocksizeAto){
        DebugStop();
    }
#ifdef PZDEBUG
    int64_t nblockMF = blockMF.NBlocks();
    int64_t nblockAto = blockAto->NBlocks();
    if(seqMF < nblockMF-1)
    {
        int64_t PosMF = blockMF.Position(seqMF);
        int64_t nextPosMF = blockMF.Position(seqMF+1);
        if(nextPosMF-PosMF != blocksizeMF) DebugStop();
    }
    if(seqAto < nblockAto-1)
    {
        int64_t PosAto = blockAto->Position(seqAto);
        int64_t nextPosAto = blockAto->Position(seqAto+1);
        if(nextPosAto-PosAto != blocksizeAto) DebugStop();
    }
#endif
    for (int idf=0; idf<blocksizeAto; idf++) {
        (*blockAto)(seqAto, 0, idf, 0) =  blockMF(seqMF, 0, idf, 0);
    }
//    Print(std::cout, cmesh);

}
void TPZMFSolutionTransfer::Match::TransferToMultiphysics(TPZCompMesh * cmesh){
    cmesh->InitializeBlock();
    TPZBlock<STATE> &blockToTransfer = cmesh->Block();
    TPZBlock<STATE>* blockAto = fblockTarget.first;
    int seqtoTrans = fblockTarget.second;
    int seqAto = fblocknumber;
    
    int blocksizeAto = blockAto->Size(seqAto);
    int blocksizeTarg = blockToTransfer.Size(seqtoTrans);
    if(blocksizeAto!=blocksizeTarg){
        DebugStop();
    }
    for (int idf=0; idf<blocksizeAto; idf++) {
        blockToTransfer.Put(seqtoTrans, idf, 0, blockAto->Get(seqAto, idf, 0));
    }
}

void TPZMFSolutionTransfer::MeshTransferData::TransferToMultiphysics(){
   
    int nmatch = fconnecttransfer.size();
    for (int imatch=0; imatch<nmatch; imatch++) {
        fconnecttransfer[imatch].TransferToMultiphysics(fcmesh_ori);
    }
}
void TPZMFSolutionTransfer::MeshTransferData::TransferFromMultiphysics(){
    int nmatch = fconnecttransfer.size();
    for (int imatch=0; imatch<nmatch; imatch++) {
        fconnecttransfer[imatch].TransferFromMultiphysics(fcmesh_ori);
    }
}

void TPZMFSolutionTransfer::MeshTransferData::BuildTransferData(TPZCompMesh* cmesh){
    TPZMultiphysicsCompMesh *multcmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);
    
    TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cmesh);
    TPZCompMesh * targetmesh = 0;
    
    if(subcmesh){
        targetmesh =subcmesh;
    }
    
    if(multcmesh){
        targetmesh =multcmesh;
    }

    int nels = targetmesh->NElements();
    for (int iel =0; iel <nels; iel++){
        TPZCompEl *celtarget = cmesh->Element(iel);
        TPZMultiphysicsElement *mul = dynamic_cast<TPZMultiphysicsElement *>(celtarget);
        TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>(celtarget);
        TPZElementGroup *elgroup;
        
        if(cond){
            mul =dynamic_cast<TPZMultiphysicsElement *>(cond->ReferenceCompEl());
            
            elgroup = dynamic_cast<TPZElementGroup *>(cond->ReferenceCompEl());
            
        }
        
        bool condition1 = !mul && !cond;
        bool condition2 = cond && (!mul && !elgroup);
        
        if((condition1||condition2) )  {
            continue;
        }
        int nsubel =0;
        TPZVec<TPZCompEl*> subels;
        
        if(((elgroup) && (!mul)) || mul){
            
            if((elgroup) && (!mul)){
                subels = elgroup->GetElGroup();
            }
            else{
                subels.resize(1);
                subels[0] = mul;
            }
            nsubel = subels.size();
            for(int isub =0; isub<nsubel; isub++){
                TPZCompEl *subcel = subels[isub];
                mul =dynamic_cast<TPZMultiphysicsElement *>(subcel);
                TPZVec<int> act_spacesEl = mul->GetActiveApproxSpaces();
                int initialtest =0;
                for (int iespacetes=0; iespacetes<act_spacesEl.size(); iespacetes++) {
                    if (act_spacesEl[iespacetes]==0) {
//                        initialtest += mul->Element(iespacetes)->NConnects();
                        continue;
                    }
                    TPZCompEl *celfrom = mul->Element(iespacetes);
                    if (!celfrom) {
                        continue;
                    }
                    
                    int nconnects = celfrom->NConnects();
                    for (int icon = 0 ; icon < nconnects; icon++){
                        TPZConnect &conectFrom = celfrom->Connect(icon);
                        int seqnumberAto = conectFrom.SequenceNumber();
                        TPZConnect &conectTarget = mul->Connect(initialtest + icon);
                        int seqnumberMF = conectTarget.SequenceNumber();
                        Match currentmatch;
                        currentmatch.fblocknumber = seqnumberAto;
                        TPZBlock<STATE> *blockAto = &celfrom->Mesh()->Block();
                        std::pair<TPZBlock<STATE> *, int> target = std::make_pair(blockAto,seqnumberMF);
                        currentmatch.fblockTarget = target;
                        fconnecttransfer.push_back(currentmatch);
                    }
                    initialtest += nconnects;
                }
            }
        }
    }
}


void TPZMFSolutionTransfer::TransferFromMultiphysics(){
    
    int nsolutionstransfers = fmeshTransfers.size();
    for (int isoltrans=0; isoltrans<nsolutionstransfers ; isoltrans++) {
        fmeshTransfers[isoltrans].TransferFromMultiphysics();
    }
    
}
void TPZMFSolutionTransfer::TransferToMultiphysics(){

    int nsolutionstransfers = fmeshTransfers.size();
    for (int isoltrans=0; isoltrans<nsolutionstransfers ; isoltrans++) {
        fmeshTransfers[isoltrans].TransferToMultiphysics();
    }
}

void TPZMFSolutionTransfer::BuildTransferData(TPZCompMesh* cmesh){
    
    MeshTransferData transdata;
    transdata.fcmesh_ori =cmesh;
    transdata.BuildTransferData(cmesh);
    fmeshTransfers.push_back(transdata);
    int nels = cmesh->NElements();
    for (int iel=0; iel<nels; iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        TPZSubCompMesh * subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!subcmesh) {
            continue;
        }
        MeshTransferData transdatasub;
        transdatasub.fcmesh_ori =subcmesh;
        transdatasub.BuildTransferData(subcmesh);
        fmeshTransfers.push_back(transdatasub);
    }
}

void TPZMFSolutionTransfer::MeshTransferData::Print(std::ostream &out)
{
    out << "MeshTransferData " << (void *) fcmesh_ori;
    int64_t trsize = fconnecttransfer.size();
    out << " number of connect transfer blocks " << trsize << std::endl;
    for (int64_t i=0; i<trsize; i++) {
        fconnecttransfer[i].Print(out, fcmesh_ori);
    }
}
void TPZMFSolutionTransfer::Print(std::ostream &out)
{
    out << __PRETTY_FUNCTION__;
    int64_t nmesh = fmeshTransfers.size();
    out << " Number of meshes " << nmesh << std::endl;
    for (int64_t i=0; i<nmesh; i++) {
        fmeshTransfers[i].Print(out);
    }
}

void TPZMFSolutionTransfer::Match::Print(std::ostream &out, TPZCompMesh *cmesh)
{
    out << "MeshA " << (void *) &cmesh->Block() << " block " << fblockTarget.second
    << " MeshB " << fblockTarget.first << " block " << fblocknumber << std::endl;
    int blnumA = fblockTarget.second;
    int blnumB = fblocknumber;
    TPZBlock<STATE> *blockA = &cmesh->Block();
    TPZBlock<STATE> *blockB = fblockTarget.first;
    int64_t nblA = cmesh->Block().NBlocks();
    int64_t nblB = fblockTarget.first->NBlocks();
    out << "BlockA pos " << blockA->Position(fblockTarget.second)
    << " size " << blockA->Size(fblockTarget.second);
    if(blnumA < nblA-1) out << " next pos " << blockA->Position(fblockTarget.second+1);
    out << std::endl;
    out << "BlockB pos " << blockB->Position(fblocknumber)
    << " size " << blockB->Size(fblocknumber);
    if(blnumB < nblB-1) out << " next pos " << blockB->Position(fblocknumber+1);
    out << std::endl;
    int blsizeA = cmesh->Block().Size(fblockTarget.second);
    int blsizeB = fblockTarget.first->Size(fblocknumber);
    for (int i=0; i<blsizeA; i++) {
        out << cmesh->Block()(fblockTarget.second,0,i,0) << " ";
    }
    out << std::endl;
    for (int i=0; i<blsizeB; i++) {
        out << (*(fblockTarget.first))(fblocknumber,0,i,0) << " ";
    }
    out << std::endl;

}
