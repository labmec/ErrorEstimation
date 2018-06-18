//
//  TPZPostProcessError.cpp
//  PZ
//
//  Created by Philippe Devloo on 6/30/16.
//
//

#include "TPZPostProcessError.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "tpzcompmeshreferred.h"


#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZFrontStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZMaterial.h"
#include "pzbndcond.h"
#include "TPZVecL2.h"
#include "pzl2projection.h"
#include "TPZMixedErrorEstimate.h"
#include "mixedpoisson.h"

TPZPostProcessError::TPZPostProcessError(TPZCompMesh * origin) : fMeshVector(5,0)
{
    fMeshVector[0] = origin;
    CreateAuxiliaryMeshes();
}


TPZPostProcessError::TPZPostProcessError(TPZVec<TPZCompMesh *> &meshvec)
{
    fMeshVector = meshvec;
    TPZCompMesh *multiphysics = meshvec[1];
    this->fSolution = multiphysics->Solution();
    this->fBlock = multiphysics->Block();
    int64_t ncon = multiphysics->NConnects();
    fConnectSeqNumbers.resize(ncon);
    fConnectSizes.resize(ncon);
    for (int64_t i=0; i<ncon; i++) {
        fConnectSeqNumbers[i] = multiphysics->ConnectVec()[i].SequenceNumber();
        int64_t seqnum = fConnectSeqNumbers[i];
        TPZConnect &c = multiphysics->ConnectVec()[i];
        c.SetSequenceNumber(seqnum);
        int blsize = this->fBlock.Size(seqnum);
        int neq = c.NShape()*c.NState();
        fConnectSizes[i] = neq;
        if (neq != blsize) {
            DebugStop();
        }

    }
    BuildPatchStructures();
    
#ifdef PZDEBUG
    {
        std::ofstream out("../patchinfo.txt");
        PrintPatchInformation(out);
        std::ofstream out2("../mphysics.txt");
        multiphysics->Print(out2);
    }
#endif
    
}

TPZPatch TPZPostProcessError::BuildPatch(TPZCompElSide &seed)
{
    // seed an element/side of multiphysics mesh
    
    // connected : all elements that will compose the patch
    TPZStack<TPZCompElSide> connected;
    // build the set of elements which contain the node
    TPZCompMesh *cmesh = fMeshVector[1];
    TPZGeoMesh *gmesh = fMeshVector[1]->Reference();
    int meshdim = gmesh->Dimension();
    
    // the geometric mesh needs to point to the multiphysics mesh
    if (gmesh->Reference() != fMeshVector[1]) {
        DebugStop();
    }
    TPZGeoElSide gelside = seed.Reference();
    if (gelside.Dimension() != 0) {
        DebugStop();
    }
    // all the elements that form the patch
    std::set<TPZCompEl *> patchelelements;
    std::set<int64_t> connectset;
    std::set<int64_t> internalconnectset, boundaryconnectset;
    // all geometric elements that are neighbour of the node and have a computational element associated
    gelside.ConnectedCompElementList(connected,0,0);
    connected.Push(gelside.Reference());
    while (connected.size()) {
        TPZCompElSide complocside = connected.Pop();
        patchelelements.insert(complocside.Element());
        TPZCompEl *cel = complocside.Element();
        int nc = cel->NConnects();
        for (int ic=0; ic<nc; ic++) {
            connectset.insert(cel->ConnectIndex(ic));
        }
        TPZStack<TPZGeoElSide> highsides;
        TPZGeoElSide geolocside = complocside.Reference();
        geolocside.Element()->AllHigherDimensionSides(geolocside.Side(),meshdim-1,highsides);
        int nsides = highsides.size();
        for (int is = 0; is<nsides; is++) {
            // this is typical for HDiv approximation spaces
            TPZGeoEl *geolocel = geolocside.Element();
            if (geolocel->SideDimension(highsides[is].Side()) != meshdim-1) {
                DebugStop();
            }
            // the element sides of dimension meshdim-1 starting from the corner node
            TPZGeoElSide geoloclocside(geolocside.Element(),highsides[is].Side());
            // now we should include the elements connected to that side
            geoloclocside.HigherLevelCompElementList2(connected,0,0);
        }
    }
    // the connects to be included are those who receive all contributions
    // build a nelconnected data structure for the elements in the patch
    std::map<int64_t,int> nelconnected;
    for (std::set<TPZCompEl *>::iterator it = patchelelements.begin(); it != patchelelements.end(); it++)
    {
        TPZStack<int64_t> nodelist;
        (*it)->BuildConnectList(nodelist);
        int nc = nodelist.size();
        for (int ic=0; ic<nc; ic++) {
            nelconnected[nodelist[ic]]++;
        }
    }
    // distinguish between internal nodes and external nodes
    for (std::map<int64_t,int>::iterator it = nelconnected.begin(); it != nelconnected.end(); it++) {
        int64_t cindex = it->first;
        int nelc = it->second;
        if (nelc == cmesh->ConnectVec()[cindex].NElConnected()) {
            internalconnectset.insert(cindex);
        }
        else
        {
            boundaryconnectset.insert(cindex);
        }
    }
    TPZPatch result;
    result.fConnectIndices.Resize(internalconnectset.size(), -1);
    result.fElIndices.Resize(patchelelements.size(), -1);
    result.fBoundaryConnectIndices.Resize(boundaryconnectset.size(), -1);
    int count = 0;
    for (std::set<int64_t>::iterator it = internalconnectset.begin(); it != internalconnectset.end(); it++) {
        result.fConnectIndices[count++] = *it;
    }
    count = 0;
    for (std::set<TPZCompEl *>::iterator it = patchelelements.begin(); it != patchelelements.end(); it++) {
        result.fElIndices[count++] = (*it)->Index();
    }
    count = 0;
    for (std::set<int64_t>::iterator it = boundaryconnectset.begin(); it != boundaryconnectset.end(); it++) {
        result.fBoundaryConnectIndices[count++] = *it;
    }
    
    return result;
}

// build vector of patches of a same color
void TPZPostProcessError::BuildPatchStructures()
{
    // vector indicating which connect indices of the H1 mesh have been processed
    TPZVec<int> connectprocessed(fMeshVector[4]->NConnects(),0);
    bool connectfailed = true;
    
    // load the references of all elements of the mixed mesh
    fMeshVector[1]->Reference()->ResetReference();
    fMeshVector[1]->LoadReferences();
    
    // connectfailed will be true if there is an element patch that could not be inserted
    while (connectfailed) {
        connectfailed = false;
        // we are creating a new color
        int64_t numvecpatch = fVecVecPatches.size();
        fVecVecPatches.resize(numvecpatch+1);
        // fillin is a vector which contains 0 if the connect has not been touched by any patch
        TPZManVector<int> fillin(fMeshVector[1]->NConnects(),0);
        // loop over all elements of the patch mesh
        for (int64_t el=0; el<fMeshVector[4]->NElements(); el++) {
            // take an element of the H1 mesh
            TPZCompEl *cel = fMeshVector[4]->Element(el);
            if(!cel) continue;
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            TPZGeoEl *gel = cel->Reference();
            int ncorner = gel->NCornerNodes();
            // patches are associated with corner nodes
            for (int i=0; i<ncorner; i++) {
                bool vertexfailed = false;
                // this doesnt make sense because the connectindex of the patch mesh has
                // no relation with the connectprocessed
                int64_t locconnectindex = intel->ConnectIndex(i);
                TPZConnect &c = intel->Connect(i);
                if (c.HasDependency()) {
                    continue;
                }
                if (connectprocessed[locconnectindex] == 1) {
                    continue;
                }
                TPZGeoElSide gelside(gel,i);
                TPZCompElSide celside(gelside.Reference());
                // celside belongs to the multiphysics mesh
                TPZPatch locpatch = BuildPatch(celside);

                for (int64_t ic = 0; ic < locpatch.fConnectIndices.size(); ic++) {
                    if (fillin[locpatch.fConnectIndices[ic]] != 0) {
                        vertexfailed = true;
                        connectfailed = true;
                        break;
                    }
                }
                for (int64_t ic = 0; ic < locpatch.fBoundaryConnectIndices.size(); ic++) {
                    if (fillin[locpatch.fBoundaryConnectIndices[ic]] != 0) {
                        vertexfailed = true;
                        connectfailed = true;
                        break;
                    }
                }
                if (vertexfailed == false) {
                    // all systems are go !!
                    locpatch.fPartitionConnectIndex = intel->ConnectIndex(i);
                    TPZManVector<REAL,3> co(3);
                    gel->Node(i).GetCoordinates(co);
                    locpatch.fCo = co;
                    fVecVecPatches[numvecpatch].Push(locpatch);
                    for (int64_t ic = 0; ic < locpatch.fConnectIndices.size(); ic++) {
                        fillin[locpatch.fConnectIndices[ic]] = 1;
                    }
                    for (int64_t ic = 0; ic < locpatch.fBoundaryConnectIndices.size(); ic++) {
                        fillin[locpatch.fBoundaryConnectIndices[ic]] = 1;
                    }
                    connectprocessed[locconnectindex] = 1;
                }
                if (vertexfailed == true) {
                    break;
                }
            }
        }
        PrintPartitionDiagnostics(numvecpatch, std::cout);
    }
}

// print the relevant information of the patches
void TPZPostProcessError::PrintPatchInformation(std::ostream &out)
{
    out << "Number of colors = " << fVecVecPatches.size() << std::endl;
    for (int64_t color = 0; color < fVecVecPatches.size(); color++)
    {
        int64_t numpatches = fVecVecPatches[color].size();
        out << "color number " << color << std::endl;
        for (int64_t p = 0; p < numpatches; p++) {
            out << "patch number " << p << std::endl;
            fVecVecPatches[color][p].Print(out);
        }
    }
}

// compute the estimated H1 seminorm errors
void TPZPostProcessError::ComputeHDivSolution()
{
    TPZCompMesh *meshmixed = fMeshVector[1];
    
    int64_t nval = fMeshVector[4]->Solution().Rows();
    for (int64_t i=0; i<nval; i++) {
        fMeshVector[4]->Solution()(i,0) = 1.;
    }
    
    int ModelDimension = meshmixed->Dimension();
    TPZAnalysis an(meshmixed);
#ifdef PZDEBUG
    int numthreads = 0;
#else
    int numthreads = 8;
#endif
#ifdef USING_MKL2
    TPZSymetricSpStructMatrix strmat(meshmixed);
    strmat.SetNumThreads(numthreads);
    an.SetStructuralMatrix(strmat);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(meshmixed);
    strmat.SetNumThreads(numthreads);
    strmat.SetDecomposeType(ELDLt);
    //		TPZSkylineStructMatrix strmat3(cmesh);
    //        strmat3.SetNumThreads(8);
#endif
    
    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;

    an.Run();
    
    /** Variable names for post processing */
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("POrder");
    scalnames.Push("Pressure");
    vecnames.Push("Flux");
    
    {
        std::stringstream sout;
        sout << "../" << "Poisson" << ModelDimension << "HDiv" << ".vtk";
        an.DefineGraphMesh(ModelDimension,scalnames,vecnames,sout.str());
    }

    
    an.PostProcess(1,ModelDimension);

}

// compute the estimated H1 seminorm errors
void TPZPostProcessError::ComputeElementErrors(TPZVec<STATE> &elementerrors)
{
    
    TPZCompMesh *meshmixed = fMeshVector[1];
    {
        std::ofstream out5("../CMeshH1.txt");
        fMeshVector[0]->Print(out5);
        std::ofstream out("../CMeshError.txt");
        meshmixed->Print(out);
        std::ofstream out2("../PressureCMesh.txt");
        fMeshVector[3]->Print(out2);
        std::ofstream out3("../FluxCMesh.txt");
        fMeshVector[2]->Print(out3);
        std::ofstream out4("../PartitionCMesh.txt");
        fMeshVector[4]->Print(out4);
        
    }
    TPZCompMesh *meshweight = fMeshVector[4];
    TPZVec<TPZCompEl *> elpointers(meshmixed->NElements());
    int64_t nelem = meshmixed->NElements();
    for (int64_t i=0; i< nelem; i++) {
        elpointers[i] = meshmixed->Element(i);
    }
    int64_t ncon = meshmixed->NConnects();
    int64_t nblocks = meshmixed->Block().NBlocks();
    TPZManVector<int64_t> origseqnum(ncon);
    for (int64_t ic=0; ic<ncon; ic++) {
        origseqnum[ic] = meshmixed->ConnectVec()[ic].SequenceNumber();
    }
    
    int ncolors = fVecVecPatches.size();
    for (int color = 0; color < ncolors; color++)
    {
        meshmixed->Solution().Zero();
        for (int64_t el=0; el<nelem; el++) {
            meshmixed->ElementVec()[el] = elpointers[el];
        }

        TPZVec<TPZCompEl *> activel(nelem,0);
        TPZManVector<int64_t> permute(nblocks,-1);
        
        meshweight->Solution().Zero();
        int npatch = fVecVecPatches[color].size();
        int64_t seqcount = 0;
        int64_t nintequations = 0;
        int64_t totalequations = 0;
        int64_t numintconnects = 0;
        for (int64_t patch=0; patch<npatch; patch++) {
            for(auto cindex:fVecVecPatches[color][patch].fConnectIndices)
            {
                TPZConnect &c = meshmixed->ConnectVec()[cindex];
                if(!c.HasDependency()) numintconnects++;
            }
        }
        int64_t extseqcount = numintconnects;
        for (int64_t patch=0; patch<npatch; patch++) {
            {
                int64_t partitionindex = fVecVecPatches[color][patch].fPartitionConnectIndex;
                int64_t seqnum = meshweight->ConnectVec()[partitionindex].SequenceNumber();
                meshweight->Block()(seqnum,0,0,0) = 1.;
            }
            {
                int64_t nel = fVecVecPatches[color][patch].fElIndices.size();
                for (int64_t el=0; el<nel; el++) {
                    int64_t elindex = fVecVecPatches[color][patch].fElIndices[el];
                    activel[elindex] = elpointers[elindex];
                }
            }
            {
                int64_t ncon = fVecVecPatches[color][patch].fConnectIndices.size();
                for (int64_t ic=0; ic<ncon; ic++) {
                    int64_t cindex = fVecVecPatches[color][patch].fConnectIndices[ic];
                    TPZConnect &c = meshmixed->ConnectVec()[cindex];
                    int64_t seqnum = c.SequenceNumber();
                    if (seqnum != -1)
                    {
                        int64_t consize = meshmixed->ConnectVec()[cindex].NShape() * meshmixed->ConnectVec()[cindex].NState();
                        nintequations += consize;
                        if(!c.HasDependency())
                        {
                            if(permute[seqnum] != -1) DebugStop();
                            
                            permute[seqnum] = seqcount++;
                            totalequations += consize;
                        }
                    }
                }
            }
            meshweight->LoadSolution(meshweight->Solution());
            {
                int64_t ncon = fVecVecPatches[color][patch].fBoundaryConnectIndices.size();
                for (int64_t ic=0; ic<ncon; ic++) {
                    int64_t cindex = fVecVecPatches[color][patch].fBoundaryConnectIndices[ic];
                    TPZConnect &c = meshmixed->ConnectVec()[cindex];
                    int64_t seqnum = c.SequenceNumber();
                    if (seqnum != -1)
                    {
                        if (permute[seqnum] != -1) {
                            DebugStop();
                        }
                        int64_t consize = meshmixed->ConnectVec()[cindex].NShape() * meshmixed->ConnectVec()[cindex].NState();
                        if(!c.HasDependency())
                        {
                            permute[seqnum] = extseqcount++;
                            totalequations += consize;
                        }
                    }
                }
            }

        }
#ifdef PZDEBUG
        {
            std::set<int64_t> permval;
            for (auto it:permute) {
                if (it == -1) {
                    continue;
                }
                if (permval.find(it) != permval.end()) {
                    DebugStop();
                }
                permval.insert(it);
            }
            if(permval.size() != extseqcount) DebugStop();
        }
#endif
        for (int64_t ic = 0; ic<nblocks; ic++) {
            if (permute[ic] == -1) {
                permute[ic] = extseqcount++;
            }
        }
        if(extseqcount != nblocks) DebugStop();
        meshmixed->Permute(permute);
        int64_t nactiveel = 0;
        for (int64_t el=0; el<nelem; el++) {
            meshmixed->ElementVec()[el] = activel[el];
            if (activel[el]) {
                nactiveel++;
            }
        }
        meshmixed->ComputeNodElCon();
        int64_t nequations = meshmixed->NEquations();
        if (nequations != totalequations) {
            DebugStop();
        }
        // Saddle permute would put the pressure equations after the boundary equations (and break the code)
//        meshmixed->SaddlePermute();
        meshmixed->ExpandSolution();

        std::cout << "Color " << color << " Number of active elements " << nactiveel << " Number of equations " << nequations
        << "\nNumber of internal equations " << nintequations << std::endl;
//        PrintPartitionDiagnostics(color, std::cout);
        nequations = meshmixed->NEquations();
        TPZAnalysis an(meshmixed,false);

        
        TPZSkylineStructMatrix strmat(meshmixed);
        int numthreads = 0;
        strmat.SetNumThreads(numthreads);
//        strmat.SetEquationRange(0, nequations);
        
        strmat.SetEquationRange(0, nintequations);
        an.SetStructuralMatrix(strmat);
        //		TPZSkylineStructMatrix strmat3(cmesh);
        //        strmat3.SetNumThreads(8);
        
        TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
        direct->SetDirect(ELDLt);
        an.SetSolver(*direct);
        delete direct;
        direct = 0;
        
        an.Assemble();
        
        TPZAutoPointer<TPZMatrix<STATE> > globmat = an.Solver().Matrix();
        for (int64_t p=0; p<npatch; p++) {
            TPZPatch &patch = fVecVecPatches[color][p];
            if (!PatchHasBoundary(patch))
            {
                int64_t firstlagrangeequation = patch.FirstLagrangeEquation(meshmixed);
                STATE diag = globmat->GetVal(firstlagrangeequation, firstlagrangeequation);
                diag += 1.;
                globmat->Put(firstlagrangeequation, firstlagrangeequation, diag);
            }
        }

        an.Solve();
        
        TPZStepSolver<STATE> &step = dynamic_cast<TPZStepSolver<STATE> &>(an.Solver());
        std::list<int64_t> singular = step.Singular();
        
        if (singular.size())
        {
            std::cout << "The following equations were flagged as singular\n";
            for (std::list<int64_t>::iterator it=singular.begin(); it != singular.end(); it++) {
                std::cout << *it << " ";
            }
            std::cout << std::endl;
        }
        else
        {
            std::cout << "The matrix has no singularity\n";
        }
        an.LoadSolution();
        {
            std::ofstream out("../meshmixed.txt");
            meshmixed->Print(out);
        }

        // now we have a partial solution
        /** Variable names for post processing */
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("POrder");
        scalnames.Push("Pressure");
        vecnames.Push("Flux");
        an.SetStep(color);

        {
            int ModelDimension = meshmixed->Dimension();
            std::stringstream sout;
            sout << "../" << "Poisson" << ModelDimension << "HDiv" << ".vtk";
            an.DefineGraphMesh(ModelDimension,scalnames,vecnames,sout.str());
        }
        
        an.PostProcess(1,meshmixed->Dimension());

        for (int64_t el=0; el<nelem; el++) {
            meshmixed->ElementVec()[el] = elpointers[el];
        }
        meshmixed->ComputeNodElCon();
        TransferAndSumSolution(meshmixed);
        ResetState();

        an.SetStep(color);

        // now draw the full mesh with the summed flux
        {
            int ModelDimension = meshmixed->Dimension();
            std::stringstream sout;
            sout << "../" << "FullPoisson" << ModelDimension << "HDiv" << ".vtk";
            an.DefineGraphMesh(ModelDimension,scalnames,vecnames,sout.str());
        }
        
        an.PostProcess(1,meshmixed->Dimension());
        fMeshVector[2]->Solution().Zero();
        fMeshVector[3]->Solution().Zero();
        an.Solution().Zero();

    }
    TPZAnalysis an(meshmixed,false);
    an.Solution() = fSolution;
    an.LoadSolution();
    TPZManVector<TPZCompMesh *,2> mixed(2);
    mixed[0] = fMeshVector[2];
    mixed[1] = fMeshVector[3];
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mixed, meshmixed);

    TPZManVector<REAL,3> errors(3,0.);
    {
        int64_t nels = meshmixed->ElementVec().NElements();
        meshmixed->ElementSolution().Redim(nels, 3);
    }
    an.PostProcessError(errors);
    std::cout << "Estimated error " << errors << std::endl;
}

// print partition diagnostics
void TPZPostProcessError::PrintPartitionDiagnostics(int color, std::ostream &out) const
{
    TPZCompMesh *meshmixed = fMeshVector[1];
    if (color < 0 || color >= fVecVecPatches.size()) {
        DebugStop();
    }
    TPZVec<TPZPatch> &vecpatch = fVecVecPatches[color];
    int64_t numpatch = vecpatch.size();
    // determine if the patch is a boundary patch or not
    TPZVec<int> IsInternalPatch(numpatch,0);
    for (int64_t p = 0; p<numpatch; p++) {
        IsInternalPatch[p] = !PatchHasBoundary(vecpatch[p]);
    }
    out << "Number of patches " << numpatch << std::endl;
    for (int64_t p = 0; p<numpatch; p++) {
        TPZPatch &patch = vecpatch[p];
        out << "Diagnostics for patch number " << p << " of color " << color << std::endl;
        out << "Location of patch " << patch.fCo << std::endl;
        out << "Element indices " << patch.fElIndices << std::endl;
        if(IsInternalPatch[p])
        {
            out << "Patch is internal\n";
        }
        else
        {
            out << "Patch has boundary sides\n";
        }
        
        out << "Internal connects sequence numbers ";
        for (int64_t ci=0; ci < patch.fConnectIndices.size(); ci++) {
            TPZConnect &c = meshmixed->ConnectVec()[patch.fConnectIndices[ci]];
            if (c.HasDependency()) {
                out << "*";
            }
            if (c.LagrangeMultiplier() != 0) {
                out << "L";
            }
            out << c.SequenceNumber() << " ";
        }
        out << std::endl;
        out << "Boundary connects sequence numbers ";
        for (int64_t ci=0; ci < patch.fBoundaryConnectIndices.size(); ci++) {
            TPZConnect &c = meshmixed->ConnectVec()[patch.fBoundaryConnectIndices[ci]];
            if (c.HasDependency()) {
                out << "*ERROR*";
            }
            out << c.SequenceNumber() << " ";
        }
        out << std::endl;
    }
}

// determine if a given patch is boundary or not
bool TPZPostProcessError::PatchHasBoundary(TPZPatch &patch) const
{
    TPZCompMesh *meshmixed = fMeshVector[1];
    int meshdim = meshmixed->Dimension();
    int64_t numel = patch.fElIndices.size();
    bool HasBoundary = false;
    for (int64_t el=0; el<numel; el++) {
        int64_t elindex = patch.fElIndices[el];
        TPZCompEl *cel = meshmixed->Element(elindex);
        if (!cel) {
            DebugStop();
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            DebugStop();
        }
        if (gel->Dimension() != meshdim) {
            HasBoundary = true;
            break;
        }
    }
    return HasBoundary;

}

// return the first equation associated with a lagrange multiplier
int64_t TPZPatch::FirstLagrangeEquation(TPZCompMesh *cmesh) const
{
    int64_t nconnect = fConnectIndices.size();
    for (int64_t ic=0; ic<nconnect; ic++) {
        int64_t cindex = fConnectIndices[ic];
        TPZConnect &c = cmesh->ConnectVec()[cindex];
        if (c.SequenceNumber() == -1 || c.NDof() == 0 || c.LagrangeMultiplier() == 0) {
            continue;
        }
        int64_t seqnum = c.SequenceNumber();
        int64_t eq = cmesh->Block().Position(seqnum);
        return eq;
    }
    DebugStop();
    return -1;
}

// Sum the solution stored in fSolution of the second mesh to the fSolution vector
void TPZPostProcessError::TransferAndSumSolution(TPZCompMesh *cmesh)
{
    int64_t nconnect = cmesh->NConnects();
    for (int64_t ic=0; ic<nconnect; ic++) {
        TPZConnect &c = cmesh->ConnectVec()[ic];
        if (c.SequenceNumber() == -1) {
            continue;
        }
        if (c.NElConnected() == 0) {
            DebugStop();
        }
        int64_t seqnum = c.SequenceNumber();
        int neq = c.NDof();
        int64_t pos = cmesh->Block().Position(seqnum);
        
        int64_t targetseqnum = fConnectSeqNumbers[ic];
        int64_t targetpos = fBlock.Position(targetseqnum);
        
        
        for (int eq=0; eq<neq; eq++) {
            // tototototo
            fSolution(targetpos+eq,0) += cmesh->Solution()(pos+eq,0);
        }
    }
    cmesh->Solution().Zero();
}

// Reset the state of the HDiv mesh to its original structure
void TPZPostProcessError::ResetState()
{
    TPZCompMesh *multiphysics = fMeshVector[1];
    multiphysics->Block() = this->fBlock;
    multiphysics->Solution() = this->fSolution;
    for (int64_t i=0; i<fConnectSeqNumbers.size(); i++) {
        int64_t seqnum = fConnectSeqNumbers[i];
        TPZConnect &c = multiphysics->ConnectVec()[i];
        c.SetSequenceNumber(seqnum);
        int blsize = this->fBlock.Size(seqnum);
        int neq = c.NShape()*c.NState();
        if (neq != blsize) {
            DebugStop();
        }
    }
    TPZManVector<TPZCompMesh *> meshvec(2);
    meshvec[0] = fMeshVector[2];
    meshvec[1] = fMeshVector[3];
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, multiphysics);

}

// check whether the connectsizes have changed
void TPZPostProcessError::CheckConnectSizes()
{
    TPZCompMesh *multiphysics = fMeshVector[1];
    int64_t ncon = multiphysics->NConnects();
    for (int64_t i=0; i<ncon; i++) {
        TPZConnect &c = multiphysics->ConnectVec()[i];
        int neq = c.NShape()*c.NState();
        int blsize = fConnectSizes[i];
        if (neq != blsize) {
            DebugStop();
        }
    }

}

/// create a fluxmesh based on the original H1 mesh
// the flux mesh will be put in the second position of the mesh vector
void TPZPostProcessError::CreateFluxMesh()
{
    int matId = 1;
    TPZCompMesh *cmeshroot = fMeshVector[0];
    int dim = cmeshroot->Dimension();
    TPZMaterial *rootmat = cmeshroot->FindMaterial(matId);
    int nstate = rootmat->NStateVariables();
    /// criar materiais
    //  TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
    //    material->SetForcingFunction(force1);
    TPZGeoMesh *gmesh = cmeshroot->Reference();
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    for (auto it:cmeshroot->MaterialVec()) {
        int matid = it.first;
        TPZMaterial *mat = it.second;
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        if (!bnd) {
            int matId = mat->Id();
            int nstate = mat->NStateVariables();
            TPZVecL2 *material = new TPZVecL2(matId);
            material->SetDimension(dim);
            material->SetNStateVariables(nstate);
            cmesh->InsertMaterialObject(material);
        }
    }
    for (auto it:cmeshroot->MaterialVec()) {
        TPZMaterial *mat = it.second;
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        if(bnd)
        {
            TPZMaterial *matorig = bnd->Material();
            int matid = matorig->Id();
            TPZMaterial *matl2 = cmesh->FindMaterial(matid);
            TPZBndCond *bc = matl2->CreateBC(matl2, bnd->Id(), bnd->Type(), bnd->Val1(), bnd->Val2());
            cmesh->InsertMaterialObject(bc);
        }
    }
    
    
    cmesh->SetDimModel(dim);
    
    cmesh->SetAllCreateFunctionsHDiv();
    
    cmesh->SetDefaultOrder(cmeshroot->GetDefaultOrder());
    
    
    
    
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    
    /// adjust the order of the elements
    int64_t nel = cmeshroot->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmeshroot->Element(el);
        if(!cel) continue;
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        int nc = intel->NConnects();
        int order = intel->Connect(nc-1).Order();
        TPZGeoEl *gel = cel->Reference();
        TPZCompEl *celL2 = gel->Reference();
        TPZInterpolationSpace *intelHDiv = dynamic_cast<TPZInterpolationSpace *>(celL2);
        intelHDiv->SetPreferredOrder(order);
    }
    cmesh->ExpandSolution();

    //#ifdef LOG4CXX
    //    if(logdata->isDebugEnabled())
    //    {
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Computacional_1 Fluxo\n ";
    //        cmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //    }
    //#endif

    fMeshVector[2] = cmesh;
}

/// create the lagrange mesh corresponding to the flux mesh
void TPZPostProcessError::CreatePressureMesh()
{
    TPZCompMesh *fluxmesh = fMeshVector[2];
    TPZGeoMesh *gmesh = fluxmesh->Reference();
    int dim = fluxmesh->Dimension();
    
    TPZCompMeshReferred *cmesh = new TPZCompMeshReferred(gmesh);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuousReferred();
    for (auto it:fluxmesh->MaterialVec()) {
        TPZMaterial *mat = it.second;
        int matdim = mat->Dimension();
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        if (!bnd && matdim == dim) {
            int matId = mat->Id();
            int nstate = mat->NStateVariables();
            TPZVec<STATE> sol(nstate,0.);
            TPZL2Projection *material = new TPZL2Projection(matId,dim,nstate,sol);
            cmesh->InsertMaterialObject(material);
        }
    }

    gmesh->ResetReference();
    
    
    int64_t nel = fluxmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *flux = fluxmesh->Element(el);
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(flux);
        TPZGeoEl *gel = 0;
        if(intel) gel = intel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        int nconnects = intel->NConnects();
        TPZConnect &c = intel->Connect(nconnects-1);
        cmesh->SetDefaultOrder(c.Order());
        int64_t index;
        TPZCompEl *cel = cmesh->ApproxSpace().CreateCompEl(gel, *cmesh, index);
        cel->Reference()->ResetReference();
    }
    cmesh->LoadReferred(fMeshVector[0]);
    cmesh->InitializeBlock();
    for (auto &c:cmesh->ConnectVec()) {
        c.SetLagrangeMultiplier(1);
    }
    
    fMeshVector[3] = cmesh;
}

/// create the partition of unity mesh
void TPZPostProcessError::CreatePartitionofUnityMesh()
{
    TPZCompMeshReferred *pressuremesh = dynamic_cast<TPZCompMeshReferred *> (fMeshVector[3]);
    TPZGeoMesh *gmesh = pressuremesh->Reference();
    int dim = pressuremesh->Dimension();
    
    TPZCompMeshReferred *cmesh = new TPZCompMeshReferred(gmesh);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuousReferred();
    cmesh->SetDefaultOrder(1);
    
    for (auto it:pressuremesh->MaterialVec()) {
        TPZMaterial *mat = it.second;
        int matdim = mat->Dimension();
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        if (!bnd && matdim == dim) {
            int matId = mat->Id();
            int nstate = 1;
            TPZVec<STATE> sol(nstate,0.);
            TPZL2Projection *material = new TPZL2Projection(matId,dim,nstate,sol);
            cmesh->InsertMaterialObject(material);
        }
    }
    
    gmesh->ResetReference();
    
    
    int64_t nel = pressuremesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *pressure = pressuremesh->Element(el);
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(pressure);
        TPZGeoEl *gel = 0;
        if(intel) gel = intel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        int64_t index;
        cmesh->ApproxSpace().CreateCompEl(gel, *cmesh, index);
    }
    cmesh->LoadReferred(fMeshVector[0]);
    pressuremesh->LoadReferred(cmesh);
    cmesh->InitializeBlock();
    fMeshVector[4] = cmesh;

}


/// create the multiphysics mesh that will compute the projection matrix
void TPZPostProcessError::CreateMixedMesh()
{
    // the pressure mesh is the rootmesh
    TPZCompMesh *cmeshroot = fMeshVector[0];
    TPZGeoMesh *gmesh = fMeshVector[0]->Reference();
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim = gmesh->Dimension();
    
    for (auto it:cmeshroot->MaterialVec()) {
        TPZMaterial *mat = it.second;
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        if (!bnd) {
            int matId = mat->Id();
            int nstate = mat->NStateVariables();
            TPZMaterial *material = 0;
            if (nstate == 1) {
                TPZMixedErrorEstimate<TPZMixedPoisson> *locmat = new TPZMixedErrorEstimate<TPZMixedPoisson>(matId,dim);
                locmat->SetSignConvention(-1);
                material = locmat;
                locmat->SetForcingFunction(mat->ForcingFunction());
                //incluindo os dados do problema
                TPZFNMatrix<9,REAL> PermTensor(3,3,0.);
                TPZFNMatrix<9,REAL> InvPermTensor(3,3,0.);
                PermTensor.Identity();
                InvPermTensor.Identity();
                
                locmat->SetPermeabilityTensor(PermTensor, InvPermTensor);
            }
            mphysics->InsertMaterialObject(material);
        }
    }
    for (auto it:cmeshroot->MaterialVec()) {
        TPZMaterial *mat = it.second;
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        if(bnd)
        {
            TPZMaterial *matorig = bnd->Material();
            int matid = matorig->Id();
            TPZMaterial *matmixed = mphysics->FindMaterial(matid);
            TPZBndCond *bc = matmixed->CreateBC(matmixed, bnd->Id(), bnd->Type(), bnd->Val1(), bnd->Val2());
            mphysics->InsertMaterialObject(bc);
        }
    }

    TPZManVector<TPZCompMesh *> meshvec(2,0);
    meshvec[0] = fMeshVector[2];
    meshvec[1] = fMeshVector[3];
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    //Fazendo auto build
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    mphysics->CleanUpUnconnectedNodes();
    //------- Create and add group elements -------
    fMeshVector[1] = mphysics;
}


// create the meshes that allow us to compute the error estimate
void TPZPostProcessError::CreateAuxiliaryMeshes()
{
    CreateFluxMesh();
    CreatePressureMesh();
    CreateMixedMesh();
    CreatePartitionofUnityMesh();
    
    TPZCompMesh *multiphysics = fMeshVector[1];
    this->fSolution = multiphysics->Solution();
    this->fBlock = multiphysics->Block();
    int64_t ncon = multiphysics->NConnects();
    fConnectSeqNumbers.resize(ncon);
    fConnectSizes.resize(ncon);
    for (int64_t i=0; i<ncon; i++) {
        fConnectSeqNumbers[i] = multiphysics->ConnectVec()[i].SequenceNumber();
        int64_t seqnum = fConnectSeqNumbers[i];
        TPZConnect &c = multiphysics->ConnectVec()[i];
        c.SetSequenceNumber(seqnum);
        int blsize = this->fBlock.Size(seqnum);
        int neq = c.NShape()*c.NState();
        fConnectSizes[i] = neq;
        if (neq != blsize) {
            DebugStop();
        }
        
    }
    BuildPatchStructures();
    
#ifdef PZDEBUG
    {
        std::ofstream out("../patchinfo.txt");
        PrintPatchInformation(out);
        std::ofstream out2("../mphysics.txt");
        multiphysics->Print(out2);
    }
#endif

}

