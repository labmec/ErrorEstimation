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
//#include "tpzcompmeshreferred.h"


#include "TPZLinearAnalysis.h"
#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZFrontStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontSym.h"

#include "pzvisualmatrix.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZMaterial.h"
#include "TPZBndCondT.h"
#include "Projection/TPZL2Projection.h"
#include "TPZMixedErrorEstimate.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZNullMaterial.h"

#define ERRORESTIMATION_DEBUG

TPZPostProcessError::TPZPostProcessError(TPZCompMesh * origin,ProblemConfig &config) : fMeshVector(5,0)
{
    fMeshVector[Eorigin] = origin;
    fExact = &(*(config.exact));
    CreateAuxiliaryMeshes();
}


TPZPostProcessError::TPZPostProcessError(TPZVec<TPZCompMesh *> &meshvec)
{
    fMeshVector = meshvec;
    TPZCompMesh *multiphysics = meshvec[Eflux];
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
    
#ifdef ERRORESTIMATION_DEBUG
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
    TPZCompMesh *cmesh = fMeshVector[Emulti];
    TPZGeoMesh *gmesh = fMeshVector[Emulti]->Reference();
    int meshdim = gmesh->Dimension();
    
    // the geometric mesh needs to point to the multiphysics mesh
    if (gmesh->Reference() != fMeshVector[Emulti]) {
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
    
    result.fPatchIsBoundary = PatchHasBoundary(result);
    return result;
}

// build vector of patches of a same color
void TPZPostProcessError::BuildPatchStructures()
{
    // vector indicating which connect indices of the H1 mesh have been processed
    TPZVec<int64_t> connectprocessed(fMeshVector[Epatch]->NConnects(),0);
    bool connectfailed = true;
    
    // load the references of all elements of the mixed mesh
    fMeshVector[Emulti]->Reference()->ResetReference();
    fMeshVector[Emulti]->LoadReferences();
    
    // connectfailed will be true if there is an element patch that could not be inserted
    while (connectfailed) {
        connectfailed = false;
        // we are creating a new color
        int64_t numvecpatch = fVecVecPatches.size();
        fVecVecPatches.resize(numvecpatch+1);
        // fillin is a vector which contains 0 if the connect has not been touched by any patch
        TPZManVector<int> fillin(fMeshVector[Emulti]->NConnects(),0);
        // loop over all elements of the patch mesh
        for (int64_t el=0; el<fMeshVector[Epatch]->NElements(); el++) {
            // take an element of the H1 mesh
            TPZCompEl *cel = fMeshVector[Epatch]->Element(el);
            if(!cel) continue;
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            TPZGeoEl *gel = cel->Reference();
            int ncorner = gel->NCornerNodes();
            // patches are associated with corner nodes
            for (int i=0; i<ncorner; i++) {
                bool vertexfailed = false;
                // this doesnt make sense because the connectindex of the patch
                // mesh has no relation with the connectprocessed
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

void TPZPostProcessError::BuildPatchStructures2()
{
    // vector indicating which connect indices of the H1 mesh have been processed
    TPZVec<int64_t> connectprocessed(fMeshVector[Epatch]->NConnects(),0);
    
    bool connectfailed = true;
    
    // load the references of all elements of the mixed mesh
    fMeshVector[Emulti]->Reference()->ResetReference();
    fMeshVector[Emulti]->LoadReferences();
    
    // connectfailed will be true if there is an element patch that could not be inserted
    while (connectfailed) {
        connectfailed = false;

        int64_t numvecpatch = fVecVecPatches.size();
        
        // loop over all elements of the patch mesh
        for (int64_t el=0; el<fMeshVector[Epatch]->NElements(); el++) {
            // take an element of the H1 mesh
            TPZCompEl *cel = fMeshVector[Epatch]->Element(el);
            if(!cel) continue;
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            TPZGeoEl *gel = cel->Reference();
            int ncorner = gel->NCornerNodes();
            // patches are associated with corner nodes
            for (int i=0; i<ncorner; i++) {
                bool vertexfailed = false;
                // this doesnt make sense because the connectindex of the patch
                // mesh has no relation with the connectprocessed
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
                

                if (connectprocessed[locconnectindex] == 0) {
                    // all systems are go !!
                    locpatch.fPartitionConnectIndex = intel->ConnectIndex(i);
                    TPZManVector<REAL,3> co(3);
                    gel->Node(i).GetCoordinates(co);
                    locpatch.fCo = co;
                    numvecpatch = fVecVecPatches.size();
                    fVecVecPatches.resize(numvecpatch+1);
                    fVecVecPatches[numvecpatch].Push(locpatch);

                    connectprocessed[locconnectindex] = 1;
                    PrintPartitionDiagnostics(numvecpatch, std::cout);
                }
            }
        }
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
    TPZCompMesh *meshmixed = fMeshVector[Eflux];
    
    int64_t nval = fMeshVector[Eorigin]->Solution().Rows();
    TPZFMatrix<STATE> &mesh4sol = fMeshVector[Eorigin]->Solution();
    for (int64_t i=0; i<nval; i++) {
        mesh4sol(i,0) = 1.;
    }
    
    int ModelDimension = meshmixed->Dimension();
    TPZLinearAnalysis an(meshmixed);
#ifdef ERRORESTIMATION_DEBUG
    int numthreads = 0;
#else
    int numthreads = 8;
#endif
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> strmat(meshmixed);
    strmat.SetNumThreads(numthreads);
    an.SetStructuralMatrix(strmat);
#else
    TPZSSpStructMatrix<STATE> strmat(meshmixed);
    strmat.SetNumThreads(numthreads);
    an.SetStructuralMatrix(strmat);
//    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(meshmixed);
//    strmat.SetNumThreads(numthreads);
//    strmat.SetDecomposeType(ELDLt);
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

using namespace std;
// compute the estimated H1 seminorm element errors
void TPZPostProcessError::ComputeElementErrors(TPZVec<STATE> &elementerrors)
{//elementerrors isn't being filled
    
    TPZMultiphysicsCompMesh *multiphysicsmesh = dynamic_cast<TPZMultiphysicsCompMesh*>(fMeshVector[Emulti]);
    if(1){
        std::ofstream out0("../CMeshMixed.txt");
        multiphysicsmesh->Print(out0);
        std::ofstream out1("../FluxCMesh.txt");
        fMeshVector[Eflux]->Print(out1);
        std::ofstream out2("../PressureCMesh.txt");
        fMeshVector[Epressure]->Print(out2);
        std::ofstream out3("../PartitionCMesh.txt");
        fMeshVector[Epatch]->Print(out3);
        std::ofstream out4("../CMeshH1.txt");
        fMeshVector[Eorigin]->Print(out4);
        std::ofstream out5("../AveragePressureCMesh.txt");
        fMeshVector[Epressureaverage]->Print(out5);
    }
    
    int64_t averagepressureconnindex = fMeshVector[Eflux]->NConnects()+fMeshVector[Epressure]->NConnects();
    
    TPZCompMesh *meshpatch = fMeshVector[Epatch];
    
    //Vector of pointers to the computational elements of multiphysicsmesh will be used to activate elements
    TPZVec<TPZCompEl *> elpointers(multiphysicsmesh->NElements());
    int64_t nelem = multiphysicsmesh->NElements();
    for (int64_t i = 0; i < nelem; i++) {
        elpointers[i] = multiphysicsmesh->Element(i);
    }
    
    //Take the sequence number of each connect (origseqnum is not used)
    int64_t ncon = multiphysicsmesh->NConnects();
    int64_t nblocks = multiphysicsmesh->Block().NBlocks();
    TPZManVector<int64_t> origseqnum(ncon);
    for (int64_t ic = 0; ic < ncon; ic++) {
        origseqnum[ic] = multiphysicsmesh->ConnectVec()[ic].SequenceNumber();
    }
    
    //For each color compute the contribution to equilibrate flux reconstruction
    int64_t ncolors = fVecVecPatches.size();
    for (int color = 0; color < ncolors; color++){
        multiphysicsmesh->Solution().Zero();
        for (int64_t el = 0; el < nelem; el++) {
            multiphysicsmesh->ElementVec()[el] = elpointers[el];
        }
        
        TPZVec<TPZCompEl *> activel(nelem,0); //to activate/deactivate elements
        TPZManVector<int64_t> permute(nblocks,-1);
        
        meshpatch->Solution().Zero();
        int npatch = fVecVecPatches[color].size();// Number of patches associated to the color
        int64_t seqcount = 0;
        int64_t nintequations = 0; // Number of equations associated to a color
        int64_t totalequations = 0;
        int64_t numintconnects = 0; // Number of connects associated to a color
        
        for (int64_t patch = 0; patch < npatch; patch++){
            //volumetric elements of the patch should share a single average pressure connect
            int ncontocondensed = 0;
            int64_t nel = fVecVecPatches[color][patch].fElIndices.size();
            for (int64_t el = 0; el < nel; el++) {
                int64_t elindex = fVecVecPatches[color][patch].fElIndices[el];
                TPZCompEl* cel = multiphysicsmesh->Element(elindex);
                if(cel->Material()->Id()==1){
                    int64_t ncon = cel->NConnects();
                    cel->SetConnectIndex(ncon-1, averagepressureconnindex+patch);
                    ncontocondensed++;
                    std::cout << "pressure average connect index: " << cel->ConnectIndex(ncon-1) << std::endl;
                }
            }
            TPZPatch& Patch = fVecVecPatches[color][patch];
            std::cout << "fpatchIsBoundary: " << Patch.fPatchIsBoundary << std::endl;
            if(Patch.fPatchIsBoundary){
                std::set<int64_t> averagepressure;
                int64_t nel = fVecVecPatches[color][patch].fElIndices.size();
                for (int64_t el = 0; el < nel; el++) {
                    int64_t elindex = fVecVecPatches[color][patch].fElIndices[el];
                    TPZCompEl* cel = multiphysicsmesh->Element(elindex);
                    if(cel->Material()->Id()==1){
                        
                        int64_t ncon = cel->NConnects();
                        averagepressure.insert(cel->ConnectIndex(ncon-1));
                        TPZConnect &con=cel->Connect(ncon-1);
                        con.SetCondensed(true);
                        
                    }
                }
                if(averagepressure.size()!=1) DebugStop();
            }
            else{
                std::set<int64_t> averagepressure;
                int64_t nel = fVecVecPatches[color][patch].fElIndices.size();
                for (int64_t el = 0; el < nel; el++) {
                    int64_t elindex = fVecVecPatches[color][patch].fElIndices[el];
                    TPZCompEl* cel = multiphysicsmesh->Element(elindex);
                    if(cel->Material()->Id()==1){
                        
                        int64_t ncon = cel->NConnects();
                        averagepressure.insert(cel->ConnectIndex(ncon-1));
                        TPZConnect &con=cel->Connect(ncon-1);
                        con.SetCondensed(false);
                        
                    }
                }
                if(averagepressure.size()!=1) DebugStop();
            }
        }
        
        // Counts the number of internal connects associated with a color
        for (int64_t patch = 0; patch < npatch; patch++) {
            for(auto cindex:fVecVecPatches[color][patch].fConnectIndices)
            {
                TPZConnect &c = multiphysicsmesh->ConnectVec()[cindex];
                if(!c.HasDependency() && !c.IsCondensed()){
                    numintconnects++;
                }
            }
        }
        
        int64_t extseqcount = numintconnects;
        TPZFMatrix<STATE> &weightsol = meshpatch->Solution();
        std::set<int64_t> boundaryconnects;
        for (int64_t patch = 0; patch < npatch; patch++) {
            {//Modifies the solution coeficients corresponding to hat functions of a color. They satisfy the unity partition property
                int64_t partitionindex = fVecVecPatches[color][patch].fPartitionConnectIndex;
                TPZConnect& c =  meshpatch->ConnectVec()[partitionindex];
                int64_t seqnum = c.SequenceNumber(); 
                weightsol.at(meshpatch->Block().at(seqnum,0,0,0)) = 1.;
            }
            
            
            
            {//Define the activ computational elements of the multiphysicsmesh
                int64_t nel = fVecVecPatches[color][patch].fElIndices.size();
                for (int64_t el = 0; el < nel; el++) {
                    int64_t elindex = fVecVecPatches[color][patch].fElIndices[el];
                    activel[elindex] = elpointers[elindex];
                }
            }
            
            {//Loop through the internal connects of a patch to count the number of equations associated to the color
                int64_t ncon = fVecVecPatches[color][patch].fConnectIndices.size();
                for (int64_t ic = 0; ic < ncon; ic++) {
                    int64_t cindex = fVecVecPatches[color][patch].fConnectIndices[ic];
                    TPZConnect &c = multiphysicsmesh->ConnectVec()[cindex];
                    int64_t seqnum = c.SequenceNumber();
                    if (seqnum != -1)
                    {
                        int64_t consize = multiphysicsmesh->ConnectVec()[cindex].NShape() * multiphysicsmesh->ConnectVec()[cindex].NState();
                        nintequations += consize;
                        if(!c.HasDependency() && !c.IsCondensed())
                        {
                            if(permute[seqnum] != -1) DebugStop();
                            
                            permute[seqnum] = seqcount++;
                            totalequations += consize;
                            
                        }
                    }
                }
            }
            
            meshpatch->LoadSolution(meshpatch->Solution());// Necessary to expand solution to hanging nodes
            
            {//Loop over boundary connects of a patch
                int64_t ncon = fVecVecPatches[color][patch].fBoundaryConnectIndices.size();
                for (int64_t ic = 0; ic < ncon; ic++) {
                    
                    int64_t cindex = fVecVecPatches[color][patch].fBoundaryConnectIndices[ic];
                    
                    TPZConnect &c = multiphysicsmesh->ConnectVec()[cindex];
                    c.SetCondensed(true);
                    int64_t seqnum = c.SequenceNumber();
                    if (seqnum != -1)
                    {
                        if (permute[seqnum] != -1) {
                            DebugStop();
                        }
                        int64_t consize = multiphysicsmesh->ConnectVec()[cindex].NShape() * multiphysicsmesh->ConnectVec()[cindex].NState();
                        if(!c.HasDependency())
                        {
                            permute[seqnum] = extseqcount++;
                            boundaryconnects.insert(cindex);
                            //totalequations += consize;
                        }
                    }
                    else{
                        DebugStop();
                    }
                }
            }

            if(0){//MARK: Deactivate average pressure space on boundary patch
                TPZPatch& Patch = fVecVecPatches[color][patch];
                std::cout << "fpatchIsBoundary: " << Patch.fPatchIsBoundary << std::endl;
                if(Patch.fPatchIsBoundary){
                    TPZManVector<int,5> activ2 = {1,1,0,0,0};
                    int nel = Patch.fElIndices.size();
                    for(int i=0; i<nel; i++){
                        int64_t index = Patch.fElIndices[i];
                        TPZCompEl* cel = multiphysicsmesh->ElementVec()[index];
                        TPZMultiphysicsElement* multicel = dynamic_cast<TPZMultiphysicsElement*>(cel);
                        multicel->SetActiveApproxSpaces(activ2);
                    }
                }
            }
            
        }//patch loop
        
        for(auto it : boundaryconnects){
            TPZConnect &c=multiphysicsmesh->ConnectVec()[it];
            if(!c.IsCondensed()) DebugStop();
        }
        
#ifdef ERRORESTIMATION_DEBUG
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
        
        for (int64_t ic = 0; ic < nblocks; ic++) {
            if (permute[ic] == -1) {
                permute[ic] = extseqcount++;
            }
        }
        
        if(extseqcount != nblocks) DebugStop();
        
        multiphysicsmesh->Permute(permute);//Renumbering the equations
        
        int64_t nactiveel = 0;
        //Activates multiphysicsmesh comp. elements of a color, and counts them
        for (int64_t el = 0; el < nelem; el++) {
            multiphysicsmesh->ElementVec()[el] = activel[el];
            if (activel[el]) {
                nactiveel++;
            }
        }
        //Computes the number of elements connected to each connect
        multiphysicsmesh->ComputeNodElCon();
        //multiphysicsmesh->CleanUpUnconnectedNodes();
        //multiphysicsmesh->SaddlePermute();
        
        if(0){
            std::ofstream out("../multiphysicsmeshcolor.txt");
            multiphysicsmesh->Print(out);
        }
        
        int64_t nequations = multiphysicsmesh->NEquations();
        std::cout << "nequations: " << nequations << std::endl;
        std::cout << "totalequations: " << totalequations << std::endl;

        if (nequations != totalequations) {
            DebugStop();
        }
        
        // Saddle permute would put the pressure equations after the boundary equations (and break the code)
        //        multiphysicsmesh->SaddlePermute();
        
        multiphysicsmesh->ExpandSolution();
        
        std::cout << "Color " << color << " Number of active elements " << nactiveel << " Number of equations " << nequations
        << "\nNumber of internal equations " << nintequations << std::endl;
        // PrintPartitionDiagnostics(color, std::cout);
        
        
        nequations = multiphysicsmesh->NEquations();
        
        std::cout<<"Solving mixed problem"<<std::endl;
        TPZLinearAnalysis an(multiphysicsmesh,RenumType::ENone);
         
        //TPZSkylineStructMatrix<STATE> strmat(multiphysicsmesh);
        TPZFStructMatrix<STATE> strmat(multiphysicsmesh);
        //TPZSSpStructMatrix<STATE> strmat(multiphysicsmesh);
        
        int numthreads = 0;
        strmat.SetNumThreads(numthreads);
        
        strmat.SetEquationRange(0, nequations);
        an.SetStructuralMatrix(strmat);
        
        TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
        direct->SetDirect(ELDLt);
        an.SetSolver(*direct);
        delete direct;
        direct = 0;
        
//        if(0){// Iterative solver
//            //TPZMatrixSolver<double> *precond{nullptr};
//            TPZStepSolver<STATE> *iterative = new TPZStepSolver<STATE>;
//            iterative->SetSOR(40, 100.0, 0.0000001, 0);
//            an.SetSolver(*iterative);
//            delete iterative;
//            iterative = 0;
//        }
        
        an.Assemble();
        
        TPZMatrixSolver<STATE> &matsolver = an.MatrixSolver<STATE>();
        TPZAutoPointer<TPZMatrix<STATE>> globmat = matsolver.Matrix();

        if(1){
            std::stringstream matname, rhsname;
            matname << "colormatrix_" << color << ".nb";
            std::ofstream outmat(matname.str());
            globmat->Print("GK=",outmat,EMathematicaInput);
            
            TPZFMatrix<STATE> globalrhs = an.Rhs();
            rhsname << "colorrhs_" << color << ".nb";
            std::ofstream outrhs(rhsname.str());
            globalrhs.Print("F=",outrhs,EMathematicaInput);
            
//            { //Print matrix for VTK
//                TPZFMatrix<REAL> mat(100,100);
//                multiphysicsmesh->ComputeFillIn(100, mat);
//                VisualMatrix(mat, "MatrixColor.vtk");
//            }
        }
        
        for (int64_t p = 0; p < npatch; p++) {
            TPZPatch &patch = fVecVecPatches[color][p];
            if (!PatchHasBoundary(patch))
            {
                int64_t firstlagrangeequation = patch.FirstLagrangeEquation(multiphysicsmesh);
                STATE diag = globmat->GetVal(firstlagrangeequation, firstlagrangeequation);
                diag += 1.;
                globmat->Put(firstlagrangeequation, firstlagrangeequation, diag);
            }
        }
        
        an.Solve();
        
        TPZStepSolver<STATE> &step = dynamic_cast<TPZStepSolver<STATE> &>(an.MatrixSolver<STATE>());
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
        
        if(0){
            std::stringstream outt;
            outt << "compMeshColor" << color << ".txt";
            auto* fluxmesh = multiphysicsmesh->MeshVector()[0];
            std::ofstream out(outt.str());
            fluxmesh->Print(out);
            
            std::stringstream outtt;
            outtt << "CMeshPatch_color_" << color << ".txt";
            std::ofstream outt2(outtt.str());
            multiphysicsmesh->Print(outt2);
            
        }
        
        if(0){
            // now we have a partial solution (only in one color of the colors loop)
            /** Variable names for post processing */
            TPZStack<std::string> scalnames, vecnames;
            scalnames.Push("POrder");
            scalnames.Push("Pressure");
            scalnames.Push("PartialErrorColor");
            vecnames.Push("Flux");
            an.SetStep(color); //
            
            TPZVec<REAL> errorscolor(6,0);
            int n = multiphysicsmesh->NElements();
            multiphysicsmesh->ElementSolution().Resize(n,6);
            bool storeerrors = true;
            an.PostProcessError(errorscolor,storeerrors);
            int ModelDimension = multiphysicsmesh->Dimension();
            std::stringstream sout;
            sout << "../" << "Poisson" << ModelDimension << "HDiv" << ".vtk";
            an.DefineGraphMesh(ModelDimension,scalnames,vecnames,sout.str());
            an.PostProcess(0,multiphysicsmesh->Dimension());
            TPZFMatrix<STATE> sol;
            sol=multiphysicsmesh->Solution();
            std::stringstream name;
            name << "solution_" << color << ".nb";
            std::ofstream outt(name.str());
            sol.Print("pzsol=",outt,EMathematicaInput);
        }
        
        if(0){ //just to plot hat functions of one color
            if(color == 6 and meshpatch->NElements() > 12){
                //meshpatch->LoadSolution(meshpatch->Solution());// Necessary to expand solution to hanging nodes
                
                TPZStack<std::string> scalnames2, vecnames2;
                scalnames2.Push("Solution");
                std::stringstream sout;
                sout << "Hatfunction" << ".vtk";
                TPZLinearAnalysis an2(meshpatch);
                an2.DefineGraphMesh(meshpatch->Dimension(), scalnames2, vecnames2, sout.str());
                int resolution = 0;
                an2.PostProcess(resolution, meshpatch->Dimension());
            }
        }
        
        //Recovery the initial comp. elements
        for (int64_t el = 0; el < nelem; el++) {
            multiphysicsmesh->ElementVec()[el] = elpointers[el];
        }
        
        multiphysicsmesh->ComputeNodElCon();
        
        TransferAndSumSolution(multiphysicsmesh);
        
        if(1){
            auto* fluxmesh = multiphysicsmesh->MeshVector()[0];
            std::ofstream out("hdivsol.txt");
            fluxmesh->Print(out);
        }
        
        for(auto it : boundaryconnects){
            TPZConnect &c=multiphysicsmesh->ConnectVec()[it];
            c.SetCondensed(false);
        }
        
        if(0){//MARK: Activate average pressure space
            TPZManVector<int,5> activ={1,1,0,0,1};
            for(auto it:activel){
                if(!it) continue;
                TPZMultiphysicsElement* multicel = dynamic_cast<TPZMultiphysicsElement*>(it);
                multicel->SetActiveApproxSpaces(activ);
            }
        }
        
        ResetState();
        an.SetStep(color);
        
         //Now draw the full mesh with the partial sum flux of the colors loop
        if(0){
            TPZStack<std::string> scalnames, vecnames;
            scalnames.Push("POrder");
            scalnames.Push("Pressure");
            vecnames.Push("Flux");
            int ModelDimension = multiphysicsmesh->Dimension();
            std::stringstream sout;
            sout << "../" << "FullPoisson" << ModelDimension << "HDiv" << ".vtk";
            an.DefineGraphMesh(ModelDimension,scalnames,vecnames,sout.str());
            an.PostProcess(3,multiphysicsmesh->Dimension());
        }
        
        fMeshVector[Epressure]->Solution().Zero();
        fMeshVector[Epatch]->Solution().Zero();
        //multiphysicsmesh->Solution().Print("solu1");
        an.Solution().Zero();
        //multiphysicsmesh->Solution().Print("solu2");
    }// colors loop
    
    if(0){// Just to check the property of the partition of unity
        TPZFMatrix<STATE> &weightsol = meshpatch->Solution();
        for (int color_ = 0; color_ < ncolors; color_++){
            int npatch = fVecVecPatches[color_].size();
            for (int patch_ = 0; patch_ < npatch; patch_++){
                int64_t partitionindex = fVecVecPatches[color_][patch_].fPartitionConnectIndex;
                TPZConnect& c =  meshpatch->ConnectVec()[partitionindex];
                int64_t seqnum = c.SequenceNumber();
                weightsol.at(meshpatch->Block().at(seqnum,0,0,0)) = 1.;
            }
        }
        meshpatch->LoadSolution(meshpatch->Solution());// Necessary to expand solution to hanging nodes
        TPZStack<std::string> scalnames2, vecnames2;
        scalnames2.Push("Solution");
        std::stringstream sout;
        sout << "SumHatfunctions" << ".vtk";
        TPZLinearAnalysis an2(meshpatch);
        an2.DefineGraphMesh(meshpatch->Dimension(), scalnames2, vecnames2, sout.str());
        int resolution = 0;
        an2.PostProcess(resolution, meshpatch->Dimension());
    }
    
    TPZLinearAnalysis an(multiphysicsmesh,RenumType::ENone);
    an.Solution() = fSolution;
    an.LoadSolution();
    
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("POrder");
    scalnames.Push("Pressure");
    vecnames.Push("Flux");
    std::stringstream sout;
    sout << "../" << "Reconstructed Flux" << multiphysicsmesh->Dimension() << "HDiv" << ".vtk";
    an.DefineGraphMesh(multiphysicsmesh->Dimension(),scalnames,vecnames,sout.str());
    an.PostProcess(1,multiphysicsmesh->Dimension());
    
    
    TPZManVector<TPZCompMesh *,2> mixed(2);
    mixed[0] = fMeshVector[Epressure];
    mixed[1] = fMeshVector[Epatch];
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mixed, multiphysicsmesh);

    std::ofstream outfmesh("Fluxmesh.txt");
    fMeshVector[Eflux]->Print(outfmesh);
    
    {
        int64_t nels = multiphysicsmesh->ElementVec().NElements();
        multiphysicsmesh->ElementSolution().Redim(nels, 6);
    }
    
    TPZManVector<REAL,6> errors(6,0.);
    an.PostProcessError(errors);
    cout << "Global flux estimated error " << errors[2] << std::endl;
    cout << "Global residual estimate error " << errors[3] << std::endl;

    ofstream myfile;
    myfile.open("ArquivosErros_estimate.txt", ios::app);
    myfile << "\n\n Estimator errors  \n";
    myfile << "\n-------------------------------------------------- \n";
    myfile << "Order = " << multiphysicsmesh->GetDefaultOrder() << "\n";
    myfile << "DOF Total = " << multiphysicsmesh->NEquations() << "\n";
    myfile << "Global flux estimated error = " << errors[2] << "\n";
    myfile << "Global residual estimate error " << errors[3] << "\n";

    myfile.close();
    
}

// print partition diagnostics
void TPZPostProcessError::PrintPartitionDiagnostics(int64_t color, std::ostream &out) const
{
    TPZCompMesh *multiphysicsmesh = fMeshVector[Emulti];
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
            TPZConnect &c = multiphysicsmesh->ConnectVec()[patch.fConnectIndices[ci]];
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
            TPZConnect &c = multiphysicsmesh->ConnectVec()[patch.fBoundaryConnectIndices[ci]];
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
    TPZCompMesh *meshmixed = fMeshVector[Eflux];
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

// Sum the solution stored in fSolution of the multiphysics mesh to the fSolution vector
void TPZPostProcessError::TransferAndSumSolution(TPZCompMesh *cmesh)
{
    int64_t nconnect = cmesh->NConnects();
    for (int64_t ic=0; ic<nconnect; ic++) {
        TPZConnect &c = cmesh->ConnectVec()[ic];
        if (c.SequenceNumber() == -1) {
            continue;
        }
//        if (c.NElConnected() == 0) {
//            DebugStop();
//        }
        int64_t seqnum = c.SequenceNumber();
        int neq = c.NDof();
        int64_t pos = cmesh->Block().Position(seqnum);
        
        int64_t targetseqnum = fConnectSeqNumbers[ic];
        int64_t targetpos = fBlock.Position(targetseqnum);
        
        
        for (int eq=0; eq<neq; eq++) {
            // tototototo
            TPZFMatrix<STATE> &sol = cmesh->Solution();
            fSolution(targetpos+eq,0) += sol(pos+eq,0);
            //fSolution.Print(std::cout);
        }
    }
    cmesh->Solution().Zero();
}

// Reset the state of the HDiv mesh to its original structure
void TPZPostProcessError::ResetState()
{
    TPZCompMesh *multiphysics = fMeshVector[Emulti];
    multiphysics->Block() = this->fBlock;
    multiphysics->Solution() = this->fSolution;
    for (int64_t i=0; i<fConnectSeqNumbers.size(); i++) {
        int64_t seqnum = fConnectSeqNumbers[i];
        TPZConnect &c = multiphysics->ConnectVec()[i];
        int64_t prevseqnumber = c.SequenceNumber();
        int prevblocksize = this->fBlock.Size(prevseqnumber);
        c.SetSequenceNumber(seqnum);
        int blsize = this->fBlock.Size(seqnum);
        int neq = c.NShape()*c.NState();
        if (neq != blsize) {
            DebugStop();
        }
    }
    TPZManVector<TPZCompMesh *> meshvec(2);
    meshvec[0] = fMeshVector[Epressure];
    meshvec[1] = fMeshVector[Epatch];
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, multiphysics);
    
}

// check whether the connectsizes have changed
void TPZPostProcessError::CheckConnectSizes()
{
    TPZCompMesh *multiphysics = fMeshVector[Eflux];
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
    TPZCompMesh *cmeshroot = fMeshVector[Eorigin];
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
        TPZBndCondT<STATE> *bnd = dynamic_cast<TPZBndCondT<STATE> *>(mat);
        if (!bnd) {
            int matId = mat->Id();
            int nstate = mat->NStateVariables();
            TPZNullMaterial<STATE> *material = new TPZNullMaterial<STATE>(matId);
            material->SetDimension(dim);
            material->SetNStateVariables(nstate);
            cmesh->InsertMaterialObject(material);
        }
    }
    for (auto it:cmeshroot->MaterialVec()) {
        TPZMaterial *mat = it.second;
        TPZBndCondT<STATE> *bnd = dynamic_cast<TPZBndCondT<STATE> *>(mat);
        if(bnd)
        {
            TPZMaterialT<STATE> *matorig = dynamic_cast<TPZMaterialT<STATE> *>(bnd->Material());
            int matid = matorig->Id();
            TPZMaterialT<STATE> *matl2 = dynamic_cast<TPZMaterialT<STATE> *>(cmesh->FindMaterial(matid));
            TPZBndCondT<STATE> *bc = matl2->CreateBC(matl2, bnd->Id(), bnd->Type(), bnd->Val1(), bnd->Val2());
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
        intelHDiv->PRefine(order+1);
    }
    
    {
        int meshdim = cmesh->Dimension();
        int64_t nel = cmesh->NElements();
        for(int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if(gel->Dimension() < meshdim) continue;
            int nsides = gel->NSides();
            int ncon = cel->NConnects();
            int connorder = cel->Connect(ncon-1).Order();
            int maxorder = connorder;
            for(int ic = 0; ic<ncon-1; ic++) {
                int cord = cel->Connect(ic).Order();
                maxorder = cord > maxorder ? cord : maxorder;
            }
            if(maxorder > connorder) {
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                intel->ForceSideOrder(nsides-1, maxorder);
            }
        }
        
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
    
    fMeshVector[Eflux] = cmesh;
}

/// create the lagrange mesh corresponding to the flux mesh
void TPZPostProcessError::CreatePressureMesh()
{
    TPZCompMesh *fluxmesh = fMeshVector[Eflux];
    TPZGeoMesh *gmesh = fluxmesh->Reference();
    int dim = fluxmesh->Dimension();
    
    // MARK: go stepwise here, this is a big change
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    for (auto it:fluxmesh->MaterialVec()) {
        TPZMaterial *mat = it.second;
        int matdim = mat->Dimension();
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        if (!bnd && matdim == dim) {
            int matId = mat->Id();
            int nstate = mat->NStateVariables();
            TPZVec<STATE> sol(nstate,0.);
            TPZNullMaterial<> *material = new TPZNullMaterial<>(matId,dim,nstate);
            //            TPZL2Projection *material = new TPZL2Projection(matId,dim,nstate,sol);
            cmesh->InsertMaterialObject(material);
        }
    }
    
//    cmesh->AutoBuild();
//    cmesh->CleanUpUnconnectedNodes();
//    cmesh->ExpandSolution();
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
        //cmesh->SetDefaultOrder(c.Order());
        TPZCompEl *cel = cmesh->ApproxSpace().CreateCompEl(gel, *cmesh);
        TPZInterpolatedElement *intel2 = dynamic_cast<TPZInterpolatedElement*>(cel);
        if(!intel) continue;
        intel2->PRefine(c.Order());
        cel->Reference()->ResetReference();
    }
    //    cmesh->LoadReferred(fMeshVector[0]);
    cmesh->InitializeBlock();
    for (auto &c:cmesh->ConnectVec()) {
        c.SetLagrangeMultiplier(1);
    }
    
    fMeshVector[Epressure] = cmesh;
}

void TPZPostProcessError::CreateAveragePressureMesh(){
    TPZCompMesh *pressuremesh = fMeshVector[Epressure];
    TPZGeoMesh *gmesh = pressuremesh->Reference();
    int dim = pressuremesh->Dimension();
    TPZCompMesh *averagepressuremesh = new TPZCompMesh(gmesh);
    int nstate = 1;
    {
        for (auto it : pressuremesh->MaterialVec()) {
            TPZMaterial* mat = it.second;
            int matid = mat->Id();
            TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(matid,dim,nstate);
            averagepressuremesh->SetDimModel(dim);
            nullmat->SetNStateVariables(1);
            averagepressuremesh->InsertMaterialObject(nullmat);
        }
    }
    //averagepressuremesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    averagepressuremesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    averagepressuremesh->SetDefaultOrder(0);
    averagepressuremesh->AutoBuild();
    
    int64_t nconnects = averagepressuremesh->NConnects();
    for (int ic = 0; ic<nconnects; ic++) {            averagepressuremesh->ConnectVec()[ic].SetLagrangeMultiplier(3);
    }
    
    int64_t nel = averagepressuremesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(averagepressuremesh->Element(el));
        if(disc) disc->SetFalseUseQsiEta();
    }
    
    fMeshVector.Resize(6);
    fMeshVector[Epressureaverage] = averagepressuremesh;
}


/// create the partition of unity mesh
void TPZPostProcessError::CreatePartitionofUnityMesh()
{
    //    TPZCompMeshReferred *pressuremesh = dynamic_cast<TPZCompMeshReferred *> (fMeshVector[3]);
    TPZCompMesh *pressuremesh = fMeshVector[Epressure];
    TPZGeoMesh *gmesh = pressuremesh->Reference();
    int dim = pressuremesh->Dimension();
    
    //    TPZCompMeshReferred *cmesh = new TPZCompMeshReferred(gmesh);
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetAllCreateFunctionsContinuous();
    //    cmesh->ApproxSpace().SetAllCreateFunctionsContinuousReferred();
    cmesh->SetDefaultOrder(1);
    
    for (auto it:pressuremesh->MaterialVec()) {
        TPZMaterial *mat = it.second;
        int matdim = mat->Dimension();
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        if (!bnd && matdim == dim) {
            int matId = mat->Id();
            int nstate = 1;
            TPZNullMaterial<> *material = new TPZNullMaterial<>(matId,dim,nstate);
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
        cmesh->ApproxSpace().CreateCompEl(gel, *cmesh);
    }
    //    cmesh->LoadReferred(fMeshVector[0]);
    //    pressuremesh->LoadReferred(cmesh);
    cmesh->InitializeBlock();
    fMeshVector[Epatch] = cmesh;
    
}

/// create the multiphysics mesh that will compute the projection matrix
void TPZPostProcessError::CreateMixedMesh()
{
    // the H1 mesh is the rootmesh
    TPZCompMesh *cmeshroot = fMeshVector[Eorigin];
    TPZGeoMesh *gmesh = fMeshVector[Eorigin]->Reference();
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZMultiphysicsCompMesh *mphysics = new TPZMultiphysicsCompMesh(gmesh);
    
    //create material
    int dim = gmesh->Dimension();
    typedef TPZMixedDarcyFlow TPZMixedPoisson;
    
    for (auto it:cmeshroot->MaterialVec()) {
        TPZMaterialT<STATE> *mat = dynamic_cast<TPZMaterialT<STATE> *>(it.second);
        TPZBndCondT<STATE> *bnd = dynamic_cast<TPZBndCondT<STATE>*>(mat);
        if (!bnd) {
            int matId = mat->Id();
            int nstate = mat->NStateVariables();
            TPZMaterialT<STATE> *material = 0;
            if (nstate == 1) {
                TPZMixedErrorEstimate<TPZMixedPoisson> *locmat = new TPZMixedErrorEstimate<TPZMixedPoisson>(matId,dim);
                locmat->SetSignConvention(1);
                material = locmat;
                locmat->SetForcingFunction(mat->ForcingFunction(),mat->ForcingFunctionPOrder());
                
                if(fExact){ //compel setintegrationrule new int...
                    locmat->SetExactSol(fExact->ExactSolution(), mat->ForcingFunctionPOrder());
                    
                    if(0){
                        auto exactsol= locmat->ExactSol();
                        TPZVec<STATE> x(3,0);
                        x[0]=0.5;
                        x[1]=0.5;
                        TPZFMatrix<STATE> du(3,1,0);
                        TPZVec<STATE> u(3,0);
                        
                        fExact->ExactSolution()(x,u,du);
                        std::cout << "u[0]=" << u[0] << std::endl;
                        exactsol(x,u,du);
                    }
                }
                
                //incluindo os dados do problema
                locmat->SetConstantPermeability(1.);
            }
            
            mphysics->InsertMaterialObject(material);
        }
    }
    for (auto it:cmeshroot->MaterialVec()) {
        TPZMaterialT<STATE> *mat = dynamic_cast<TPZMaterialT<STATE>*>(it.second);
        TPZBndCondT<STATE> *bnd = dynamic_cast<TPZBndCondT<STATE> *>(mat);
        if(bnd)
        {
            TPZMaterialT<STATE> *matorig = dynamic_cast<TPZMaterialT<STATE> *>(bnd->Material());
            int matid = matorig->Id();
            TPZMaterialT<STATE> *matmixed = dynamic_cast<TPZMaterialT<STATE> *>(mphysics->FindMaterial(matid));
            //matmixed->SetForcingFunction(matorig->ForcingFunction(),10);
            //bnd->IntegrationRuleOrder(10);
            TPZBndCondT<STATE> *bc = matmixed->CreateBC(matmixed, bnd->Id(), bnd->Type(), bnd->Val1(), bnd->Val2());
            //bc->IntegrationRuleOrder(15);
            //bc->SetForcingFunctionBC(matorig->ForcingFunction(),10);
            mphysics->InsertMaterialObject(bc);
        }
    }
    
    TPZManVector<TPZCompMesh *> meshvec(5,0);
    meshvec[0] = fMeshVector[Eflux];
    meshvec[1] = fMeshVector[Epressure];
    meshvec[2] = fMeshVector[Epatch];
    meshvec[3] = fMeshVector[Eorigin];
    meshvec[4] = fMeshVector[Epressureaverage];
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    
    //Fazendo auto build
    TPZManVector<int> activ = {1,1,0,0,1};
    mphysics->BuildMultiphysicsSpace(activ, meshvec);
    
    mphysics->CleanUpUnconnectedNodes();
    
    //------- Create and add group elements -------
    fMeshVector[Emulti] = mphysics;
}

// Create the meshes that allow us to compute the error estimate
void TPZPostProcessError::CreateAuxiliaryMeshes()
{
    CreateFluxMesh();
    CreatePressureMesh();
    CreatePartitionofUnityMesh();
    CreateAveragePressureMesh();
    CreateMixedMesh();
    
    TPZCompMesh* cmeshmulti = fMeshVector[Emulti];
    this->fSolution = cmeshmulti->Solution();
    this->fBlock = cmeshmulti->Block();
    int64_t ncon = cmeshmulti->NConnects();
    fConnectSeqNumbers.resize(ncon); //Resize the vector of the original connect sequence numbers
    fConnectSizes.resize(ncon); //Resize the vector of connects sizes in the multiphysics mesh
    
    for (int64_t i = 0; i < ncon; i++) {
        fConnectSeqNumbers[i] = cmeshmulti->ConnectVec()[i].SequenceNumber();
        int64_t seqnum = fConnectSeqNumbers[i];
        TPZConnect &con = cmeshmulti->ConnectVec()[i];
        con.SetSequenceNumber(seqnum);
        int blsize = this->fBlock.Size(seqnum);
        int neq = con.NShape()*con.NState();
        fConnectSizes[i] = neq;
        if (neq != blsize) {
            DebugStop();
        }
    }
    
    BuildPatchStructures2(); // Use BuildPatchStructures2() for one patch by color
    
#ifdef ERRORESTIMATION_DEBUG
    {
        std::ofstream out("../patchinfo.txt");
        PrintPatchInformation(out);
        std::ofstream outp("../partitionmesh.txt");
        fMeshVector[Epatch]->Print(outp);
        std::ofstream out2("../mphysics.txt");
        cmeshmulti->Print(out2);
    }
#endif
    
}

