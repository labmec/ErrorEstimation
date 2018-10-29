//
//  TPZHybridHDivErrorEstimator.cpp
//  ErrorEstimateHDiv
//
//  Created by Philippe Devloo on 10/06/18.
//

#include "TPZHybridHDivErrorEstimator.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "mixedpoisson.h"
#include "TPZMixedHdivErrorEstimate.h"
#include "pzbuildmultiphysicsmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzanalysis.h"

#include "TPZVTKGeoMesh.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("HDivErrorEstimator"));
#endif

TPZHybridHDivErrorEstimator::~TPZHybridHDivErrorEstimator()
{
    delete fPostProcMesh[1];
    delete fPostProcMesh[2];
    delete fPostProcMesh[0];
}

/// compute the element errors comparing the reconstructed solution based on average pressures
/// with the original solution
void TPZHybridHDivErrorEstimator::ComputeErrors(TPZVec<REAL> &elementerrors, bool store)
{
    /// I havent tested to compute post processing more than a single time
    // Please test me!
    if (fPostProcMesh[0]) {
        DebugStop();
    }
    CreatePostProcessingMesh();
    ComputeElementStiffnesses();
    if(fProblemConfig.makepressurecontinuous)
    {
        ComputeAveragePressures();
    }
    ComputeNodalAverages();
    // fPostProcMesh[0] is the multiphysics mesh
    fPostProcMesh[0]->LoadSolution(fPostProcMesh[0]->Solution());
    TPZManVector<TPZCompMesh *,2> meshvec(2);
    // fPostProcMesh[1] is the flux mesh
    // fPostProcMesh[2] is the pressure mesh
    meshvec[0] = fPostProcMesh[1];
    meshvec[1] = fPostProcMesh[2];
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, fPostProcMesh[0]);
#ifdef PZDEBUG
    {
        std::ofstream out("postprocmesh.txt");
        fPostProcMesh[0]->Print(out);
        std::ofstream out1("fluxmesh.txt");
        fPostProcMesh[1]->Print(out1);
        std::ofstream out2("pressuremesh.txt");
        fPostProcMesh[2]->Print(out2);
    }
#endif
    {
        TPZAnalysis an(fPostProcMesh[0],false);
        if (fExact) {
            an.SetExact(fExact->ExactSolution());
        }
        TPZManVector<REAL> errorvec(6,0.);
        int64_t nelem = fPostProcMesh[0]->NElements();
        fPostProcMesh[0]->ElementSolution().Redim(nelem, 4);
        an.PostProcessError(errorvec);
        ComputeEffectivityIndices();
        std::cout << "Computed errors " << errorvec << std::endl;
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("PressureFem");
        scalnames.Push("PressureReconstructed");
        scalnames.Push("PressureExact");
        scalnames.Push("PressureErrorExact");
        scalnames.Push("PressureErrorEstimate");
        scalnames.Push("EnergyErrorExact");
        scalnames.Push("EnergyErrorEstimate");
        scalnames.Push("PressureEffectivityIndex");
        scalnames.Push("EnergyEffectivityIndex");
        vecnames.Push("FluxFem");
        vecnames.Push("FluxReconstructed");
        vecnames.Push("FluxExact");
        int dim = 2;
        std::string plotname;
        {
            std::stringstream out;
            out << "HybridPostProcessed_P" << fProblemConfig.porder << "_" << dim << "D_" << fProblemConfig.problemname <<
            "_NEL" << nelem <<".vtk";
            plotname = out.str();
        }
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(0,dim);
        an.SetStep(1);
        an.PostProcess(2,dim);
    }
    {
        TPZAnalysis an(fPostProcMesh[2],false);
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("State");
        int dim = this->fOriginal[0]->Reference()->Dimension()-1;
        std::string plotname;
        {
            std::stringstream out;
            out << "HybridPostProcessed_P" << fProblemConfig.porder << "_" << dim << "D_" << fProblemConfig.problemname << ".vtk";
            plotname = out.str();
        }
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(2,dim);
    }

}



/// create the post processed multiphysics mesh (which is necessarily hybridized)
void TPZHybridHDivErrorEstimator::CreatePostProcessingMesh()
{
    int nmeshes = fPostProcMesh.size();
    TPZManVector<TPZCompMesh *, 4> meshvec_Hybrid(nmeshes-1, 0);
    CloneMeshVec();
#ifdef PZDEBUG
    {
        std::ofstream out("CloneFluxMesh.txt");
        fPostProcMesh[1]->Print(out);
        std::ofstream outp("ClonePressureMesh.txt");
        fPostProcMesh[2]->Print(outp);
    }
#endif
    IncreaseSideOrders(fPostProcMesh[1]);
    for(int i=1; i<nmeshes; i++)
    {
        meshvec_Hybrid[i-1] = fPostProcMesh[i];
    }
#ifdef PZDEBUG
    {
        std::ofstream out("CloneFluxMeshBefore.txt");
        fPostProcMesh[1]->Print(out);
    }
#endif

    fHybridizer.ComputePeriferalMaterialIds(meshvec_Hybrid);
    fHybridizer.ComputeNState(meshvec_Hybrid);
    /// insert the material objects for HDivWrap and LagrangeInterface
    fHybridizer.InsertPeriferalMaterialObjects(meshvec_Hybrid);
    fHybridizer.HybridizeInternalSides(meshvec_Hybrid);
    TPZCompMeshReferred *RefFluxMesh = dynamic_cast<TPZCompMeshReferred *>(fPostProcMesh[1]);
    RefFluxMesh->LoadReferred(fOriginal[1]);
#ifdef PZDEBUG
    {
        std::ofstream out("CloneFluxMesh.txt");
        RefFluxMesh->Print(out);
    }
#endif
    TPZCompMeshReferred *RefPressureMesh = dynamic_cast<TPZCompMeshReferred *>(fPostProcMesh[2]);
    RefPressureMesh->LoadReferred(fOriginal[2]);
    CreateMultiphysicsMesh();
    TPZCompMesh *cmesh_Hybrid = fPostProcMesh[0];
    fHybridizer.CreateInterfaceElements(cmesh_Hybrid, meshvec_Hybrid);
#ifdef PZDEBUG
    {
        std::ofstream out("multiphysics.txt");
        cmesh_Hybrid->Print(out);
    }
#endif
    fHybridizer.GroupElements(cmesh_Hybrid);
    cmesh_Hybrid->CleanUpUnconnectedNodes();
#ifdef PZDEBUG
    {
        std::ofstream out("multiphysicsgrouped.txt");
        cmesh_Hybrid->Print(out);
        std::ofstream outvtk("multiphysics.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh_Hybrid,outvtk);
        std::ofstream outgvtk("postprocessgmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(cmesh_Hybrid->Reference(),outgvtk);
    }
#endif

}

/// computing the element stifnesses will "automatically" compute the condensed form of the matrices
void TPZHybridHDivErrorEstimator::ComputeElementStiffnesses()
{
    for (auto cel:fPostProcMesh[0]->ElementVec()) {
        if(!cel) continue;
        TPZElementMatrix ek, ef;
        cel->CalcStiff(ek, ef);
    }
}

/// increase the side orders of the post processing mesh
void TPZHybridHDivErrorEstimator::IncreaseSideOrders(TPZCompMesh *fluxmesh)
{
    int64_t nel = fluxmesh->NElements();
    int dim = fluxmesh->Dimension();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fluxmesh->Element(el);
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != dim) {
            continue;
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        int nc = cel->NConnects();
        int order = cel->Connect(nc - 1).Order();
        int nsides = gel->NSides();
        int ncorner = gel->NCornerNodes();
        intel->SetPreferredOrder(order);
        for (int side = ncorner; side < nsides - 1; side++) {
            if (intel->NSideConnects(side)) {
                intel->SetSideOrder(side, order);
            }
        }
        //        intel->Print();
    }
    fluxmesh->InitializeBlock();

}

/// compute the average pressures of the hybridized form of the H(div) mesh
void TPZHybridHDivErrorEstimator::ComputeAveragePressures()
{
    TPZCompMesh *pressure = fOriginal[2];
    TPZCompMesh *pressureHybrid = fPostProcMesh[2];
    int InterfaceMatid = fHybridizer.LagrangeInterface;
    TPZGeoMesh *gmesh = pressure->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    pressure->LoadReferences();
    int64_t nel = pressureHybrid->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = pressureHybrid->Element(el);
        if(!cel || !cel->Reference() || cel->Reference()->Dimension() != dim-1)
        {
            continue;
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        TPZGeoEl *gel = cel->Reference();
        if (gel->MaterialId() != InterfaceMatid) {
            continue;
        }
        if (!intel || gel->Dimension() != dim-1) {
            DebugStop();
        }
        int nc = cel->NConnects();
        int order = cel->Connect(nc-1).Order();
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 1, 0);
        TPZManVector<TPZTransform<REAL> ,2> tr(2);
        tr[0] = gelside.NeighbourSideTransform(celstack[0].Reference());
        {
            TPZGeoEl *right = celstack[0].Element()->Reference();
            TPZTransform<REAL> tmp = right->SideToSideTransform(celstack[0].Side(), right->NSides()-1);
            tr[0] = tmp.Multiply(tr[0]);
        }
        if (celstack.size() == 1) {
            TPZCompElSide lowlevel = gelside.LowerLevelCompElementList2(1);
            if (!lowlevel) {
                DebugStop();
            }
            celstack.Push(lowlevel);
            tr[1] = TPZTransform<REAL>(gelside.Dimension());
            gelside.SideTransform3(lowlevel.Reference(), tr[1]);
        }
        else if(celstack.size() == 2)
        {
            tr[1] = gelside.NeighbourSideTransform(celstack[1].Reference());
        }
        else
        {
            DebugStop();
        }
        {
            TPZGeoEl *right = celstack[1].Element()->Reference();
            TPZTransform<REAL> tmp = right->SideToSideTransform(celstack[1].Side(), right->NSides()-1);
            tr[1] = tmp.Multiply(tr[1]);
        }
        
        std::unique_ptr<TPZIntPoints> intp( gel->CreateSideIntegrationRule(gel->NSides()-1, 2*order));
        int nshape = intel->NShapeF();
        TPZFNMatrix<20,REAL> L2Mat(nshape,nshape,0.), L2Rhs(nshape,1,0.);
        TPZFNMatrix<220,REAL> phi(nshape,1,0.), dshape(dim,nshape);
        int64_t npoints = intp->NPoints();
        for (int64_t ip=0; ip<npoints; ip++) {
            TPZManVector<REAL,3> pt(dim-1,0.),pt1(dim,0.), pt2(dim,0.),sol1(1),sol2(1);
            REAL weight;
            intp->Point(ip, pt, weight);
            intel->Shape(pt, phi, dshape);
            tr[0].Apply(pt, pt1);
            tr[1].Apply(pt, pt2);
            celstack[0].Element()->Solution(pt1, 0, sol1);
            celstack[1].Element()->Solution(pt2, 0, sol2);
            //            std::cout << "Values " << sol1 << " " << sol2 << std::endl;
            for (int ishape=0; ishape<nshape; ishape++) {
                L2Rhs(ishape,0) += weight*phi(ishape,0)*(sol1[0]+sol2[0])/2.;
                for (int jshape = 0; jshape<nshape; jshape++) {
                    L2Mat(ishape,jshape) += weight*phi(ishape,0)*phi(jshape,0);
                }
            }
        }
        L2Mat.SolveDirect(L2Rhs, ECholesky);
        //        L2Rhs.Print("Average pressure");
        int count = 0;
        for (int ic=0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            int64_t seqnum = c.SequenceNumber();
            int64_t pos = pressureHybrid->Block().Position(seqnum);
            int ndof = c.NShape()*c.NState();
            for (int idf = 0; idf<ndof; idf++) {
                pressureHybrid->Solution()(pos+idf,0) = L2Rhs(count++);
            }
        }
    }
    TPZManVector<TPZCompMesh *,2> meshvec(2);
    meshvec[0] = fPostProcMesh[1];
    meshvec[1] = fPostProcMesh[2];
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, fPostProcMesh[0]);
}

/// set the cornernode values equal to the averages
void TPZHybridHDivErrorEstimator::ComputeNodalAverages()
{
    TPZCompMesh *pressureHybrid = fPostProcMesh[2];
    int InterfaceMatid = fHybridizer.LagrangeInterface;
    TPZGeoMesh *gmesh = pressureHybrid->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    pressureHybrid->LoadReferences();
    int lagrangematid = fHybridizer.LagrangeInterface;
    TPZMaterial *mat = pressureHybrid->FindMaterial(lagrangematid);
    if(!mat) DebugStop();
    int nstate = mat->NStateVariables();
    if(dim != 2)
    {
        DebugStop();
    }
    int64_t nel = pressureHybrid->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = pressureHybrid->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        if(gel->Dimension() != 1 || gel->MaterialId() != InterfaceMatid)
        {
            continue;
        }
        for (int side=0; side<2; side++) {
            TPZGeoElSide gelside(gel,side);
            TPZStack<TPZCompElSide> celstack;
            gelside.ConnectedCompElementList(celstack, 1, 0);
            celstack.Push(gelside.Reference());
            TPZManVector<STATE,3> averageval(nstate,0.);
            std::set<int64_t> connects;
            for (int elc=0; elc<celstack.size(); elc++) {
                TPZCompElSide celside = celstack[elc];
                TPZGeoElSide gelside0 = celside.Reference();
                if (gelside0.Element()->Dimension() != 1) {
                    continue;
                }
                TPZInterpolatedElement *intel1 = dynamic_cast<TPZInterpolatedElement *>(celside.Element());
                if(!intel1 || intel1->NConnects() != 3) DebugStop();
                int64_t conindex = intel1->ConnectIndex(celside.Side());
                connects.insert(conindex);
                TPZConnect &c = intel1->Connect(celside.Side());
                int64_t seqnum = c.SequenceNumber();
                if(c.NState() != nstate || c.NShape() != 1) DebugStop();
                for (int istate = 0; istate<nstate; istate++) {
                    averageval[istate] += pressureHybrid->Block().Get(seqnum, 0, istate, 0);
                }
            }
            auto ncontr = connects.size();
            for (int istate = 0; istate<nstate; istate++) {
                averageval[istate] /= ncontr;
            }
            for(auto conindex : connects)
            {
                TPZConnect &c = pressureHybrid->ConnectVec()[conindex];
                int64_t seqnum = c.SequenceNumber();
                if(c.NState() != nstate || c.NShape() != 1) DebugStop();
                for (int istate = 0; istate<nstate; istate++) {
#ifdef LOG4CXX
                    if (logger->isDebugEnabled())
                    {
                        std::stringstream sout;
                        sout << "value before " << pressureHybrid->Block()(seqnum, 0, istate, 0) <<
                        " value after " << averageval[istate] << " diff " << pressureHybrid->Block()(seqnum, 0, istate, 0)-averageval[istate] << " ncontr " << ncontr;
                            //            res2.Print("Residual",sout);
                        LOGPZ_DEBUG(logger,sout.str())
                    }
#endif

                    pressureHybrid->Block()(seqnum, 0, istate, 0) = averageval[istate];
                }
            }
        }
    }
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = pressureHybrid->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        if(gel->Dimension() != 1 || gel->MaterialId() != InterfaceMatid)
        {
            continue;
        }
        for (int side=0; side<2; side++) {
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (gelside != neighbour) {
                if(IsDirichletCondition(neighbour))
                {
                    TPZConnect &c = intel->Connect(side);
                    int nstate = c.NState();
                    TPZManVector<STATE,3> vals(nstate,0.);
                    GetDirichletValue(neighbour, vals);
                    int64_t seqnum = c.SequenceNumber();
                    for (int ist = 0; ist<nstate; ist++) {
#ifdef LOG4CXX
                        if (logger->isDebugEnabled())
                        {
                            std::stringstream sout;
                            sout << "value before " << pressureHybrid->Block()(seqnum, 0, ist, 0) <<
                            " value after " << vals[ist];
                            //            res2.Print("Residual",sout);
                            LOGPZ_DEBUG(logger,sout.str())
                        }
#endif
                        pressureHybrid->Block()(seqnum,0,ist,0) = vals[ist];
                    }
                }
                neighbour = neighbour.Neighbour();
            }
        }
    }
    TPZManVector<TPZCompMesh *,2> meshvec(2);
    meshvec[0] = fPostProcMesh[1];
    meshvec[1] = fPostProcMesh[2];
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, fPostProcMesh[0]);
}



/// clone the meshes into the post processing mesh
void TPZHybridHDivErrorEstimator::CloneMeshVec()
{
    CreateFluxMesh();
    CreatePressureMesh();
//    for (int i = 3; i < fOriginal.size(); i++) {
//        fPostProcMesh[i] = fOriginal[i]->Clone();
//    }

}

/// create the multiphysics mesh using TPZCompMeshReferred and TPZMixedErrorEstimate material
void TPZHybridHDivErrorEstimator::CreateMultiphysicsMesh()
{
    // the pressure mesh is the rootmesh
    TPZCompMesh *cmeshroot = fOriginal[0];
    TPZGeoMesh *gmesh = cmeshroot->Reference();
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    fHybridizer.InsertPeriferalMaterialObjects(mphysics);

    //criando material
    int dim = gmesh->Dimension();
    
    for (auto it:cmeshroot->MaterialVec()) {
        TPZMaterial *mat = it.second;
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        if (!bnd) {
            int matId = mat->Id();
            if (fHybridizer.IsPeriferalMaterialId(matId)) {
                continue;
            }
            int nstate = mat->NStateVariables();
            TPZMaterial *material = 0;
            if (nstate == 1) {
                TPZMixedPoisson *mix = dynamic_cast<TPZMixedPoisson *>(mat);
                if(!mix)
                {
                    DebugStop();
                }
                TPZMixedHDivErrorEstimate<TPZMixedPoisson> *locmat = new TPZMixedHDivErrorEstimate<TPZMixedPoisson>(matId,dim);
                material = locmat;
                locmat->SetForcingFunction(mat->ForcingFunction());
                if(fExact)
                {
                    locmat->SetForcingFunctionExact(fExact->Exact());
                }
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
            if (bnd->ForcingFunction()) {
                bc->TPZDiscontinuousGalerkin::SetForcingFunction(bnd->ForcingFunction());
            }
            mphysics->InsertMaterialObject(bc);
        }
    }
    
    TPZManVector<TPZCompMesh *> meshvec(2,0);
    meshvec[0] = fPostProcMesh[1];
    meshvec[1] = fPostProcMesh[2];
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    //Fazendo auto build
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    mphysics->CleanUpUnconnectedNodes();
#ifdef PZDEBUG
    {
        std::ofstream out("multiphysics.txt");
        mphysics->Print(out);
    }
#endif
    //------- Create and add group elements -------
    fPostProcMesh[0] = mphysics;

}

/// clone the fluxmesh but using TPZCompElReferred objects
void TPZHybridHDivErrorEstimator::CreateFluxMesh()
{
    TPZGeoMesh *gmesh = fOriginal[0]->Reference();
    TPZCompMeshReferred *fluxmesh = new TPZCompMeshReferred(gmesh);
    TPZCompMesh *origflux = fOriginal[1];
    origflux->CopyMaterials(*fluxmesh);
    gmesh->ResetReference();
    int meshdim = gmesh->Dimension();
    fluxmesh->ApproxSpace().SetAllCreateFunctionsHDivReferred(gmesh->Dimension());
    for (auto celorig:origflux->ElementVec()) {
        if(!celorig) continue;
        TPZGeoEl *gel = celorig->Reference();
        int64_t index;
        int geldim = gel->Dimension();
        fluxmesh->ApproxSpace().CreateCompEl(gel, *fluxmesh, index);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(fluxmesh->Element(index));
        for(int side=0; side<gel->NSides(); side++)
        {
            
            if(gel->SideDimension(side) < meshdim-1) continue;
            int conindex = intel->SideConnectLocId(0,side);
            TPZConnect &corig = celorig->Connect(conindex);
            TPZConnect &newcon = intel->Connect(conindex);
            if (newcon.Order() != corig.Order()) {
                intel->SetSideOrder(side, corig.Order());
            }
            if (side == gel->NSides()-1) {
                intel->SetPreferredOrder(newcon.Order());
            }
        }
    }
    if(fUpliftPostProcessMesh)
    {
        int64_t nel = fluxmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = fluxmesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if (gel->Dimension() == meshdim) {
                int side = gel->NSides() - 1;
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
                int conindex = intel->SideConnectLocId(0, side);
                TPZConnect &c = intel->Connect(conindex);
                int porder = c.Order();
                intel->SetSideOrder(side, porder + fUpliftPostProcessMesh);
                intel->SetPreferredOrder(porder+fUpliftPostProcessMesh);
            }
        }

    }
    fluxmesh->LoadReferred(origflux);
    fluxmesh->InitializeBlock();
    fPostProcMesh[1] = fluxmesh;
}

/// clone the pressure mesh using TPZCompElReferred objects
void TPZHybridHDivErrorEstimator::CreatePressureMesh()
{
    TPZGeoMesh *gmesh = fOriginal[0]->Reference();
    TPZCompMeshReferred *pressuremesh = new TPZCompMeshReferred(gmesh);
    TPZCompMesh *origpressure = fOriginal[2];
    origpressure->CopyMaterials(*pressuremesh);
    gmesh->ResetReference();
    int meshdim = gmesh->Dimension();
    pressuremesh->ApproxSpace().SetAllCreateFunctionsContinuousReferred();
    pressuremesh->ApproxSpace().CreateDisconnectedElements(true);
    for (auto celorig:origpressure->ElementVec()) {
        if(!celorig) continue;
        TPZGeoEl *gel = celorig->Reference();
        TPZInterpolatedElement *intelorig = dynamic_cast<TPZInterpolatedElement *>(celorig);
        int order = intelorig->GetPreferredOrder();
        int64_t index;
        TPZCompEl *cel = pressuremesh->ApproxSpace().CreateCompEl(gel, *pressuremesh, index);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        intel->PRefine(order);
        cel->Reference()->ResetReference();
    }
    if(fUpliftPostProcessMesh)
    {
        int64_t nel = pressuremesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = pressuremesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if (gel->Dimension() == meshdim) {
//                int side = gel->NSides() - 1;
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
                int conindex = intel->NConnects()-1;
                TPZConnect &c = intel->Connect(conindex);
                int porder = c.Order();
                intel->PRefine(porder+fUpliftPostProcessMesh);
                intel->SetPreferredOrder(porder+fUpliftPostProcessMesh);
            }
        }
    }

    pressuremesh->LoadReferred(origpressure);
    pressuremesh->InitializeBlock();
    fPostProcMesh[2] = pressuremesh;
}

/// compute the effectivity indices of the pressure error and flux error and store in the element solution
void TPZHybridHDivErrorEstimator::ComputeEffectivityIndices()
{
    TPZCompMesh *cmesh = fPostProcMesh[0];
    if (cmesh->ElementSolution().Cols() != 4) {
        DebugStop();
    }
    int64_t nrows = cmesh->ElementSolution().Rows();
    cmesh->ElementSolution().Resize(nrows, 6);
    for (int64_t el=0; el<nrows; el++) {
        for(int i=0; i<3; i+=2)
        {
            if(!IsZero(cmesh->ElementSolution()(el,i)))
            {
                cmesh->ElementSolution()(el,4+i/2) = cmesh->ElementSolution()(el,i+1)/cmesh->ElementSolution()(el,i);
            }
            else
            {
                cmesh->ElementSolution()(el,4+i/2) = 1.;
            }
        }
    }
}

/// returns true if the material associated with the element is a boundary condition
/// and if the boundary condition is dirichlet type
bool TPZHybridHDivErrorEstimator::IsDirichletCondition(TPZGeoElSide gelside)
{
    TPZGeoEl *gel = gelside.Element();
    int matid = gel->MaterialId();
    TPZMaterial *mat = fPostProcMesh[0]->FindMaterial(matid);
    TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
    if(!bc) return false;
    int typ = bc->Type();
    if(typ == 0) return true;
    return false;
}

/// return the value of the Dirichlet condition
void TPZHybridHDivErrorEstimator::GetDirichletValue(TPZGeoElSide gelside, TPZVec<STATE> &vals)
{
    TPZGeoEl *gel = gelside.Element();
    int matid = gel->MaterialId();
    TPZMaterial *mat = fPostProcMesh[0]->FindMaterial(matid);
    TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
    if(!bc) DebugStop();
    int typ = bc->Type();
    if(typ != 0) DebugStop();
    TPZManVector<REAL,3> xco(3,0.);
    if (bc->HasForcingFunction()) {
        gel->NodePtr(gelside.Side())->GetCoordinates(xco);
        bc->ForcingFunction()->Execute(xco, vals);
    }
    else
    {
        int nv = vals.size();
        for (int iv=0; iv<nv; iv++) {
            vals[iv] = bc->Val2()(iv,0);
        }
    }
}
