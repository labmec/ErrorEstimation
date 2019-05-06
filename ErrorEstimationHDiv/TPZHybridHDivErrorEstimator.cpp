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
#include "pzmat1dlin.h"
#include "TPZVecL2.h"
#include "TPZNullMaterial.h"

#include "TPZParFrontStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzskylstrmatrix.h"

#include "TPZMultiphysicsCompMesh.h"
#include "pzmultiphysicscompel.h"

#include "TPZVTKGeoMesh.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("HDivErrorEstimator"));
#endif

TPZHybridHDivErrorEstimator::~TPZHybridHDivErrorEstimator()
{
    TPZVec<TPZCompMesh *> meshvec = fPostProcMesh.MeshVector();
    for (int i = 0; i<2; i++) {
        delete meshvec[i];
    }
}

/// compute the element errors comparing the reconstructed solution based on average pressures
/// with the original solution
void TPZHybridHDivErrorEstimator::ComputeErrors(TPZVec<REAL> &elementerrors, bool store)
{
        TPZAnalysis an(&fPostProcMesh,false);
    
        if (fExact) {
            an.SetExact(fExact->ExactSolution());
        }


        TPZManVector<REAL> errorvec(6,0.);
        int64_t nelem = fPostProcMesh.NElements();
        fPostProcMesh.LoadSolution(fPostProcMesh.Solution());
        fPostProcMesh.ExpandSolution();
        fPostProcMesh.ElementSolution().Redim(nelem, 4);
    
        an.PostProcessError(errorvec);//calculo do erro com sol exata e aprox

        std::cout << "Computed errors " << errorvec << std::endl;
        
        ComputeEffectivityIndices();
    
        PostProcessing(an);
        
    
    
}
void TPZHybridHDivErrorEstimator::PostProcessing(TPZAnalysis &an){
    
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
    int dim = fPostProcMesh.Reference()->Dimension();
    std::string plotname;
    {
        std::stringstream out;
        out << "HybridPostProcessed_P" << fProblemConfig.porder << "_" << dim << "D_" << fProblemConfig.problemname << ".vtk";
        plotname = out.str();
    }
    an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
    an.PostProcess(2,dim);
    //        an.SetStep(1);
    //        an.PostProcess(0,dim);
    {
        TPZAnalysis an(fPostProcMesh.MeshVector()[1],false);
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("State");
        int dim = this->fOriginal->Reference()->Dimension()-1;
        std::string plotname;
        {
            std::stringstream out;
            out << "LagrangeMultiplierPostProces _" << fProblemConfig.problemname << ".vtk";
            plotname = out.str();
        }
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(2,dim);
    }
    
}

void TPZHybridHDivErrorEstimator::PostProcessingHybridMesh(){
    
    int nmeshes = fPostProcMesh.MeshVector().size();
    TPZManVector<TPZCompMesh *, 4> meshvec_Hybrid(nmeshes, 0);
    
    CloneMeshVec();//copia as malhas de fluxo e pressao
    
    
//    fPostProcMesh[1]=fOriginal[3];
//    fPostProcMesh[2]=fOriginal[4];
    
    
    //fPostProcMesh[0] ainda esta vazio
    
#ifdef PZDEBUG
    {
        std::ofstream out("NewFluxMesh.txt");
        fPostProcMesh.MeshVector()[0]->Print(out);
        std::ofstream outp("NewPressureMesh.txt");
        fPostProcMesh.MeshVector()[1]->Print(outp);
    }
#endif
    
    
    
    
    
//#ifdef PZDEBUG
//    {
//        std::ofstream out("FluxBorder.txt");
//        fPostProcMesh[1]->Print(out);
//        std::ofstream out2("PressureBorder.txt");
//        fPostProcMesh[2]->Print(out2);
//    }
//#endif
    
    //aumenta a ordem do fluxo de borda e pressao de borda
    //IncreaseSideOrders(fPostProcMesh[1]);
   // IncreasePressureSideOrders(fPostProcMesh[2]);

    
    for(int i=0; i<nmeshes; i++)
    {
        meshvec_Hybrid[i] = fPostProcMesh.MeshVector()[i];
    }
//#ifdef PZDEBUG
//    {
//        std::ofstream out("EnrichedFluxBorder.txt");
//        fPostProcMesh[1]->Print(out);
//        std::ofstream out2("EnrichedPressureBorder.txt");
//        fPostProcMesh[2]->Print(out2);
//    }
//#endif
    
    

    
    CreateMultiphysicsMesh();
    
    TPZCompMesh *cmesh_Hybrid = &fPostProcMesh;
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

/// create the post processed multiphysics mesh (which is necessarily hybridized)
void TPZHybridHDivErrorEstimator::CreatePostProcessingMesh()
{
    // initialize the post processing mesh
    fPostProcMesh.SetReference(fOriginal->Reference());
    fOriginal->CopyMaterials(fPostProcMesh);
    // switch the material from mixed to TPZMixedHdivErrorEstimate...
    SwitchMaterialObjects();

    TPZManVector<TPZCompMesh *,4> mesh_vectors(4,0);
    mesh_vectors[2] = fOriginal->MeshVector()[0];
    mesh_vectors[3] = fOriginal->MeshVector()[1];
    mesh_vectors[0] = mesh_vectors[2]->Clone();
    mesh_vectors[1] = mesh_vectors[3]->Clone();
    // clone the flux and pressure mesh
    if (!fOriginalIsHybridized)
    {
        fHybridizer.ComputePeriferalMaterialIds(mesh_vectors);
        fHybridizer.ComputeNState(mesh_vectors);
        fHybridizer.HybridizeInternalSides(mesh_vectors);
        
//        #ifdef PZDEBUG
//            {
//                std::ofstream out("OriginalFlux.txt");
//                fOriginal[3]->Print(out);
//                std::ofstream out2("OriginalPotential.txt");
//                fOriginal[4]->Print(out2);
//                std::ofstream out3("OriginalMeshHybrid.txt");
//                fOriginal[0]->Print(out3);
//            }
//        #endif
    }
    else
    {
        IdentifyPeripheralMaterialIds();
    }

   
    IncreasePressureSideOrders(mesh_vectors[1]);//malha da pressao
    IncreaseSideOrders(fPostProcMesh.MeshVector()[0]);//malha do fluxo


    
    TPZManVector<int> active(4,0);
    active[0] = 1;
    active[1] = 1;
    fPostProcMesh.BuildMultiphysicsSpace(active, mesh_vectors);
    // construction of the multiphysics mesh
        //cria elementos de interface
    fHybridizer.CreateInterfaceElements(&fPostProcMesh);
    fHybridizer.GroupElements(&fPostProcMesh);
    fPostProcMesh.CleanUpUnconnectedNodes();
//#ifdef PZDEBUG
//        {
//            std::ofstream out("multiphysicsgrouped.txt");
//            cmesh_Hybrid->Print(out);
//            std::ofstream outvtk("multiphysics.vtk");
//            TPZVTKGeoMesh::PrintCMeshVTK(cmesh_Hybrid,outvtk);
//            std::ofstream outgvtk("postprocessgmesh.vtk");
//            TPZVTKGeoMesh::PrintGMeshVTK(cmesh_Hybrid->Reference(),outgvtk);
//        }
//#endif
    
    
//#ifdef PZDEBUG
//    {
//        std::ofstream out("EnrichedFluxBorder.txt");
//        fPostProcMesh[1]->Print(out);
//    }
//#endif
    
    


}

/// computing the element stifnesses will "automatically" compute the condensed form of the matrices
void TPZHybridHDivErrorEstimator::ComputeElementStiffnesses()
{
    {
    std::ofstream out ("MeshToComputeStiff.txt");
    fPostProcMesh.Print(out);
    }
    
    TPZCompMesh *hybrid=&fPostProcMesh;

    for (auto cel:fPostProcMesh.ElementVec()) {
        if(!cel) continue;
            TPZElementMatrix ek, ef;
            cel->CalcStiff(ek, ef);
        
    }
    
    
}

/// increase the side orders of the post processing flux mesh
void TPZHybridHDivErrorEstimator::IncreaseSideOrders(TPZCompMesh *mesh)
{
    int64_t nel = mesh->NElements();
    int dim = mesh->Dimension();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = mesh->Element(el);
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
                TPZConnect &c = intel->SideConnect(0, side);
                if(c.Order() != order)
                {
                    intel->SetSideOrder(side, order);
                }
            }
        }
        //        intel->Print();
    }
    mesh->InitializeBlock();

}

void TPZHybridHDivErrorEstimator::IncreasePressureSideOrders(TPZCompMesh *mesh)
{
    /// esta ordem eh mal calculado!!!
    int OringOrder = mesh->GetDefaultOrder();
    
    int64_t nel = mesh->NElements();
    int dim = mesh->Dimension();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = mesh->Element(el);
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension()== dim){
            continue;
        }
        
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 1, 0);
        int nhighside=celstack.NElements();
        //TPZStack<STATE> vecOrder;
        int maxOrder=0;
        for (int side=0;side<nhighside;side++){
        TPZInterpolatedElement *intelS = dynamic_cast<TPZInterpolatedElement*>(celstack[side].Element());
            int orderEl=intelS->GetPreferredOrder();
            
            maxOrder = orderEl < maxOrder ? maxOrder : orderEl;
           // vecOrder.Push(orderEl);
        }

        
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
//        int nc = cel->NConnects();
//        int order = cel->Connect(nc-1).Order()+1;
        int order=maxOrder;//OringOrder;
        int nsides = gel->NSides();
        int ncorner = gel->NCornerNodes();
        intel->SetPreferredOrder(order);
        for (int side = ncorner; side < nsides; side++) {
            if (intel->NSideConnects(side)) {
                intel->SetSideOrder(side, maxOrder/*order*/);
            }
        }
        //        intel->Print();
    }
    mesh->InitializeBlock();
    
}

/// compute the average pressures of the hybridized form of the H(div) mesh
void TPZHybridHDivErrorEstimator::ComputeAveragePressures()
{
    TPZCompMesh *pressure = fOriginal->MeshVector()[1];
    TPZCompMesh *pressureHybrid = fPostProcMesh.MeshVector()[1];
    int fInterfaceMatid = fHybridizer.fLagrangeInterface;
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
        if (gel->MaterialId() != fInterfaceMatid) {
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
            celstack[0].Element()->Solution(pt1, 0, sol1);//solucao a esquerda
            celstack[1].Element()->Solution(pt2, 0, sol2);//solucao a direita
             //           std::cout << "Values " << sol1 << " " << sol2 << std::endl;
            //projecao L2 da media das soluceos no espaco Lh, do esqueleto da malha
            for (int ishape=0; ishape<nshape; ishape++) {
                L2Rhs(ishape,0) += weight*phi(ishape,0)*(sol1[0]+sol2[0])/2.;
                for (int jshape = 0; jshape<nshape; jshape++) {
                    L2Mat(ishape,jshape) += weight*phi(ishape,0)*phi(jshape,0);
                }
            }
        }
        L2Mat.SolveDirect(L2Rhs, ECholesky);
        //apos este passo temos uma pressao que é continua ao longo das interfaces dos elementos, nos esqueletos. Falta suavizar nos vértices
       // L2Rhs.Print("Average pressure");
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
    meshvec[0] = fPostProcMesh.MeshVector()[0];
    meshvec[1] = fPostProcMesh.MeshVector()[1];
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, &fPostProcMesh);
}

/// set the cornernode values equal to the averages
void TPZHybridHDivErrorEstimator::ComputeNodalAverages()
{
    TPZCompMesh *pressureHybrid = fPostProcMesh.MeshVector()[1];
    int fInterfaceMatid = fHybridizer.fLagrangeInterface;
    TPZGeoMesh *gmesh = pressureHybrid->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    pressureHybrid->LoadReferences();
    int lagrangematid = fHybridizer.fLagrangeInterface;
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
        if(gel->Dimension() != 1 || gel->MaterialId() != fInterfaceMatid)
        {
            continue;
        }
        //percorre cada no do elemento de inteface
        for (int side=0; side<2; side++) {
            TPZGeoElSide gelside(gel,side);
            TPZStack<TPZCompElSide> celstack;
            gelside.ConnectedCompElementList(celstack, 1, 0);
            celstack.Push(gelside.Reference());
            TPZManVector<STATE,3> averageval(nstate,0.);
            std::set<int64_t> connects;
            //para cada elemento que tem este no procede como segue
            for (int elc=0; elc<celstack.size(); elc++) {
                TPZCompElSide celside = celstack[elc];
                TPZGeoElSide gelside0 = celside.Reference();
                if (gelside0.Element()->Dimension() != 1) {
                    continue;
                }
                TPZInterpolatedElement *intel1 = dynamic_cast<TPZInterpolatedElement *>(celside.Element());
                if(!intel1 || intel1->NConnects() != 3) DebugStop();
                int64_t conindex = intel1->ConnectIndex(celside.Side());
                connects.insert(conindex);//insere os conects associado a este no
                TPZConnect &c = intel1->Connect(celside.Side());
                int64_t seqnum = c.SequenceNumber();
                if(c.NState() != nstate || c.NShape() != 1) DebugStop();
                //soma a soluçao dos elementos que possuem este no
                for (int istate = 0; istate<nstate; istate++) {
                    averageval[istate] += pressureHybrid->Block().Get(seqnum, 0, istate, 0);
                    
                }//somou todas as medias de pressao nas interfaces??
            }
            
            //por fim divide pelo numero de conects
            auto ncontr = connects.size();
            for (int istate = 0; istate<nstate; istate++) {
                averageval[istate] /= ncontr;
             //   std::cout<<"valor nos vertices "<<averageval[istate]<<std::endl;
            }//constroi o operador averaging
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
        if(gel->Dimension() != 1 || gel->MaterialId() != fInterfaceMatid)
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
    meshvec[0] = fPostProcMesh.MeshVector()[0];
    meshvec[1] = fPostProcMesh.MeshVector()[1];
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, &fPostProcMesh);
}



/// clone the meshes into the post processing mesh
void TPZHybridHDivErrorEstimator::CloneMeshVec()
{
//    CreateFluxMesh();
//    CreatePressureMesh();
    
    
    
    for (int i = 0; i < fOriginal->MeshVector().size(); i++) {
        fPostProcMesh.MeshVector()[i] = fOriginal->MeshVector()[i]->Clone();
    }
    
    

}
void TPZHybridHDivErrorEstimator::CreateMultiphysicsHybridMesh(){
    // the pressure mesh is the rootmesh
    TPZCompMesh *cmeshroot = fOriginal->MeshVector()[1];
    TPZGeoMesh *gmesh = cmeshroot->Reference();
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    
    //criando material
    int dim = gmesh->Dimension();
//    int nmat=cmeshroot->MaterialVec().size();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    for (auto it:cmeshroot->MaterialVec()) {
        TPZMaterial *mat = it.second;
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        if (bnd){
            TPZMaterial *matorig = bnd->Material();
            int matid = matorig->Id();
            TPZMaterial *matmixed = mphysics->FindMaterial(matid);
            TPZBndCond *bc = matmixed->CreateBC(matmixed, bnd->Id(), bnd->Type(), bnd->Val1(), bnd->Val2());
            if (bnd->ForcingFunction()) {
                bc->TPZDiscontinuousGalerkin::SetForcingFunction(bnd->ForcingFunction());
            }
            mphysics->InsertMaterialObject(bc);
        }
        else{
            int matId = mat->Id();
            int nstate = mat->NStateVariables();
            TPZMaterial *material = 0;
            if (nstate == 1  ) {
                TPZMixedPoisson *mix = dynamic_cast<TPZMixedPoisson *>(mat);
                if(mix){
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
                mphysics->InsertMaterialObject(material);
                }
                else{
                    
                    int lagrangematid = fHybridizer.fLagrangeInterface;
                    TPZMat1dLin *linear=new TPZMat1dLin(lagrangematid/*matId*/);//dynamic_cast<TPZMat1dLin *>(mat);
                    
                    TPZFNMatrix<1,STATE> xk(nstate,nstate,0.),xb(nstate,nstate,0.),xc(nstate,nstate,0.),xf(nstate,1,0.);
                   // material = locmat;
                   // mphysics->InsertMaterialObject(material);
                    mphysics->InsertMaterialObject(linear);
                    linear->SetMaterial(xk, xc, xb, xf);
                    
                    
                }
                
            }
            
        }
       
        
        
    }
    
    TPZManVector<TPZCompMesh *> meshvec(2,0);
    meshvec[0] = fPostProcMesh.MeshVector()[0];
    meshvec[1] = fPostProcMesh.MeshVector()[1];
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    //Fazendo auto build
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    mphysics->CleanUpUnconnectedNodes();
//#ifdef PZDEBUG
//    {
//        std::ofstream out("multiphysics.txt");
//        mphysics->Print(out);
//    }
//#endif
    
    //------- Create and add group elements -------
    DebugStop();

  
}
/// create the multiphysics mesh using TPZCompMeshReferred and TPZMixedErrorEstimate material
void TPZHybridHDivErrorEstimator::CreateMultiphysicsMesh()
{
    DebugStop();
    // the pressure mesh is the rootmesh
    TPZCompMesh *cmeshroot = fOriginal->MeshVector()[1];//fOriginal[0];
    
    TPZGeoMesh *gmesh = cmeshroot->Reference();
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();

    //por enquanto definir na mao os Id do esqueleto
    
    fHybridizer.fHDivWrapMatid=2;
    fHybridizer.fLagrangeInterface=3;
    fHybridizer.fInterfaceMatid=4;
    

 //   TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    
    if(fOriginalIsHybridized){
        
        //here we construct multiphysics mesh with active space as de postproces mesh
        TPZMultiphysicsCompMesh *mphysics = new TPZMultiphysicsCompMesh(gmesh);
        
        std::ofstream outgeo("geometria.txt");
        mphysics->Reference()->Print(outgeo);
        
        //Have to include the materials. Here we just did a copy of previous materials
        
       
        int dim = fOriginal->Reference()->Dimension();
        //
        for (auto it:fOriginal->MaterialVec()) {
            TPZMaterial *mat = it.second;
                int nstate = mat->NStateVariables();
                TPZMaterial *material = 0;
                int matId = mat->Id();
                if (nstate == 1  ) {
                    TPZMixedPoisson *mix = dynamic_cast<TPZMixedPoisson *>(mat);
                    if(mix){
                    TPZMixedHDivErrorEstimate<TPZMixedPoisson> *locmat = new TPZMixedHDivErrorEstimate<TPZMixedPoisson>(matId,dim);
                        material=locmat;
                        locmat->SetForcingFunction(mat->ForcingFunction());
                        if(fExact)
                        {
                            locmat->SetForcingFunctionExact(fExact->Exact());
                        }
                        TPZFNMatrix<9,REAL> PermTensor(3,3,0.);
                        TPZFNMatrix<9,REAL> InvPermTensor(3,3,0.);
                        PermTensor.Identity();
                        InvPermTensor.Identity();
                        locmat->SetPermeabilityTensor(PermTensor, InvPermTensor);
                        locmat->SetDimension(dim);
                        mphysics->InsertMaterialObject(material);
                        
                    }
                }
        }
        
 
        //
        
        fOriginal->CopyMaterials(*mphysics);
        
        TPZManVector<TPZCompMesh *,6> mp_meshes_vec(4);
        //flux and pressure reconstructed
        mp_meshes_vec[0] = fPostProcMesh.MeshVector()[0];
        mp_meshes_vec[1] = fPostProcMesh.MeshVector()[1];
        //flux and pressure original
        mp_meshes_vec[2] = fOriginal->MeshVector()[0];
        mp_meshes_vec[3] = fOriginal->MeshVector()[1];
        
        mphysics->SetDimModel(gmesh->Dimension());
        
        //the active space will be the flux an pressure reconstructed
        TPZManVector<int,5>  active_approx_spaces(4,0);
        active_approx_spaces[0]=1;
        active_approx_spaces[1]=1;
        
        //definir quais materias deverao estar no contibute...melhorar isso
 
        mphysics->BuildMultiphysicsSpace( active_approx_spaces, mp_meshes_vec);

        
//        TPZManVector<TPZCompMesh *> meshvec(2,0);
//        meshvec[0] = fPostProcMesh[1];
//        meshvec[1] = fPostProcMesh[2];
        
//        mphysics->SetAllCreateFunctionsMultiphysicElem();
////        mphysics->AutoBuild();
 //       mphysics->CleanUpUnconnectedNodes();
//
//        TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
//
//        TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
 //       TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
        
        
      //  fOriginal[0]->CopyMaterials(*mphysics);
        
    }
    else{
        TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
        
        fHybridizer.InsertPeriferalMaterialObjects(mphysics);

        //criando material
        int dim = gmesh->Dimension();
        {
            int n=fHybridizer.fLagrangeInterface;
            int n1=fHybridizer.fHDivWrapMatid;
            int n2=fHybridizer.fInterfaceMatid;
            std::cout<<"---PerifericalMaterialId --- "<<std::endl;
            std::cout <<" LagrangeInterface = "<<n<<std::endl;
            std::cout <<" HDivWrapMatid = "<<n1<<std::endl;
            std::cout <<" InterfaceMatid = "<<n2<<std::endl;
        }
        
        
        for (auto it:cmeshroot->MaterialVec()) {
            TPZMaterial *mat = it.second;
            TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
            if (!bnd) {
                int matId = mat->Id();
                std::cout<<" PeriferalMaterialId "<<fHybridizer.IsPeriferalMaterialId(matId)<<" Id "<<matId<<std::endl;
                
                if (fHybridizer.IsPeriferalMaterialId(matId)) {
                    continue;
                }
                int nstate = mat->NStateVariables();
                TPZMaterial *material = 0;
                if (nstate == 1  ) {
                    TPZMixedPoisson *mix = dynamic_cast<TPZMixedPoisson *>(mat);
                    if(!mix){
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
                    locmat->SetDimension(dim);
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
        meshvec[0] = fPostProcMesh.MeshVector()[0];
        meshvec[1] = fPostProcMesh.MeshVector()[1];
        
        mphysics->SetAllCreateFunctionsMultiphysicElem();
        //Fazendo auto build
        mphysics->AutoBuild();
        mphysics->CleanUpUnconnectedNodes();
        
        TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
        
        TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
        
        mphysics->CleanUpUnconnectedNodes();
//        fPostProcMesh[0] = mphysics;
    
    }
    
//    TPZManVector<TPZCompMesh *> meshvec(2,0);
//    meshvec[0] = fPostProcMesh[1];
//    meshvec[1] = fPostProcMesh[2];
//
//    mphysics->SetAllCreateFunctionsMultiphysicElem();
//    //Fazendo auto build
//    mphysics->AutoBuild();
//    mphysics->CleanUpUnconnectedNodes();
//
//    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
//
//    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
//    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
//
//    mphysics->CleanUpUnconnectedNodes();
//    fPostProcMesh[0] = mphysics;

}

void TPZHybridHDivErrorEstimator::NewCreateFluxMesh(){
    //malha original do fluxo
    DebugStop();
    TPZGeoMesh *gmesh = fOriginal->MeshVector()[0]->Reference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    int dim = gmesh->Dimension();
    TPZCompMesh *fluxoriginal=fOriginal->MeshVector()[0];
    
    for (auto it : fOriginal->MaterialVec()) {
        TPZMaterial *mat = it.second;
        int matid = mat->Id();
        TPZBndCond *bc=dynamic_cast<TPZBndCond *>(mat);
            if(bc){
                TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 1.);
                TPZBndCond *bc2 = mat->CreateBC(mat, matid, 0, val1, val2);
                cmesh->InsertMaterialObject(bc2);
            }
        
            else{
                TPZVecL2 *mix=dynamic_cast<TPZVecL2 *>(mat);
                if(mix){
                    TPZVecL2 *mix2 = new TPZVecL2(matid);
                    mix2->SetDimension(dim);
                    cmesh->InsertMaterialObject(mix2);
                }
              
            }

    }
 
   fluxoriginal->CopyMaterials(*cmesh);
    cmesh->SetDefaultOrder(cmesh->GetDefaultOrder());
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    
    
//    fPostProcMesh[1] = cmesh;
    
}

/// clone the fluxmesh but using TPZCompElReferred objects
void TPZHybridHDivErrorEstimator::CreateFluxMesh()
{
    DebugStop();
    TPZGeoMesh *gmesh = fOriginal->Reference();
    TPZCompMeshReferred *fluxmesh = new TPZCompMeshReferred(gmesh);
    TPZCompMesh *origflux = fOriginal->MeshVector()[0];
    origflux->CopyMaterials(*fluxmesh);
    gmesh->ResetReference();
    int meshdim = gmesh->Dimension();
    fluxmesh->ApproxSpace().SetAllCreateFunctionsHDivReferred(gmesh->Dimension());
    for (auto celorig:origflux->ElementVec()) {
        if(!celorig) continue;
        TPZGeoEl *gel = celorig->Reference();
        int64_t index;
      //  int geldim = gel->Dimension();
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
    
//    fPostProcMesh[1]=fluxmesh;
    
    
}


void TPZHybridHDivErrorEstimator::NewCreatePressureMesh(){
    DebugStop();
    TPZGeoMesh *gmesh = fOriginal->Reference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    TPZCompMesh *orinalpressure=fOriginal->MeshVector()[1];
    for (auto it : fOriginal->MeshVector()[1]->MaterialVec()) {
        TPZMaterial *mat=it.second;
        int maId = mat->Id();
        TPZMixedPoisson *mix = new TPZMixedPoisson(maId, cmesh->Dimension());
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    int porder=cmesh->GetDefaultOrder();
    cmesh->SetDefaultOrder(porder);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->AutoBuild();
    int64_t n_connects = cmesh->NConnects();
    for (int64_t i = 0; i < n_connects; ++i) {
        cmesh->ConnectVec()[i].SetLagrangeMultiplier(1);
    }
    
  //  orinalpressure->CopyMaterials(*cmesh);
    cmesh->InitializeBlock();
//    fPostProcMesh[2] = cmesh;
    
    
    
}
/// clone the pressure mesh using TPZCompElReferred objects
void TPZHybridHDivErrorEstimator::CreatePressureMesh()
{
    DebugStop();
    TPZGeoMesh *gmesh = fOriginal->Reference();
    TPZCompMeshReferred *pressuremesh = new TPZCompMeshReferred(gmesh);
    TPZCompMesh *origpressure = fOriginal->MeshVector()[1];
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
//    fPostProcMesh[2] =pressuremesh;
   
 
}

/// compute the effectivity indices of the pressure error and flux error and store in the element solution
void TPZHybridHDivErrorEstimator::ComputeEffectivityIndices()
{
    /**The  ElementSolution() is a matrix with 4 cols,
     col 0: pressure estimate error
     col 1: pressure exact error
     col 2: flux estimate error
     col 3: flux exact error
     Is increased 2 cols on ElementSolution() to store the efectivity index for pressure and flux
     **/
    TPZCompMesh *cmesh = &fPostProcMesh;
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
//                REAL ErrorEstimate=cmesh->ElementSolution()(el,i+1);
//                REAL ErrorExact=cmesh->ElementSolution()(el,i);
                cmesh->ElementSolution()(el,4+i/2) =cmesh->ElementSolution()(el,i+1)/cmesh->ElementSolution()(el,i);
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
    TPZMaterial *mat = fPostProcMesh.FindMaterial(matid);
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
    TPZMaterial *mat = fPostProcMesh.FindMaterial(matid);
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


void TPZHybridHDivErrorEstimator::PotentialReconstruction(){
    /// I havent tested to compute post processing more than a single time
    // Please test me!
    if (fPostProcMesh.MeshVector().size()) {
        DebugStop();
    }
    
    //Create the post processing mesh (hybridized H(div) mesh) with increased approximation order
    // for the border fluxes
    // in the future we can opt to create an H(1) post processing mesh
    CreatePostProcessingMesh();
    
    //calculando media das pressoes internas e valor nos vertices
    
    if(fProblemConfig.makepressurecontinuous)
    {
        ComputeAveragePressures();
    }
    
    ComputeNodalAverages();
    
    PlotLagrangeMultiplier("AfterNodalAverage");
#ifdef PZDEBUG
    {
        std::ofstream out("MeshWhithSmoothPressure.txt");
        fPostProcMesh.Print(out);
        std::ofstream out2("PressureMeshSmooth.txt");
        fPostProcMesh.MeshVector()[1]->Print(out2);
        std::ofstream out3("FluxMeshSmooth.txt");
        fPostProcMesh.MeshVector()[0]->Print(out3);
    }
#endif
    

    //Resolver problema local com potencial continuo como condicao de Dirichlet
    
    ComputeElementStiffnesses();
    
    fPostProcMesh.LoadSolution(fPostProcMesh.Solution());
    
#ifdef PZDEBUG
    {
        std::ofstream out("MeshAposLoadSol.txt");
        fPostProcMesh.Print(out);

    }
#endif
 
    
    TPZManVector<TPZCompMesh *,2> meshvec(2);
    // fPostProcMesh[1] is the flux mesh
    // fPostProcMesh[2] is the pressure mesh
    meshvec[0] = fPostProcMesh.MeshVector()[0];
    meshvec[1] = fPostProcMesh.MeshVector()[1];
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, &fPostProcMesh);
#ifdef PZDEBUG
    {
        std::ofstream out("PressureAposLoadSol.txt");
        fPostProcMesh.MeshVector()[1]->Print(out);
        std::ofstream out2("FluxAposLoadSol.txt");
        fPostProcMesh.MeshVector()[0]->Print(out2);

    }
#endif
    

}

void TPZHybridHDivErrorEstimator::PlotLagrangeMultiplier(const std::string &filename, bool reconstructed)
{
    TPZCompMesh *pressure = fPostProcMesh.MeshVector()[1];
    
    if(reconstructed == false) pressure = fOriginal->MeshVector()[1];
    
    TPZAnalysis an(pressure,false);
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("State");
    
    int dim = pressure->Reference()->Dimension()-1;
    std::string plotname;
    {
        std::stringstream out;
        out << filename << ".vtk";
        plotname = out.str();
    }
    an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
    an.PostProcess(2,dim);
    

}

/// identify the peripheral material objects and store the information in fHybridizer
void TPZHybridHDivErrorEstimator::IdentifyPeripheralMaterialIds()
{
    int dim = fOriginal->Dimension();
    // identify the material id for interface elements
    int64_t nel = fOriginal->NElements();
    for(int64_t el = 0; el < nel; el++)
    {
        TPZCompEl *cel = fOriginal->Element(el);
        TPZMultiphysicsInterfaceElement *interface = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
        if (interface) {
            int matid = interface->Reference()->MaterialId();
            fHybridizer.fInterfaceMatid = matid;
            break;
        }
    }
    /// identify the material id of the pressure
    TPZCompMesh *pressure_mesh = fOriginal->MeshVector()[1];
    nel = pressure_mesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = pressure_mesh->Element(el);
        if(cel)
        {
            TPZGeoEl *gel = cel->Reference();
            if(gel && gel->Dimension() == dim-1)
            {
                fHybridizer.fLagrangeInterface = gel->MaterialId();
                break;
            }
        }
    }
    /// identify the material id of boundary elements
    TPZCompMesh *fluxmesh = fOriginal->MeshVector()[0];
    for (auto map_pair : fluxmesh->MaterialVec()) {
        TPZMaterial *mat = map_pair.second;
        TPZNullMaterial *nullmat = dynamic_cast<TPZNullMaterial *>(mat);
        if(nullmat)
        {
            fHybridizer.fHDivWrapMatid = nullmat->Id();
            break;
        }
    }
}

/// switch material object from mixed poisson to TPZMixedHdivErrorEstimate
void TPZHybridHDivErrorEstimator::SwitchMaterialObjects()
{
    for(auto matid : fPostProcMesh.MaterialVec())
    {
        TPZMixedPoisson *mixpoisson = dynamic_cast<TPZMixedPoisson *> (matid.second);
        if(mixpoisson)
        {
            TPZMixedHDivErrorEstimate<TPZMixedPoisson> *newmat = new TPZMixedHDivErrorEstimate<TPZMixedPoisson>(*mixpoisson);
            
            if(fExact)
            {
                newmat->SetForcingFunctionExact(fExact->Exact());
            }

            for (auto bcmat : fPostProcMesh.MaterialVec()) {
                TPZBndCond *bc = dynamic_cast<TPZBndCond *>(bcmat.second);
                if (bc) {
                    bc->SetMaterial(newmat);
                }
            }
            fPostProcMesh.MaterialVec()[newmat->Id()] = newmat;
            delete mixpoisson;
        }
    }
}

