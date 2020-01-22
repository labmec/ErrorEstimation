//
// Created by gustavo on 30/05/19.
//

#include "TPZHDivErrorEstimatorH1.h"
#include "mixedpoisson.h"
#include "pzbndcond.h"
#include "TPZHDivErrorEstimateMaterial.h"
#include "TPZVTKGeoMesh.h"
#include "pzintel.h"
#include "pzcondensedcompel.h"
#include "pzbuildmultiphysicsmesh.h"

/// create the post processed multiphysics mesh
void TPZHDivErrorEstimatorH1::CreatePostProcessingMesh()
{
    
    std::cout<<"Creating post processed meshes"<<std::endl;
    // initialize the post processing mesh
    fPostProcMesh.SetReference(fOriginal->Reference());
    int dim = fOriginal->Dimension();
    fOriginal->CopyMaterials(fPostProcMesh);
    // switch the material from mixed to TPZMixedHdivErrorEstimate...
    SwitchMaterialObjects();
    
    TPZManVector<TPZCompMesh *,4> mesh_vectors(4,0);
    mesh_vectors[2] = fOriginal->MeshVector()[0];//flux
    mesh_vectors[3] = fOriginal->MeshVector()[1];//potential
    // create a copy of the pressure mesh
    mesh_vectors[0] = fOriginal->MeshVector()[1]->Clone();//H1 mesh
    mesh_vectors[1] = CreateDiscontinuousMesh(mesh_vectors[0]);//L2 mesh
    
    if(!fOriginalIsHybridized)
    {
        // if the original mesh is not hybridized, it would be much more complicated to generate the
        // pressure interface elements etc
        DebugStop();
    }
    
    IncreasePressureOrders(mesh_vectors[0]);
    IdentifyPeripheralMaterialIds();
   
    // increase the order of the dim-1 elements to the maximum of both neighbouring elements
    IncreasePressureSideOrders(mesh_vectors[0]);//malha da pressao

    if(dim == 3)
    {
        CreateEdgeSkeletonMesh(mesh_vectors[0]);
    }
#ifdef PZDEBUG2
    {
        std::ofstream out("EnrichedPressure.txt");
        mesh_vectors[0]->Print(out);
        
        std::ofstream out2("DiscontinuousMesh.txt");
        mesh_vectors[1]->Print(out2);
    }
#endif
    
    
    
    TPZManVector<int> active(4,0);
    active[0] = 1;
    active[1] = 1;
    fPostProcMesh.BuildMultiphysicsSpace(active, mesh_vectors);
//    {
//        std::ofstream out("multiphysicsWithnoInterface.txt");
//        fPostProcMesh.Print(out);
//    }
    
    if(fperformUplift){
        
        std::cout<<"Solving local Neumann problem"<<std::endl;
        // compute a higher order pressure solution: compute the Local Neumann problem
        UpliftPressure();
        {
            TPZAnalysis an(&fPostProcMesh,false);
            
            TPZStack<std::string> scalnames, vecnames;
            scalnames.Push("PressureReconstructed");
            
            int dim = 2;
            std::string plotname("PotentialUplifting.vtk");
            an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
            an.PostProcess(2, dim);
            
            
            
        }
        
    }
    
    

    

    // transfer the continuous pressures to the multiphysics space
    {
        TPZManVector<TPZCompMesh *,2> meshvec(2);
        meshvec[0] = fPostProcMesh.MeshVector()[0];
        meshvec[1] = fPostProcMesh.MeshVector()[1];
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, &fPostProcMesh);
    }

    
    //post processing for local problem
    {
#ifdef PZDEBUG
        {
            std::ofstream out("MeshPosNeumann.txt");
            fPostProcMesh.Print(out);
            
        }
#endif
        
        TPZAnalysis an(fPostProcMesh.MeshVector()[0],false);

        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("State");

        int dim = 2;
        std::string plotname("LocalNeumannProblem.vtk");
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(2, dim);
        
        
    }
    
    PreparePostProcElements();
    

}



/// return a pointer to the pressure mesh
TPZCompMesh *TPZHDivErrorEstimatorH1::PressureMesh()
{
    return fPostProcMesh.MeshVector()[0];
}

/// prepare the elements of postprocmesh to compute the pressures with increased accuracy
void TPZHDivErrorEstimatorH1::PreparePostProcElements()
{
    
    /// we increase the connects of the borders of the pressure elements so that
    // the condensed compel does not condense these equations
    // then if we compute the stiffness matrix and load the solution the internal solution
    // is updated
    fPostProcMesh.ComputeNodElCon();
    int64_t nel = fPostProcMesh.NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fPostProcMesh.Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        if(gel->Dimension() != fPostProcMesh.Dimension()) continue;
        TPZMultiphysicsElement *mcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if(!mcel) DebugStop();
        TPZCompEl *subcel = mcel->Element(0);
        TPZInterpolatedElement *subintel = dynamic_cast<TPZInterpolatedElement *>(subcel);
        if(!subintel) DebugStop();
        int nsides = gel->NSides();
        for (int side = 0; side<nsides-1; side++) {
            if(subintel->NSideConnects(side) == 0) continue;
            int connectindex = subintel->MidSideConnectLocId(side);
            TPZConnect &connect = cel->Connect(connectindex);
            connect.IncrementElConnected();
        }
    }
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fPostProcMesh.Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        if(gel->Dimension() != fPostProcMesh.Dimension()) continue;
        TPZCondensedCompEl *condense = new TPZCondensedCompEl(cel,false);
    }
    for (auto matit : fPostProcMesh.MaterialVec()) {
        TPZMaterial *mat = matit.second;
        TPZHDivErrorEstimateMaterial *errormat = dynamic_cast<TPZHDivErrorEstimateMaterial *>(mat);
        if(errormat)
        {
            errormat->fNeumannLocalProblem = false;
        }
    }
    fPostProcMesh.CleanUpUnconnectedNodes();
}


/// compute a more precise approximation for the pressure
void TPZHDivErrorEstimatorH1::IncreasePressureOrders(TPZCompMesh *pressuremesh)
{
//    Increase the order of the pressure elements
    int64_t nel = pressuremesh->NElements();
    pressuremesh->Reference()->ResetReference();
    for(int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = pressuremesh->Element(el);
        if(!cel ||cel->Reference()->Dimension() != pressuremesh->Dimension()) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if(!intel) DebugStop();
        int nc = intel->NConnects();
        int order = intel->Connect(nc-1).Order();
        intel->PRefine(order+this->fUpliftOrder);
    }
    pressuremesh->ExpandSolution();
}

/// create a constant pressure mesh used for uplifting the pressure
TPZCompMesh *TPZHDivErrorEstimatorH1::CreateDiscontinuousMesh(const TPZCompMesh *pressuremesh)
{
    TPZCompMesh *cmesh = new TPZCompMesh(pressuremesh->Reference());
    for (auto matptr: pressuremesh->MaterialVec()) {
        if(matptr.second->Dimension() == pressuremesh->Dimension())
        {
            TPZMaterial *newmat = matptr.second->NewMaterial();
            cmesh->InsertMaterialObject(newmat);
        }
    }
    pressuremesh->Reference()->ResetReference();
    int64_t nel = pressuremesh->NElements();
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->SetDefaultOrder(0);//space of constant functions
    for (int64_t el=0; el<nel; el++) {
        const TPZCompEl *cel = pressuremesh->Element(el);
        if(!cel || cel->Reference()->Dimension() != pressuremesh->Dimension())
        {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        int64_t index;
        cmesh->CreateCompEl(gel, index);
    }
    cmesh->ExpandSolution();
    cmesh->ComputeNodElCon();
    return cmesh;
}

/// switch material object from mixed poisson to TPZMixedHdivErrorEstimate
void TPZHDivErrorEstimatorH1::SwitchMaterialObjects()
{
    // switch the material of the HDiv approximation to a material for an H1 approximation
    for(auto matid : fPostProcMesh.MaterialVec())
    {
        TPZMixedPoisson *mixpoisson = dynamic_cast<TPZMixedPoisson *> (matid.second);
        if(mixpoisson)
        {
            int dim = mixpoisson->Dimension();
            int matid = mixpoisson->Id();
            
            TPZHDivErrorEstimateMaterial *newmat = new TPZHDivErrorEstimateMaterial(*mixpoisson);
            
            if(fExact)
            {
                newmat->SetForcingFunction(fExact->Exact());
                newmat->SetForcingFunction(fExact->ForcingFunction());
    
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

/// Compute an uplifted solution for the pressure
void TPZHDivErrorEstimatorH1::UpliftPressure()
{
    if(fUpliftOrder < 0) return;
    fPostProcMesh.ComputeNodElCon();
    int64_t nel = fPostProcMesh.NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fPostProcMesh.Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        if(gel->Dimension() != fPostProcMesh.Dimension()) continue;
        TPZCondensedCompEl *cond = new TPZCondensedCompEl(cel,false);
#ifdef PZDEBUG
        int nc = cond->NConnects();
        for (int ic=0; ic<nc; ic++) {
            if(cond->Connect(ic).IsCondensed() == false)
            {
                DebugStop();
            }
        }
#endif
    }
    fPostProcMesh.ExpandSolution();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fPostProcMesh.Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        if(gel->Dimension() != fPostProcMesh.Dimension()) continue;
        TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>(cel);
        if(!cond)
        {
            DebugStop();
        }
        TPZElementMatrix ek,ef;
        //stifness matrix for local Neumann problem
        cond->CalcStiff(ek,ef);
        //solve the local Neumann problem
        cond->LoadSolution();
    }
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fPostProcMesh.Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        if(gel->Dimension() != fPostProcMesh.Dimension()) continue;
        TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>(cel);
        if(!cond)
        {
            DebugStop();
        }
        cond->Unwrap();
    }
}

/// copy the solution from the neighbouring skeleton elements
void TPZHDivErrorEstimatorH1::CopySolutionFromSkeleton()
{
    TPZCompMesh *pressuremesh = PressureMesh();
    {
        std::ofstream out("pressuremesh.txt");
        pressuremesh->Print(out);
    }
    pressuremesh->Reference()->ResetReference();
    pressuremesh->LoadReferences();
    int dim = pressuremesh->Dimension();
    int64_t nel = pressuremesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = pressuremesh->Element(el);
        if(!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if(!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if(gel->Dimension() != dim) continue;
        int nsides = gel->NSides();
        for (int is=0; is<nsides; is++) {
            TPZGeoElSide gelside(gel,is);
            if(IsDirichletCondition(gelside)){
                std::cout<<"bc side em copyskeleton \n";
                continue;
                
            }
            
            
            TPZConnect &c = intel->Connect(is);
            int64_t c_seqnum = c.SequenceNumber();
            int c_blocksize = c.NShape()*c.NState();
            //TPZGeoElSide gelside(gel,is);
            TPZStack<TPZCompElSide> celstack;
            gelside.EqualLevelCompElementList(celstack, 1, 0);
            int nst = celstack.NElements();
            for(int ist=0; ist<nst; ist++)
            {
                TPZCompElSide cneigh = celstack[ist];
                TPZGeoElSide gneigh = cneigh.Reference();
                if (gneigh.Element()->MaterialId() == this->fHybridizer.fLagrangeInterface) {
                    TPZInterpolatedElement *intelneigh = dynamic_cast<TPZInterpolatedElement *>(cneigh.Element());
                    if(!intelneigh) DebugStop();
                    TPZConnect &con_neigh = intelneigh->Connect(cneigh.Side());
                    int64_t con_seqnum = con_neigh.SequenceNumber();
                    int con_size = con_neigh.NState()*con_neigh.NShape();
                    if(con_size != c_blocksize) DebugStop();
                    for (int ibl = 0; ibl<con_size; ibl++) {
                        pressuremesh->Block()(c_seqnum,0,ibl,0) = pressuremesh->Block()(con_seqnum,0,ibl,0);
                    }
                    break;
                }
                // all elements must have at least one neighbour of type skeleton
                if(ist == nst-1)
                {
                    std::cout << "Connect " << is << " from element el " << el << " was not updated \n";
                }
            }
        }
    }
	{
		std::ofstream out("pressuremesh2.txt");
		pressuremesh->Print(out);
	}
}

