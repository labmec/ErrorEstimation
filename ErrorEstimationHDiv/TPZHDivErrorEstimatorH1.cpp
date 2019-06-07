//
// Created by gustavo on 30/05/19.
//

#include "TPZHDivErrorEstimatorH1.h"

/// create the post processed multiphysics mesh (which is necessarily hybridized)
void TPZHDivErrorEstimatorH1::CreatePostProcessingMesh()
{
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
    mesh_vectors[0] = fOriginal->MeshVector()[1]->Clone();
    mesh_vectors[1] = CreateDiscontinuousMesh(mesh_vectors[0]);
    
    if(!fOriginalIsHybridized)
    {
        // if the original mesh is not hybridized, it would be much more complicated to generate the
        // pressure interface elements etc
        DebugStop();
    }
    IdentifyPeripheralMaterialIds();
    int lastmatid = fPostProcMesh.MaterialVec().rbegin()->first;
    fSkeletonMatId = lastmatid+1;
    // increase the order of the dim-1 elements to the maximum of both neighbouring elements
    IncreasePressureSideOrders(mesh_vectors[0]);//malha da pressao

    if(dim == 3)
    {
        CreateEdgeSkeletonMesh(mesh_vectors[0]);
    }
#ifdef PZDEBUG
    {
        std::ofstream out("EnrichedFluxBorder.txt");
        mesh_vectors[0]->Print(out);
        
        std::ofstream out2("EnrichedPressure.txt");
        mesh_vectors[1]->Print(out2);
    }
#endif
    
    
    
    TPZManVector<int> active(4,0);
    active[0] = 1;
    active[1] = 1;
    fPostProcMesh.BuildMultiphysicsSpace(active, mesh_vectors);
    {
        std::ofstream out("multiphysicsWithnoInterface.txt");
        fPostProcMesh.Print(out);
    }

    DebugStop();
}

/// compute a more precise approximation for the pressure
void TPZHDivErrorEstimatorH1::UpliftPressure()
{
    // fPostProcMesh has been created
    // for each element of dimension dim - compute the stiffness matrix, invert it and load the solution
}

/// create a constant pressure mesh used for uplifting the pressure
TPZCompMesh *TPZHDivErrorEstimatorH1::CreateDiscontinuousMesh(const TPZCompMesh *pressuremesh)
{
    DebugStop();
}

/// switch material object from mixed poisson to TPZMixedHdivErrorEstimate
void TPZHDivErrorEstimatorH1::SwitchMaterialObjects()
{
    // switch the material of the HDiv approximation to a material for an H1 approximation
    DebugStop();
}

