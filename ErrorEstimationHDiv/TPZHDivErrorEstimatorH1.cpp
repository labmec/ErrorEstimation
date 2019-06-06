//
// Created by gustavo on 30/05/19.
//

#include "TPZHDivErrorEstimatorH1.h"

#include "TPZVTKGeoMesh.h"

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
    mesh_vectors[0] = mesh_vectors[3]->Clone();//potential reconstructed

    mesh_vectors[1] = mesh_vectors[3]->Clone();// discontinous space
    mesh_vectors[1]->SetDefaultOrder(0); // Constant function
    DebugStop();

    if (!fOriginalIsHybridized)
    {
        fHybridizer.ComputePeriferalMaterialIds(mesh_vectors);
        fHybridizer.ComputeNState(mesh_vectors);
        fHybridizer.HybridizeInternalSides(mesh_vectors);
        int lastmatid = fPostProcMesh.MaterialVec().rbegin()->first;
        fSkeletonMatId = lastmatid+1;
#ifdef PZDEBUG
        {
            std::ofstream out("OriginalFlux.txt");
            fOriginal[3].Print(out);
            std::ofstream out2("OriginalPotential.txt");
            fOriginal[4].Print(out2);
            std::ofstream out3("OriginalMeshHybrid.txt");
            fOriginal[0].Print(out3);
        }
#endif
    }
    else
    {
        IdentifyPeripheralMaterialIds();
        int lastmatid = fPostProcMesh.MaterialVec().rbegin()->first;
        fSkeletonMatId = lastmatid+1;
    }

    // increase the order of the dim-1 elements to the maximum of both neighbouring elements
    // malha da pressao
    IncreasePressureSideOrders(mesh_vectors[0]);//malha da pressao
    // TODO isso é o uplift pressure?

    //IncreaseSideOrders(mesh_vectors[0]);//malha do fluxo

    if(dim == 3)
    {
        CreateEdgeSkeletonMesh(mesh_vectors[1]);
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

    // construction of the multiphysics mesh
    //cria elementos de interface
    fHybridizer.CreateInterfaceElements(&fPostProcMesh);
    fHybridizer.GroupElements(&fPostProcMesh);
    fPostProcMesh.CleanUpUnconnectedNodes();

#ifdef PZDEBUG
    {
        std::ofstream out("multiphysicsgrouped.txt");
        fPostProcMesh.Print(out);
//            std::ofstream outvtk("multiphysics.vtk");
//            TPZVTKGeoMesh::PrintCMeshVTK(cmesh_Hybrid,outvtk);
        std::ofstream outgvtk("postprocessgmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fPostProcMesh.Reference(),outgvtk);
    }
#endif


}

/// compute a more precise approximation for the pressure
void TPZHDivErrorEstimatorH1::UpliftPressure()
{

}
