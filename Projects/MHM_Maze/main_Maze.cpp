#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Maze_common.h"
int main(){
    TPZLogger::InitializePZLOG();
    
    ConfigCasesMaze ConfCasesMeze;
    ConfCasesMeze.SetImageName("⁨../Mazes/maze128x128.png");
    ConfCasesMeze.SetImperviousMatPermeability(1);//pouco permeavel
    ConfCasesMeze.SetPermeableMatPermeability(1000000);//dentro do labirinto
    ConfCasesMeze.SetFluxOrder(1);
    ConfCasesMeze.SetPressureOrder(1);
    ConfCasesMeze.SetCCPressureIn(100);//pressao na entrada
    ConfCasesMeze.SetCCPressureOut(1);//pressao na saida
    ConfCasesMeze.SetMHMOpenChannel(false);
    ConfCasesMeze.SetVTKName("maze128x128.vtk");

    MHMTest(ConfCasesMeze);
    return 0;
}


int MHMTest(ConfigCasesMaze &Conf){

    TRunConfig Configuration;

    TPZGeoMesh *gmeshcoarse = GenerateGeoMesh(Conf.GetImageName(), 32, 32);
    {
        std::ofstream file(Conf.GetVTKName());
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshcoarse, file);
    }

    int interface_mat_id = 600;
    Conf.SetMHMOpenChannel(true);
    bool OpenChannel = Conf.GetMHMOpenChannel();

    TPZAutoPointer<TPZMHMixedMeshChannelControl> MHMixed;

    {
        TPZAutoPointer<TPZGeoMesh> gmeshauto = new TPZGeoMesh(*gmeshcoarse);
        {
            std::ofstream out("gmeshauto.txt");
            gmeshauto->Print(out);
        }
        TPZMHMixedMeshChannelControl *mhm = new TPZMHMixedMeshChannelControl(gmeshauto);
        TPZVec<int64_t> coarseindices;
        ComputeCoarseIndices(gmeshauto.operator->(), coarseindices);
        gmeshauto->AddInterfaceMaterial(1, 2, interface_mat_id);
        gmeshauto->AddInterfaceMaterial(2, 1, interface_mat_id);


        // criam-se apenas elementos geometricos
        mhm->DefinePartitionbyCoarseIndices(coarseindices);
        MHMixed = mhm;

        TPZMHMixedMeshChannelControl &meshcontrol = *mhm;
        {
            std::set<int> matids;
            matids.insert(1);
            matids.insert(2);
            mhm->fMaterialIds = matids;
            matids.clear();
            matids.insert(-1);
            matids.insert(-2);
            matids.insert(-3);
            matids.insert(-4);
            matids.insert(-5);
            matids.insert(-6);
            mhm->fMaterialBCIds = matids;
        }

        InsertMaterialObjects(*mhm);

        meshcontrol.SetInternalPOrder(1);
        meshcontrol.SetSkeletonPOrder(1);

//        meshcontrol.DivideSkeletonElements(2);
        meshcontrol.DivideBoundarySkeletonElements();

        bool substructure = true;
        std::map<int, std::pair<TPZGeoElSide, TPZGeoElSide>> test;
        if (OpenChannel) {
            TPZCompMesh *flux_temp = MixedTest(Conf,32,32);
            std::cout << "flux_temp norm of solution " << Norm(flux_temp->Solution()) << std::endl;
            test = IdentifyChanel(flux_temp);
            flux_temp->Reference()->ResetReference();
//            delete flux_temp;
        }

        meshcontrol.BuildComputationalMesh(substructure, OpenChannel, test);

#ifdef ERRORESTIMATION_DEBUG
        if (1) {
            std::ofstream file("GMeshControlHDiv.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(meshcontrol.GMesh().operator->(), file);
        }
#endif

        std::cout << "MHM Hdiv Computational meshes created\n";

        std::cout << "Number of equations MHMixed " << MHMixed->CMesh()->NEquations() << std::endl;

    }
    
    TPZCompMesh *MixedMesh = MHMixed->CMesh().operator->();
    

    SolveProblem(MHMixed->CMesh(), MHMixed->GetMeshes(), Conf.GetExactSolution(),  Conf.GetVTKName(), Configuration);
    ProblemConfig config;
    config.dimension = 2;
    config.exact = nullptr;
    config.problemname = "MazeHdiv128x128";
    config.dir_name = "Results128x128";
    config.porder = 3;
    config.hdivmais = 3;
    config.materialids = {1, 2};
    config.bcmaterialids = {-1, -2, -3, -4, -5, -6};
    config.makepressurecontinuous = true;
    config.ndivisions = 0;
    config.gmesh = MixedMesh->Reference();

    TPZMultiphysicsCompMesh *originalMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(MHMixed->CMesh().operator->());
    bool postProcWithHdiv = false;
//    TPZDarcyMHMHDivErrorEstimator ErrorEstimator(*originalMesh, MHMixed.operator->(), postProcWithHdiv);
//    EstimateError(ErrorEstimator, config);
    //LocateElementsToAdapt(ErrorEstimator, config);

    return 0;
}

