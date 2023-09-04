//
//  ToolsMHM.cpp
//  AdaptivityTest
//
//  Created by Denise De Siqueira on 01/11/19.
//

#include "ToolsMHM.h"
#include <iostream>
#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzreal.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoelside.h"
#include "tpzgeoelrefpattern.h"
#include "tpzautopointer.h"
#include "TPZLinearAnalysis.h"

#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZStructMatrixT.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"


#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"

#include "TPZVTKGeoMesh.h"

#include "TPZHybridizeHDiv.h"
#include "TPZGenGrid2D.h"

#include "TPZNullMaterial.h"

#include <string>
#include <cmath>
#include <set>
#include <TPZMFSolutionTransfer.h>
#include "TPZPersistenceManager.h"

#include "TPZMHMHDivErrorEstimator.h"

TPZCompMesh *CMeshFlux(TPZGeoMesh * gmesh,int pOrder){
    
    int impervious_mat = 1;
    int permeable_mat = 2;
    int dim = gmesh->Dimension();
    
    TPZVec<REAL> convdir(dim,0.0);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    
    //Definition of the approximation space:
    
    TPZNullMaterial<> *mat_0 = new TPZNullMaterial<>(impervious_mat,dim);
    TPZNullMaterial<> *mat_1 = new TPZNullMaterial<>(permeable_mat,dim);
    
    //  inserting volumetric materials objects
    cmesh->InsertMaterialObject(mat_0);
    cmesh->InsertMaterialObject(mat_1);
    
    cmesh->SetAllCreateFunctionsHDiv(); //Creating H(div) functions
    
    
    int type_D = 0;
    int type_N = 1;
    TPZFMatrix<STATE> val1(1, 1, 0.);
    TPZManVector<REAL,1> val2(1, 0.);
    
    // Insert boundary conditions
    //Neumann boundary conditions (flux = 0)
    int right_bc_id = -2;
    TPZBndCondT<STATE> * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_N, val1, val2);
    cmesh->InsertMaterialObject(right_bc);
    
    int left_bc_id = -4;
    TPZBndCondT<STATE> * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_N, val1, val2);
    cmesh->InsertMaterialObject(left_bc);
    
    int bottom_bc_1id = -1;
    TPZBndCondT<STATE>* bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_N, val1, val2);
    cmesh->InsertMaterialObject(bottom_bc_1);
    
    int top_bc_1id = -3;
    TPZBndCondT<STATE> * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_N, val1, val2);
    cmesh->InsertMaterialObject(top_bc_1);
    
    
    //Dirichlet Conditions (p=1 in, p=0 out)
    int bottom_bc_id = -5;
    TPZBndCondT<STATE> * bottom_bc = mat_0->CreateBC(mat_0, bottom_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(bottom_bc);
    
    int top_bc_id = -6;
    TPZBndCondT<STATE> * top_bc = mat_0->CreateBC(mat_0, top_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(top_bc);
    
    cmesh->SetName("LaberintoTest");
    cmesh->AutoBuild();
    
#ifdef ERRORESTIMATION_DEBUG
    std::ofstream file("cmesh_flux.txt");
    cmesh->Print(file);
#endif
    
    return cmesh;
    
}

// compute the coarse indices of the geometric mesh
void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices)
{
    //    {
    //        std::ofstream out("gmeshref.txt");
    //        gmesh->Print(out);
    //    }
    coarseindices.Resize(gmesh->NElements());
    int count = 0;
    for (int64_t el=0; el<gmesh->NElements(); el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel || gel->Dimension() != gmesh->Dimension()) continue;
        if(gel->Father()) continue;
        coarseindices[count] = el;
        count++;
    }
    coarseindices.Resize(count);
}


std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>> IdentifyChanel(TPZCompMesh *cmesh){
    
    {
        std::ofstream out("fluxmesh.txt");
        cmesh->Print(out);
    }
    cmesh->LoadReferences();
    int nelements = cmesh->NElements();
    TPZGeoElSide first_gelside;
    for(int iel=0; iel<=nelements; iel++){
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel){continue;}
        TPZGeoEl *gel = cel->Reference();
        if(gel->MaterialId() == -5)
        {
            TPZGeoElSide gelside(gel);
            first_gelside = gelside.Neighbour();
            break;
        }
    }
    if(!first_gelside) DebugStop();
    int count=0;
    int count_chain=0;
    double Flux_Max = 0.0;
    std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>> chain;
    bool exit=false;
    
    while(exit==false){
        
        // if their is a neighbour along side 6 with material id -6, we found the exit
        TPZGeoElSide exit_test = first_gelside;
        exit_test.SetSide(6);
        TPZGeoElSide candidate_exist =exit_test.Neighbour();
        if(candidate_exist && candidate_exist.Element()->MaterialId()==-6){
            std::cout<<"Cadena encontrada con exito";
            exit = true;
            break;
        }
        
        
        //    while(candidate.Element()->Dimension()!=2){
        //
        //            candidate=candidate.Neighbour();
        //
        //    }
        int side_in = first_gelside.Side();
        // loop over the one dimensional sides
        for(int ican=4; ican<8; ican++){
            // if the one-d side is the entry side, do not consider
            if(side_in == ican ){continue;}

            first_gelside.SetSide(ican);
            // find a neighbour of dimension 2
            TPZGeoElSide candidate = first_gelside.Neighbour();
            while(candidate.Element()->Dimension()!=2){
                
                candidate=candidate.Neighbour();
            }
            if(candidate.Element() == first_gelside.Element()) DebugStop();
            
            //compute the flux values at the center of the 3 neighbours
            TPZVec<REAL> qsi(2);
            qsi[0]=0;
            qsi[1]=0;
            int var=1;
            TPZVec<STATE> sol;
            TPZCompEl *cel = candidate.Element()->Reference();
            cel->Solution(qsi, var,sol);
            double Flux_can_Mag = sqrt(sol[0]*sol[0] + sol[1]*sol[1]);
            
            if(Flux_can_Mag > Flux_Max){
                Flux_Max =Flux_can_Mag;
                first_gelside.SetSide(ican);
                // first is the exit side
                // second is the entry side
                chain[count].first= first_gelside;
                chain[count].second =candidate;
            }
            if(Flux_Max == 0.)
            {
                std::cout << "Mesh solution norm " << Norm(cmesh->Solution()) << std::endl;
                int nc = cel->NConnects();
                for (int ic=0; ic<nc; ic++) {
                    TPZConnect &c = cel->Connect(ic);
                    c.Print(*cmesh);
                }
                DebugStop();
            }
        }
        Flux_Max=0.0;
        first_gelside =chain[count].second;
        if(!first_gelside) DebugStop();
        count++;
    }
    
    //here
    
    cmesh->LoadReferences();
    TPZGeoMesh *gmesh =cmesh->Reference();
//    for(auto it:gmesh->ElementVec()){
//        int father_index = it->FatherIndex();
//        int element = it->Index();
//        std::cout<<"Element: "<<element<<" father: "<<father_index<<std::endl;
//
//    }
    
    
    int n_el_chain = chain.size();
    // same type as chain
    std::map<int, std::pair<TPZGeoElSide, TPZGeoElSide>> skelchanel;
    // the skel channel will filter the elements of chain that
    // have different father index
    // create a boundary element of matid 10 on the interfaces between macro
    // elements
    int count_skel_chanel=0;
    int matId_skel_chanel = 10;
    count =0;
    for(auto it : chain){
        TPZGeoEl *first_element = it.second.first.Element();
        TPZGeoEl *second_element = it.second.second.Element();
        
        int first_father_index = first_element->LowestFather()->Index();
        int second_father_index = second_element->LowestFather()->Index();
        
        if(0)
        {
            TPZManVector<REAL> xic1(2),xic2(2),xc1(3),xc2(3);
            first_element->CenterPoint(8, xic1);
            first_element->X(xic1, xc1);
            second_element->CenterPoint(8, xic2);
            second_element->X(xic2, xc2);
            std::cout << "fathers " << first_father_index << ' ' << second_father_index << " x1 " << xc1 << " x2 " << xc2 << std::endl;
        }
        
        if(first_father_index!=second_father_index){
            TPZGeoElBC(it.second.first, matId_skel_chanel);
            skelchanel[count_skel_chanel] = it.second;
            count_skel_chanel++;
        }
        
    }
    std::ofstream out("TestMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return skelchanel;
    
}

TPZGeoMesh *CreateCircleGeoMesh() {

    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);

    TPZVec<REAL> coord(3, 0.);

    // Inserts node at origin
    gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[0].Initialize(coord, *gmesh);

    // Inserts circumference nodes
    for (int64_t i = 0; i < 16; i++) {
        const REAL step = M_PI / 8;
        coord[0] = cos(i * step);
        coord[1] = sin(i * step);
        std::cout<<"{"<<coord[0]<<"," <<coord[1]<<"},"<<std::endl;



        const int64_t newID = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[newID].Initialize(coord, *gmesh);
    }

    int matIdTriangle = 1, matIdArc = 2;

    // Inserts triangle elements
    TPZManVector<int64_t, 3> nodesIdVec(3);
    for (int64_t i = 0; i < 7; i++) {
        nodesIdVec[0] = 0;
        nodesIdVec[1] = 1 + 2 * i;
        nodesIdVec[2] = 3 + 2 * i;
        std::cout<<"{"<<nodesIdVec[0] <<"," <<nodesIdVec[1] <<", "<<nodesIdVec[2] <<"},"<<std::endl;
        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle>>(nodesIdVec, matIdTriangle, *gmesh);
    }
    nodesIdVec[0] = 0;
    nodesIdVec[1] = 15;
    nodesIdVec[2] = 1;

    new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle>>(nodesIdVec, matIdTriangle, *gmesh);
    // Inserts arc elements
    for (int64_t i = 0; i < 7; i++) {
        nodesIdVec[0] = 1 + 2 * i;
        nodesIdVec[1] = 3 + 2 * i;
        nodesIdVec[2] = 2 + 2 * i;
        new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, matIdArc, *gmesh);
    }

    nodesIdVec[0] = 15;
    nodesIdVec[1] = 1;
    nodesIdVec[2] = 16;
    //para o broblema do douglas matIdArc=3
    new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, 3, *gmesh);
    //    // Finally, inserts line elements to complete boundary
    //    nodesIdVec.Resize(2);
    //    nodesIdVec[0] = 0;
    //    nodesIdVec[1] = 1;
    //    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, matIdArc, *gmesh);
    //
    //    nodesIdVec[0] = 0;
    //    nodesIdVec[1] = 14;
    //    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, matIdArc, *gmesh);

    gmesh->BuildConnectivity();

    return gmesh;
}

TPZGeoMesh *CreateLMHMMesh(int nDiv, TPZVec<int64_t>& coarseIndexes) {
    
    int factor = (int)(pow(2, nDiv) + 0.5);
    
    int xElements = factor * 6;
    int yElements = factor * 4;
    
    TPZManVector<int> nx(2);
    nx[0] = xElements;
    nx[1] = yElements;
    
    TPZManVector<REAL> x0(3,0.), x1(3,1.);
    x1[2] = 0.;
    
    
    TPZGenGrid2D gen(nx, x0, x1, 1, 0);
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gen.Read(gmesh);
    gen.SetBC(gmesh, 4, -1);
    gen.SetBC(gmesh, 5, -1);
    gen.SetBC(gmesh, 6, -1);
    gen.SetBC(gmesh, 7, -1);
    
    gmesh->SetDimension(2);
    gmesh->BuildConnectivity();
    
    // Assigns matIDs to create L elements
    {
        int coarseIndex = 0;
        
        // Get number of 2D elements
        int64_t nelem = xElements * yElements;
        
        coarseIndexes.Resize(nelem, -1);
        for (int64_t elem = 0; elem < nelem; elem++) {
            TPZGeoEl *gel = gmesh->ElementVec()[elem];
            if (gel->Dimension() != 2) DebugStop();
            
            int lineInPattern = elem / nx[0] % 4;
            int colInPattern = elem % nx[0] % 6;
            
            // IDs of elements in the neighbourhood to which a coarse index has been already assigned
            int leftEl = elem - 1;
            int bottomEl = elem - nx[0];
            
            if (lineInPattern == 0) {
                if (colInPattern % 2 == 0) {
                    coarseIndexes[elem] = coarseIndex;
                    coarseIndex++;
                } else {
                    coarseIndexes[elem] = coarseIndexes[leftEl];
                }
            } else if (lineInPattern == 1) {
                if (colInPattern == 0 || colInPattern == 2 || colInPattern == 5) {
                    coarseIndexes[elem] = coarseIndexes[bottomEl];
                } else if (colInPattern == 4) {
                    coarseIndexes[elem] = coarseIndexes[leftEl];
                } else {
                    coarseIndexes[elem] = coarseIndex;
                    coarseIndex++;
                }
            } else if (lineInPattern == 2) {
                if (colInPattern == 1 || colInPattern == 4) {
                    coarseIndexes[elem] = coarseIndexes[bottomEl];
                } else if (colInPattern == 2) {
                    coarseIndexes[elem] = coarseIndexes[leftEl];
                } else {
                    coarseIndexes[elem] = coarseIndex;
                    coarseIndex++;
                }
            } else if (lineInPattern == 3) {
                if (colInPattern == 0) {
                    coarseIndexes[elem] = coarseIndexes[bottomEl];
                } else if (colInPattern % 2 == 1) {
                    coarseIndexes[elem] = coarseIndexes[leftEl];
                } else {
                    coarseIndexes[elem] = coarseIndexes[leftEl] + 1;
                }
            }
        }
    }
    
//    for(int ilinha =0 ; ilinha <coarseIndexes.size() ;ilinha++){
//        std::cout<<"coarseIndex "<<coarseIndexes[ilinha]<<std::endl;
//    }
    
    return gmesh;
}

void SolveProblem(TPZAutoPointer<TPZCompMesh> cmesh, const TPZVec<TPZAutoPointer<TPZCompMesh>> &compmeshes,
                  TPZAnalyticSolution *analytic, const std::string &prefix, TRunConfig config) {
    //calculo solution
    bool shouldrenumber = true;
    TPZLinearAnalysis an(cmesh,shouldrenumber);
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<> strmat(cmesh.operator->());
    strmat.SetNumThreads(0/*config.n_threads*/);
#else
    TPZSkylineStructMatrix<STATE> strmat(cmesh.operator->());
    strmat.SetNumThreads(config.n_threads);
#endif


#ifdef ERRORESTIMATION_DEBUG
    if(0)
    {
        std::ofstream file("MeshToSolveProblem.txt");
        cmesh->Print(file);
    }
#endif


    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "Assembling\n";
    an.Assemble();

    std::cout << "Solving\n";
    an.Solve();
    std::cout << "Finished\n";
    an.LoadSolution(); // compute internal dofs


    TPZMFSolutionTransfer transfer;
    transfer.BuildTransferData(cmesh.operator->());
//    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(compmeshes, cmesh);
    transfer.TransferFromMultiphysics();

    if(0)
    {
        std::ofstream out1("mfmesh.txt");
        cmesh->Print(out1);
        std::ofstream out2("flux.txt");
        compmeshes[0]->Print(out2);
        std::ofstream out3("pressure.txt");
        compmeshes[1]->Print(out3);
        std::ofstream out4("transfer.txt");
        transfer.Print(out4);
    }

    TPZStack<std::string> scalnames,vecnames;
    TPZMaterial *mat = cmesh->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }
    if (analytic)
    {
        an.SetExact(analytic->ExactSolution());
        scalnames.Push("ExactPressure");
        vecnames.Push("ExactFlux");
    }
    {
        scalnames.Push("Pressure");
        //  scalnames.Push("Permeability");
        vecnames.Push("Flux");
        //   vecnames.Push("Derivative");
    }


    std::string plotname;
    {
        std::stringstream out;
        out << "MHMHdiv"<<"_kin" <<config.pOrderInternal << "ksk_"<<config.pOrderSkeleton << "hsk_" <<config.numDivSkeleton<<"hin_"<< config.numHDivisions<< ".vtk";
        plotname = out.str();

    }
    int resolution=0;
    an.DefineGraphMesh(cmesh->Dimension() , scalnames, vecnames, plotname);
    an.PostProcess(resolution,cmesh->Dimension() );


    if(analytic)
    {
        TPZManVector<REAL> errors(4,0.);
        an.SetThreadsForError(config.n_threads);
        an.SetExact(analytic->ExactSolution());
        an.PostProcessError(errors,false);

        //Erro

        std::ofstream myfile;
        myfile.open("ArquivosErrosMHM.txt", std::ios::app);
        myfile << "\n\n Error for MHM formulation " ;
        myfile << "\n-------------------------------------------------- \n";
        myfile << "Ndiv = " << config.numHDivisions << " Order Internal= " << config.pOrderInternal <<" Order Skeleton= " << config.pOrderSkeleton <<"\n";
        myfile << "DOF Total = " << cmesh->NEquations() << "\n";
        myfile << "Energy norm = " << errors[0] << "\n";//norma energia
        myfile << "error norm L2 = " << errors[1] << "\n";//norma L2
        myfile << "Semi norm H1 = " << errors[2] << "\n";//norma L2
        myfile.close();

    }

}
