//
// Created by victor on 16/03/2021.
//

#include "MeshInit.h"
#include "DataStructure.h"
#include <TPZMultiphysicsCompMesh.h>
#include "TPZMatLaplacianHybrid.h"
#include "TPZNullMaterial.h"
#include "TPZGenGrid2D.h"
#include "TPZCompMeshTools.h"
#include "TPZCompElDisc.h"
#include "TPZHybMixDiffMaterial.h"
#include "ForcingFunction.h"
#include "DarcyFlow/TPZHybridDarcyFlow.h"
#include "InputTreatment.h"

void InsertMaterialMixed_MultiK(TPZMultiphysicsCompMesh *cmesh_mixed, ProblemConfig &config, PreConfig &pConfig){
    int matID_Q1 = 2;
    int matID_Q2 = 3;
    int dim = config.gmesh->Dimension();
    int dirichlet = 0;
    int neumann = 1;

    TLaplaceExample1 *mat1 = new TLaplaceExample1,*mat2 = new TLaplaceExample1;
    mat1->fExact = ChooseAnaliticSolution(pConfig);
    mat2->fExact = ChooseAnaliticSolution(pConfig);

    TPZFNMatrix<9,REAL> K, invK;
    K.Resize(3,3); invK.Resize(3,3);
    K.Identity();invK.Identity();
    K(0,0) = K(1,1) = pConfig.perm_Q1;
    invK(0,0) = invK(1,1) = 1./pConfig.perm_Q1;
    mat1->setPermeabilyTensor(K,invK);
    K(0,0) = K(1,1) = pConfig.perm_Q2;
    invK(0,0) = invK(1,1) =  1./pConfig.perm_Q2;
    mat2->setPermeabilyTensor(K,invK);

    TPZMixedPoisson *material_Q1 = new TPZMixedPoisson(matID_Q1,dim); //Using standard PermealityTensor = Identity.
    TPZMixedPoisson *material_Q2 = new TPZMixedPoisson(matID_Q2,dim);
    int pOrder =5;
    material_Q1->SetForcingFunction(config.exact->ForceFunc(), pOrder);
    material_Q1->SetExactSol(config.exact->ExactSolution(), pOrder);
    material_Q2->SetForcingFunction(config.exact->ForceFunc(), pOrder);
    material_Q2->SetExactSol(config.exact->ExactSolution(), pOrder);

    material_Q1->SetConstantPermeability(pConfig.perm_Q1);
    material_Q2->SetConstantPermeability(pConfig.perm_Q2);

    cmesh_mixed->InsertMaterialObject(material_Q1);
    cmesh_mixed->InsertMaterialObject(material_Q2);

    //Boundary Conditions
    TPZFMatrix<STATE> val1(2,2,0.);
    TPZManVector<STATE, 2> val2(2, 0.);

    TPZBndCondT<STATE>* BCond0_Q1 = material_Q1->CreateBC(material_Q1, -5, dirichlet, val1, val2);
    BCond0_Q1->SetForcingFunctionBC(config.exact.operator*().ExactSolution(),pOrder);

    auto *BCond0_Q2 = material_Q2->CreateBC(material_Q2, -6, dirichlet, val1, val2);
    BCond0_Q2->SetForcingFunctionBC(config.exact.operator*().ExactSolution(), pOrder);

    //auto *BCond1_Q1 = material_Q1->CreateBC(material_Q1, -8, neumann, val1, val2);
    //auto *BCond1_Q2 = material_Q1->CreateBC(material_Q2, -9, neumann, val1, val2);

    cmesh_mixed->InsertMaterialObject(BCond0_Q1);
    cmesh_mixed->InsertMaterialObject(BCond0_Q2);
    //cmesh_mixed->InsertMaterialObject(BCond1_Q1);
    //cmesh_mixed->InsertMaterialObject(BCond1_Q2);
}

void InsertMaterialMixed_ELaplace(TPZMultiphysicsCompMesh *cmesh_mixed, ProblemConfig &config, PreConfig &pConfig){
    int dim = config.gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    TPZMixedPoisson *material = new TPZMixedPoisson(matID, dim); //Using standard PermealityTensor = Identity.
    if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
        material->SetExactSol(config.exact->ExactSolution(), material->PolynomialOrderExact());
    }
    cmesh_mixed->InsertMaterialObject(material);

    //Boundary Conditions
    TPZFMatrix<STATE> val1(2, 2, 0.);
    TPZManVector<STATE, 2> val2(2, 0.);

    auto *BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2);
    auto *BCond1 = material->CreateBC(material, -2, neumann, val1, val2);
    auto *BCond2 = material->CreateBC(material, -3, neumann, val1, val2);
    if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
        int pOrder=5;
        BCond0->SetForcingFunctionBC(config.exact.operator*().ExactSolution(),pOrder=5);
    }
    int pOrder=5;
    DebugStop();
    //    BCond2->SetForcingFunctionBC(LinearFunc,pOrder);

    cmesh_mixed->InsertMaterialObject(BCond0);
    cmesh_mixed->InsertMaterialObject(BCond1);
    cmesh_mixed->InsertMaterialObject(BCond2);
}

void InsertMaterialHybrid_MultiK(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, ProblemConfig &config, PreConfig &pConfig){

    TPZGeoMesh* gmesh = cmesh_H1Hybrid->Reference();
    int matID_Q1 = 2;
    int matID_Q2 = 3;
    int dim = gmesh->Dimension();
    int dirichlet = 0;
    int neumann = 1;

    TLaplaceExample1 *mat1 = new TLaplaceExample1,*mat2 = new TLaplaceExample1;
    mat1->fExact = ChooseAnaliticSolution(pConfig);
    mat2->fExact = ChooseAnaliticSolution(pConfig);

    TPZFNMatrix<9,REAL> K, invK;
    K.Resize(3,3); invK.Resize(3,3);
    K.Identity();invK.Identity();
    K(0,0) = K(1,1) = pConfig.perm_Q1;
    invK(0,0) = invK(1,1) = 1./pConfig.perm_Q1;
    mat1->setPermeabilyTensor(K,invK);
    K(0,0) = K(1,1) = pConfig.perm_Q2;
    invK(0,0) = invK(1,1) =  1./pConfig.perm_Q2;
    mat2->setPermeabilyTensor(K,invK);

    TPZMatLaplacianHybrid *material_Q1 = new TPZMatLaplacianHybrid(matID_Q1, dim);
    TPZMatLaplacianHybrid *material_Q2 = new TPZMatLaplacianHybrid(matID_Q2, dim);

    material_Q1->SetConstantPermeability(pConfig.perm_Q1);
    material_Q2->SetConstantPermeability(pConfig.perm_Q2);

    cmesh_H1Hybrid->InsertMaterialObject(material_Q1);
    cmesh_H1Hybrid->InsertMaterialObject(material_Q2);

    int porder = 3;
    if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
        //material_Q1->SetForcingFunction(mat1->ForceFunc(), porder);
        //TPZMatErrorCombinedSpaces<STATE> *err1 = material_Q1;
        //err1->SetExactSol(mat1->ExactSolution(), porder);
        material_Q1->SetExactSol(config.exact->ExactSolution(), porder);

        //material_Q2->SetForcingFunction(mat2->ForceFunc(), porder);
        //TPZMatErrorCombinedSpaces<STATE> *err2 = material_Q2;
        //err2->SetExactSol(mat2->ExactSolution(), porder);
        material_Q1->SetExactSol(config.exact->ExactSolution(), porder);
    }

    // Inserts boundary conditions
    TPZFMatrix<STATE> val1(1, 1, 0.);
    TPZManVector<STATE, 1> val2(1, 0.);
    auto *BCond0_Q1 = material_Q1->CreateBC(material_Q1, -5, dirichlet, val1, val2);
    auto *BCond0_Q2 = material_Q2->CreateBC(material_Q2, -6, dirichlet, val1, val2);
    if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
        int pOrder=5;
        BCond0_Q1->SetForcingFunctionBC(mat1->ExactSolution(),pOrder);
        BCond0_Q2->SetForcingFunctionBC(mat2->ExactSolution(),pOrder);
    }

    auto *BCond1_Q1 = material_Q1->CreateBC(material_Q1, -8, neumann, val1, val2);
    auto *BCond1_Q2 = material_Q1->CreateBC(material_Q1, -9, neumann, val1, val2);

    cmesh_H1Hybrid->InsertMaterialObject(BCond0_Q1);
    cmesh_H1Hybrid->InsertMaterialObject(BCond0_Q2);
    cmesh_H1Hybrid->InsertMaterialObject(BCond1_Q1);
    cmesh_H1Hybrid->InsertMaterialObject(BCond1_Q2);
}

void InsertMaterialHybrid_ELaplace(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,ProblemConfig &config, PreConfig &pConfig){
    TPZGeoMesh *gmesh = cmesh_H1Hybrid->Reference();
    int dim = gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    TPZMatLaplacianHybrid *material = new TPZMatLaplacianHybrid(matID, dim);
    material->SetConstantPermeability(1.);
    cmesh_H1Hybrid->InsertMaterialObject(material);

    int porder = 3;
    TPZMatErrorCombinedSpaces<STATE> *err = material;
    err->SetExactSol(config.exact.operator*().ExactSolution(), porder);

    // Inserts boundary conditions
    TPZFMatrix<STATE> val1(1, 1, 0.);
    TPZManVector<STATE, 1> val2(1, 0.);
    auto *BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2);
    auto *BCond1 = material->CreateBC(material, -2, neumann, val1, val2);
    auto *BCond2 = material->CreateBC(material, -3, neumann, val1, val2);

    int pOrder=5;
    if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
        BCond0->SetForcingFunctionBC(config.exact.operator*().ExactSolution(),pOrder);
    }
    DebugStop();
    //BCond2->SetForcingFunctionBC(LinearFunc,pOrder);

    cmesh_H1Hybrid->InsertMaterialObject(BCond0);
    cmesh_H1Hybrid->InsertMaterialObject(BCond1);
    cmesh_H1Hybrid->InsertMaterialObject(BCond2);
}

void SetFExact(TLaplaceExample1 *mat1, TLaplaceExample1 *mat2,PreConfig &pConfig){
    switch(pConfig.type) {
        case 0:
            mat1->fExact = TLaplaceExample1::ESinSin;
            mat2->fExact = TLaplaceExample1::ESinSin;
            break;
        case 1:
            mat1->fExact = TLaplaceExample1::EArcTan;
            mat2->fExact = TLaplaceExample1::EArcTan;
            break;
        case 2:
            mat1->fExact = TLaplaceExample1::ESteklovNonConst;
            mat2->fExact = TLaplaceExample1::ESteklovNonConst;
            break;
        default:
            DebugStop();
            break;
    }
    mat1->fSignConvention = 1;
    mat2->fSignConvention = 1;
}

void SetMultiPermeMaterials(TPZGeoMesh* gmesh){

    TPZAdmChunkVector<TPZGeoEl *> &elvec = gmesh->ElementVec();
    int numEl = elvec.NElements();
    TPZVec<int64_t> nodeInd;
    int numberNodes;
    TPZManVector<REAL,2> elCenterCoord = {0.,0.,0.};
    int meshDim = gmesh ->Dimension();

    if (meshDim == 2) {
        //Getting center coordinates of the element
        for (int ind = 0; ind < numEl; ind++) {
            TPZGeoEl *gel = elvec[ind];
            gel->GetNodeIndices(nodeInd);
            numberNodes = gel->NNodes();
            elCenterCoord[0] = elCenterCoord[1] = 0.;
            for (int nodeIter = 0; nodeIter < numberNodes; nodeIter++) {
                TPZGeoNode nodo = gmesh->NodeVec()[nodeInd[nodeIter]];
                elCenterCoord[0] += nodo.Coord(0);
                elCenterCoord[1] += nodo.Coord(1);
            }
            elCenterCoord[0] /= numberNodes;
            elCenterCoord[1] /= numberNodes;

            // Changing the permeability
            //std::cout <<"el: " << ind  << "x_ave: " << elCenterCoord[0] << "y_ave: " << elCenterCoord[1] << std::endl;
            if ((elCenterCoord[0] >= 0. && elCenterCoord[1] >= 0.) ||
                (elCenterCoord[0] < 0. && elCenterCoord[1] <= 0.)) { // first/third quadrants
                if (gel->Dimension() == gmesh->Dimension()) { /*if(gel->MaterialId() == 1)*/ gel->SetMaterialId(2); }
                else {/*if (gel->MaterialId() < 0 && gel->MaterialId() > -5)*/ gel->SetMaterialId(-5);
                } // Assumes dirichlet boundary condition
            } else {
                if (gel->Dimension() == gmesh->Dimension()) { /*if(gel->MaterialId() == 1)*/ gel->SetMaterialId(3); }
                else {/*if (gel->MaterialId() < 0 && gel->MaterialId() > -5)*/ gel->SetMaterialId(-6);
                }// Assumes dirichlet boundary condition
            }
        }
    }
    if (meshDim == 3) {
        //Getting center coordinates of the element
        for (int ind = 0; ind < numEl; ind++) {
            TPZGeoEl *gel = elvec[ind];
            gel->GetNodeIndices(nodeInd);
            numberNodes = gel->NNodes();
            elCenterCoord[0] = elCenterCoord[1] = 0., elCenterCoord[2] = 0.;
            for (int nodeIter = 0; nodeIter < numberNodes; nodeIter++) {
                TPZGeoNode nodo = gmesh->NodeVec()[nodeInd[nodeIter]];
                elCenterCoord[0] += nodo.Coord(0);
                elCenterCoord[1] += nodo.Coord(1);
                elCenterCoord[2] += nodo.Coord(2);
            }
            elCenterCoord[0] /= numberNodes;
            elCenterCoord[1] /= numberNodes;
            elCenterCoord[2] /= numberNodes;

            // Changing the permeability
            //std::cout <<"el: " << ind  << "x_ave: " << elCenterCoord[0] << "y_ave: " << elCenterCoord[1] << std::endl;
            if ((elCenterCoord[0] > 0. && elCenterCoord[1] > 0. &&  elCenterCoord[2] > 0.)  || //x,y,z > 0
                (elCenterCoord[0] > 0. && elCenterCoord[1] < 0. &&  elCenterCoord[2] < 0.)  || //y,z < 0, x > 0
                (elCenterCoord[0] < 0. && elCenterCoord[1] < 0. &&  elCenterCoord[2] > 0.)  || //x,y < 0, z > 0
                (elCenterCoord[0] < 0. && elCenterCoord[1] > 0. &&  elCenterCoord[2] < 0.))    //x,z < 0, y > 0
            {
                if (gel->Dimension() == meshDim) gel->SetMaterialId(2);
                else gel->SetMaterialId(-5);
            }
            else
            {
                if (gel->Dimension() == meshDim) gel->SetMaterialId(3);
                else gel->SetMaterialId(-6);
            }
        }
    }
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}

void BuildFluxMesh(TPZCompMesh *cmesh_flux, ProblemConfig &config, PreConfig &pConfig){

    int dim = config.gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    cmesh_flux->SetDefaultOrder(config.k);
    cmesh_flux->SetDimModel(config.gmesh->Dimension());

    cmesh_flux->SetAllCreateFunctionsHDiv();

    auto *material = new TPZNullMaterial(matID);
    material->SetDimension(dim);
    cmesh_flux->InsertMaterialObject(material);

    //(EE)Create Boundary conditions
    TPZFMatrix<STATE> val1(2,2,0.);
    TPZManVector<STATE, 2> val2(2,0.);

    //(PROPRIETARY)
    auto *BCond0 = material->CreateBC(material,-1,dirichlet,val1,val2);
    cmesh_flux->InsertMaterialObject(BCond0);

    auto *BCond1 = material->CreateBC(material,-2,neumann,val1,val2);
    cmesh_flux->InsertMaterialObject(BCond1);

    if (pConfig.type == 4) {
        auto *BCond2 = material->CreateBC(material, -3, neumann, val1, val2);
        cmesh_flux->InsertMaterialObject(BCond2);
    }
    //This part is for multi-K problemas
    if (pConfig.type == 2) BuildFluxMesh_MultiK(cmesh_flux,config);

    cmesh_flux->AutoBuild();
    cmesh_flux->InitializeBlock();
}

void CreateMixedAtomicMeshes(TPZVec<TPZCompMesh *> &meshvec, PreConfig &eData, ProblemConfig &config){
    //Flux mesh creation
    TPZCompMesh *cmesh_flux = new TPZCompMesh(config.gmesh);
    BuildFluxMesh(cmesh_flux,config,eData);
    {
        int64_t nelem = cmesh_flux->NElements();
        int dimgrid = cmesh_flux->Dimension();
        for (int64_t el = 0; el<nelem; el++) {
            TPZCompEl *cel = cmesh_flux->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            int dimgel = gel->Dimension();
            int nconnects = cel->NConnects();
            int nSides = gel->NSides(cmesh_flux->Dimension()-1);
            //if(dimgrid != 2) DebugStop();
            if(dimgel == dimgrid){
                if(nconnects != nSides+1){
                    DebugStop();
                }
                for(int ic = 0 ; ic < nSides; ic ++){
                    cel->Connect(ic).SetLagrangeMultiplier(3);
                }
            }
            if(dimgel == dimgrid-1){
                if(nconnects != nSides){
                    DebugStop();
                }
                cel->Connect(0).SetLagrangeMultiplier(3);
            }
            if (dimgel > dimgrid || dimgel < dimgrid - 1){
                DebugStop();
            }
        }
    }

    //Potential mesh creation
    TPZCompMesh *cmesh_p = new TPZCompMesh(config.gmesh);
    BuildPotentialMesh(cmesh_p,config,eData);
    {
        int64_t nconnects = cmesh_p->NConnects();
        for (int ic=0; ic<nconnects; ic++) {
            cmesh_p->ConnectVec()[ic].SetLagrangeMultiplier(1);
        }
    }


    TPZCompMesh *gspace = new TPZCompMesh(config.gmesh);
    {
        InsertNullSpaceMaterialIds(gspace,config);
        gspace->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
        gspace->SetDefaultOrder(0);//sao espacos de pressao media
        gspace->AutoBuild();
        int64_t nconnects = gspace->NConnects();
        for (int ic = 0; ic<nconnects; ic++) {
            gspace->ConnectVec()[ic].SetLagrangeMultiplier(2);
        }
        int64_t nel = gspace->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(gspace->Element(el));
            if(disc) disc->SetFalseUseQsiEta();
        }
    }
    TPZCompMesh *average = new TPZCompMesh(config.gmesh);
    {
        InsertNullSpaceMaterialIds(average,config);
        average->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
        average->SetDefaultOrder(0);
        average->AutoBuild();
        int64_t nconnects = average->NConnects();
        for (int ic = 0; ic<nconnects; ic++) {
            average->ConnectVec()[ic].SetLagrangeMultiplier(4);
        }
        int64_t nel = average->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(average->Element(el));
            if(disc) disc->SetFalseUseQsiEta();
        }
    }

    meshvec[0] = cmesh_flux;
    meshvec[1] = cmesh_p;
    meshvec[2] = gspace;
    meshvec[3] = average;
    TPZCompMeshTools::AdjustFluxPolynomialOrders(meshvec[0], config.n); //Increases internal flux order by "hdivmais"
    TPZCompMeshTools::SetPressureOrders(meshvec[0], meshvec[1]);//Set the pressure order the same as the internal flux
}

void InsertNullSpaceMaterialIds(TPZCompMesh *nullspace, ProblemConfig &config)
{
    for (auto matid:config.materialids) {
        auto *nullmat = new TPZNullMaterial(matid);
        nullmat->SetDimension(config.gmesh->Dimension());
        nullmat->SetNStateVariables(1);
        nullspace->InsertMaterialObject(nullmat);
    }
}


void BuildFluxMesh_MultiK(TPZCompMesh *cmesh_flux, ProblemConfig &config) {

    int dim = config.gmesh->Dimension();
    int matID_Q1 = 2;
    int matID_Q2 = 3;
    int dirichlet = 0;
    int neumann = 1;

    auto *material_Q1 = new TPZNullMaterial<STATE>(matID_Q1);
    auto *material_Q2 = new TPZNullMaterial<STATE>(matID_Q2);
    material_Q1->SetDimension(dim);
    material_Q2->SetDimension(dim);
    cmesh_flux->InsertMaterialObject(material_Q1);
    cmesh_flux->InsertMaterialObject(material_Q2);

    //(EE)Create Boundary conditions
    TPZFMatrix<STATE> val1(2, 2, 0.);
    TPZManVector<STATE, 2> val2(2, 0.);

    //(PROPRIETARY)
    auto *BCond0_Q1 = material_Q1->CreateBC(material_Q1, -5, dirichlet, val1, val2);
    auto *BCond0_Q2 = material_Q1->CreateBC(material_Q2, -6, dirichlet, val1, val2);
    cmesh_flux->InsertMaterialObject(BCond0_Q1);
    cmesh_flux->InsertMaterialObject(BCond0_Q2);

    auto *BCond1_Q1 = material_Q1->CreateBC(material_Q1, -8, neumann, val1, val2);
    auto *BCond1_Q2 = material_Q2->CreateBC(material_Q2, -9, neumann, val1, val2);
    cmesh_flux->InsertMaterialObject(BCond1_Q1);
    cmesh_flux->InsertMaterialObject(BCond1_Q2);
}

void BuildPotentialMesh(TPZCompMesh *cmesh_p, ProblemConfig &config, PreConfig &pConfig){
    int matID = 1;
    int dim = config.gmesh->Dimension();

    int potential_order = config.k + config.n;
    cmesh_p->SetDefaultOrder(potential_order);
    cmesh_p->SetDimModel(dim);

    cmesh_p->SetAllCreateFunctionsContinuous(); //H1 functions
    cmesh_p->ApproxSpace().CreateDisconnectedElements(true);

    auto *material = new TPZNullMaterial<STATE>(matID);
    material->SetDimension(dim);
    cmesh_p->InsertMaterialObject(material);

    if (pConfig.type == 2) BuildPotentialMesh_MultiK(cmesh_p, config);

    cmesh_p->AutoBuild();
    cmesh_p->ExpandSolution();

    TPZAdmChunkVector<TPZConnect> &nodeIt = cmesh_p->ConnectVec();
    for(auto &nodo : nodeIt){
        nodo.SetLagrangeMultiplier(1);
    }
}

void BuildPotentialMesh_MultiK(TPZCompMesh *cmesh_p, ProblemConfig &config){
    int matID_Q1 = 2;
    int matID_Q2 = 3;
    int dim = config.gmesh->Dimension();

    auto *material_Q1 = new TPZNullMaterial(matID_Q1); material_Q1->SetDimension(dim);
    cmesh_p->InsertMaterialObject(material_Q1);

    auto *material_Q2 = new TPZNullMaterial(matID_Q2); material_Q2->SetDimension(dim);
    cmesh_p->InsertMaterialObject(material_Q2);
}

void InsertMaterialMixed(TPZMultiphysicsCompMesh *cmesh_mixed, ProblemConfig config, PreConfig &pConfig){

    int dim = config.gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    cmesh_mixed->SetDefaultOrder(config.k + config.n);
    cmesh_mixed->SetDimModel(dim);
    cmesh_mixed->SetAllCreateFunctionsMultiphysicElem();

    switch(pConfig.type){
        case 2:
            InsertMaterialMixed_MultiK(cmesh_mixed, config, pConfig);
            break;
        case 4:
            InsertMaterialMixed_ELaplace(cmesh_mixed, config, pConfig);
            break;
        default:
            TPZMixedPoisson *material = new TPZMixedPoisson(matID, dim); //Using standard PermealityTensor = Identity.
            if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
                material->SetForcingFunction(config.exact->ForceFunc(), material->ForcingFunctionPOrder());
                material->SetExactSol(config.exact->ExactSolution(), material->PolynomialOrderExact());
            }
            cmesh_mixed->InsertMaterialObject(material);

            //Boundary Conditions
            TPZFMatrix<STATE> val1(2, 2, 0.);
            TPZManVector<STATE, 2> val2(2, 0.);

            auto *BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2);
            auto *BCond1 = material->CreateBC(material, -2, neumann, val1, val2);
            int pOrder=5;
            if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
                BCond0->SetForcingFunctionBC(config.exact.operator*().ExactSolution(),pOrder);
                BCond1->SetForcingFunctionBC(config.exact.operator*().ExactSolution(),pOrder);
            }

            cmesh_mixed->InsertMaterialObject(BCond0);
            cmesh_mixed->InsertMaterialObject(BCond1);
            break;
    }
}

void InsertMaterialMixHyb(TPZMultiphysicsCompMesh *multMesh, PreConfig &pConfig, ProblemConfig &config){

    int dim = pConfig.dim;
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    if(pConfig.type != 2) {
        multMesh->SetDefaultOrder(pConfig.k + pConfig.n);
        multMesh->SetDimModel(dim);
        multMesh->SetAllCreateFunctionsMultiphysicElem();

        TPZHybMixDiffMaterial *material = new TPZHybMixDiffMaterial(matID, dim); //Using standard PermealityTensor = Identity.
        if(pConfig.debugger) {
            material->SetForcingFunction(config.exact->ForceFunc(), material->ForcingFunctionPOrder());
            material->SetExactSol(config.exact->ExactSolution(), material->PolynomialOrderExact());
        } else {
            //material->SetInternalFlux(1.);
        }
        multMesh->InsertMaterialObject(material);

        //Boundary Conditions
        TPZFMatrix<STATE> val1(2, 2, 0.);
        TPZManVector<STATE, 2> val2(21, 0.);

        int pOrder=5;
        auto *BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2);
        if(pConfig.debugger)
        {
            BCond0->SetForcingFunctionBC(config.exact->ExactSolution(),pOrder);
        }

        auto *BCond1 = material->CreateBC(material, -2, neumann, val1, val2);

        multMesh->InsertMaterialObject(BCond0);
        multMesh->InsertMaterialObject(BCond1);
    }

    else {
        DebugStop();
    }
}

void InsertMaterialHybrid(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, ProblemConfig &config, PreConfig &pConfig)
{
    TPZGeoMesh *gmesh = cmesh_H1Hybrid->Reference();
    int dim = gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    switch(pConfig.type){
        case 2:
            InsertMaterialHybrid_MultiK(cmesh_H1Hybrid, config, pConfig);
            break;
        case 4:
            InsertMaterialHybrid_ELaplace(cmesh_H1Hybrid, config, pConfig);
            break;
        default:
            auto *material = /*new TPZHybridDarcyFlow(matID,dim);*/new TPZMatLaplacianHybrid(matID, dim);
            TPZMatErrorCombinedSpaces<STATE> *err = material;
            material->SetConstantPermeability(1.);
            cmesh_H1Hybrid->InsertMaterialObject(material);

            if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
                material->SetForcingFunction(config.exact->ForceFunc(), pConfig.integrationorder);
                err->SetExactSol(config.exact->ExactSolution(), pConfig.integrationorder);
            }

            // Inserts boundary conditions
            TPZFMatrix<STATE> val1(1, 1, 0.);
            TPZManVector<STATE, 2> val2(11, 0.);
            auto *BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2);
            auto *BCond1 = material->CreateBC(material, -2, neumann, val1, val2);

            if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
                BCond0->SetForcingFunctionBC(config.exact->ExactSolution(),pConfig.integrationorder);
                BCond1->SetForcingFunctionBC(config.exact->ExactSolution(),pConfig.integrationorder);
            }

            cmesh_H1Hybrid->InsertMaterialObject(BCond0);
            cmesh_H1Hybrid->InsertMaterialObject(BCond1);
            break;
    }
}

TPZCompMesh* InsertCMeshH1(ProblemConfig &config, PreConfig &pConfig) {

    TPZCompMesh* cmesh = new TPZCompMesh(config.gmesh);
    TPZDarcyFlow* mat = 0;
    int matID_Q1 = 2;
    int matID_Q2 = 3;
    int dirichlet = 0;
    int neumann = 1;
    int dim = config.gmesh->Dimension();

    if(pConfig.type != 2) {
        for (auto matid : config.materialids) {
            TPZDarcyFlow *mix = new TPZDarcyFlow(matid, cmesh->Dimension());
//            mix->SetExactSol(config.exact->ExactSolution(), mix->PolynomialOrderExact());
//            mix->SetForcingFunction(config.exact->ForceFunc(), mix->ForcingFunctionPOrder());
            
            mix->SetExactSol(config.exact->ExactSolution(), pConfig.integrationorder);
            mix->SetForcingFunction(config.exact->ForceFunc(), pConfig.integrationorder);

            if (!mat) mat = mix;
            cmesh->InsertMaterialObject(mix);
        }

        //int pOrder=18;
        for (auto matid : config.bcmaterialids) {
            TPZFNMatrix<1, REAL> val1(1, 1, 0.);
            TPZManVector<STATE, 2> val2(11, 0.);
            int bctype = 0;
            auto *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
            bc->SetForcingFunctionBC(config.exact->ExactSolution(),pConfig.integrationorder);
            cmesh->InsertMaterialObject(bc);
        }
    }
    else{

        TLaplaceExample1 *mat1 = new TLaplaceExample1,*mat2 = new TLaplaceExample1;
        SetFExact(mat1,mat2,pConfig);

        TPZFNMatrix<9,REAL> K, invK;
        K.Resize(3,3); invK.Resize(3,3);
        K.Identity();invK.Identity();
        K(0,0) = K(1,1) = pConfig.perm_Q1;
        invK(0,0) = invK(1,1) = 1./pConfig.perm_Q1;
        mat1->setPermeabilyTensor(K,invK);
        K(0,0) = K(1,1) = pConfig.perm_Q2;
        invK(0,0) = invK(1,1) =  1./pConfig.perm_Q2;
        mat2->setPermeabilyTensor(K,invK);

        TPZDarcyFlow *material_Q1 = new TPZDarcyFlow(matID_Q1, dim);
        TPZDarcyFlow *material_Q2 = new TPZDarcyFlow(matID_Q2, dim);

        material_Q1->SetConstantPermeability(pConfig.perm_Q1);
        material_Q2->SetConstantPermeability(pConfig.perm_Q2);

        if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
            material_Q1->SetForcingFunction(mat1->ForceFunc(), material_Q1->ForcingFunctionPOrder());
            material_Q1->SetExactSol(mat1->ExactSolution(), material_Q2->PolynomialOrderExact());

            material_Q2->SetForcingFunction(mat2->ForceFunc(), material_Q2->ForcingFunctionPOrder());
            material_Q2->SetExactSol(mat2->ExactSolution(), material_Q2->PolynomialOrderExact());
        }

        cmesh->InsertMaterialObject(material_Q1);
        cmesh->InsertMaterialObject(material_Q2);

        // Inserts boundary conditions
        TPZFMatrix<STATE> val1(1, 1, 0.);
        TPZManVector<STATE, 2> val2(11, 1.);
        auto *BCond0_Q1 = material_Q1->CreateBC(material_Q1, -5, dirichlet, val1, val2);
        auto *BCond0_Q2 = material_Q2->CreateBC(material_Q2, -6, dirichlet, val1, val2);
        int pOrder=5;
        if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
            BCond0_Q1->SetForcingFunctionBC(mat1->ExactSolution(),pOrder);
            BCond0_Q2->SetForcingFunctionBC(mat2->ExactSolution(),pOrder);
        }
        auto *BCond1_Q1 = material_Q1->CreateBC(material_Q1, -8, neumann, val1, val2);
        auto *BCond1_Q2 = material_Q1->CreateBC(material_Q1, -9, neumann, val1, val2);

        cmesh->InsertMaterialObject(BCond0_Q1);
        cmesh->InsertMaterialObject(BCond0_Q2);
        cmesh->InsertMaterialObject(BCond1_Q1);
        cmesh->InsertMaterialObject(BCond1_Q2);
    }

    cmesh->SetDefaultOrder(config.k);//ordem
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();

    cmesh->AutoBuild();


    return cmesh;
}

void FluxErrorInsertMaterial(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, ProblemConfig &config, PreConfig &pConfig)
{
    TPZGeoMesh *gmesh = cmesh_H1Hybrid->Reference();
    int dim = gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    // Creates Poisson material
    TPZMatLaplacianHybrid *material = new TPZMatLaplacianHybrid(matID, dim);
    material->SetConstantPermeability(1.);
    cmesh_H1Hybrid->InsertMaterialObject(material);

    int porder = 3;
    material->SetForcingFunction(config.exact->ForceFunc(), porder);
    TPZMatErrorCombinedSpaces<STATE> *err = material;
    err->SetExactSol(config.exact->ExactSolution(), porder);

    // Inserts boundary conditions
    TPZFMatrix<STATE> val1(1, 1, 0.);
    TPZManVector<STATE, 2> val2(11, 1.);
    auto *BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2);

    if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
        int pOrder=5;
        BCond0->SetForcingFunctionBC(config.exact->ExactSolution(),pOrder=5);
    }
    auto *BCond1 = material->CreateBC(material, -2, neumann, val1, val2);

    cmesh_H1Hybrid->InsertMaterialObject(BCond0);
    cmesh_H1Hybrid->InsertMaterialObject(BCond1);

}
