//
//  New_main.cpp
//  ErrorEstimateHDiv
//
//  Created by Denise De Siqueira on 01/04/19.
//

#include "pzlog.h"
#include "TPZRefPatternDataBase.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"

#include "tpzquadraticquad.h"
#include "tpzquadratictrig.h"
#include "tpzgeoelrefpattern.h"
#include "tpzarc3d.h"


#include "ProblemConfig.h"

#include "mixedpoisson.h"
#include "TPZVecL2.h"
#include "pzbndcond.h"

#include "pzintel.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZHybridizeHDiv.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZMultiphysicsCompMesh.h"

#include "TPZHybridHDivErrorEstimator.h"
#include "Tools.h"

#include "TPZBFileStream.h"
#include <tuple>
#include <memory>

void PlotLagrangreMultiplier(TPZCompMesh *cmesh);
void SolveHybridProblem(TPZCompMesh *Hybridmesh,int n2);
TPZManVector<REAL,3> ParametricCircle(REAL radius,REAL theta);
TPZGeoMesh *MakeCircle( int ndiv);

bool IsgmeshReader=true;


int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    
    ProblemConfig config;
    config.porder = 1;
    config.hdivmais = 1;
    
    config.ndivisions=3;
    config.prefine=false;
    config.makepressurecontinuous = true;
    config.exact.fExact = TLaplaceExample1::ESinMark;//ESinSinDirNonHom;//ESinSin;//EArcTanSingular;//ESinMark;//EArcTan;//
    config.problemname = "ESinMark";//"ESinSinDirNonHom";//"ESinSin";////"EArcTanSingular_PRef";//""ArcTang";//
    
 //   FunctionTest();
    
    int dim=2;
    
    //malha geometrica
    TPZGeoMesh *gmesh = nullptr;
    
    if(IsgmeshReader){

        std::string meshfilename = "../LCircle.msh";
        TPZGmshReader gmsh;
        gmsh.GetDimNamePhysical()[1]["dirichlet"] =2;
        gmsh.GetDimNamePhysical()[2]["domain"] = 1;
        
        config.materialids.insert(1);
        config.bcmaterialids.insert(2);
        
        gmsh.SetFormatVersion("4.1");
        gmesh = gmsh.GeometricGmshMesh(meshfilename);
        gmsh.PrintPartitionSummary(std::cout);
        gmesh->SetDimension(dim);
        config.gmesh = gmesh;
    }
    
    else{
    
    
//    int nDiv=6;
//    int nelem = 1<<nDiv;
    gmesh = CreateGeoMesh(1);
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.bcmaterialids.insert(-2);
    config.gmesh = gmesh;
    }
    
    UniformRefinement(config.ndivisions, gmesh);
    
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmeshInitial.txt");
        gmesh->Print(out2);
        
    }
    
    TPZManVector<TPZCompMesh*, 2> meshvec_HDiv(2, 0);
    TPZMultiphysicsCompMesh *cmesh_HDiv = CreateHDivMesh(config);//Hdiv x L2
    cmesh_HDiv->InitializeBlock();
    meshvec_HDiv = cmesh_HDiv->MeshVector();
    
    //cria malha hibrida
    
    TPZHybridizeHDiv hybrid;
    auto HybridMesh = hybrid.Hybridize(cmesh_HDiv);
    HybridMesh->CleanUpUnconnectedNodes();//enumerar adequadamente os connects
    HybridMesh->AdjustBoundaryElements();
    delete cmesh_HDiv;
    delete meshvec_HDiv[0];
    delete meshvec_HDiv[1];
    
    int n=hybrid.fLagrangeInterface;
    int n1=hybrid.fHDivWrapMatid;
    int n2=hybrid.fInterfaceMatid;
    
    std::cout<<"---Original PerifericalMaterialId --- "<<std::endl;
    std::cout <<" LagrangeInterface = "<<n<<std::endl;
    std::cout <<" HDivWrapMatid = "<<n1<<std::endl;
    std::cout <<" InterfaceMatid = "<<n2<<std::endl;
    
    
    cmesh_HDiv=(HybridMesh);//malha hribrida
    meshvec_HDiv[0] = (HybridMesh)->MeshVector()[0];//malha Hdiv
    meshvec_HDiv[1] = (HybridMesh)->MeshVector()[1];//malha L2
    
    
    SolveHybridProblem(cmesh_HDiv,n2);
    
    
    {
        
//                std::ofstream outgeo("HrybridGeometria.txt");
//                std::get<0>(HybridMesh)->Reference()->Print(outgeo);
//                std::ofstream out("OriginalHybridMesh.txt");
//                std::get<0>(HybridMesh)->Print(out);
        //
        std::ofstream out2("OriginalFluxMesh.txt");
        meshvec_HDiv[0]->Print(out2);
        
        std::ofstream out3("OriginalPotentialMesh.txt");
        meshvec_HDiv[1]->Print(out3);
        
        
    }
    
    PlotLagrangreMultiplier(meshvec_HDiv[1]);
    
    
    //reconstroi potencial e calcula o erro
    {
        TPZHybridHDivErrorEstimator HDivEstimate(*cmesh_HDiv);
        
        HDivEstimate.fProblemConfig = config;
        HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
        HDivEstimate.SetAnalyticSolution(config.exact);
        
        HDivEstimate.PotentialReconstruction();
        
        TPZManVector<REAL> elementerrors;
        HDivEstimate.ComputeErrors(elementerrors);
        
    }
    
    delete cmesh_HDiv;
    delete meshvec_HDiv[0];
    delete meshvec_HDiv[1];
    return 0;
    
    
}


void SolveHybridProblem(TPZCompMesh *Hybridmesh,int n2){
    
    TPZAnalysis an(Hybridmesh);
    
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(Hybridmesh);
    strmat.SetNumThreads(0);
    //        strmat.SetDecomposeType(ELDLt);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(Hybridmesh);
    strmat.SetNumThreads(2);
    //        TPZSkylineStructMatrix strmat3(cmesh_HDiv);
    //        strmat3.SetNumThreads(8);
#endif
    
    
    std::set<int> matIds;
    matIds.insert(1);
    matIds.insert(-1);
    matIds.insert(-2);
   
   // matIds.insert(4);
    matIds.insert(n2);
    strmat.SetMaterialIds(matIds);

    an.SetStructuralMatrix(strmat);
    
    
    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    scalnames.Push("ExactPressure");
    vecnames.Push("Flux");
    vecnames.Push("ExactFlux");
   // scalnames.Push("Divergence");
    an.DefineGraphMesh(2, scalnames, vecnames, "OriginalHybrid_Problem.vtk");
    //        meshvec_Hybrid[1]->Solution().Print("Press");
    // Post processing
    an.PostProcess(2,2);
    
}
void PlotLagrangreMultiplier(TPZCompMesh *cmesh){
    
    TPZAnalysis an(cmesh,false);
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("State");
    
    int dim = cmesh->Reference()->Dimension()-1;
    std::string plotname;
    {
        std::stringstream out;
        out << "OriginalLagrangeMultiplier" << ".vtk";
        plotname = out.str();
    }
    an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
    an.PostProcess(2,dim);
    
}

TPZGeoMesh *MakeCircle( int ndiv)
{
    int fbc1=-1;
    int dim=2;
    int matId=1;
    
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    geomesh->SetDimension(dim);
    
    int nodes =  17 + 8;
    REAL radius = 1.0;
    REAL innerradius = radius/2.0;
    geomesh->SetMaxNodeId(nodes-1);
    geomesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,7> Node(nodes);
    
    TPZManVector<int64_t,8> TopolQQuadrilateral(8);
    TPZManVector<int64_t,8> TopolQuadrilateral(4);
    TPZManVector<int64_t,6> TopolQTriangle(6);
    TPZManVector<int64_t,2> TopolLine(2);
    TPZManVector<int64_t,3> TopolArc(3);
    TPZManVector<REAL,3> coord(3,0.);
    TPZVec<REAL> xc(3,0.);
    
    
    int nodeindex = 0;
    //8 nos na circunferencia de raio 1/2
    for (int inode = 0; inode < 8 ; inode++) {
        // i node
        coord = ParametricCircle(innerradius, inode * M_PI/4.0);
        geomesh->NodeVec()[nodeindex].SetCoord(coord);
        geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
        nodeindex++;
    }
    
    //8 nos na circunferencia de raio 1
    
    for (int inode = 0; inode < 8 ; inode++) {
        // i node
        coord = ParametricCircle(radius, inode * M_PI/4.0);
        geomesh->NodeVec()[nodeindex].SetCoord(coord);
        geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
        nodeindex++;
    }
    

    
    // 24 node id at circle center
    coord = ParametricCircle(0.0,0.0);
    geomesh->NodeVec()[nodeindex].SetCoord(coord);
    geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    int elementid = 0;
    
    TopolQTriangle[0] = 16;
    TopolQTriangle[1] = 8;
    TopolQTriangle[2] = 10;
    TopolQTriangle[3] = 0;
    TopolQTriangle[4] = 9;
    TopolQTriangle[5] = 2;
    new TPZGeoElRefPattern<pzgeom::TPZQuadraticTrig  > (elementid,TopolQTriangle, matId,*geomesh);
    elementid++;

    TopolQTriangle[0] = 16;
    TopolQTriangle[1] = 10;
    TopolQTriangle[2] = 12;
    TopolQTriangle[3] = 2;
    TopolQTriangle[4] = 11;
    TopolQTriangle[5] = 4;
    new TPZGeoElRefPattern<  pzgeom::TPZQuadraticTrig  > (elementid,TopolQTriangle, matId,*geomesh);
    elementid++;
    
  
    TopolQTriangle[0] = 16;
    TopolQTriangle[1] = 12;
    TopolQTriangle[2] = 14;
    TopolQTriangle[3] = 4;
    TopolQTriangle[4] = 13;
    TopolQTriangle[5] = 6;
    new TPZGeoElRefPattern<  pzgeom::TPZQuadraticTrig  > (elementid,TopolQTriangle, matId,*geomesh);
    elementid++;
 
    
    
    
    // outer arcs bc's
    
    TopolArc[0] = 8;
    TopolArc[1] = 10;
    TopolArc[2] = 9;
    //    TopolArc[0] = 8;
    //    TopolArc[1] = 9;
    //    TopolArc[2] = 10;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fbc1,*geomesh);
    elementid++;
    
    TopolArc[0] = 10;
    TopolArc[1] = 12;
    TopolArc[2] = 11;
    //    TopolArc[0] = 10;
    //    TopolArc[1] = 11;
    //    TopolArc[2] = 12;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fbc1,*geomesh);
    elementid++;
    
    TopolArc[0] = 12;
    TopolArc[1] = 14;
    TopolArc[2] = 13;
    //    TopolArc[0] = 12;
    //    TopolArc[1] = 13;
    //    TopolArc[2] = 14;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fbc1,*geomesh);
    
    
    
    
    geomesh->BuildConnectivity();
    
    int nref = ndiv;
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < nref; iref++) {
        int nel = geomesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = geomesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }
            gel->Divide(sons);
        }
    }
    
    
    std::ofstream out("CircleMixed.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, out, true);
    
    return geomesh;
}

TPZManVector<REAL,3> ParametricCircle(REAL radius,REAL theta)
{
    TPZManVector<REAL,3> xcoor(3,0.0);
    xcoor[0] = radius * cos(theta);
    xcoor[1] = radius * sin(theta);
    xcoor[2] = 0.0 ;
    return xcoor;
}
