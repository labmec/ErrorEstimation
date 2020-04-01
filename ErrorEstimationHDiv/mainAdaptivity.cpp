//
// Created by Gustavo Batistela on 12/07/19.
//

#include "TPZRefPatternDataBase.h"
#include "TPZGmshReader.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "TPZGeoLinear.h"

#include <cmath>
#include <tuple>

#include "ProblemConfig.h"
#include "Tools.h"
#include "TPZHybridHDivErrorEstimator.h"
#include "TPZHDivErrorEstimatorH1.h"

//#include "pzelchdiv.h"

bool readGeoMeshFromFile = false;
bool postProcessWithHDiv = false;
int refinementSteps = 13;

void TracingTriangleBug(TPZMultiphysicsCompMesh* multiphysics);

int main(int argc, char* argv[]) {

#ifdef LOG4CXX
    InitializePZLOG();
#endif

    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);

    // Creates geometric mesh

    TPZGeoMesh *gmeshOriginal =
        nullptr; // CreateLShapeMesh(bcIDs);//CreateLShapeMesh(bcIDs);

    ProblemConfig config;

    config.porder = 1;
    config.hdivmais = 1;

    config.dimension = 2;
    config.makepressurecontinuous = true;

    config.exact = new TLaplaceExample1;
    config.exact.operator*().fExact = TLaplaceExample1::ESinSin;

    config.dir_name = "TriangularLShapeMesh";
    config.problemname = "ESinSinMark_UniRef";

    std::string command = "mkdir " + config.dir_name;
    system(command.c_str());

    if (readGeoMeshFromFile) {
        std::string meshfilename = "../LMesh.msh"; //"../LMesh3.msh";

        TPZGmshReader gmsh;
        gmsh.GetDimNamePhysical()[1]["dirichlet"] = 2;
        gmsh.GetDimNamePhysical()[2]["domain"] = 1;
        gmsh.SetFormatVersion("4.1");
        gmeshOriginal = gmsh.GeometricGmshMesh(meshfilename);
        gmsh.PrintPartitionSummary(std::cout);
        config.materialids.insert(1);
        config.bcmaterialids.insert(2);

    }
else

    {
        TPZManVector<int, 4> bcIDs(8, -1);
        gmeshOriginal = CreateLShapeMesh(bcIDs);//CreateGeoMesh(1, bcIDs);//
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        //UniformRefinement(3, gmeshOriginal);
        
    }
    
    gmeshOriginal->SetDimension(config.dimension);
    config.gmesh = gmeshOriginal;
    
    
        
        
    TPZManVector<TPZCompMesh*, 2> meshvec_HDiv(2, 0);
        
    TPZMultiphysicsCompMesh* cmesh_HDiv = CreateHDivMesh(config); //Hdiv x L2
    cmesh_HDiv->InitializeBlock();
    SolveMixedProblem(cmesh_HDiv, config);
    
    meshvec_HDiv = cmesh_HDiv->MeshVector();
    
    
    for (int iSteps = 1; iSteps < refinementSteps; iSteps++) {
    
    
        config.adaptivityStep = iSteps;
        
        UniformRefinement(iSteps, gmeshOriginal);
        
    
        

              TPZManVector<TPZCompMesh*, 2> meshvec_HDiv(2, 0);
              
              TPZMultiphysicsCompMesh* cmesh_HDiv = nullptr;
              
              
              cmesh_HDiv = CreateHDivMesh(config);//Hdiv x L2
              cmesh_HDiv->InitializeBlock();
               #ifdef PZDEBUG2
              {
                  std::ofstream out("MultiPhysicsMesh.txt");
                  cmesh_HDiv->Print(out);
                  std::ofstream outvtk("MultiPhysicsMesh.vtk");
                  
                  TPZVTKGeoMesh::PrintCMeshVTK(cmesh_HDiv,outvtk);
        
                  
              }
              #endif
              
              meshvec_HDiv = cmesh_HDiv->MeshVector();
              
              //cria malha hibrida
              TPZHybridizeHDiv hybrid;
              auto HybridMesh = hybrid.Hybridize(cmesh_HDiv);
              HybridMesh->CleanUpUnconnectedNodes();
              HybridMesh->AdjustBoundaryElements();
              
              delete cmesh_HDiv;
              delete meshvec_HDiv[0];
              delete meshvec_HDiv[1];
              
              cmesh_HDiv = (HybridMesh);//malha hribrida
              meshvec_HDiv[0] = (HybridMesh)->MeshVector()[0];//malha Hdiv
              meshvec_HDiv[1] = (HybridMesh)->MeshVector()[1];//malha L2
              
              SolveHybridProblem(cmesh_HDiv, hybrid.fInterfaceMatid, config,false);
    
    
   
        //reconstroi potencial e calcula o erro
        {
            TPZHybridHDivErrorEstimator HDivEstimate(*cmesh_HDiv);
            HDivEstimate.fProblemConfig = config;
            HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
            HDivEstimate.SetAnalyticSolution(config.exact);
            HDivEstimate.fUpliftPostProcessMesh = config.hdivmais;
            HDivEstimate.fPostProcesswithHDiv = postProcessWithHDiv;
        
            HDivEstimate.PotentialReconstruction();
            TPZManVector<REAL> elementerrors;
            HDivEstimate.ComputeErrors(elementerrors);
            hAdaptivity(&HDivEstimate.fPostProcMesh, gmeshOriginal, config);
    }
    
   // return 0;
        
    }
}

        


// TODO: Turn this into a Unit Test
    void TracingTriangleBug(TPZMultiphysicsCompMesh* multiphysics) {
        
        // Gets H(div) mesh
        TPZCompMesh* fluxMesh = multiphysics->MeshVector()[0];
        
        int64_t nel = fluxMesh->NElements();
        for (int64_t i = 0; i < nel; i++) {
            TPZCompEl* cel = fluxMesh->Element(i);
            
            if (cel->Material()->Id() != 1) continue;
            if (cel->Dimension() != 2) continue;
            
            const auto fluxEl = dynamic_cast<TPZInterpolatedElement*> (cel);
            if (!fluxEl) DebugStop();
            
            // Gets flux element
            const auto gel = fluxEl->Reference();
            
            // Initialize material requirements
            TPZMaterialData elData;
            fluxEl->InitMaterialData(elData);
            
            // Gets last connect, which contains internal shape functions
            TPZConnect& con = cel->Connect(fluxEl->NConnects() - 1);
            const int nInternalPhi = con.NShape();
            
            // Iterates element edges (d-1 sides)
            const int nNodes = gel->NCornerNodes();
            
            // Creates integration rule on edge
            const int pOrderIntRule = 3;
            TPZIntPoints* intRule = gel->CreateSideIntegrationRule(gel->NSides() - 1, pOrderIntRule);
            
            TPZManVector<REAL, 3> xi(2, 0);
            REAL w;
            
            // Stores results of the integration
            int nShapeF = fluxEl->NShapeF();
            TPZFMatrix<REAL> numericalIntegration(nShapeF, 1, 0);
            
            TPZManVector<REAL, 2> gradPhi7(2, 0);
            const int npts = intRule->NPoints();
            for (auto ipt = 0; ipt < npts; ipt++) {
                intRule->Point(ipt, xi, w);
                
                fluxEl->ComputeRequiredData(elData, xi);
                elData.ComputeFunctionDivergence();
                
                if (ipt == 0) {
                    elData.Print(std::cout);
                }
                
                TPZManVector<REAL, 3> phi6(3, 0);
                TPZManVector<REAL, 3> phi7(3, 0);
                REAL div6 = elData.divphi(6, 0);
                REAL div7 = elData.divphi(7, 0);
                
                int vec6Id = elData.fVecShapeIndex[6].first;
                int vec7Id = elData.fVecShapeIndex[7].first;
                
                int phi6Id = elData.fVecShapeIndex[6].second;
                int phi7Id = elData.fVecShapeIndex[7].second;
                
                for (int i = 0; i < 3; i++) {
                    phi6[i] = elData.fMasterDirections(i, vec6Id) * elData.phi(phi6Id, 0);
                    phi7[i] = elData.fMasterDirections(i, vec7Id) * elData.phi(phi7Id, 0);
                }
                
                std::cout << "xi: " << xi << " phi6 = " << phi6 << " phi7 = " << phi7 << " div6 = " << div6
                          << " div7 = " << div7 << std::endl;
                
                gradPhi7[0] += w * elData.dphi(0, phi7Id);
                gradPhi7[1] += w * elData.dphi(1, phi7Id);
                
                std::cout << "gradPhi7: " << elData.dphi(0, phi7Id) << " " << elData.dphi(1, phi7Id) << '\n';
                std::cout << "gradPhi6: " << elData.dphi(0, phi6Id) << " " << elData.dphi(1, phi6Id) << '\n';
                
                for (int iPhi = 0; iPhi < nShapeF; iPhi++) {
                    numericalIntegration(iPhi, 0) += w * elData.divphi(iPhi, 0);
                }
                
            }
            
            for (int i = 0; i < numericalIntegration.Rows(); i++) {
                std::cout << numericalIntegration(i, 0) << '\n';
            }
            
            std::cout << "gradPhiX = " << gradPhi7[0] << " gradPhiY = " << gradPhi7[1] << '\n';
            
            std::cout << '\n';
            
        }
    }
