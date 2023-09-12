
#include "common_files.h"
#include "TPZLinearAnalysis.h"
#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "pzcmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "DarcyFlow/TPZHybridDarcyFlow.h"
#include "TPZAnalyticSolution.h"
#include "TPZMaterial.h"
#include "TPZMHMeshControl.h"


void SolveProblem(TPZMultiphysicsCompMesh &cmesh, std::stringstream &rootname)
{
    int numthreads = 1;
    bool optimizeBandwidth = true;
    bool plotting = true;
    TPZLinearAnalysis an(&cmesh, optimizeBandwidth); //Creates the object that will manage the analysis of the problem
#ifdef PZ_USING_MKL
    TPZSymetricSpStructMatrix matskl(cmesh);
#else
    TPZSkylineStructMatrix<STATE> matskl(&cmesh); // asymmetric case ***
#endif
    matskl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    
    std::cout << "Assemble matrix with NDoF = " << cmesh.NEquations() << "." << std::endl;
    an.Assemble(); //Assembles the global stiffness matrix (and load vector)
    std::cout << "Assemble finished." << std::endl;
    
    //            an.PostProcessError(Errors,std::cout);
    
#ifdef PZDEBUG
    //Imprimir Matriz de rigidez Global:
    if (false) {
        std::ofstream filestiff("stiffness.nb");
        an.MatrixSolver<STATE>().Matrix()->Print("K1 = ", filestiff, EMathematicaInput);
        
        std::ofstream filerhs("rhs.nb");
        an.Rhs().Print("R = ", filerhs, EMathematicaInput);
    }
#endif
    
    std::cout << "Solving." << std::endl;
    an.Solve();
    std::cout << "Solved." << std::endl;
    
    
    {
        TPZStepSolver<STATE> solver;
        an.SetSolver(solver);
    }
#ifdef PZDEBUG
    if (0) {
        std::ofstream file("file.txt");
        an.Solution().Print("sol=", file, EMathematicaInput);
        
    }
#endif
    
    
    if (plotting) {
        std::string plotfile;
        {
            std::stringstream sout;
            sout << rootname.str() << ".vtk";
            plotfile = sout.str();
        }
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("SigmaX");
        scalnames.Push("SigmaY");
        scalnames.Push("TauXY");
        vecnames.Push("displacement");
        static int count = 0;
        an.SetStep(count);
        an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
        an.PostProcess(2);
        count++;
    }
    
#ifdef PZDEBUG
    //Imprimindo vetor solução:
    {
        TPZFMatrix<STATE> solucao = cmesh.Solution(); //Pegando o vetor de solução, alphaj
        std::ofstream solout("sol.nb");
        solucao.Print("Sol", solout, EMathematicaInput); //Imprime na formatação do Mathematica
        
        std::ofstream fileAlpha("alpha.nb");
        an.Solution().Print("Alpha = ", fileAlpha, EMathematicaInput);
    }
#endif
    
}
//   matids.clear();
//   matids.insert(-1);
//   TPZManVector<STATE,3> result;
//  result = cmesh_m_HDiv->Integrate("state",matids);
//  std::cout << "Sigma Y"  << result << std::endl;


//    //Calculo do erro
//    std::cout << "Computing Error " << std::endl;
void ComputeError(TPZCompMesh &cmesh, std::stringstream &rootname, int href, int pref, TPZMHMeshControl &mhm)
{
    int numthreads = 1;
    int nelx = 1 << href;
    TPZMaterial *mat = cmesh.FindMaterial(matID);
    TPZMatErrorCombinedSpaces<STATE> *materr = dynamic_cast<TPZMatErrorCombinedSpaces<STATE> *>(mat);
    TPZManVector<REAL, 6> Errors(materr->NEvalErrors());

    std::stringstream sout;
    sout << rootname.str();
    sout << "_" << mhm.fpOrderSkeleton << "_Error.nb";
    std::ofstream ErroOut(sout.str(), std::ios::app);
    ErroOut << "(* Type of simulation " << rootname.str() << " *)\n";
    ErroOut << "(* Number of elements " << nelx << " *)" << std::endl;
    ErroOut << "(* Type of Element ";

    ErroOut << " *)\n";
    ErroOut << "(* Number of Condensed equations " << mhm.fGlobalSystemWithLocalCondensationSize << " *)" << std::endl;
    ErroOut << "(* Number of equations before condensation " << mhm.fGlobalSystemSize << " *)" << std::endl;
    ErroOut << "(*\n";
    TPZLinearAnalysis an(&cmesh,false);
    an.SetThreadsForError(numthreads);
    bool store_errors = true;
    cmesh.ExpandElementSolution(Errors.size());
    std::cout << "Computing errors." << std::endl;
    an.PostProcessError(Errors, store_errors, ErroOut);
    std::cout << "Computed errors." << std::endl;
    ErroOut << "nelx ribporder internalporder n_condensed - n_total - error_sigma - error_energy - error_div_sigma - error_u - error_r - error_as\n";
    ErroOut << "*)\n";
    TPZManVector<STATE, 10> output(Errors.size() + 5, 0);
    output[0] = nelx;
    output[1] = mhm.fpOrderSkeleton;
    output[2] = mhm.fpOrderInternal;
    output[3] = cmesh.NEquations();
    output[4] = mhm.fGlobalSystemSize;
    for (int i = 0; i < Errors.size(); i++) {
        output[5 + i] = Errors[i];
    }
    ErroOut << "Error[[" << href + 1 << "," << pref + 1 << "]] = {" << output << "};\n";
    
    std::cout << "Errors = " << Errors << std::endl;
}

