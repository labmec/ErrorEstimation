//
//  MHM_Performance.cpp
//  Evaluates the performance of the MHM mesh builder.
//
//  Created by Gustavo A. Batistela on 17/06/2020.
//

#include "ToolsMHM.h"
#include "boost/date_time/posix_time/posix_time.hpp"

using namespace boost::posix_time;
using time_duration = ptime::time_system_type::time_duration_type;

template <typename... Functions> void TimeMeasurer(const std::string &testName, Functions &&... fs);

void TestMHMPerformance(int nCoarseDivisions, int nInternalRefinements, std::stringstream &results);

int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    std::stringstream results;
    results << ":: MHM Performance Test Results ::\n\n";
    results << "Coarse Grid    Internal Refinements    Elapsed Time\n";

    for (int nCoarseDiv = 1; nCoarseDiv < 8; nCoarseDiv++) {
        for (int nFineRefs = 1; nFineRefs < 4; nFineRefs++) {
            TestMHMPerformance(nCoarseDiv, nFineRefs, results);
        }
        results << '\n';
    }

    std::string outFileName = "MHMPerformanceTest.txt";
    std::ofstream output;
    output.open(outFileName);
    output << results.str();
    output.close();

    return 0;
}

template <typename... Functions> void TimeMeasurer(const std::string &testName, Functions &&... fs) {
    ptime start_time = microsec_clock::local_time();

    auto list = {(fs(), 0)...};

    ptime end_time = microsec_clock::local_time();
    time_duration t = end_time - start_time;
    std::cout << testName << " test duration: " << t << '\n';
}

void TestMHMPerformance(int nCoarseDivisions, int nInternalRefinements, std::stringstream &results) {

    ProblemConfig config;
    config.porder = 2;
    config.hdivmais = 1;
    config.ndivisions = nCoarseDivisions;

    TPZManVector<int, 4> bcids(4, -1);
    TPZGeoMesh *gmesh = CreateGeoMesh(nCoarseDivisions, bcids);

    config.gmesh = gmesh;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);

    {
        std::stringstream fileName;
        fileName << "MHM_Coarse_" << nCoarseDivisions << "x" << nCoarseDivisions << "_InternalRef_"
                 << nInternalRefinements << ".vtk";
        std::ofstream out(fileName.str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }

    UniformRefinement(nInternalRefinements, gmesh);
    DivideLowerDimensionalElements(gmesh);

    {
        std::stringstream fileName;
        fileName << "MHM_Fine_" << nCoarseDivisions << "x" << nCoarseDivisions << "_InternalRef_"
                 << nInternalRefinements << ".vtk";
        std::ofstream out(fileName.str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }

    TPZMultiphysicsCompMesh *cmeshHDiv = CreateHDivMesh(config);
    cmeshHDiv->InitializeBlock(); // TODO: do I need this?

    {
        ptime start_time = microsec_clock::local_time();

        TPZAutoPointer<TPZGeoMesh> gmeshauto(config.gmesh);
        auto *mhm = new TPZMHMixedMeshControl(gmeshauto);

        // This function marks as subdomains every element which does not have a father, i.e., the elements which are
        // not created by a refinement.
        // It works in this case because the geometric mesh is created as a 'n x n' quadrilateral grid.
        // If it was created as a single quadrilateral and then refined 'n' times, the resulting MHM mesh would have
        // only one subdomain.
        TPZVec<int64_t> elementIndexes;
        ComputeCoarseIndices(gmeshauto.operator->(), elementIndexes);

        mhm->DefinePartitionbyCoarseIndices(elementIndexes);

        // Indicate material indices to the MHM control structure
        mhm->fMaterialIds = config.materialids;
        mhm->fMaterialBCIds = config.bcmaterialids;

        // Insert the material objects in the multiphysics mesh
        InsertMaterialObjects(*mhm, config);

        // General approximation order settings
        mhm->SetInternalPOrder(config.porder);
        mhm->SetSkeletonPOrder(config.porder);
        mhm->SetHdivmaismaisPOrder(config.hdivmais);

        // Refine skeleton elements
        mhm->DivideSkeletonElements(nInternalRefinements);
        mhm->DivideBoundarySkeletonElements();

        // Creates MHM mesh
        bool substructure = true;
        mhm->BuildComputationalMesh(substructure);

        ptime end_time = microsec_clock::local_time();
        time_duration t = end_time - start_time;

        results << "        " << nCoarseDivisions << "x" << nCoarseDivisions
        << "                       " << nInternalRefinements << "    " << t << '\n';
    }
}
