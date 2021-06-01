//
//  ToolsMHM.hpp
//  AdaptivityTest
//
//  Created by Denise De Siqueira on 01/11/19.
//

#ifndef ToolsMHM_hpp
#define ToolsMHM_hpp

#include <stdio.h>

#endif /* ToolsMHM_hpp */


#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <fstream>
#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"
#include "pzreal.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzautopointer.h"

#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "pzfstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "TPZGeoLinear.h"

//#include "TPZMatLaplacianHybrid.h"
#include "TPZLagrangeMultiplier.h"

#include "pzbuildmultiphysicsmesh.h"
#include "pzelementgroup.h"
#include "TPZCompMeshTools.h"
#include "pzcondensedcompel.h"
#include "pzfunction.h"
#include "pzgraphmesh.h"
#include "pzfmatrix.h"

#include "pzlog.h"

#include "TPZVTKGeoMesh.h"
#include "pzvisualmatrix.h"

//#include "pzgengrid.h"

#include "TPZGenGrid2D.h"
#include "TPZExtendGridDimension.h"
#include "pzcheckgeom.h"

#include "TPZMHMeshControl.h"
#include "TPZMHMixedMeshControl.h"
#include "TPZMHMixedMeshChannelControl.h"

#include "TPZMHMixedHybridMeshControl.h"
#include "TPZHybridizeHDiv.h"
//#include "ConfigCasesMaze.h"
#include <iostream>
#include <string>
#include <cmath>
#include <set>

#include "TPZPersistenceManager.h"

#include "TPZMHMHDivErrorEstimator.h"
#include "Tools.h"
#include "TPZGenGrid2D.h"

struct TRunConfig
{
    int nelxcoarse = -1;
    int nelycoarse = -1;
    int numHDivisions = 0;
    int pOrderInternal = 1;
    int hdivmaismais = 1;
    int numDivSkeleton = 0;
    int pOrderSkeleton = 1;
    int Hybridize = 0;
    int Condensed = 1;
    int LagrangeMult = 0;
    int newline = 0;
    int n_threads = 0;
    bool MHM_HDiv_Elast = false;

    /// number of equations when not condensing anything
    int64_t fGlobalSystemSize = -1;
    /// number of equations considering local condensation
    int64_t fGlobalSystemWithLocalCondensationSize = -1;
    /// number of equations of the global system
    int64_t fNumeq = -1;

    REAL fDeltaT = 1.;

    /// number of timesteps
    int64_t nTimeSteps = 10;

    std::ostream &InlinePrint(std::ostream &out)
    {
//        out << "nelxCoarse " << nelxcoarse << " nelyCoarse " << nelycoarse << " numHDiv " << numHDivisions << " porderInternal " << pOrderInternal << " numDivSkeleton " << numDivSkeleton
//        << " porderSkeleton " << pOrderSkeleton << " Hybridize " << Hybridize << " Condensed " << Condensed << " LagrangeMult " << LagrangeMult
//        << " sysnocondense " << fGlobalSystemSize << " syslocalcondense " << fGlobalSystemWithLocalCondensationSize << " neq " << fNumeq;
        // <inputs>: n_dx n_dy k_skeleton m_div k_subelement
        out << "n_dx " << nelxcoarse << " n_dy " << nelycoarse  << " k_skeleton " << pOrderSkeleton << " m_div " << numHDivisions << " k_subelement " << pOrderInternal << " upscaling_dof " << fNumeq << " total_dof " << fGlobalSystemSize;
        return out;
    }
    std::ostream &MathematicaInlinePrint(std::ostream &out)
    {
        out << "nelxCoarse, " << nelxcoarse << ", nelyCoarse, " << nelycoarse << " ,numHDiv, " << numHDivisions << " ,porderInternal, " << pOrderInternal << " ,numDivSkeleton, " << numDivSkeleton
            << " ,porderSkeleton, " << pOrderSkeleton << " ,Hybridize, " << Hybridize << " ,Condensed, " << Condensed << " ,LagrangeMult, " << LagrangeMult
            << " ,sysnocondense, " << fGlobalSystemSize << " ,syslocalcondense, " << fGlobalSystemWithLocalCondensationSize << " ,neq, " << fNumeq;
        return out;
    }

    std::ostream &ConfigPrint(std::ostream &out)
    {
        out << nelxcoarse << "x" << nelycoarse << "_HSkel" << numDivSkeleton << "_pSkel" << pOrderSkeleton << "_HDiv" << numHDivisions << "_pInt" << pOrderInternal;
        return out;
    }
};
// Creating the computational flux mesh
TPZCompMesh *CMeshFlux(TPZGeoMesh * gmesh,int pOrder);

// Create a geometric mesh with the given parameters, nx and ny are the coarse elements number. The total number of elements are defined by the image read.
TPZGeoMesh *GenerateGeoMesh(std::string name, int nx, int ny);

// Create a geoElSide map to be open because it has flux.
// cmesh is the flux mesh obtained from a Hdiv previus simulation.
std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>> IdentifyChanel (TPZCompMesh *cmesh);

// Compute the geometric mesh coarse indices
void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices);

TPZGeoMesh *CreateLMHMMesh(int nDiv, TPZVec<int64_t>& coarseIndexes);

/// Solve the problem composed of a multiphysics mesh composed of compmeshes - applies to MHM and MHM-H(div)
void SolveProblem(const TPZAutoPointer<TPZCompMesh>& cmesh, const TPZVec<TPZAutoPointer<TPZCompMesh> >& compmeshes, TPZAnalyticSolution *analytic, const std::string& prefix, TRunConfig config);
