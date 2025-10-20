//
//  TPZGeoPatch.cpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 22/09/25.
//

#include "TPZGeoPatch.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZVTKGenerator.h"

    // Default constructor
    TPZGeoPatch::TPZGeoPatch() {
        // Initialization code here
    }

    // Copy constructor
    TPZGeoPatch::TPZGeoPatch(const TPZGeoPatch &other) : fElementIndexes(other.fElementIndexes),
                                                       fConnectIndexes(other.fConnectIndexes),
                                                       fHatFunctionValues(other.fHatFunctionValues) {
        // Copy data from other to this
    }

    // Move constructor
    TPZGeoPatch::TPZGeoPatch(TPZGeoPatch &&other) noexcept : fElementIndexes(std::move(other.fElementIndexes)),
                                               fConnectIndexes(std::move(other.fConnectIndexes)),
                                               fHatFunctionValues(std::move(other.fHatFunctionValues)) {
        // Move data from other to this
    }

    // Destructor
    TPZGeoPatch::~TPZGeoPatch() {
        // Cleanup code here
    }

    // Assignment operator
    TPZGeoPatch &TPZGeoPatch::operator=(const TPZGeoPatch &other) {
        if (this != &other) {
            // Copy data from other to this
        }
        fElementIndexes = other.fElementIndexes;
        fConnectIndexes = other.fConnectIndexes;
        fHatFunctionValues = other.fHatFunctionValues;

        return *this;
    }

    // Move assignment operator
    TPZGeoPatch &TPZGeoPatch::operator=(TPZGeoPatch &&other) noexcept {
        if (this != &other) {
            // Move data from other to this
        }
        fElementIndexes = std::move(other.fElementIndexes);
        fConnectIndexes = std::move(other.fConnectIndexes);
        fHatFunctionValues = std::move(other.fHatFunctionValues);
        return *this;
    }

void TPZGeoPatch::Print(std::ostream &out) const {
    out << "Contained element indices " << fElementIndexes << std::endl;
    int64_t ncon = fConnectIndexes.size();
    for(int64_t ic = 0; ic<ncon; ic++) {
        out << "con " << fConnectIndexes[ic] << " mult " << fHatFunctionValues[ic] << std::endl;
    }
}

void TPZGeoPatch::PlotPatch(const std::string &plotname, int step, TPZMultiphysicsCompMesh &cmesh) {
    int dim = cmesh.Dimension();
    auto &meshvec = cmesh.MeshVector();
    TPZCompMesh *distflux = meshvec[4];
    {
        TPZFMatrix<STATE> &sol = distflux->Solution();
        sol.Zero();
    }
    TPZGeoMesh *gmesh = distflux->Reference();
    // this is to allow for "loading" the unit value in the computational mesh
    distflux->LoadReferences();
    TPZCompMesh *patch = meshvec[2];
    {
        TPZFMatrix<STATE> &sol = patch->Solution();
        sol.Zero();
    }
    int64_t ncon = fConnectIndexes.size();
    for(int ic = 0; ic<ncon; ic++) {
        int64_t cindex = fConnectIndexes[ic];
        REAL val = fHatFunctionValues[ic];
        int64_t seqnum = patch->ConnectVec()[cindex].SequenceNumber();
        int64_t pos = patch->Block().Position(seqnum);
        TPZFMatrix<STATE> &sol = patch->Solution();
        sol(pos,0) = val;
    }
    // this part is to visualize the elements associated with the patch
    // set value 1 in the discontinuous mesh
    for(auto el : fElementIndexes) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->Dimension() != dim) {
//            std::cout << "gel index " << gel->Index() << " has dimension " << gel->Dimension() << std::endl;
            continue;
        }
        TPZCompEl *cel = gel->Reference();
        if(!cel) {
            std::cout << "gel index " << gel->Index() << " has no computational element\n";
            continue;
        }
        TPZConnect &c = cel->Connect(0);
        int64_t seqnum = c.SequenceNumber();
        int64_t pos = distflux->Block().Position(seqnum);
        TPZFMatrix<STATE> &sol = distflux->Solution();
        sol(pos,0) = 1.;

    }
    TPZStack<std::string> fields;
    fields.Push("Partition");
    fields.Push("DistFlux");
    TPZVTKGenerator vtk(&cmesh, fields, plotname, 0);
    vtk.SetStep(step);
    vtk.Do();
}

// load the patch values in the patch mesh
void TPZGeoPatch::LoadPatchVelues(TPZCompMesh *patch) {
    
    int64_t ncon = fConnectIndexes.size();
    TPZFMatrix<STATE> &sol = patch->Solution();
    for(int ic = 0; ic<ncon; ic++) {
        int64_t cindex = fConnectIndexes[ic];
        REAL val = fHatFunctionValues[ic];
        int64_t seqnum = patch->ConnectVec()[cindex].SequenceNumber();
        int64_t pos = patch->Block().Position(seqnum);
        sol(pos,0) = val;
    }
}

// zero the patch values in the patch mesh
void TPZGeoPatch::ZeroPatchValues(TPZCompMesh *patch) {
    int64_t ncon = fConnectIndexes.size();
    TPZFMatrix<STATE> &sol = patch->Solution();
    for(int ic = 0; ic<ncon; ic++) {
        int64_t cindex = fConnectIndexes[ic];
        REAL val = fHatFunctionValues[ic];
        int64_t seqnum = patch->ConnectVec()[cindex].SequenceNumber();
        int64_t pos = patch->Block().Position(seqnum);
        sol(pos,0) = 0.;
    }
}

/// @brief this method will expand the patch datastructure to include the data
void TPZGeoPatch::AddConnect(int64_t connectindex, STATE value, int64_t geonodeindex) {
    int64_t size = fConnectIndexes.size();
    fConnectIndexes.Resize(size+1, connectindex);
    fHatFunctionValues.Resize(size+1, value);
#ifdef PZDEBUG
    if(fGeoNodes.find(geonodeindex) != fGeoNodes.end())
    {
        DebugStop();
    }
#endif
    fGeoNodes.insert(geonodeindex);
}

