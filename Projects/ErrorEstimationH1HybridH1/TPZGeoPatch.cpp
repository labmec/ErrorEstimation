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
            fConnectIndexValues(other.fConnectIndexValues), fGeoNodes(other.fGeoNodes), fArea(other.fArea), fCenter(other.fCenter) {
        // Copy data from other to this
    }

    // Move constructor
    TPZGeoPatch::TPZGeoPatch(TPZGeoPatch &&other) noexcept : fElementIndexes(std::move(other.fElementIndexes)),
            fConnectIndexValues(std::move(other.fConnectIndexValues)),
            fGeoNodes(std::move(other.fGeoNodes)),fArea(other.fArea),
            fCenter(std::move(other.fCenter))
{
        // Move data from other to this
                fArea = other.fArea;
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
        fConnectIndexValues = other.fConnectIndexValues;
        fGeoNodes = other.fGeoNodes;
        fArea = other.fArea;
        fCenter = other.fCenter;

        return *this;
    }

    // Move assignment operator
    TPZGeoPatch &TPZGeoPatch::operator=(TPZGeoPatch &&other) noexcept {
        if (this != &other) {
            // Move data from other to this
        }
        fElementIndexes = std::move(other.fElementIndexes);
        fConnectIndexValues = std::move(other.fConnectIndexValues);
        fGeoNodes = std::move(other.fGeoNodes);
        fCenter = std::move(other.fCenter);
        fArea = other.fArea;
        return *this;
    }

void TPZGeoPatch::Print(std::ostream &out) const {
    out << "Contained element indices ";
    for(int64_t it : fElementIndexes) out << it << " ";out << std::endl;
    out << "Connect index and values \n";
    for (auto it : fConnectIndexValues) {
        out << it.first << "|" << it.second << " ";
    }
    out << std::endl;
    out << "Geometric nodes "; for(int64_t node : fGeoNodes) out << node << " ";
    out << std::endl;
    out << "Area " << fArea << " Center " << fCenter << std::endl;
}

void TPZGeoPatch::PlotPatch(const std::string &plotname, int step, TPZMultiphysicsCompMesh &cmesh) {
    int dim = cmesh.Dimension();
    auto &meshvec = cmesh.MeshVector();
    // mesh 4 has the distributed flux value
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
    for(auto it : fConnectIndexValues) {
        auto [cindex, val] = it;
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
    
    TPZFMatrix<STATE> &sol = patch->Solution();
    for(auto it : fConnectIndexValues) {
        auto [cindex, val] = it;
        int64_t seqnum = patch->ConnectVec()[cindex].SequenceNumber();
        int64_t pos = patch->Block().Position(seqnum);
        sol(pos,0) = val;
    }
}

// zero the patch values in the patch mesh
void TPZGeoPatch::ZeroPatchValues(TPZCompMesh *patch) {
    TPZFMatrix<STATE> &sol = patch->Solution();
    for(auto it : fConnectIndexValues) {
        auto [cindex, val] = it;
        int64_t seqnum = patch->ConnectVec()[cindex].SequenceNumber();
        int64_t pos = patch->Block().Position(seqnum);
        sol(pos,0) = 0.;
    }
}

/// @brief this method will expand the patch datastructure to include the data
void TPZGeoPatch::AddConnect(int64_t connectindex, STATE value, int64_t geonodeindex) {
    
#ifdef PZDEBUG
    if(fGeoNodes.find(geonodeindex) != fGeoNodes.end() ||
       fConnectIndexValues.find(connectindex) != fConnectIndexValues.end())
    {
        DebugStop();
    }
#endif
    fConnectIndexValues[connectindex] = value;
    fGeoNodes.insert(geonodeindex);
}

/// Sum two patches
TPZGeoPatch TPZGeoPatch::operator+(const TPZGeoPatch &other) const {
    TPZGeoPatch result(*this);
    result.fElementIndexes.insert(other.fElementIndexes.begin(),other.fElementIndexes.end());
    for(auto itother : other.fConnectIndexValues) {
        auto itfind = result.fConnectIndexValues.find(itother.first);
        if(itfind == result.fConnectIndexValues.end()) {
            result.fConnectIndexValues[itother.first] = itother.second;
        } else {
            itfind->second += itother.second;
        }
    }
    result.fGeoNodes.insert(other.fGeoNodes.begin(),other.fGeoNodes.end());
    return result;
}

TPZGeoPatch &TPZGeoPatch::operator+=(const TPZGeoPatch &other)  {
    fElementIndexes.insert(other.fElementIndexes.begin(),other.fElementIndexes.end());
    for(auto itother : other.fConnectIndexValues) {
        auto itfind = fConnectIndexValues.find(itother.first);
        if(itfind == fConnectIndexValues.end()) {
            fConnectIndexValues[itother.first] = itother.second;
        } else {
            itfind->second += itother.second;
        }
    }
    fGeoNodes.insert(other.fGeoNodes.begin(),other.fGeoNodes.end());
    return *this;
}

void TPZGeoPatch::ComputeAreaCenter(TPZGeoMesh *gmesh, const std::set<int64_t> &elindices, TPZVec<REAL> &center, REAL &area) {
    int dim = gmesh->Dimension();
    // order of the integration rule. SBFem elements have non constant jacobian
    int order = 3;
    area = 0.;
    center.Fill(0.);
    // temporary variables
    TPZManVector<REAL,3> ksi(dim,0.), x(3);
    REAL weight = 0.;
    TPZFNMatrix<9,REAL> jac(dim,dim),jacinv(dim,dim), axes(dim,3);
    REAL detjac;
    
    for(auto it : elindices) {
        TPZGeoEl *gel = gmesh->Element(it);
        if(gel->Dimension() != dim) continue;
        auto intrule = gel->CreateSideIntegrationRule(gel->NSides()-1, order);
        int np = intrule->NPoints();
        for(int ip =0; ip<np; ip++) {
            intrule->Point(ip, ksi, weight);
            gel->Jacobian(ksi, jac, axes, detjac, jacinv);
            gel->X(ksi,x);
            area += fabs(detjac)*weight;
            for(int i=0; i<3; i++) center[i] += x[i]*fabs(detjac)*weight;
        }
        delete intrule;
    }
    if(area != 0.) {
        for(int i=0; i<3; i++) center[i] /= area;
    }
}
