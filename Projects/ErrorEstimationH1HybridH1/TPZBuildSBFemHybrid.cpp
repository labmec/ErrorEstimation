//
//  TPZBuildSBFemHybrid.cpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on [date].
//

#include "TPZBuildSBFemHybrid.h"
#include "TPZMaterial.h"
#include "pzgeoel.h"
#include "pzgeoelbc.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "TPZSBFemVolume.h"
#include "TPZSBFemElementGroup.h"
#include "pzcondensedcompel.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("TPZBuildSBFemHybrid");
#endif

    /// @brief Duplicate the skeleton elements and associate them with scaling centers
    void DuplicateSkeletonElements();

        // Destructor
    TPZBuildSBFemHybrid::~TPZBuildSBFemHybrid() {

    }


void TPZBuildSBFemHybrid::DuplicateSkeletonElements() {
    std::set<int> volmatids;
    for(auto it : fMatIdTranslation) {
        volmatids.insert(it.second);
    }
    std::set<int> skelgeneration = fBoundaryMatIds;
    skelgeneration.insert(fSkeletonMatId);
    // Implementation for duplicating skeleton elements
    int64_t nel = fGMesh->NElements();
    if(fElementPartition.size() < fGMesh->NElements()+10) {
        fElementPartition.Resize(fGMesh->NElements()+50,-1);
    }
    int dim = fGMesh->Dimension();
    for (int64_t i = 0; i < nel; i++) {
        TPZGeoEl *gel = fGMesh->Element(i);
        int matid = gel->MaterialId();
        if (gel && skelgeneration.find(matid) != skelgeneration.end()) {
            if(gel->HasSubElement() || gel->Dimension() != dim-1) {
                continue;
            }
            // gel is a skeleton element
            /// look for neighbouring volumetric elements
            std::list<TPZGeoElSide> volumetric_neighbours;
            TPZGeoElSide skelside(gel);
            for(TPZGeoElSide neighbour = skelside.Neighbour(); neighbour != skelside; neighbour = neighbour.Neighbour()) {
                int neighmatid = neighbour.Element()->MaterialId();
                if(volmatids.find(neighmatid) == volmatids.end()) {
                    continue;
                }
                volumetric_neighbours.push_front(neighbour);
            }
            if(volumetric_neighbours.size() == 0 || volumetric_neighbours.size() > 2) {
                DebugStop();
            }
            auto firstneigh = *volumetric_neighbours.begin();
            int64_t elpartition = fElementPartition[firstneigh.Element()->Index()];
            if(elpartition == -1) {
                std::cout << "For skel" << gel->Index() << " volumetric neighbour " << *volumetric_neighbours.begin() << " has no partition\n";
                DebugStop();
            }
            if(gel->MaterialId() != fSkeletonMatId) {
                // we put the hybrid sequence into the same element partition
                // this implies that when we create the element group and condensed element
                // the boundary will be included
                // we identified a boundary element. Insert a skeleton element
                fElementPartition[gel->Index()] = elpartition;
                if(volumetric_neighbours.size() > 1) DebugStop();
                TPZGeoElBC gbc(firstneigh,fSkeletonMatId);
                fElementPartition[gbc.CreatedElement()->Index()] = elpartition;
                TPZGeoElSide locskelside(gbc.CreatedElement());
                TPZGeoElBC gbc2(locskelside,fInterfaceMaterialIds.first);
                 fElementPartition[gbc2.CreatedElement()->Index()] = elpartition;

            } else {
                // make sure the skeleton element is the first neighbour of the volumetric element
                if(firstneigh.Neighbour() != skelside) {
 //                   std::cout << "Adjusting the neighbouring sequence ";
                    skelside.RemoveConnectivity();
                    firstneigh.InsertConnectivity(skelside);
                    if(firstneigh.Neighbour() != skelside) DebugStop();
                }
                // put the skeleton element into the partition of first volumetric neighbour
                fElementPartition[gel->Index()] = elpartition;
                TPZGeoElBC gbc(skelside,fInterfaceMaterialIds.first);
                 fElementPartition[gbc.CreatedElement()->Index()] = elpartition;
            }
            int64_t gmeshnel = fGMesh->NElements();
            int64_t partsize = fElementPartition.size();
            if(partsize < gmeshnel+10) {
                fElementPartition.Resize(gmeshnel+50,-1);
            }
            if(volumetric_neighbours.size() == 2) {
                // create a flux element
                TPZGeoElSide intfaceneigh = skelside.Neighbour();
                if(intfaceneigh.Element()->MaterialId() != fInterfaceMaterialIds.first) DebugStop();
                // this is where we create the flux element
                TPZGeoElBC flux(intfaceneigh, fFluxMaterialId);
                TPZGeoElSide second_neigh = *volumetric_neighbours.rbegin();
                int64_t second_partition = fElementPartition[second_neigh.Element()->Index()];
                if(second_partition == -1) {
                    DebugStop();
                }
                
                // Duplicate the skeleton element and associate it with a scaling center
                TPZGeoElBC skelclone(second_neigh,fSkeletonMatId);
                TPZGeoEl *newgel = skelclone.CreatedElement();
                fElementPartition[newgel->Index()] = second_partition;
               TPZGeoElBC gbc_new(newgel,fInterfaceMaterialIds.second);
               fElementPartition[gbc_new.CreatedElement()->Index()] = second_partition;
            }
        }
    }
    fElementPartition.Resize(fGMesh->NElements(), -1);
}

/// @brief Create the skeleton approximation space, one scaling center at a time
void TPZBuildSBFemHybrid::CreateSkeletonApproximationSpace(TPZCompMesh &cmesh) {
    // Implementation for creating the skeleton approximation space
    cmesh.SetAllCreateFunctionsContinuous();
    cmesh.SetDefaultOrder(fSkeletonPOrder);
    int dim = fGMesh->Dimension();
    
    std::map<int64_t,std::set<int64_t>> center_to_boundary;
    int64_t nel = fGMesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (gel && gel->Dimension() == dim-1) {
            int64_t elpartition = fElementPartition[el];
            int gelmatid = gel->MaterialId();
            if(elpartition != -1 && gelmatid == fSkeletonMatId) {
                center_to_boundary[elpartition].insert(gel->Index());
            }
        }
    }

    cmesh.Reference()->ResetReference();
    for(auto it : center_to_boundary) {
        int64_t center = it.first;
        std::set<int64_t> &boundary = it.second;
        for(auto index : boundary) {
            cmesh.CreateCompEl(fGMesh->Element(index));
        }
        for(auto index : boundary) {
            fGMesh->Element(index)->ResetReference();
        }
    }
    // create the HDiv lagrange multipliers
    cmesh.SetAllCreateFunctionsHDiv();
    cmesh.SetDefaultOrder(fLagrangePOrder);
    std::set<int> lagrangematids(fBoundaryMatIds);
    lagrangematids.insert(fFluxMaterialId);
    cmesh.AutoBuild(lagrangematids);
}

void TPZBuildSBFemHybrid::CreateInterfaceElements(TPZCompMesh &cmesh) {
    // Implementation for creating interface elements
    std::set<int> fluxandbound(fBoundaryMatIds);
    fluxandbound.insert(fFluxMaterialId);
    int dim = cmesh.Dimension();
    TPZGeoMesh *gmesh = cmesh.Reference();
    gmesh->ResetReference();
    cmesh.LoadReferences();
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || gel->Dimension() != dim-1 || gel->MaterialId() != fSkeletonMatId || gel->HasSubElement()) {
            continue;
        }
        int64_t elpartition = fElementPartition[el];
        if (elpartition == -1) {
            DebugStop();
        }
        TPZGeoElSide skelside(gel);
        TPZGeoElSide interface;
        for(TPZGeoElSide neighbour = skelside.Neighbour(); neighbour != skelside; neighbour = neighbour.Neighbour()) {
            int neighmatid = neighbour.Element()->MaterialId();
            if(neighmatid != fInterfaceMaterialIds.first && neighmatid != fInterfaceMaterialIds.second) {
                continue;
            }
            int64_t neighpartition = fElementPartition[neighbour.Element()->Index()];
            if(neighpartition == elpartition) {
                interface = neighbour;
                break;
            }
        }
        if (!interface) {
            DebugStop();
        }
        TPZGeoElSide fluxside = skelside.HasNeighbour(fluxandbound);
        if (!fluxside) {
            DebugStop();
        }
        auto cskelside = skelside.Reference();
        auto cfluxside = fluxside.Reference();
        auto cel = new TPZInterfaceElement(cmesh, interface.Element(), cskelside, cfluxside);

    }
}

/// create geometric volumetric elements
// the lower dimensional elements already exist (e.g. all connects have been created
void TPZBuildSBFemHybrid::CreateVolumetricElements(TPZCompMesh &cmesh) {
    TPZGeoMesh *gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    cmesh.LoadReferences();
    int dim = gmesh->Dimension();
        // all computational elements have been loaded
    std::set<int> matidstarget;
    for (std::map<int,int>::iterator it = fMatIdTranslation.begin(); it!= fMatIdTranslation.end(); it++) {
        int64_t mat = it->second;
        if (cmesh.FindMaterial(mat)) {
            matidstarget.insert(it->second);
        }
    }
    cmesh.ApproxSpace().SetAllCreateFunctionsSBFem(dim);
    cmesh.AutoBuild(matidstarget);
    PutInSBFemGroup(cmesh);
}

/// @put the SBFem elements into an SBFemElementGroup
void TPZBuildSBFemHybrid::PutInSBFemGroup(TPZCompMesh &cmesh) {
    int64_t numgroups = fPartitionCenterNode.size();
    int64_t groupelementindices(numgroups);

    if(fPOrderBubbleFunctions > 0) {
        TPZSBFemElementGroup::SetDefaultPolynomialOrder(fPOrderBubbleFunctions);
        // TPZSBFemElementGroup::gPolynomialShapeFunctions = true;
    }
    TPZVec<int64_t> elementgroupindices(numgroups);
    for (int64_t el=0; el<numgroups; el++) {
        TPZSBFemElementGroup* cel = new TPZSBFemElementGroup(cmesh);
        elementgroupindices[el] = cel->Index();
    }
    int64_t nel = cmesh.NElements();
    int dim = cmesh.Dimension();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh.Element(el);
        if (!cel) {
            continue;
        }
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
        if (sbfem) {
            TPZGeoEl *sbgel = sbfem->Reference();
            int geldim = sbgel->Dimension();
            int64_t partition = fElementPartition[sbgel->Index()];
            if(partition == -1) DebugStop();
            int side = sbgel->FirstSide(geldim-1);
            TPZGeoElSide gelside(sbgel,side);
            TPZGeoElSide skelside = gelside.HasNeighbour(fSkeletonMatId);
            if(!skelside) DebugStop();
            int64_t skelpartition = fElementPartition[skelside.Element()->Index()];
            if(skelpartition != partition) {
                skelside++;
                skelside = skelside.HasNeighbour(fSkeletonMatId);
                skelpartition = fElementPartition[skelside.Element()->Index()];
                if(skelpartition != partition) DebugStop();
            }
            TPZCompEl *cskel = skelside.Element()->Reference();
            sbfem->SetSkeleton(cskel->Index());
            int64_t celgroupindex = elementgroupindices[partition];
            TPZCompEl *celgr = cmesh.Element(celgroupindex);
            TPZSBFemElementGroup *sbfemgr = dynamic_cast<TPZSBFemElementGroup *>(celgr);
            if (!sbfemgr) {
                DebugStop();
            }
            sbfemgr->AddElement(sbfem);
        }
    }
    for (int64_t el=0; el<numgroups; el++) {
        int64_t index;
        
        index = elementgroupindices[el];
        TPZCompEl *cel = cmesh.Element(index);
        TPZSBFemElementGroup *sbfemgroup = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (!sbfemgroup) {
            DebugStop();
        }
        const TPZVec<TPZCompEl *> &subgr = sbfemgroup->GetElGroup();
        int64_t nsub = subgr.NElements();
        for (int64_t is=0; is<nsub; is++) {
            TPZCompEl *cel = subgr[is];
            TPZSBFemVolume *femvol = dynamic_cast<TPZSBFemVolume *>(cel);
            if (!femvol) {
                DebugStop();
            }
            femvol->SetElementGroupIndex(index);
        }
        if (nsub == 0)
        {
            delete sbfemgroup;
        }
    }
    if(fPOrderBubbleFunctions > 0) {
        for (int64_t el=0; el<numgroups; el++) {
            int64_t index;
            
            index = elementgroupindices[el];
            TPZCompEl *cel = cmesh.Element(index);
            TPZSBFemElementGroup *sbfemgroup = dynamic_cast<TPZSBFemElementGroup *>(cel);
            if (!sbfemgroup) {
                DebugStop();
            }
            sbfemgroup->InitializeInternalConnect();
        }
    }

}

/// @brief Initializet he Lagrange levels of the connects
void TPZBuildSBFemHybrid::InitializeLagrangeLevels(TPZCompMesh &cmesh) {
    // for each element, all displacement lagrange levels will be 1 except one corner node which will be 2
    // the flux dofs will have level 1

    int64_t npartitions = fElementPartition.size();
    TPZVec<int64_t> partitionconnect(npartitions,-1);
    std::set<int> matids = fBoundaryMatIds;
    matids.insert(fFluxMaterialId);
    int64_t nel = cmesh.NElements();
    for (int64_t el = 0; el < nel; el++)
    {
        TPZCompEl *cel = cmesh.Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) {
            TPZSBFemElementGroup *sbfem = dynamic_cast<TPZSBFemElementGroup *>(cel);
            if(!sbfem) {
                DebugStop();
            }
            int nconnect = sbfem->NConnects();
            for (int ic = 0; ic < nconnect; ic++) {
                TPZConnect &c = sbfem->Connect(ic);
                c.SetLagrangeMultiplier(0);
            }
            sbfem->Connect(0).SetLagrangeMultiplier(2);
        } else {
            int gelmatid = gel->MaterialId();
            if(matids.find(gelmatid) == matids.end()) {
                continue;
            }
            int nconnect = cel->NConnects();
            for (int ic = 0; ic < nconnect; ic++) {
                TPZConnect &c = cel->Connect(ic);
                c.SetLagrangeMultiplier(1);
            }
        }
    }
}

/// @brief group and condense the elements
void TPZBuildSBFemHybrid::GroupAndCondenseElements(TPZCompMesh &cmesh) {
    // Implementation of the grouping and condensation of elements
    int64_t nelgroup = this->fPartitionCenterNode.size();
    TPZVec<TPZElementGroup *> groupvec(nelgroup);
    for(int64_t ig = 0; ig < nelgroup; ig++) {
        TPZElementGroup *elgr = new TPZElementGroup(cmesh);
        groupvec[ig] = elgr;
    }
    int64_t nel = cmesh.NElements();
    for (int64_t el = 0; el < nel; el++) {
        int partition = -1;
        TPZCompEl *cel = cmesh.Element(el);
        if (!cel) continue;
        // Group elements based on some criteria
        TPZSBFemElementGroup *sbfem = dynamic_cast<TPZSBFemElementGroup *>(cel);
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
        if (sbfem) {
            int64_t partition = GetPartition(sbfem);
            groupvec[partition]->AddElement(sbfem);
            continue;
        } 
        if (elgr) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            DebugStop();
        }
        int64_t groupindex = fElementPartition[gel->Index()];
        if(groupindex != -1) {
            if(groupindex < 0 || groupindex >= groupvec.size()) {
                std::cout << "gel index " << gel->Index() << " has group index " << groupindex << std::endl;
                std::cout << fElementPartition << std::endl;
            }
            groupvec[groupindex]->AddElement(cel);
        }
    }
    // increment nelconnected of connect with lagrange level 2
    cmesh.ComputeNodElCon();
    int64_t nconnects = cmesh.NConnects();
    for (int64_t ic = 0; ic < nconnects; ic++) {
        TPZConnect &c = cmesh.ConnectVec()[ic];
        if (c.LagrangeMultiplier() == 2) {
            c.IncrementElConnected();
        }
    }
    // Condense each group
    for (int64_t ig = 0; ig < nelgroup; ig++) {
        TPZElementGroup *elgr = groupvec[ig];
        if (!elgr) DebugStop();
        TPZCondensedCompElT<STATE> *condensed = new TPZCondensedCompElT<STATE>(elgr);
    }
    cmesh.ComputeNodElCon();
    cmesh.CleanUpUnconnectedNodes();
}

