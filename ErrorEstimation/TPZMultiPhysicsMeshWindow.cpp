#include "TPZMultiPhysicsMeshWindow.h"

TPZMultiPhysicsMeshWindow::TPZMultiPhysicsMeshWindow() : TPZMultiphysicsCompMesh() {
    // Default constructor implementation
}

TPZMultiPhysicsMeshWindow::TPZMultiPhysicsMeshWindow(const TPZMultiPhysicsMeshWindow &other) : TPZMultiphysicsCompMesh(other), m_connect_correspondence(other.m_connect_correspondence),
    m_Referred(other.m_Referred)
{
    // Copy constructor implementation
}

TPZMultiPhysicsMeshWindow& TPZMultiPhysicsMeshWindow::operator=(const TPZMultiPhysicsMeshWindow &other) {
    if (this != &other) {
        TPZCompMesh::operator=(other);
        // Copy any additional members here if needed
    }
    m_connect_correspondence = other.m_connect_correspondence;
    m_Referred = other.m_Referred;
    return *this;
}

TPZMultiPhysicsMeshWindow::~TPZMultiPhysicsMeshWindow() {
    // Destructor implementation
}

/// @brief Build the multiphysics space for the corresponding geometric elements
void TPZMultiPhysicsMeshWindow::BuildMultiphysicsSpace(const TPZVec<int64_t> &gelindexes) {
    // Implementation of building the multiphysics space
    SetAllCreateFunctionsMultiphysicElem();
    /// this will create a multiphysics element for each geometric element in gelindexes for which a material is defined
    fCreate.BuildMesh(*this, gelindexes);
    AddElements();
    AddConnects();
}
/// @brief Add the connects from the atomic meshes
void TPZMultiPhysicsMeshWindow::AddConnects() {
    // Implementation of adding connects from atomic meshes
    int64_t NMPhysicsConnects = this->NConnects();
    if(NMPhysicsConnects != 0) {
        std::cerr << "Warning: The multiphysics mesh already has connects. They should be cleared before adding new connects." << std::endl;
        DebugStop();
    }
    int64_t NMPhysicsElements = this->NElements();
    if(NMPhysicsElements == 0) {
        std::cerr << "Warning: The multiphysics mesh has no elements." << std::endl;
        DebugStop();
    }
    TPZGeoMesh *gmesh = Reference();
    if (!gmesh) {
        std::cerr << "Error: No geometric mesh associated with the multiphysics mesh." << std::endl;
        DebugStop();
    }
    // gmesh->ResetReference();
    auto meshvec = this->MeshVector();
    auto active_spaces = this->GetActiveApproximationSpaces();
    TPZManVector<std::set<int64_t>, 7> connectIndicesPerSpace(meshvec.size());
    for(int imesh = 0; imesh < meshvec.size(); imesh++) {
        if(active_spaces[imesh] == 0) continue;
        TPZCompMesh *cmesh = meshvec[imesh];
        if (!cmesh) continue;
        for (int64_t el = 0; el < NMPhysicsElements; el++) {
            TPZCompEl *cel = this->Element(el);
            TPZMultiphysicsElement *mcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if (!mcel) continue;
            TPZCompEl *celatomic = mcel->Element(imesh);
            if(!celatomic) continue;
            std::set<int64_t> connectlist;
            celatomic->BuildConnectList(connectlist);
            for (auto cindex : connectlist) {
                connectIndicesPerSpace[imesh].insert(cindex);
            }
            
        }
    }
    TPZManVector<std::map<int64_t,int64_t>, 7> connectMapPerSpace(meshvec.size());
    m_connect_correspondence.Resize(meshvec.size());
    for(int imesh = 0; imesh < meshvec.size(); imesh++) {
        if(active_spaces[imesh] == 0) continue;
        int64_t count = 0;
        int64_t nconnects = connectIndicesPerSpace[imesh].size();
        m_connect_correspondence[imesh].resize(nconnects);
        for (auto it = connectIndicesPerSpace[imesh].begin(); it != connectIndicesPerSpace[imesh].end(); it++) {
            m_connect_correspondence[imesh][count] = *it;
            count++;
        }
        BuildConnectMap(imesh, connectMapPerSpace[imesh]);
    }
    /// Create the connects and copy their data structure
    int64_t seqnum = 0;
    for (int i_as = 0; i_as < meshvec.size(); i_as++)
    {
        if (active_spaces[i_as] == 0) {
            continue;
        }
        TPZCompMesh *cmesh = meshvec[i_as];
        for(auto it : m_connect_correspondence[i_as]) {
            TPZConnect &c = cmesh->ConnectVec()[it];
            int64_t new_index = this->AllocateNewConnect(c);
            TPZConnect &cnew = this->ConnectVec()[new_index];
            cnew.CopyFrom(c, connectMapPerSpace[i_as]);
            cnew.SetCondensed(false);
            cnew.SetSequenceNumber(seqnum++);
//            std::cout << "Copying connect " << it << " from mesh " << i_as << " to new connect " << new_index << std::endl;
        }
    }
    for (int64_t el = 0; el < NMPhysicsElements; el++) {
        TPZCompEl *cel = this->Element(el);
        TPZMultiphysicsElement *mcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mcel) continue;
        int nc = 0;
        for(int imesh = 0; imesh < meshvec.size(); imesh++) {
            if(active_spaces[imesh] == 0) continue;
            TPZCompEl *celatomic = mcel->Element(imesh);
            if(!celatomic) continue;
            int ncon = celatomic->NConnects();
            nc += ncon;
        }
        TPZManVector<int64_t> newindices(nc);
        int count = 0;
        for(int imesh = 0; imesh < meshvec.size(); imesh++) {
            if(active_spaces[imesh] == 0) continue;
            TPZCompEl *celatomic = mcel->Element(imesh);
            if(!celatomic) continue;
            int ncon = celatomic->NConnects();
            for (int icon = 0; icon < ncon; icon++) {
                int64_t cindex = celatomic->ConnectIndex(icon);
                if(connectMapPerSpace[imesh].find(cindex) == connectMapPerSpace[imesh].end()) DebugStop();
                int64_t newindex = connectMapPerSpace[imesh][cindex];
                newindices[count++] = newindex;
            }
        }
        mcel->SetConnectIndexes(newindices);
    }

    InitializeBlock();
}


/// @brief build the map between connect indices of the atomic mesh and the multiphysics mesh
void TPZMultiPhysicsMeshWindow::BuildConnectMap(int imesh, std::map<int64_t, int64_t> &connectmap) {
    connectmap.clear();
    if(imesh < 0 || imesh >= m_connect_correspondence.size()) {
        std::cerr << "Error: Invalid mesh index in BuildConnectMap." << std::endl;
        DebugStop();
    }
    int64_t firstconnect = 0;
    for (int im = 0; im<imesh; im++) firstconnect += m_connect_correspondence[im].size();
    for (int64_t i = 0; i < m_connect_correspondence[imesh].size(); i++) {
        int64_t atomic_connect_index = m_connect_correspondence[imesh][i];
        connectmap[atomic_connect_index] = i+firstconnect;
    }
}

    /// @brief  Add the atomic elements to the multiphysics elements
void TPZMultiPhysicsMeshWindow::AddElements() {
    TPZGeoMesh * geometry = Reference();
//    geometry->ResetReference();
    int64_t n_cels = NElements();
    // for each geometric element, the computational multiphysics element
    int64_t n_gels = geometry->NElements();
//    TPZVec<TPZCompEl *> Referred(n_gels);
    auto meshvec = this->MeshVector();
    auto m_active_approx_spaces = this->GetActiveApproximationSpaces();
    int n_approx_spaces = meshvec.size();
    if(m_Referred.size() != n_approx_spaces) {
        if(m_Referred.size()) DebugStop();
        m_Referred.Resize(n_approx_spaces);
        for(int iappr = 0; iappr < n_approx_spaces; iappr++) {
            m_Referred[iappr].Resize(n_gels, 0);
            if(meshvec[iappr]) LoadReferred(meshvec[iappr], m_Referred[iappr]);
        }
    }
    for(int i_as = 0; i_as < n_approx_spaces; i_as++)
    {
        /// for a given atomic space, load the references
        TPZCompMesh *atom = meshvec[i_as];
        if(!atom) continue;
//        Referred.Fill(0);
//        LoadReferred(atom, Referred);
        // atom->LoadReferences(Referred);
        int64_t icel;
        // loop over the multiphysics elements
        for(icel=0; icel < n_cels; icel++)
        {
            TPZCompEl * cel = ElementVec()[icel];
            TPZMultiphysicsElement * mfcel = dynamic_cast<TPZMultiphysicsElement *> (cel);
            if(mfcel)
            {
                int64_t found = 0;
                int64_t gelindex = mfcel->ReferenceIndex();
                TPZCompEl *celatom = m_Referred[i_as][gelindex];
                
                if (celatom) {
                    mfcel->AddElement(celatom, i_as);
                    continue;
                }
                else
                {
                    // look for a multiphysics element in the ancestral tree
                    TPZGeoEl *gel = geometry->Element(gelindex);
                    TPZGeoEl *gelF = gel;
                    while(gelF->Father())
                    {
                        gelF = gelF->Father();
                        int gelFindex = gelF->Index();
                        if (m_Referred[i_as][gelFindex]) {
#ifdef PZDEBUG
                            if (gelF->MaterialId() != gel->MaterialId()) {
                                DebugStop();
                            }
#endif
                            mfcel->AddElement(m_Referred[i_as][gelFindex], i_as);
                            found = true;
                            break;
                        }
                    }
                }
                if (!found) {
                    mfcel->AddElement(0, i_as);
                }
            }
            else {
                DebugStop();
            }
        }
    }
    
    for (int64_t icel = 0; icel < n_cels; icel++) {
        TPZCompEl *cel = Element(icel);
        TPZMultiphysicsElement *mfel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfel) {
            continue;
        }
#ifdef PZDEBUG
        {
            int ncontained = 0;
            for(int i=0; i<n_approx_spaces; i++)
            {
                TPZCompEl *atcel = mfel->Element(i);
                if(atcel) ncontained++;
            }
            if(ncontained == 0)
            {
                TPZGeoEl *gel = cel->Reference();
                std::cout << "Multiphysics element " << icel << " with matid " <<
                    gel->MaterialId() << " does not refer to any elements "
                << " geometric index " << gel->Index() << std::endl;
                TPZCompMesh *flux = meshvec[0];
                fReference->ResetReference();
                flux->LoadReferences();
                std::cout << "Flux reference " << (void *) gel->Reference() << std::endl;
                DebugStop();
            }
        }
#endif
        mfel->SetActiveApproxSpaces(m_active_approx_spaces);
        mfel->InitializeIntegrationRule();
    }
    
}
