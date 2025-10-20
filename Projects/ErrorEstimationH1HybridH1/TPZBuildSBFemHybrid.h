#ifndef TPZBUILDSBFEMHYBRID_H
#define TPZBUILDSBFEMHYBRID_H

#include "TPZBuildSBFem.h"

class TPZBuildSBFemHybrid : public TPZBuildSBFem
{
public:

    /// @brief  The default constructor
    TPZBuildSBFemHybrid() = default;

    /// simple constructor
    TPZBuildSBFemHybrid(TPZAutoPointer<TPZGeoMesh> & gmesh, int skeletonmatid, std::map<int,int> &matidtranslation) : TPZBuildSBFem(gmesh,skeletonmatid,matidtranslation),
        fInterfaceMaterialIds(0,0) {

    }

    TPZBuildSBFemHybrid(TPZAutoPointer<TPZGeoMesh> & gmesh, int skeleton = 0) : TPZBuildSBFem(gmesh,skeleton), fInterfaceMaterialIds(0,0)
    {
    }

    TPZBuildSBFemHybrid(const TPZBuildSBFemHybrid &copy) : TPZBuildSBFem(copy),
        fInterfaceMaterialIds(copy.fInterfaceMaterialIds),
        fFluxMaterialId(copy.fFluxMaterialId),
        fLagrangePOrder(copy.fLagrangePOrder)
    {
    }

    TPZBuildSBFemHybrid(const TPZBuildSBFem &copy) : TPZBuildSBFem(copy),
        fInterfaceMaterialIds(0,0),
        fFluxMaterialId(0),
        fLagrangePOrder(0)
    {
    }
    // Destructor
    virtual ~TPZBuildSBFemHybrid() override;

    /// @brief Set the interface material IDs
    void SetInterfaceMaterialIds(int matid1, int matid2) {
        fInterfaceMaterialIds = std::make_pair(matid1, matid2);
    }

    /// @brief Set the flux material ID
    void SetFluxMaterialId(int matid) {
        fFluxMaterialId = matid;
    }

    /// @brief Get the flux material ID
    int GetFluxMaterialId() const {
        return fFluxMaterialId;
    }

    /// @brief Get the interface material IDs
    std::pair<int,int> GetInterfaceMaterialIds() const {
        return fInterfaceMaterialIds;
    }


    /// @brief Get the Lagrange multiplier polynomial order
    int GetLagrangePOrder() const {
        return fLagrangePOrder;
    }


    /// @brief Set the Lagrange multiplier polynomial order
    void SetLagrangePOrder(int porder) {
        fLagrangePOrder = porder;
    }

    /// @brief Duplicate the skeleton elements and associate them with scaling centers
    /// this method should be called after the volumetric SBFem elements are created
    void DuplicateSkeletonElements();

    /// @brief Create the skeleton approximation space, one scaling center at a time
    void CreateSkeletonApproximationSpace(TPZCompMesh &cmesh);

    /// create geometric volumetric elements
// the lower dimensional elements already exist (e.g. all connects have been created
    virtual void CreateVolumetricElements(TPZCompMesh &cmesh) override;


    /// @brief Create interface elements
    void CreateInterfaceElements(TPZCompMesh &cmesh);

protected:

    /// @put the SBFem elements into an SBFemElementGroup
    void PutInSBFemGroup(TPZCompMesh &cmesh);
public:
    /// @brief Initialize the Lagrange levels of the connects
    void InitializeLagrangeLevels(TPZCompMesh &cmesh);

    /// @brief group and condense the elements
    void GroupAndCondenseElements(TPZCompMesh &cmesh);

private:
    // Private members specific to TPZBuildSBFemHybrid
    protected:

    /// @brief  The interface material IDs
    std::pair<int,int> fInterfaceMaterialIds = {0,0};

    /// @brief the flux material ID
    int fFluxMaterialId = 0;

    /// @brief The Lagrange multiplier polynomial order
    int fLagrangePOrder = 0;
};

#endif // TPZBUILDSBFEMHYBRID_H
