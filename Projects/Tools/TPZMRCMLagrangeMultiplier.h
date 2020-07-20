//
//  TPZMRCMLagrangeMultiplier.h
//  PZ
//
//  Created by Philippe Devloo on 5/2/14.
//
//

#ifndef __PZ__TPZMRCMLagrangeMultiplier__
#define __PZ__TPZMRCMLagrangeMultiplier__

#include <iostream>
#include "pzdiscgal.h"

/// Material which implements a Lagrange Multiplier
class TPZMRCMLagrangeMultiplier : public TPZDiscontinuousGalerkin
{
    
    /// Number of state variables
    int fNStateVariables;
    
    /// Dimensiona associated with the material
    int fDimension;
    
    STATE fMultiplier;
    
    /// weight factor for the flux penalty
    STATE fBeta;
    
    public :
	/** @brief Simple constructor */
	TPZMRCMLagrangeMultiplier() : TPZRegisterClassId(&TPZMRCMLagrangeMultiplier::ClassId),
    TPZDiscontinuousGalerkin()
    {
        
    }
	/** @brief Constructor with the index of the material object within the vector */
	TPZMRCMLagrangeMultiplier(int nummat, int dimension, int nstate) : TPZRegisterClassId(&TPZMRCMLagrangeMultiplier::ClassId),
    TPZDiscontinuousGalerkin(nummat), fNStateVariables(nstate), fDimension(dimension), fMultiplier(1.)
    {
        
    }
	
	/** @brief Copy constructor */
	TPZMRCMLagrangeMultiplier(const TPZMRCMLagrangeMultiplier &copy) : TPZRegisterClassId(&TPZMRCMLagrangeMultiplier::ClassId),
    TPZDiscontinuousGalerkin(copy), fNStateVariables(copy.fNStateVariables), fDimension(copy.fDimension), fMultiplier(copy.fMultiplier)
    {
        
    }
    
    TPZMRCMLagrangeMultiplier &operator=(const TPZMRCMLagrangeMultiplier &copy)
    {
        TPZDiscontinuousGalerkin::operator=(copy);
        fNStateVariables = copy.fNStateVariables;
        fDimension = copy.fDimension;
        fMultiplier = copy.fMultiplier;
        return *this;
    }
    
    TPZMaterial *NewMaterial() override
    {
        return new TPZMRCMLagrangeMultiplier(*this);
    }
    
	/** @brief Destructor */
	virtual ~TPZMRCMLagrangeMultiplier()
    {
        
    }
    
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const override
    {
        return fDimension;
    }
    
    virtual void SetMultiplier(STATE mult)
    {
        fMultiplier = mult;
    }
    
	
	virtual std::string Name() override
    {
        return "TPZMRCMLagrangeMultiplier";
    }
	
	/**
	 * @brief Fill material data parameter with necessary requirements for the ContributeInterface method.
     * @since April 10, 2007
	 */
	/**
	 * Here, in base class, all requirements are considered as necessary. \n
	 * Each derived class may optimize performance by selecting only the necessary data.
	 */
	virtual void FillDataRequirementsInterface(TPZMaterialData &data) override
    {
        data.SetAllRequirements(false);
        data.fNeedsHSize = true;
    }
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of interface elements */
    virtual void FillDataRequirementsInterface(TPZMaterialData &data, TPZVec<TPZMaterialData > &datavec_left, TPZVec<TPZMaterialData > &datavec_right) override
    {
        data.SetAllRequirements(false);
//        data.fNeedsNormal = true;
        int nref_left = datavec_left.size();
        for(int iref = 0; iref<nref_left; iref++){
            datavec_left[iref].SetAllRequirements(false);
//            datavec_left[iref].fNeedsSol    = true;
//            datavec_left[iref].fNeedsNormal = true;
        }
        int nref_right = datavec_right.size();
        for(int iref = 0; iref<nref_right; iref++){
            datavec_right[iref].SetAllRequirements(false);
//            datavec_right[iref].fNeedsSol    = true;
//            datavec_right[iref].fNeedsNormal = true;
        }
        
    }
	
    /**
     * @{
     * @name Contribute methods
     * @}
     */
    
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override
    {
        DebugStop();
    }
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override
    {
        DebugStop();
    }
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)  override {
        DebugStop();
    }
    
	virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef)  override {
		DebugStop();
	}
	virtual void ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {
		DebugStop();
	}
    
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    
//    virtual void ContributeInterface(TPZVec<TPZMaterialData> &datavec, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec,
//                                     REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
//    {
//        DebugStop();
//    }
	
	/**
	 * @brief Computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since June 5, 2012
	 */
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
	
	/**
	 * @brief It computes a contribution to residual vector at one integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef) override;
	
    
	/**
	 * @brief Computes a contribution to residual vector at one integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ef [out] is the load vector
	 * @since June 5, 2012
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ef) override
    {
        ContributeInterface(data, dataleft[0], dataright[0], weight, ef);
    }
	
    
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override
    {
        DebugStop();
    }
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point to multiphysics simulation
	 * @param data [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since February 21, 2013
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override
    {
        DebugStop();
    }
    
	
	/**
	 * @brief It computes a contribution to residual vector at one BC integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc) override
    {
        DebugStop();
    }
	
    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    void SolutionDisc(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, int var, TPZVec<STATE> &Solout)
    {
        std::cout << __PRETTY_FUNCTION__ << " should never be called\n";
        DebugStop();
    }
	
    /** @} */
	
	/**
	 * @brief Dicontinuous galerkin materials implement contribution of discontinuous elements and interfaces.
	 * @since Feb 05, 2004
	 */
	/**
	 * @brief Computes interface jump = leftu - rightu
	 * @since Feb 14, 2006
	 */
	virtual void InterfaceJump(TPZVec<REAL> &x, TPZSolVec &leftu,TPZSolVec &rightu,TPZSolVec &jump) override
    {
        DebugStop();
    }
	
	
	
	virtual int NStateVariables() const override
    {
        return fNStateVariables;
    }
	
	
	virtual void ContributeInterfaceErrors(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
										   REAL weight,
										   TPZVec<STATE> &nkL,
										   TPZVec<STATE> &nkR,
										   int &errorid)  override {
		PZError << "Method not implemented\n";
	}
	
	virtual void ContributeInterfaceBCErrors(TPZMaterialData &data, TPZMaterialData &dataleft,
											 REAL weight,
											 TPZVec<STATE> &nk,
											 TPZBndCond &bc,
											 int &errorid)  override {
		PZError << "Method not implemented\n";
	}
	
    /** @{
     * @name Save and Load methods
     */
    
	/** @brief Unique identifier for serialization purposes */
	public:
virtual int ClassId() const override;

	
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context) override;
	
    /**
     * @}
     */
};

#endif /* defined(__PZ__TPZMRCMLagrangeMultiplier__) */
