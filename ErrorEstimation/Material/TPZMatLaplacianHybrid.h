//
//  TPZMatLaplacianHybrid.hpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 14/07/19.
//

#ifndef TPZMatLaplacianHybrid_hpp
#define TPZMatLaplacianHybrid_hpp

#include <stdio.h>
#include "DarcyFlow/TPZDarcyFlow.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"


typedef TPZDarcyFlow TPZMatLaplacian;

class TPZMatLaplacianHybrid : public TPZMatCombinedSpacesT<STATE>, public TPZMatErrorCombinedSpaces<STATE>, public TPZDarcyFlow
{
    
public:
    
    TPZMatLaplacianHybrid(int matid, int dim);
        
    TPZMatLaplacianHybrid();
    
    TPZMatLaplacianHybrid(const TPZMatLaplacian &copy);
    
    virtual ~TPZMatLaplacianHybrid();
    
    TPZMatLaplacianHybrid &operator=(const TPZMatLaplacianHybrid &copy);
    
    virtual TPZMaterial *NewMaterial() const override;
    
    virtual int VariableIndex(const std::string &name) const override;
    int NSolutionVariables(int var) const override;
    
    virtual int NEvalErrors() const override {return 4;}
    
    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE> > &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE> > &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE> > &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc) override;
    
    virtual void Solution(const TPZVec<TPZMaterialDataT<STATE> > &datavec, int var, TPZVec<STATE> &Solout)override;

    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE> > &data, TPZVec<REAL> &errors) override;

    void VectorialProd(TPZVec<REAL> & ivec, TPZVec<REAL> & jvec, TPZVec<REAL> & kvec, bool unitary);

    virtual void FillBoundaryConditionDataRequirements(int type,TPZVec<TPZMaterialDataT<STATE> > &datavec) const override
    {
        datavec[0].fNeedsNormal = true;
    }
    

    virtual void ErrorsBC(TPZVec<TPZMaterialDataT<STATE> > &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors,TPZBndCondT<STATE> &bc){
    }
    
    virtual void ErrorsBC(TPZVec<TPZMaterialDataT<STATE> > &data, TPZVec<REAL> &errors,TPZBndCondT<STATE> &bc)
    {
    }


    /** @brief Gets the order of the integration rule necessary to integrate an element with polinomial order p */
    ///  HDiv simulations use an additional integration order.
    ///  In order to test if HybridH1 and HDiv are equivalents, both integration orders shall be equivalent.
    ///  virtual int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const
    virtual int IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const;

    virtual int ClassId() const override;
    
    
    virtual void Write(TPZStream &buf, int withclassid) const override;
    
    virtual void Read(TPZStream &buf, void *context) override;

    /** @brief Creates an associated boundary condition.
     @param[in] reference The volumetric material associated with the BC.
     @param[in] id Boundary condition identifier.
     @param[in] type Type of the boundary condition.
     @param[in] val1 Value to be set at the element matrix.
     @param[in] val2 Value to be set at the rhs vector.
    */
    virtual TPZBndCondT<STATE>* CreateBC(TPZMaterial *reference,
                                        int id, int type,
                                        const TPZFMatrix<STATE> &val1,
                                        const TPZVec<STATE> &val2) override
    {
        return new  TPZBndCondBase<STATE,TPZMatCombinedSpacesBC<STATE> , TPZMatErrorCombinedSpacesBC<STATE>  >
        (reference,id, type,val1,val2);
    }



};

#endif /* TPZMatLaplacianHybrid_hpp */
