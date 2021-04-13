//
//  TPZMatLaplacianHybrid.hpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 14/07/19.
//

#ifndef TPZMatLaplacianHybrid_hpp
#define TPZMatLaplacianHybrid_hpp

#include <stdio.h>
#include "TPZMatLaplacian.h"

class TPZMatLaplacianHybrid : public TPZMatLaplacian
{
    
public:
    
    TPZMatLaplacianHybrid(int matid, int dim);
    
    TPZMatLaplacianHybrid(int matid)
    : TPZRegisterClassId(&TPZMatLaplacianHybrid::ClassId), TPZMatLaplacian(matid)
    {
        
    }
    
    TPZMatLaplacianHybrid();
    
    TPZMatLaplacianHybrid(const TPZMatLaplacian &copy);
    
    virtual ~TPZMatLaplacianHybrid();
    
    TPZMatLaplacianHybrid &operator=(const TPZMatLaplacianHybrid &copy);
    
    virtual TPZMaterial *NewMaterial() override;
    
    virtual int VariableIndex(const std::string &name) override;
    int NSolutionVariables(int var)override;
    
    virtual int NEvalErrors()  override {return 4;}
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)override;

    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)override;

    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<REAL> &errors)override;

    void VectorialProd(TPZVec<REAL> & ivec, TPZVec<REAL> & jvec, TPZVec<REAL> & kvec, bool unitary);

    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec) override
    {
        datavec[0].fNeedsNormal = true;
    }
    

    virtual void ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors,TPZBndCond &bc)override{
    }
    
    virtual void ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<REAL> &errors,TPZBndCond &bc)override
    {
    }


    /** @brief Gets the order of the integration rule necessary to integrate an element with polinomial order p */
    ///  HDiv simulations use an additional integration order.
    ///  In order to test if HybridH1 and HDiv are equivalents, both integration orders shall be equivalent.
    virtual int IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const override;

    virtual int ClassId() const override;
    
    
    virtual void Write(TPZStream &buf, int withclassid) const override;
    
    virtual void Read(TPZStream &buf, void *context) override;
    

    
};

#endif /* TPZMatLaplacianHybrid_hpp */
