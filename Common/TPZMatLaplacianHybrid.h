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
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
    

    
    
    virtual int ClassId() const override;
    
    
    virtual void Write(TPZStream &buf, int withclassid) const override;
    
    virtual void Read(TPZStream &buf, void *context) override;
    

    
};

#endif /* TPZMatLaplacianHybrid_hpp */
