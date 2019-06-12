//
//  TPZHDivErrorEstimateMaterial.cpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 11/06/19.
//

#include "TPZHDivErrorEstimateMaterial.h"

TPZHDivErrorEstimateMaterial::TPZHDivErrorEstimateMaterial(int matid, int dim) : TPZMixedPoisson(matid,dim)
{
    
}

TPZHDivErrorEstimateMaterial::TPZHDivErrorEstimateMaterial() : TPZMixedPoisson()
{
    
}

TPZHDivErrorEstimateMaterial::TPZHDivErrorEstimateMaterial(const TPZMixedPoisson &copy) : TPZMixedPoisson(copy)
{
    
}

TPZHDivErrorEstimateMaterial::TPZHDivErrorEstimateMaterial(const TPZHDivErrorEstimateMaterial &copy) : TPZMixedPoisson(copy)
{
    
}

TPZHDivErrorEstimateMaterial::~TPZHDivErrorEstimateMaterial()
{
    
}

TPZHDivErrorEstimateMaterial &TPZHDivErrorEstimateMaterial::operator=(const TPZHDivErrorEstimateMaterial &copy)
{
    TPZMixedPoisson::operator=(copy);
    return *this;
}


void TPZHDivErrorEstimateMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}
