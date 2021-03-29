//
//  TPZHDivErrorEstimateMaterial.hpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 11/06/19.
//
/**
 This material implement the Ainsworth proposal
 **/

#ifndef TPZHDivErrorEstimateMaterial_hpp
#define TPZHDivErrorEstimateMaterial_hpp

#include <stdio.h>
#include "mixedpoisson.h"

class TPZHDivErrorEstimateMaterial : public TPZMixedPoisson {

public:
    TPZHDivErrorEstimateMaterial(int matid, int dim);

    TPZHDivErrorEstimateMaterial();
    
    TPZHDivErrorEstimateMaterial(const TPZHDivErrorEstimateMaterial &copy);
    
    TPZHDivErrorEstimateMaterial(const TPZMixedPoisson &copy);
    
    virtual ~TPZHDivErrorEstimateMaterial();
    
    TPZHDivErrorEstimateMaterial &operator=(const TPZHDivErrorEstimateMaterial &copy);

    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef) override;

    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;

    virtual void FillDataRequirements(TPZVec<TPZMaterialData> &datavec) override;
    virtual void FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec) override;

    bool fNeumannLocalProblem = false;

    void SetNeumannProblem(bool neumannProblem) {
        fNeumannLocalProblem = neumannProblem;
    }

    virtual int NEvalErrors() override { return 5; } // erro de oscilacao de dados tbem

    /// Compute the error and error estimate
    // error[0] - error computed with exact pressure
    // error[1] - error computed with reconstructed pressure
    // error[2] - energy error computed with exact solution
    // error[3] - energy error computed with reconstructed solution
    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<REAL> &errors) override;

    void ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact,
                  TPZVec<REAL> &errors, TPZBndCond &bc) override;

    virtual int VariableIndex(const std::string &name) override;
    virtual int NSolutionVariables(int var) override;
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override;

    // Returns the first non-null approximation space index, which will be the
    // H1 reconstruction space
    int FirstNonNullApproxSpaceIndex(TPZVec<TPZMaterialData> &datavec);
};

#endif /* TPZHDivErrorEstimateMaterial_hpp */
