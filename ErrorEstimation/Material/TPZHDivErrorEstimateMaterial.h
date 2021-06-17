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

#include <cstdio>
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZMaterialDataT.h"

class TPZHDivErrorEstimateMaterial: public virtual TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
        TPZMatErrorCombinedSpaces<STATE>, TPZDarcyFlowInterface>, public TPZMixedDarcyFlow {

public:
    TPZHDivErrorEstimateMaterial(int matid, int dim);

    TPZHDivErrorEstimateMaterial();
    
    TPZHDivErrorEstimateMaterial(const TPZHDivErrorEstimateMaterial &copy);
    
    explicit TPZHDivErrorEstimateMaterial(const TPZMixedDarcyFlow &copy);
    
    ~TPZHDivErrorEstimateMaterial() override;
    
    TPZHDivErrorEstimateMaterial &operator=(const TPZHDivErrorEstimateMaterial &copy);

    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef) override;

    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

    // TODO add doc
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    // TODO add doc
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    bool fNeumannLocalProblem = false;

    void SetNeumannProblem(bool neumannProblem) {
        fNeumannLocalProblem = neumannProblem;
    }

    int NEvalErrors() override { return 5; }

    /// Compute the error and error estimate
    // error[0] - error computed with exact pressure
    // error[1] - error computed with reconstructed pressure
    // error[2] - energy error computed with exact solution
    // error[3] - energy error computed with reconstructed solution
    void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;

    void ErrorsBC(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact,
                  TPZVec<REAL> &errors, TPZBndCond &bc);

    [[nodiscard]] int VariableIndex(const std::string &name) const override;
    [[nodiscard]] int NSolutionVariables(int var) const override;

    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) override;

    // Returns the first non-null approximation space index, which will be the
    // H1 reconstruction space
    static int FirstNonNullApproxSpaceIndex(const TPZVec<TPZMaterialDataT<STATE>> &datavec);
};

#endif /* TPZHDivErrorEstimateMaterial_hpp */
