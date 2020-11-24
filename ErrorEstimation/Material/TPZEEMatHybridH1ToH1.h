//
//  TPZHybridH1ErrorEstimateMaterial.hpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 11/06/19.
//


#ifndef ERRORESTIMATION_TPZEEMatHybridH1ToH1_H
#define ERRORESTIMATION_TPZEEMatHybridH1ToH1_H

#include <stdio.h>
#include "TPZMatLaplacianHybrid.h"
#include "TPZEEMatHybridH1ToHDiv.h"

class TPZEEMatHybridH1ToH1: public TPZMatLaplacianHybrid
{

public:


    TPZEEMatHybridH1ToH1(int matid, int dim);

    TPZEEMatHybridH1ToH1();

    TPZEEMatHybridH1ToH1(TPZEEMatHybridH1ToHDiv &hdivEE);

    TPZEEMatHybridH1ToH1(const TPZEEMatHybridH1ToH1 &copy);

    TPZEEMatHybridH1ToH1(const TPZMatLaplacianHybrid &copy);

    virtual ~TPZEEMatHybridH1ToH1();

    TPZEEMatHybridH1ToH1 &operator=(const TPZEEMatHybridH1ToH1 &copy);


    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;


    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override;
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec) override;

    virtual void UpdateBCValues(TPZVec<TPZMaterialData> &datavec);

    /** @brief Gets the order of the integration rule necessary to integrate an element with polinomial order p */
    ///  HDiv simulations use an additional integration order.
    ///  In order to test if HybridH1 and HDiv are equivalents, both integration orders shall be equivalent.
    virtual int IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const override;


    bool fNeumannLocalProblem = true;

    // TODO: Suport sigma errors (return 5)
    virtual int NEvalErrors() override {return 5;}//erro de oscilacao de dados tbem

    /// Compute the error and error estimate
    // error[0] - error computed with exact pressure
    // error[1] - error computed with reconstructed pressure
    // error[2] - energy error computed with exact solution
    // error[3] - energy error computed with reconstructed solution
    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors) override;

    void ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors,TPZBndCond &bc)override;


    virtual int VariableIndex(const std::string &name) override;
    virtual int NSolutionVariables(int var) override;
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var,
                          TPZVec<STATE> &Solout) override;

    // Returns the first non-null approximation space index, which will be the
    // H1 reconstruction space
    int FirstNonNullApproxSpaceIndex(TPZVec<TPZMaterialData> &datavec);
};
#endif //ERRORESTIMATION_TPZEEMatHybridH1ToH1_H
