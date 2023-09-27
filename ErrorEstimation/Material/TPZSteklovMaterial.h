//
//  TPZSteklovMaterial.hpp
//  SteklovProblem
//
//  Created by Denise De Siqueira on 27/11/20.
//
/**
  This material implement the variational form for Steklov problem: find (sigma, u, lambda) such that:

  sigma = -K grad u
  div(sigma) = 0
  sigma.n = lambda u

  |A B^T||sigma| = Inverse(lambda)|C 0||sigma|
  |B 0  ||u    |                  |0 0||u    |

 **/

#ifndef TPZSteklovMaterial_hpp
#define TPZSteklovMaterial_hpp

#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include <cstdio>

class TPZSteklovMaterial : public TPZMixedDarcyFlow {

public:
    TPZSteklovMaterial(int matid, int dim);

    TPZSteklovMaterial();

    TPZSteklovMaterial(const TPZSteklovMaterial &copy);

    explicit TPZSteklovMaterial(const TPZMixedDarcyFlow &copy);

    ~TPZSteklovMaterial() override;

    TPZSteklovMaterial &operator=(const TPZSteklovMaterial &copy);

    void Contribute(const TPZVec<TPZMaterialDataT<STATE> > &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE> > &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override;
};

#endif /* TPZSteklovMaterial_hpp */
