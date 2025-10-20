//
//  TPZH1ErrorHybridH1EstimateMaterial.hpp
//  HybridH1vsMixed
//
//  Created by Philippe Devloo on 01/05/25.
//

#ifndef TPZH1ErrorHybridH1EstimateMaterial_hpp
#define TPZH1ErrorHybridH1EstimateMaterial_hpp

#include <stdio.h>
#include "DarcyFlow/TPZHybridDarcyFlow.h"

class TPZH1ErrorHybridH1EstimateMaterial : public TPZHybridDarcyFlow {
    
public:
    enum MMeshPositions {Eflux = 0, Epressure = 1, Epatch = 2, Eorigin = 3, Eaveragepressure = 4};
    
    enum MErrorPositions {EH1, EHybrid, EEstimate, EResidual, EResidualWeighted };
    /**
     * @brief Default constructor
     */
    TPZH1ErrorHybridH1EstimateMaterial() : TPZHybridDarcyFlow() {
        
    }

    /**
     * @brief Class constructor
     * @param [in] id material id
     * @param [in] dim problem dimension
     */
    TPZH1ErrorHybridH1EstimateMaterial(int id, int dim) : TPZHybridDarcyFlow(id,dim) {
        
    }

    TPZH1ErrorHybridH1EstimateMaterial(const TPZH1ErrorHybridH1EstimateMaterial &copy) : TPZHybridDarcyFlow(copy) {

    }
    TPZH1ErrorHybridH1EstimateMaterial(const TPZDarcyFlow &copy) : TPZHybridDarcyFlow(copy) {

    }
    TPZH1ErrorHybridH1EstimateMaterial &operator=(const TPZH1ErrorHybridH1EstimateMaterial &copy){
        TPZHybridDarcyFlow::operator=(copy);
        return *this;
    }
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec)const override;
    
    virtual void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;
    


    /**
     * @brief Returns a 'std::string' with the name of the material
     */
    [[nodiscard]] std::string Name() const override { return "TPZH1ErrorHybridH1EstimateMaterial"; }

    virtual int NEvalErrors()  const override {return 5;}

    virtual void ErrorNames(TPZVec<std::string> &names) const override {
        names[0] = "H1_error";
        names[1] = "Hybrid_error";
        names[2] = "Estimated_error";
        names[3] = "Residual";
        names[4] = "Residual_weighted";
    }
    /** @name Contribute
        @ingroup Contribute*/
    /** @{ */
    /**
     * @brief It computes a contribution to the stiffness matrix
     * and load vector at one integration point.
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     */
    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                            REAL weight,TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef) override;
    /**@}*/
    /** @name ContributeSBFemRhs
        @ingroup Contribute*/
    /** @{ */
    /**
     * @brief It computes a contribution to the stiffness matrix
     * and load vector at one integration point.
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     */
    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                            const TPZFMatrix<CSTATE> &phi, const TPZFMatrix<CSTATE> &dphix,
                            REAL weight,
                            TPZFMatrix<CSTATE> &ef);
    /**@}*/

    /** @name ContributeBC
        @ingroup Contribute*/
    /**@{*/
    /**
     * @brief It computes a contribution to the stiffness matrix
     * and load vector at one BC integration point.
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     * @param[in] bc is the boundary condition material
     */
    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                              REAL weight, TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef,
                              TPZBndCondT<STATE> &bc) override;

    /**@}*/
    /** @brief Returns the solution associated with a given index
        based on the finite element approximation.
        @param[in] datavec Stores all the input data.
        @param[in] var Index of the queried solution
        @param[out] sol FEM Solution at the integration point
    */
    virtual void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                          int var, TPZVec<STATE> &sol) override;
    /**
     * @brief Returns an integer associated with a post-processing variable name
     * @param [in] name string containing the name of the post-processing variable. Ex: "Pressure".
     */
    [[nodiscard]] int VariableIndex(const std::string &name) const override;

    /**
     * @brief Returns an integer with the dimension of a post-processing variable
     * @param [in] var index of the post-processing variable, according to TPZDarcyFlow::VariableIndex method.
     */
    [[nodiscard]] int NSolutionVariables(int var) const override;

    //! @name Error
    /** @{*/
    /*!
      \brief Calculates the error at a given point x.
      \param[in] datavec input data
      \param[out] errors The calculated errors.
     */
    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data,
                        TPZVec<REAL> &errors) override;

    /**
     * @brief Returns an unique class identifier
     */
    [[nodiscard]] int ClassId() const override;
    
    /** @brief Writes this object to the TPZStream buffer. Include the classid if `withclassid = true` */
    virtual void Write(TPZStream &buf, int withclassid) const override
    {
        
    }
    
    /** @brief Reads an objects from the TPZStream buffer. */
    virtual void Read(TPZStream &buf, void *context) override
    {
        
    }


    /**
     * @brief Creates another material of the same type
     */
    [[nodiscard]] TPZMaterial *NewMaterial() const override {
        return new TPZH1ErrorHybridH1EstimateMaterial(*this);
    }

    /**
     * @brief Prints data associated with the material.
     */
    void Print(std::ostream & out) const override {
        TPZHybridDarcyFlow::Print(out);
    }
};

#endif /* TPZH1ErrorHybridH1EstimateMaterial_hpp */
