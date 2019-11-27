//
//  ProblemConfig.h
//  ErrorEstimation
//
//  Created by Philippe Devloo on 09/05/18.
//

#ifndef ProblemConfig_h
#define ProblemConfig_h

#include <set>
#include "TPZAnalyticSolution.h"

/// class to guide the error estimator
struct ProblemConfig
{
    /// geometric mesh on which the computational meshes are based
    TPZGeoMesh *gmesh = 0;
    /// polynomial order of the original mesh
    int porder = 1;
    /// increment in internal order of flux and pressure
    int hdivmais = 1;
    /// option to compute the error based on continuous pressures or not
    bool makepressurecontinuous = 0;
    
    /// number of uniform refinements applied to the mesh
    int ndivisions = 1;
    int adaptivityStep = 0;
    int dimension = 0;
    bool prefine = false;
    bool steklovexample = false;
    bool GalvisExample = false;
    bool TensorNonConst = false;
    bool MeshNonConvex = false;
    STATE alpha=1;
    /// directory where the files will be stored
    std::string dir_name = ".";
    /// name identifying the problem
    std::string problemname;
    /// set of materialids in the mesh
    std::set<int> materialids;
    /// set of boundary condition material ids
    std::set<int> bcmaterialids;
    /// exact solution
    TLaplaceExample1 exact;


    ProblemConfig() {};

    ProblemConfig(const ProblemConfig &cp) : gmesh(cp.gmesh),
                                             porder(cp.porder),
                                             hdivmais(cp.hdivmais),
                                             makepressurecontinuous(cp.makepressurecontinuous),
                                             problemname(cp.problemname),
                                             materialids(cp.materialids),
                                             bcmaterialids(cp.bcmaterialids),
                                             exact(cp.exact),
                                             ndivisions(cp.ndivisions),
                                             prefine(cp.prefine),
                                             alpha(cp.alpha),
                                             dir_name(cp.dir_name),
                                             steklovexample(cp.steklovexample),
                                             GalvisExample(cp.GalvisExample),
                                             dimension(cp.dimension),
                                             adaptivityStep(cp.adaptivityStep),
                                            TensorNonConst(cp.TensorNonConst),
                                            MeshNonConvex(cp.MeshNonConvex)
    {
    }

    ProblemConfig &operator=(const ProblemConfig &cp) {
        gmesh = cp.gmesh;
        porder = cp.porder;
        hdivmais = cp.hdivmais;
        makepressurecontinuous = cp.makepressurecontinuous;
        problemname = cp.problemname;
        materialids = cp.materialids;
        bcmaterialids = cp.bcmaterialids;
        exact = cp.exact;
        dimension = cp.dimension;

        ndivisions = cp.ndivisions;
        prefine = cp.prefine;
        adaptivityStep = cp.adaptivityStep;

        alpha = cp.alpha;
        dir_name = cp.dir_name;
        steklovexample = cp.steklovexample;
        GalvisExample = cp.GalvisExample;
        TensorNonConst = cp.TensorNonConst;
        MeshNonConvex = cp.MeshNonConvex;


        return *this;
    }

    // Getters and Setters
    TPZGeoMesh* getGmesh() const { return gmesh; }
    void setGmesh(TPZGeoMesh* Gmesh) { gmesh = Gmesh; }

    int getPorder() const { return porder; }
    void setPorder(int pOrder) { porder = pOrder; }

    int getHdivmais() const { return hdivmais; }
    void setHdivmais(int hDivMais) { hdivmais = hDivMais; }

    bool isMakepressurecontinuous() const { return makepressurecontinuous; }
    void setMakepressurecontinuous(bool makePressureContinuous) { makepressurecontinuous = makePressureContinuous; }

    int getNdivisions() const { return ndivisions; }
    void setNdivisions(int nDivisions) { ndivisions = nDivisions; }

    int getAdaptivityStep() const { return adaptivityStep; }
    void setAdaptivityStep(int adaptivitystep) { adaptivityStep = adaptivitystep; }

    int getDimension() const { return dimension; }
    void setDimension(int Dimension) { dimension = Dimension; }

    bool isPrefine() const { return prefine; }
    void setPrefine(bool pRefine) { prefine = pRefine; }

    bool isSteklovexample() const { return steklovexample; }
    void setSteklovexample(bool steklovExample) { steklovexample = steklovExample; }

    bool isGalvisExample() const { return GalvisExample; }
    void setGalvisExample(bool galvisExample) { GalvisExample = galvisExample; }

    bool isTensorNonConst() const { return TensorNonConst; }
    void setTensorNonConst(bool tensorNonConst) { TensorNonConst = tensorNonConst; }

    bool isMeshNonConvex() const { return MeshNonConvex; }
    void setMeshNonConvex(bool meshNonConvex) { MeshNonConvex = meshNonConvex; }

    STATE getAlpha() const { return alpha; }
    void setAlpha(STATE Alpha) { alpha = Alpha; }

    const string& getDirName() const { return dir_name; }
    void setDirName(const string& dirName) { dir_name = dirName; }

    const string& getProblemname() const { return problemname; }
    void setProblemname(const string& problemName) { problemname = problemName; }

    const set<int>& getMaterialids() const { return materialids; }
    void setMaterialids(const set<int>& materialIDs) { materialids = materialIDs; }

    const set<int>& getBcmaterialids() const { return bcmaterialids; }
    void setBcmaterialids(const set<int>& bcmaterialIDs) { bcmaterialids = bcmaterialIDs; }

    const TLaplaceExample1& getExact() const { return exact; }
    void setExact(const TLaplaceExample1& Exact) { exact = Exact; }
};

#endif /* ProblemConfig_h */
