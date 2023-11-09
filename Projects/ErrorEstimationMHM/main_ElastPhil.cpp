//
// Created by Gustavo Batistela on 3/31/21.
//

#include "pzgmesh.h"
#include <Pre/TPZGenGrid3D.h>
#include <Pre/TPZMHMixedMeshControl.h>
#include <TPZMFSolutionTransfer.h>
//#include <Tools.h>
#include <ToolsMHM.h>
#include <Util/pzlog.h>
#include "TPZLinearAnalysis.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "Elasticity/TPZMixedElasticityND.h"
#include "TPZBndCondT.h"
#include "Pre/TPZMHMeshControl.h"
#include "TPZElasticityMHMHDivErrorEstimator.h"
#include "TPZHDivApproxCreator.h"
#include "TPZSimpleTimer.h"
#include "TPZTensor.h"
#include "TPZVTKGenerator.h"

enum EMatid  {ENone, EDomain, EBoundary, EMHM};
void RunElasticityProblem(int nCoarseDiv, int nInternalRef);
void IdentifyMHMDomain(TPZGeoMesh *gmesh, TPZVec<int> &domain);
void EstimateErrorElasticity(ProblemConfig &config, TPZMHMixedMeshControl *mhm);
TPZGeoMesh *CreateQuadGeoMesh(int nCoarseDiv, int nInternalRef);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);

void AdjustSkelSideOrient(TPZCompMesh *hdivmesh);

// check if the restraints are consistent
void CheckRestraintConsistencies(TPZCompMesh *hdivmesh, TPZVec<int> &domain);

void Substructure(TPZCompMesh *cmesh, TPZVec<int> &domains);

void GroupElements(TPZCompMesh *submesh);

void CondenseElements(TPZCompMesh *submesh);
void MakeConnectsInternal(TPZSubCompMesh *submesh);
const int global_nthread = 16;
void CheckBCs(TPZGeoMesh* gmesh);
void DivideGMesh(TPZGeoMesh *gmesh, int internaldiv, int skeletondiv);


void CheckNormalFluxes(TPZMultiphysicsCompMesh* cmesh, TPZAnalyticSolution* analy);

void Project(TPZCompEl* cel, TPZAnalyticSolution* analy, TPZFMatrix<REAL>& normalmat);


int main() {
    TPZLogger::InitializePZLOG();
    gRefDBase.InitializeAllUniformRefPatterns();

    // const std::set<int> nCoarseDiv = {3, 4, 5, 6};
    const std::set<int> nCoarseDiv = {4};
    const std::set<int> nInternalRef = {0};
    // const std::set<int> nInternalRef = {0, 1, 2, 3};
    for (const auto coarse_div : nCoarseDiv) {
        for (const auto internal_ref : nInternalRef) {
            RunElasticityProblem(coarse_div, internal_ref);
        }
    }

    return 0;
}

void RunElasticityProblem(const int nCoarseDiv, const int nInternalRef) {
    ProblemConfig config; 
    config.dimension = 2;
    config.exactElast = new TElasticity2DAnalytic;
    // config.exactElast.operator*().fProblemType = TElasticity2DAnalytic::EDispy;
    config.exactElast.operator*().fProblemType = TElasticity2DAnalytic::EThiago;
    config.problemname = "Elasticity";
    config.dir_name = "Journal";
    config.porder = 1;
    config.hdivmais = 1;
    config.materialids.insert(EDomain);
    config.bcmaterialids.insert(EBoundary);
    config.makepressurecontinuous = true;
    config.problemtype = ProblemConfig::TProbType::EElasticity;

    config.ndivisions = nCoarseDiv;
    config.ninternalref = nInternalRef;
    config.gmesh = CreateQuadGeoMesh(nCoarseDiv, nInternalRef);

    auto nEL = config.gmesh->NElements();
    for (auto i = 0; i < nEL; i++){
        TPZGeoEl *gel = config.gmesh->Element(i); 
        if (gel->Dimension() != config.gmesh->Dimension()) continue;
        int nsides = gel->NSides();
        int ncorner = gel->NCornerNodes();
        for (int iside = ncorner; iside < nsides-1; iside++){
            TPZGeoElSide gelside(gel,iside);
            if (!gelside.HasNeighbour(EMHM) && !gelside.HasNeighbour(EBoundary) ){
                TPZGeoElBC gelbc(gelside,EMHM);
            }
        }
        

    }
    // config.gmesh->BuildConnectivity();
    // config.gmesh->LoadReferences();

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "-GeoMesh.vtk";
        std::ofstream file(fileName);
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, file);
    }   
    int DIM = config.gmesh->Dimension();
    // if (nrefinternal || nrefskel) {
        // if(nrefskel > nrefinternal) DebugStop();
        DivideGMesh(config.gmesh, 1, 1);
    // }
    CheckBCs(config.gmesh);
    TPZVec<int> domain(config.gmesh->NElements(),-1);
    IdentifyMHMDomain(config.gmesh, domain);
    {
        std::ofstream out("domains.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, out, domain, true);
    }
    
    TPZHDivApproxCreator hdivCreator(config.gmesh);
    hdivCreator.HdivFamily() = HDivFamily::EHDivStandard;
    hdivCreator.ProbType() = ProblemType::EElastic;
    hdivCreator.IsRigidBodySpaces() = true;
    hdivCreator.SetDefaultOrder(config.porder);
    hdivCreator.SetExtraInternalOrder(0);
    hdivCreator.SetShouldCondense(false);
    //  hdivCreator.HybridType() = HybridizationType::EStandard;
    hdivCreator.HybridType() = HybridizationType::ENone;
    
    TPZAnalyticSolution *gAnalytic = 0;
    TPZMixedElasticityND* matelastic = 0;
    TElasticity2DAnalytic *elas = new TElasticity2DAnalytic;
    elas->gE = 250.;
    elas->gPoisson = 0.;
    elas->fProblemType = config.exactElast.operator*().fProblemType;
    gAnalytic = elas;
    matelastic = new TPZMixedElasticityND(EDomain, elas->gE, elas->gPoisson, 0, 0, 0 /*planestress*/, 2);

    //Insert Materials
    matelastic->SetExactSol(gAnalytic->ExactSolution(),4);
    matelastic->SetForcingFunction(gAnalytic->ForceFunc(), 3);
    hdivCreator.InsertMaterialObject(matelastic);

    TPZFMatrix<STATE> val1(2,2,0.);
    TPZManVector<STATE> val2(2,0.);
    val2[0] = 0.;
    TPZBndCondT<STATE> *BCond1 = matelastic->CreateBC(matelastic, EBoundary, 0, val1, val2);
    // BCond1->SetExactSol(gAnalytic->ExactSolution(),4);
    BCond1->SetForcingFunctionBC(gAnalytic->ExactSolution(), 3);
    hdivCreator.InsertMaterialObject(BCond1);

    // Interface materials for MHM
    val2[0] = 0.;
    TPZBndCondT<STATE> *MHM = matelastic->CreateBC(matelastic, EMHM, 0, val1, val2);
    hdivCreator.InsertMaterialObject(MHM);



    TPZMultiphysicsCompMesh *cmeshmulti = nullptr;
    {
        std::cout << "\n---------------- Creating Space -----------------" << std::endl;
        hdivCreator.CheckSetupConsistency();
        hdivCreator.SetMeshElementType();

        int lagLevelCounter = 1;
        TPZManVector<TPZCompMesh*,7> meshvec(5);
        hdivCreator.CreateAtomicMeshes(meshvec,lagLevelCounter);
        
        // {
        //     TPZCompMesh *hdivmesh = meshvec[0];
        //     if(nInternalRef == 0) AdjustSkelSideOrient(hdivmesh);
        //     CheckRestraintConsistencies(hdivmesh, domain);
        // }

        // possibly decrease polynomial order of skeleton
        // TPZCompMesh* cmeshflux = meshvec[0];
        // config.gmesh->ResetReference();
        // cmeshflux->LoadReferences();
        // if(pord != pordskel){
        //     if(pordskel > pord) DebugStop(); // cannot happen!
        //     for (int el = 0; el < cmeshflux->NElements(); el++) {
        //         TPZCompEl* cel = cmeshflux->Element(el);
        //         if(!cel) continue;
        //         TPZGeoEl* gel = cel->Reference();
        //         if (gel->MaterialId() != EMHM) continue;
                
        //         TPZInterpolatedElement* intel = dynamic_cast<TPZInterpolatedElement*>(cel);
        //         if(!intel) DebugStop();
                
        //         intel->PRefine(pordskel);
        //     }

//             {
//                 TPZCompMesh *hdivmesh = meshvec[0];
// //                AdjustSkelSideOrient(hdivmesh);
//                 CheckRestraintConsistencies(hdivmesh, domain);
//             }
        // }
        // cmeshflux->ExpandSolution();
        hdivCreator.CreateMultiPhysicsMesh(meshvec,lagLevelCounter,cmeshmulti);
    }
    TPZMultiphysicsCompMesh *cmesh = cmeshmulti;

    Substructure(cmesh, domain);






    std::cout << "NEquations " << cmesh->NEquations() << std::endl;
    //Create analysis environment
#ifdef PZDEBUG
    {
        std::string txt = "multiphysicsmesh.txt";
        std::ofstream myfile(txt);
        cmesh->Print(myfile);
    }
#endif
    TPZLinearAnalysis an(cmesh);
    
    //Solve problem
    SolveProblemDirect(an,cmesh);
    if(0)
    {
        auto SetDisp = [](TPZCompMesh *cmesh)
        {
            cmesh->Solution().Zero();
            TPZBlock &block = cmesh->Block();
            int64_t ncon = cmesh->NConnects();
            TPZFMatrix<STATE> &sol = cmesh->Solution();
            for(int64_t ic = 0; ic<ncon; ic++) {
                TPZConnect &c = cmesh->ConnectVec()[ic];
                if(c.LagrangeMultiplier() == 4) {
                    int64_t seq = c.SequenceNumber();
                    int64_t pos = block.Position(seq);
                    sol(pos,0) = 1.;
                }
            }
            cmesh->LoadSolution(sol);
            if(0)
            {
                std::ofstream out("mesh.txt");
                cmesh->Print(out);
            }
            std::cout << "printed\n";
        };
        SetDisp(cmesh);
        if(0) {
            int64_t nel = cmesh->NElements();
            for(int64_t el = 0; el<nel; el++) {
                TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cmesh->Element(el));
                if(sub) {
                    SetDisp(sub);
                }
            }
        }
    }
    cmesh->TransferMultiphysicsSolution();
    
    // Now, with the solution, let's check if the normal fluxes at
    // the skeleton are very close to what we get from the exact solution by doing sigma n
    // CheckNormalFluxes(cmesh,gAnalytic);
//    cmesh->MeshVector()[0]->Solution().Zero();
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(cmesh->MeshVector(), cmesh);
    cmesh->LoadSolution(cmesh->Solution());
    
#ifdef PZDEBUG
    {
        std::string txt = "multiphysicsmesh.txt";
        std::ofstream myfile(txt);
        cmesh->Print(myfile);
    }
#endif
    

    // -------> Calculating error
    std::cout << "------- Starting PostProc Error -------" << std::endl;
    an.SetExact(gAnalytic->ExactSolution());
    an.SetThreadsForError(global_nthread);
    std::ofstream ErroOut("myerrors.txt", std::ios::app);
    TPZMaterial *mat = cmesh->FindMaterial(EDomain);
    TPZMatErrorCombinedSpaces<STATE> *materr = dynamic_cast<TPZMatErrorCombinedSpaces<STATE>*>(mat);
    TPZManVector<REAL, 10> Errors(materr->NEvalErrors());
    bool store_errors = true;
    Errors.Fill(0.);
    cmesh->ElementSolution().Redim(cmesh->NElements(), Errors.size());
    {
        int64_t nel = cmesh->NElements();
        for(int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *>(cel);
            if(submesh) {
                submesh->ElementSolution().Redim(submesh->NElements(), Errors.size());
            }
        }
    }
    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
    an.PostProcessError(Errors, store_errors, ErroOut);
    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time PostProc Error = " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count()/1000. << " s" << std::endl;
    
    std::cout << "Computed errors." << std::endl;
    // error_sigma - error_energy - error_div_sigma - error_u - error_r - error_as - energy_norm_exact_sol
    std::cout.precision(15);
    std::cout.setf(std::ios::scientific);
    std::cout << "Errors = ";
    std::cout << Errors << std::endl;
    std::cout << "Errors = ";
    std::cout << Errors << std::endl;
    std::cout << "ErrorsFixedPrecision = ";
    std::cout.precision(15);
    std::cout.setf(std::ios::scientific);
    std::cout << Errors << std::endl;

    //  for (int i = 0; i < Errors.size(); i++) {
    //    std::cout << Errors[i] << std::endl;
    //  }
//    cmesh->ElementSolution().Print("element error");
    {
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(cmesh->MeshVector(), cmesh);
        TPZSimpleTimer postProc("Post processing2");
        const std::string plotfile = "solution";
        constexpr int vtkRes{0};
        
        
        TPZVec<std::string> fields = {
            // "ExactDisplacement",
            // "ExactStress",
            "Displacement",
            "SigmaX",
            "SigmaY",
            "TauXY",
            "ElementSigmaError"
        };
        auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
        vtk.SetNThreads(global_nthread);
        vtk.Do();
    }





}

TPZGeoMesh *CreateQuadGeoMesh(int nCoarseDiv, int nInternalRef) {

    TPZManVector<int, 4> bcIDs(4, EBoundary);
    TPZGeoMesh *gmesh = Tools::CreateGeoMesh(nCoarseDiv, bcIDs);
    gmesh->SetDimension(2);

    Tools::UniformRefinement(nInternalRef, gmesh);
    Tools::DivideLowerDimensionalElements(gmesh);

    return gmesh;
}




void EstimateErrorElasticity(ProblemConfig &config, TPZMHMixedMeshControl *mhm) {

    std::cout << "\nError Estimation processing for MHM-Hdiv problem " << std::endl;

    // Error estimation
    TPZMultiphysicsCompMesh *originalMesh = dynamic_cast<TPZMultiphysicsCompMesh *>(mhm->CMesh().operator->());
    if (!originalMesh) DebugStop();

    bool postProcWithHDiv = false;
    TPZElasticityMHMHDivErrorEstimator ErrorEstimator(*originalMesh, mhm, postProcWithHDiv);
    ErrorEstimator.SetAnalyticSolution(config.exactElast);
    
    ErrorEstimator.PrimalReconstruction();

    std::string command = "mkdir -p " + config.dir_name;
    system(command.c_str());

    TPZManVector<REAL, 6> errors;
    TPZManVector<REAL, 6> elementerrors;
    std::stringstream outVTK;
    outVTK << config.dir_name << "/" << config.problemname << "-" << config.ndivisions << "-" << config.ninternalref
           << "-Errors.vtk";
    std::string outVTKstring = outVTK.str();
    ErrorEstimator.ComputeErrors(errors, elementerrors, outVTKstring);

    {
        std::string fileName = config.dir_name + "/" + config.problemname + "-GlobalErrors.txt";
        std::ofstream file(fileName, std::ios::app);
        Tools::PrintElasticityErrors(file, config, errors);
    }
}
static std::set<int64_t> process;

void SetSubdomain(TPZGeoEl *gel, int domain, TPZVec<int> &alldomains)
{
    int64_t index = gel->Index();
//    std::cout << "Inserting " << index << " d " << domain << std::endl;
    alldomains[index] = domain;
    if(gel->HasSubElement()) {
        int nsubel = gel->NSubElements();
        if(nsubel){
//            std::cout << "nsubelements " << nsubel << std::endl;
        }
        else {
            std::cout << "I should stop\n";
        }
        for(int is=0; is<nsubel; is++) {
            SetSubdomain(gel->SubElement(is), domain, alldomains);
        }
    }
}


void ProcessElements(TPZGeoMesh *gmesh, TPZVec<int> &domain) {
    while(process.size()) {
        int64_t index = *process.begin();
        process.erase(index);
        TPZGeoEl *gel = gmesh->Element(index);
        int64_t gelindex = gel->Index();
//        std::cout << "Processing " << gelindex << " d " << domain[gelindex] << std::endl;
        int firstside = gel->FirstSide(1);
        int lastside = gel->NSides()-1;
        if(gel->Dimension() == 1) lastside++;
        int skeletonmatid = EMHM;
        int geldomain = domain[gel->Index()];
        if(geldomain == -1) DebugStop();
        if(gel->MaterialId() == EMHM) DebugStop();
        for(int side = firstside; side < lastside; side++) {
            TPZGeoElSide gelside(gel,side);
            if(gelside.HasNeighbour(skeletonmatid)) {
//                std::cout << "Element " << gelindex << " has skeleton neighbour along side " << side << "\n";
                continue;
            }
            TPZGeoElSide neighbour = gelside.Neighbour();
            int64_t indexneigh = neighbour.Element()->Index();
            if(domain[indexneigh] == -1) {
                process.insert(indexneigh);
                SetSubdomain(neighbour.Element(), geldomain, domain);
            }
        }
    }
}
void IdentifyMHMDomain(TPZGeoMesh *gmesh, TPZVec<int> &domain)
{
    int64_t nel = gmesh->NElements();
    domain.resize(nel);
    domain.Fill(-1);
    int currentdomain = 0;
    bool found = true;
    while(found) {
        found = false;
        for (int64_t el = 0; el < nel; el++) {
            if(domain[el] != -1) continue;
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->MaterialId() == EMHM) continue;
            if(gel->Father()) continue;
            SetSubdomain(gel, currentdomain, domain);
            found = true;
            process.insert(el);
            std::cout << "Inserting " << el << " d " << currentdomain << std::endl;
            ProcessElements(gmesh, domain);
            break;
        }
        currentdomain++;
    }
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->MaterialId() == EMHM) continue;
        if(domain[el] == -1) {
            DebugStop();
        }
    }

}

// check if the restraints are consistent
void CheckRestraintConsistencies(TPZCompMesh *hdivmesh, TPZVec<int> &domain) {
    int64_t nel = hdivmesh->NElements();
    hdivmesh->Reference()->ResetReference();
    hdivmesh->LoadReferences();
    TPZFMatrix<STATE> &sol = hdivmesh->Solution();
    sol.Zero();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *celskel = hdivmesh->Element(el);
        TPZGeoEl *gelskel = celskel->Reference();
        if(gelskel->MaterialId() == EMHM) {
            TPZFNMatrix<9,REAL> gradx(3,1),jac(1,1),jacinv(1,1),axes(1,3);
            REAL detjac;
            TPZManVector<REAL,2> ksi(1);
            gelskel->CenterPoint(gelskel->NSides()-1, ksi);
            gelskel->GradX(ksi, gradx);
            gelskel->Jacobian(gradx, jac, axes, detjac, jacinv);
            TPZManVector<REAL,3> skelnormal(3);
            // for (int i = 0; i<3; i++) {
            //     int j = (i+1)%3;
            //     int k = (i+2)%3;
            //     skelnormal[i] = axes(0,j)*1.-axes(0,k)*1.;
            // }
            skelnormal[0] = -axes(0,1);
            skelnormal[1] = axes(0,0);
            TPZInterpolatedElement *celskelint = dynamic_cast<TPZInterpolatedElement *>(celskel);
            int skelsideorient = celskelint->GetSideOrient(gelskel->NSides()-1);
//            std::cout << "skel " << el << " side orient " << skelsideorient << std::endl;
            TPZGeoElSide gelside(gelskel);
            REAL skelarea = gelside.Area();
            TPZStack<TPZCompElSide> elsidevec;
            int onlyinterpolated = 1;
            int removeduplicated = 0;
            gelside.HigherLevelCompElementList2(elsidevec, onlyinterpolated, removeduplicated);
            gelside.EqualLevelCompElementList(elsidevec, onlyinterpolated, removeduplicated);
            REAL smallarea[2] = {0.};
            // separate the smaller neighbours in two sets
            std::list<TPZCompElSide> twolists[2];

            int nelside = elsidevec.NElements();
            for (int is = 0; is<nelside; is++) {
                auto celside = elsidevec[is];
                auto smallgelside = celside.Reference();
                if(smallgelside.Dimension() != 1) continue;
                TPZInterpolatedElement *smallintel = dynamic_cast<TPZInterpolatedElement *>(celside.Element());
                int sideorient = smallintel->GetSideOrient(celside.Side());
                TPZManVector<REAL,3> smallnormal(3);
                TPZManVector<REAL,3> smallksi(1);
                smallgelside.CenterPoint(smallksi);
                smallgelside.Normal(smallksi, smallnormal);
                REAL inner = 0.;
                for(int i=0; i<3; i++) inner += skelnormal[i]*smallnormal[i];
                if(fabs(fabs(inner)-1.) > 1.e-8) DebugStop();
//                std::cout << "inner " << inner << " small sideorient " << sideorient << " skel orient " << skelsideorient << std::endl;
                if(inner < 0.) {
                    twolists[0].push_back(celside);
                    smallarea[0] += smallgelside.Area();
                } else {
                    twolists[1].push_back(celside);
                    smallarea[1] += smallgelside.Area();
                }
            }
            if(fabs(smallarea[0]-skelarea) > 1.e-9 ) DebugStop();
            if(fabs(smallarea[1]-skelarea) > 1.e-9 ) DebugStop();
            std::set<int> alldom;
            for (int is = 0; is<nelside; is++) {
                auto celside = elsidevec[is];
                auto smallgelside = celside.Reference();
                if(smallgelside.Dimension() != 1) continue;
                int64_t index = smallgelside.Element()->Index();
                int dom = domain[index];
                alldom.insert(dom);
            }
            if(alldom.size() != 2) DebugStop();
            // the normal component of one of the two lists matches the current normal
            TPZConnect &skelconnect = celskel->Connect(0);
            int64_t seqnum = skelconnect.SequenceNumber();
            int blsize = hdivmesh->Block().Size(seqnum);
            int64_t pos = hdivmesh->Block().Position(seqnum);
            TPZIntPoints *intrule = gelskel->CreateSideIntegrationRule(gelskel->NSides()-1, 3);
            int npoints = intrule->NPoints();
            REAL weight;
            TPZManVector<REAL,2> point(1,0.);
            TPZManVector<REAL,9> sol3d(6), point3d(2);
            TPZManVector<REAL,3> solskel(3);
            TPZManVector<REAL,3> sol3dnormal(3);
            for(int idof = 0; idof<blsize; idof++) {
                sol(pos+idof,0) = 1.;
                // update the dependent connect values
                for(auto &it : twolists[0]) {
                    int wrong = 0;
                    it.Element()->LoadSolution();
                    TPZGeoEl *gelsmall = it.Element()->Reference();
                    TPZTransform<> trskel = gelsmall->ComputeParamTrans(gelskel, gelskel->NSides()-1, it.Side());
                    TPZTransform<> trsmall = gelsmall->SideToSideTransform(it.Side(), gelsmall->NSides()-1);
                    for(int ip = 0; ip < npoints; ip++) {
                        intrule->Point(ip, point, weight);
                        TPZManVector<REAL,2> pointskel(1);
                        trskel.Apply(point, pointskel);
                        trsmall.Apply(point, point3d);
                        celskel->Solution(pointskel, 0, solskel);
                        it.Element()->Solution(point3d, 0, sol3d);
                        sol3dnormal.Fill(0.);
                        for(int i = 0; i<2; i++) for(int j=0; j<2; j++) {
                            sol3dnormal[i] += skelnormal[j]*sol3d[i*2+j];
                        }
                        for(int i=0; i<2; i++) {
                            REAL diff = fabs(sol3dnormal[i]-solskel[i]);
                            if(fabs(diff) > 1.e-7) {
                                wrong++;
                                std::cout <<  "oposite normal sol3dnormal " << sol3dnormal << " solskel " << solskel << std::endl;
                            }
                        }
                    }
                    if(wrong) {
                        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(it.Element());
//                        std::cout << "small index " << intel->Index() << " large index " << el << std::endl;
                        TPZConnect &c = intel->TPZInterpolationSpace::SideConnect(0, it.Side());
//                        std::cout << "small side orient " << intel->GetSideOrient(it.Side()) << " large side orient " << skelsideorient << std::endl;
//                        std::cout << "connect index " << intel->SideConnectIndex(0, it.Side()) << std::endl;
                        c.Print(*hdivmesh);
                        wrong = 0;
                    }
                }
                for(auto &it : twolists[1]) {
                    int wrong = 0;
                    it.Element()->LoadSolution();
                    TPZGeoEl *gelsmall = it.Element()->Reference();
                    TPZTransform<> trskel = gelsmall->ComputeParamTrans(gelskel, gelskel->NSides()-1, it.Side());
                    TPZTransform<> trsmall = gelsmall->SideToSideTransform(it.Side(), gelsmall->NSides()-1);
                    for(int ip = 0; ip < npoints; ip++) {
                        intrule->Point(ip, point, weight);
                        TPZManVector<REAL,2> pointskel(1);
                        trskel.Apply(point, pointskel);
                        trsmall.Apply(point, point3d);
                        celskel->Solution(pointskel, 0, solskel);
                        it.Element()->Solution(point3d, 0, sol3d);
                        sol3dnormal.Fill(0.);
                        for(int i = 0; i<2; i++) for(int j=0; j<2; j++) {
                            sol3dnormal[i] += skelnormal[j]*sol3d[i*2+j];
                        }
                        for(int i=0; i<2; i++) {
                            REAL diff = fabs(sol3dnormal[i]-solskel[i]);
                            if(fabs(diff) > 1.e-7) {
                                std::cout <<  "aligned normal sol3dnormal " << sol3dnormal << " solskel " << solskel << std::endl;
                                wrong++;
                            }
                        }
                    }
                    if(wrong) {
                        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(it.Element());
                        std::cout << "small index " << intel->Index() << " large index " << el << std::endl;
                        TPZConnect &c = intel->TPZInterpolationSpace::SideConnect(0, it.Side());
                        std::cout << "connect index " << intel->SideConnectIndex(0, it.Side()) << std::endl;
                        c.Print(*hdivmesh);
                        wrong = 0;
                    }
                }
            }
            delete intrule;
        }
    }
}

void AdjustSkelSideOrient(TPZCompMesh *hdivmesh) {
    int64_t nel = hdivmesh->NElements();
    hdivmesh->Reference()->ResetReference();
    hdivmesh->LoadReferences();
    TPZFMatrix<STATE> &sol = hdivmesh->Solution();
    sol.Zero();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *celskel = hdivmesh->Element(el);
        TPZGeoEl *gelskel = celskel->Reference();
        if(gelskel->MaterialId() == EMHM) {
            TPZFNMatrix<9,REAL> gradx(3,1),jac(1,1),jacinv(1,1),axes(1,3);
            REAL detjac;
            TPZManVector<REAL,2> ksi(1);
            gelskel->CenterPoint(gelskel->NSides()-1, ksi);
            gelskel->GradX(ksi, gradx);
            gelskel->Jacobian(gradx, jac, axes, detjac, jacinv);
            TPZManVector<REAL,3> skelnormal(3);
            for (int i = 0; i<3; i++) {
                int j = (i+1)%3;
                int k = (i+2)%3;
                skelnormal[i] = axes(0,j)*1.-axes(0,k)*1.;
            }
            TPZInterpolatedElement *celskelint = dynamic_cast<TPZInterpolatedElement *>(celskel);
            int skelsideorient = celskelint->GetSideOrient(gelskel->NSides()-1);
//            std::cout << "skel " << el << " side orient " << skelsideorient << std::endl;
            TPZGeoElSide gelside(gelskel);
            REAL skelarea = gelside.Area();
            TPZStack<TPZCompElSide> elsidevec;
            int onlyinterpolated = 1;
            int removeduplicated = 0;
            gelside.HigherLevelCompElementList2(elsidevec, onlyinterpolated, removeduplicated);
            gelside.EqualLevelCompElementList(elsidevec, onlyinterpolated, removeduplicated);
            REAL smallarea[2] = {0.};
            // separate the smaller neighbours in two sets
            std::list<TPZCompElSide> twolists[2];
            
            int nelside = elsidevec.NElements();
            for (int is = 0; is<nelside; is++) {
                auto celside = elsidevec[is];
                auto smallgelside = celside.Reference();
                if(smallgelside.Dimension() != 1) continue;
                TPZInterpolatedElement *smallintel = dynamic_cast<TPZInterpolatedElement *>(celside.Element());
                int sideorient = smallintel->GetSideOrient(celside.Side());
                TPZManVector<REAL,3> smallnormal(3);
                TPZManVector<REAL,3> smallksi(1);
                smallgelside.CenterPoint(smallksi);
                smallgelside.Normal(smallksi, smallnormal);
                REAL inner = 0.;
                for(int i=0; i<3; i++) inner += skelnormal[i]*smallnormal[i];
                if(fabs(fabs(inner)-1.) > 1.e-8) DebugStop();
                
//                std::cout << "inner " << inner << " sideorient " << sideorient << std::endl;
                if(inner < 0.) {
                    if(skelsideorient != -sideorient) {
//                        std::cout << "Adjusting skel " << el << " side orientation to " << -sideorient << std::endl;
                    }
                    celskelint->SetSideOrient(gelskel->NSides()-1, -sideorient);
                } else {
                    if(skelsideorient != sideorient) {
//                        std::cout << "Adjusting skel " << el << " side orientation to " << sideorient << std::endl;
                    }
                    celskelint->SetSideOrient(gelskel->NSides()-1, sideorient);
                }
                break;
            }
        }
    }
}


void Substructure(TPZCompMesh *cmesh, TPZVec<int> &domains)
{
//    std::cout << domains << std::endl;
    std::map<int,TPZSubCompMesh *> submeshes;
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) {
            DebugStop();
        }
        int64_t gelindex = gel->Index();
        int domain = domains[gelindex];
        if(domain == -1) continue;
        if(submeshes.find(domain) == submeshes.end()) {
            TPZSubCompMesh *sub = new TPZSubCompMesh(*cmesh);
            submeshes[domain] = sub;
        }
        TPZSubCompMesh *sub = submeshes[domain];
        sub->TransferElement(cmesh, el);
        
    }
    for (auto it : submeshes) {
        it.second->LoadReferences();
        it.second->ExpandSolution();
        GroupElements(it.second);
    }
    cmesh->ComputeNodElCon();
    for (auto it : submeshes) {
        MakeConnectsInternal(it.second);
    }
    for (auto it : submeshes) {
        CondenseElements(it.second);
        {
            std::ofstream out("submesh.txt");
            it.second->Print(out);
        }
    }
    for (auto it : submeshes) {
        TPZCompMesh *cmesh = it.second;
        int64_t neq = cmesh->NEquations();
        std::cout << "Configuring submesh neq = " << neq << std::endl;
#if defined(__x86_64__) || defined(__x86_64)
        it.second->SetAnalysisSparse(global_nthread);
#elif defined(__arm__) || defined(__aarch64__)
//        it.second->SetAnalysisFStruct(0);
        it.second->SetAnalysisSkyline();
#endif
    }
    cmesh->ComputeNodElCon();
    int64_t ncon = cmesh->NConnects();
    for(int64_t ic = 0; ic<ncon; ic++)
    {
        TPZConnect &c = cmesh->ConnectVec()[ic];
        if(c.NElConnected() == 0 && c.HasDependency()) c.RemoveDepend();
    }
    cmesh->CleanUpUnconnectedNodes();
}


void GroupElements(TPZCompMesh *cmesh) {
    // look for the boundary elements and group them with their neighbour
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) {
            DebugStop();
        }
        int64_t gelindex = gel->Index();
        if(gel->Dimension() == 1) {
            if(gel->MaterialId() == EMHM) continue;
            TPZGeoElSide gelside(gel);
            TPZGeoElSide neighbour = gelside.Neighbour();
            TPZGeoEl *neigh = neighbour.Element();
            TPZCompEl *cneigh = neigh->Reference();
            if(neigh->Dimension() != 2) DebugStop();
            int firstside = neigh->FirstSide(1);
            int lastside = neigh->NSides()-1;
            std::set<TPZCompEl *> grouped;
            grouped.insert(cneigh);
            for(int side=firstside; side<lastside; side++) {
                TPZGeoElSide gelside(neigh,side);
                TPZGeoElSide neighbour = gelside.Neighbour();
                if(neighbour.Element()->Dimension() == 1 && neighbour.Element()->MaterialId() != EMHM) {
                    TPZCompEl *cel = neighbour.Element()->Reference();
                    if(!cel) DebugStop();
                    grouped.insert(cel);
                }
            }
            if(grouped.size() == 1) DebugStop();
            TPZElementGroup *celgr = new TPZElementGroup(*cmesh);
            for(auto it : grouped) {
                celgr->AddElement(it);
            }
        }
    }
}

void CondenseElements(TPZCompMesh *submesh)
{
    submesh->ComputeNodElCon();
//    {
//        std::ofstream out("submesh.txt");
//        submesh->Print(out);
//    }
    int64_t nel = submesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = submesh->Element(el);
        if(!cel) continue;
        int ncon = cel->NConnects();
        bool found = false;
        for(int ic=0; ic<ncon; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if(c.LagrangeMultiplier() == 4) {
                c.IncrementElConnected();
                found = true;
                break;
            }
        }
        if(!found) continue;
//        for(int ic=0; ic<ncon; ic++) {
//            cel->Connect(ic).Print(*submesh);
//        }
        
        TPZCondensedCompEl *condense = new TPZCondensedCompElT<STATE>(cel);
    }
    submesh->CleanUpUnconnectedNodes();
}

void MakeConnectsInternal(TPZSubCompMesh *submesh) {
    int ncon = submesh->NConnects();
    int count = 0;
    bool found = false;
    for(int ic = 0; ic<ncon; ic++) {
        TPZConnect &c = submesh->Connect(ic);
        if(c.LagrangeMultiplier() == 4) {
            c.IncrementElConnected();
            count++;
        }
        if(count == 1) {
            submesh->MakeAllInternal();
            submesh->ComputeNodElCon();
            submesh->CleanUpUnconnectedNodes();
            found = true;
            break;
        }
    }
    if(!found) {
        std::cout << "No Lagrange multiplier of level 4\n";
        DebugStop();
    }
}



void CheckBCs(TPZGeoMesh* gmesh){
    const int64_t nel = gmesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl* gel = gmesh->Element(iel);
        if(gel->Dimension() != 2) continue;
        int firsts = gel->FirstSide(1);
        for( int s = firsts; s<gel->NSides()-1; s++)
        {
            TPZGeoElSide gelside(gel,s);
            if(gelside.Neighbour() == gelside)
            {
                DebugStop();
            }
        }
    }
}

void DivideGMesh(TPZGeoMesh *gmesh, int internaldiv, int skeletondiv){
    for (int div=0; div<internaldiv; div++) {
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel || gel->MaterialId() == EMHM) continue; // == EMHM skips all interfaces between macro domains
            if(gel->HasSubElement()) continue;
            TPZManVector<TPZGeoEl *> subels;
            gel->Divide(subels);
        }
    }
    for (int div=0; div<skeletondiv; div++) {
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel || gel->MaterialId() != EMHM) continue; // != EMHM skips everyone that is not interface between macro domains
            if(gel->HasSubElement()) continue;
            TPZManVector<TPZGeoEl *> subels;
            gel->Divide(subels);
        }
    }
}


void CheckNormalFluxes(TPZMultiphysicsCompMesh* cmesh, TPZAnalyticSolution* analy) {

    REAL maxnormdiff = std::numeric_limits<REAL>::min();
    for(int i = 0 ; i < cmesh->NElements() ; i++){
        TPZCompEl* cel = cmesh->Element(i);
        if(!cel) continue;
        if(!cel->Reference()) continue; // subcompmesh
        if(cel->Reference()->MaterialId() != EMHM) continue;
        if(cel->Dimension() != 1) DebugStop();
        
        // Computing the solution at the center of the element
        TPZCompEl* celhdiv = cmesh->MeshVector()[0]->Element(i);
        TPZGeoEl* gel = cel->Reference();
        TPZManVector<STATE,1> qsi = {0};
        TPZManVector<STATE,3> x = {0,0,0};
        gel->CenterPoint(gel->NSides()-1, qsi);
        gel->X(qsi, x);
        TPZVec<STATE> sol;
        celhdiv->Solution(qsi, 0, sol);
//        std::cout << "------ x = " << x << " ------" << std::endl;
//        std::cout << "sigmanormal_approx = " << sol << std::endl;
        
        // Computing the axes of the element to get the normal
        TPZFNMatrix<9,STATE> sigma(2,2,0.);
        analy->Sigma(x, sigma);
        TPZFMatrix<STATE> gradx, jac, axes, jacinv;
        REAL detjac = -1.;
        gel->GradX(qsi, gradx);
        gel->Jacobian(gradx, jac, axes, detjac, jacinv);
        TPZManVector<REAL,3> axes1(3,0.), axes2(3,0.), normal(3,0.);
        for(int idim = 0 ; idim < 2 ; idim++){
            axes1[idim] = axes(0,idim);
            axes2[idim] = axes(1,idim);
        }
        Cross(axes1, axes2, normal);
        TPZFMatrix<REAL> normalmat(3,1,0.);
        for(int idim = 0 ; idim < 3 ; idim++) normalmat(idim,0) = normal[idim];
        
        // Computing sigma n to get the normal flux in the skeleton
        TPZFMatrix<STATE> signormal_exact(3,3,0.);
        sigma.Multiply(normalmat, signormal_exact);
//        std::cout << "sigmanormal_exact = " << signormal_exact(0,0);
//        for(int idim = 1 ; idim < 3 ; idim++) std::cout << ", " <<  signormal_exact(idim,0);
//        std::cout << std::endl;
        
        // Projecting sigmanormal_exact into the space of this element
        Project(celhdiv,analy,normalmat);
        
        // Takig the norm of the difference between signormal_exact and signormal_approx (sol)
        REAL normdiff = 0.;
        TPZManVector<REAL,3> diff(3,0.);
        for(int idim = 0 ; idim < 3 ; idim++){
            normdiff += ( signormal_exact(idim,0) - sol[idim] ) * ( signormal_exact(idim,0) - sol[idim] );
        }
        normdiff = sqrt(normdiff);
        if(normdiff > maxnormdiff) maxnormdiff = normdiff;
//        std::cout << "normdiff = " << normdiff << std::endl;
        if(normdiff > 3.0) {
//            std::cout << "Big difference!!" << std::endl;
        }
    }
    cmesh->LoadSolutionFromMeshes();
//    std::cout << "maxnormdiff = " << maxnormdiff << std::endl;
}


void Project(TPZCompEl* celhdiv, TPZAnalyticSolution* analy, TPZFMatrix<REAL>& normalmat) {
    if(!celhdiv) DebugStop();
    
    // Some checks
    TPZGeoEl* gel = celhdiv->Reference();
    if(celhdiv->NConnects() != 1) DebugStop(); // Should be hdivbound
    if(gel->Type() != EOned) DebugStop(); // our mesh should only have tets
    TPZInterpolatedElement* intel = dynamic_cast<TPZInterpolatedElement*>(celhdiv);
    if(!intel) DebugStop();
    
    // Initializing material data
    TPZMaterialDataT<STATE> data;
    intel->InitMaterialData(data);
    TPZManVector<REAL> qsi(2,0.);
    int64_t nshape = data.phi.Rows();
    
    if(nshape != 3) DebugStop(); // Should we accept, at this point, more than 3 shape functions (linear skel)?
     
    // lhs and rhs
    TPZFNMatrix<81,REAL> matL2(nshape,nshape,0.);
    TPZFNMatrix<9,REAL> L2proj(nshape,3,0.);
    
    pztopology::TPZTriangle::IntruleType intrule(5); // Sigma is function of sin, cos and exp
    int npts = intrule.NPoints();
    const int neq = celhdiv->NEquations();
    for(int ip = 0 ; ip < npts ; ip++) {
        REAL weight;
        intrule.Point(ip, qsi, weight);
        intel->ComputeRequiredData(data, qsi);
        
        // Compute sigma.normal at this intpoint
        TPZFNMatrix<9,STATE> sigma(3,3,0.);
        analy->Sigma(data.x, sigma);
        TPZFMatrix<STATE> signormal_exact(3,1,0.);
        sigma.Multiply(normalmat, signormal_exact);
//        signormal_exact.Print(std::cout);
        // Compute l2 mat
        for(int i=0; i<nshape; i++) {
            L2proj(i,0) += data.phi(i,0)* signormal_exact(0,0) *data.detjac*weight;
            L2proj(i,1) += data.phi(i,0)* signormal_exact(1,0) *data.detjac*weight;
            L2proj(i,2) += data.phi(i,0)* signormal_exact(2,0) *data.detjac*weight;
            for (int j=0; j<nshape; j++) {
                matL2(i,j) += data.phi(i,0)*data.phi(j,0)*data.detjac*weight;
            }
        }
    }
    
    // Solve system
    matL2.SolveDirect(L2proj, ELU);
//    std::cout << L2proj << std::endl;
    
    TPZCompMesh* cmeshhdiv = celhdiv->Mesh();
    TPZConnect& c = celhdiv->Connect(0); // only one connect
    int64_t seq = c.SequenceNumber();
    const int64_t firsteq = cmeshhdiv->Block().Position(seq);
    const int64_t blocksize = cmeshhdiv->Block().Size(seq);
    TPZFMatrix<STATE>& solmat = cmeshhdiv->Solution();
    int count = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            solmat(firsteq+count,0) = L2proj(i,j);
            count++;
        }
    }
    
    
    // Test
//    qsi = {0,0};
//    TPZManVector<STATE,3> x = {0,0,0};
//    gel->CenterPoint(gel->NSides()-1, qsi);
//    gel->X(qsi, x);
//    TPZVec<STATE> sol;
//    celhdiv->Solution(qsi, 0, sol);
//    std::cout << "------ x = " << x << " ------" << std::endl;
//    std::cout << "sigmanormal_approx = " << sol << std::endl;
    
}




void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    //sets number of threads to be used by the solver
    constexpr int nThreads{global_nthread};

#if defined(__x86_64__) || defined(__x86_64)
    TPZSSpStructMatrix<STATE> matskl(cmesh);
//    TPZFStructMatrix<STATE> matskl(cmesh);
#elif defined(__arm__) || defined(__aarch64__)
    TPZSkylineStructMatrix<REAL> matskl(cmesh);
//    TPZFStructMatrix<> matskl(cmesh);
//    TPZBandStructMatrix<> matskl(cmesh);
#endif
    
    matskl.SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);
    
    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    
    an.SetSolver(step);
    
    std::cout << "------- Starting Assemble -------" << std::endl;
    std::cout << "Nequations = " << an.Mesh()->NEquations() << std::endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    an.Assemble();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time Assemble = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()/1000. << " s" << std::endl;
    
    ///solves the system
    std::cout << "------- Starting Solve -------" << std::endl;
    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
    an.Solve();
    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time Solve = " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count()/1000. << " s" << std::endl;
}
