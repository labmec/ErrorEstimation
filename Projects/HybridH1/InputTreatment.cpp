//
// Created by victor on 16/03/2021.
//

#include "InputTreatment.h"
#include "DataStructure.h"
#include "MeshInit.h"
#include "Tools.h"

void Configure(ProblemConfig &config,int ndiv,PreConfig &pConfig,char *argv[]){
    ReadEntry(config, pConfig);
    config.ndivisions = ndiv;
    config.dimension = pConfig.dim;
    config.prefine = false;
    config.exact.operator*().fSignConvention = 1;
    config.exact->fDimension = config.dimension;

    bool isOriginCentered = 0; /// Wheater the domain = [0,1]x^n or [-1,1]^n
    if(pConfig.type == 2 || pConfig.type >= 7) isOriginCentered = 1;

    TPZGeoMesh *gmesh;
    TPZManVector<int, 4> bcids(4, -1);
    if(pConfig.type == 4){
        bcids[0] = -3;
        bcids[1] = bcids[3] = -2;
    }
    gmesh = Tools::CreateGeoMesh(1, bcids, config.dimension,isOriginCentered,pConfig.topologyMode);

    if(config.gmesh) delete config.gmesh;
    config.gmesh = gmesh;

    Tools::UniformRefinement(config.ndivisions, gmesh);
    
    config.ApplyDivision();
    gmesh = config.gmesh;

    // apply random refinement
    if(0)
    {
        int64_t nel = gmesh->NElements();
        int dim = gmesh->Dimension();
        int numeldiv = 5;
        int neldiv = 0;
        for (int64_t el = 0; el<nel; el++) {
            auto gel = gmesh->Element(el);
            if(gel->Dimension() != dim) continue;
            if(gel->HasSubElement()) continue;
            TPZVec<TPZGeoEl *> subels;
            gel->Divide(subels);
            neldiv++;
            if(neldiv >= numeldiv) break;
        }
        // divide the boundary elements
        nel = gmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            auto gel = gmesh->Element(el);
            if(gel->Dimension() == 0 || gel->HasSubElement()) continue;
            TPZGeoElSide gelside(gel);
            bool should_divide = false;
            TPZGeoElSide neighbour = gelside.Neighbour();
            while(neighbour != gelside) {
                if(neighbour.Element()->HasSubElement()) should_divide = true;
                neighbour = neighbour.Neighbour();
            }
            if(should_divide) {
                TPZVec<TPZGeoEl *> subels;
                gel->Divide(subels);
            }
        }
    }

    {
        std::stringstream sout;
        sout << "gmesh_divide." << config.fElIndexDivide.size() << ".vtk";
        std::ofstream out(sout.str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.bcmaterialids.insert(-2);

    if(pConfig.type == 4){
        config.bcmaterialids.insert(-3);
    }

    if (pConfig.type  == 2) {
        config.materialids.insert(2);
        config.materialids.insert(3);
        config.bcmaterialids.insert(-5);
        config.bcmaterialids.insert(-6);
        config.bcmaterialids.insert(-8);
        config.bcmaterialids.insert(-9);

        SetMultiPermeMaterials(config.gmesh);
    }

    if(pConfig.argc != 1) {
        config.k = atoi(argv[3]);
        config.n = atoi(argv[4]);
    }

    if(pConfig.debugger == true && ndiv != 0){
        Tools::DrawGeoMesh(config,pConfig);
    }
}

void DataInitialization(int argc, char *argv[],PreConfig &hybConfig,PreConfig &mixConfig){
    EvaluateEntry(argc,argv,hybConfig);
    InitializeOutstream(hybConfig,argv);

    EvaluateEntry(argc,argv,mixConfig);
    InitializeOutstream(mixConfig,argv);
}

void CopyHybSetup(PreConfig &hybConfig, PreConfig &mixConfig){
    mixConfig.k = hybConfig.k;
    mixConfig.n = hybConfig.n;
    mixConfig.problem =   hybConfig.problem;
    mixConfig.approx =    "Mixed";
    mixConfig.topology =  hybConfig.topology;
    mixConfig.maxIter = hybConfig.maxIter;

    mixConfig.refLevel = hybConfig.refLevel;
    mixConfig.debugger = hybConfig.debugger;
}

void FluxErrorConfigure(ProblemConfig &config,PreConfig &pConfig){
    //config.exact.operator*().fExact = TPZAnalyticSolution::TForce;
    config.dimension = pConfig.dim;
    config.prefine = false;
    config.exact.operator*().fSignConvention = 1;
    config.exact->fDimension = config.dimension;

    TPZGeoMesh *gmesh;
    TPZManVector<int, 4> bcids(4, -1);
    bcids[3] = -2;
    gmesh = Tools::CreateGeoMesh(1, bcids, config.dimension,0,pConfig.topologyMode);

    config.gmesh = gmesh;
    config.materialids.insert(1);
    config.bcmaterialids.insert(-1);
    config.bcmaterialids.insert(-2);

    if(pConfig.debugger == true){
        Tools::DrawGeoMesh(config,pConfig);
    }
}

void EvaluateEntry(int argc, char *argv[],PreConfig &pConfig){
    if(argc != 1 && argc != 5){
        std::cout << "Invalid entry";
        DebugStop();
    }
    if(argc == 5){
        pConfig.argc = argc;
        for(int i = 3; i < 5 ; i++)
            IsInteger(argv[i]);
        if(std::strcmp(argv[2], "H1") == 0)
            pConfig.mode = 0;
        else if(std::strcmp(argv[2], "Hybrid") == 0) {
            pConfig.mode = 1;
            if(pConfig.n < 1 ){
                std::cout << "Unstable method\n";
                DebugStop();
            }
        }
        else if(std::strcmp(argv[2], "Mixed") == 0) {
            pConfig.mode = 2;
            if(pConfig.n < 0) DebugStop();
        }
        else DebugStop();

        if(std::strcmp(argv[1], "ESinSin") == 0) {
            pConfig.type = 0;
            pConfig.problem = "ESinSin";
        }
        else if(std::strcmp(argv[1], "EArcTan") == 0) {
            pConfig.type = 1;
            if(pConfig.n < 1 ){
                std::cout << "Unstable method\n";
                DebugStop();
            }
            pConfig.problem = "EArcTan";
        }
        else if(std::strcmp(argv[1], "ESteklovNonConst") == 0) {
            pConfig.type = 2;
            if(pConfig.n < 0) DebugStop();
            pConfig.problem = "ESteklovNonConst";
        }
        else DebugStop();
    }
    else{
        if (pConfig.approx == "H1") pConfig.mode = 0;
        else if (pConfig.approx == "Hybrid")  pConfig.mode = 1;
        else if (pConfig.approx == "Mixed") pConfig.mode = 2;
        else DebugStop();

        if (pConfig.problem== "ESinSin") pConfig.type= 0;
        else if (pConfig.problem=="EArcTan")  pConfig.type = 1;
        else if (pConfig.problem == "ESteklovNonConst") pConfig.type = 2;
        else if (pConfig.problem == "EBubble2D") pConfig.type = 3;
        else if (pConfig.problem == "ELaplace") pConfig.type = 4;
        else if (pConfig.problem == "E2SinSin") pConfig.type = 5;
        else if (pConfig.problem == "E10SinSin") pConfig.type = 6;
        else if (pConfig.problem == "ESing2D") pConfig.type = 7;
        else if (pConfig.problem == "ESinMark") pConfig.type = 8;
        else if (pConfig.problem == "EProb") pConfig.type = 9;
        else DebugStop();
    }

    if (pConfig.topologyMode != -1) DebugStop();
    if (pConfig.topology == "Triangular") pConfig.topologyMode = 1;
    else if (pConfig.topology == "Quadrilateral") pConfig.topologyMode = 2;
    else if (pConfig.topology == "Tetrahedral") pConfig.topologyMode = 3;
    else if (pConfig.topology == "Hexahedral") pConfig.topologyMode = 4;
    else if (pConfig.topology == "Prism") pConfig.topologyMode = 5;
    else if (pConfig.topology == "LQuad") pConfig.topologyMode = 6;
    if (pConfig.topologyMode == -1) DebugStop();

    if(pConfig.topologyMode < 3 || pConfig.topologyMode == 6) pConfig.dim = 2;
    else pConfig.dim = 3;
}

void InitializeOutstream(PreConfig &pConfig, char *argv[]){
    //Error buffer
    if( remove( "Erro.txt" ) != 0) perror( "Error deleting file" );
    else puts( "Error log successfully deleted" );

    pConfig.Erro.open("Erro.txt",std::ofstream::app);
    pConfig.Erro << "----------COMPUTED ERRORS----------\n";

    pConfig.Log = new TPZVec<REAL>(pConfig.numErrors, -1);
    pConfig.rate = new TPZVec<REAL>(pConfig.numErrors, -1);

    ProblemConfig config;
    Configure(config,0,pConfig,argv);

    std::stringstream out;
    switch (pConfig.topologyMode) {
        case 1:
            pConfig.topologyFileName = "2D-Tri";
            break;
        case 2:
            pConfig.topologyFileName = "2D-Qua";
            break;
        case 3:
            pConfig.topologyFileName = "3D-Tetra";
            break;
        case 4:
            pConfig.topologyFileName = "3D-Hex";
            break;
        case 5:
            pConfig.topologyFileName = "3D-Prism";
            break;
        case 6:
            pConfig.topologyFileName = "2D-LQuad";
            break;
        default:
            DebugStop();
            break;
    }

    switch(pConfig.mode) {
        case 0: //H1
            out << "H1_" <<  pConfig.topologyFileName << "_" << config.problemname << "_k-"
                << config.k;
            pConfig.plotfile = out.str();
            break;
        case 1: //Hybrid
            out << "Hybrid_" <<   pConfig.topologyFileName << "_" << config.problemname  << "_k-"
                << config.k << "_n-" << config.n;
            pConfig.plotfile = out.str();
            break;
        case 2: // Mixed
            out << "Mixed_" <<  pConfig.topologyFileName << "_" << config.problemname << "_k-"
                << config.k << "_n-" << config.n;
            pConfig.plotfile = out.str();
            break;
        default:
            std::cout << "Invalid mode number";
            DebugStop();
            break;
    }
    std::string command = "mkdir -p " + pConfig.plotfile;
    system(command.c_str());

    std::string timer_name = pConfig.plotfile + "/timer.txt";

    if( remove(timer_name.c_str()) != 0)
        perror( "Error deleting file" );
    else
        puts( "Error log successfully deleted" );

    pConfig.timer.open(timer_name, std::ofstream::app);
}

void IsInteger(char *argv){
    std::istringstream ss(argv);
    int x;
    if (!(ss >> x)) {
        std::cerr << "Invalid number: " << argv << '\n';
    } else if (!ss.eof()) {
        std::cerr << "Trailing characters after number: " << argv << '\n';
    }
}

void ReadEntry(ProblemConfig &config, PreConfig &preConfig){


    config.exact = new TLaplaceExample1;
    switch(preConfig.type){
        case 0:
            config.exact.operator*().fExact = TLaplaceExample1::ESinSin;
            break;
        case 1:
            config.exact.operator*().fExact = TLaplaceExample1::EArcTan;
            break;
        case 2:
            config.exact.operator*().fExact = TLaplaceExample1::ESteklovNonConst;
            preConfig.h*=2;
            break;
        case 3:
            config.exact.operator*().fExact = TLaplaceExample1::EBubble2D;
            break;
        case 4:
            config.exact.operator*().fExact = TLaplaceExample1::ELaplace2D;
            break;
        case 5:
            config.exact.operator*().fExact = TLaplaceExample1::E2SinSin;
            break;
        case 6:
            config.exact.operator*().fExact = TLaplaceExample1::E10SinSin;
            break;
        case 7:
            DebugStop();
            //config.exact.operator*().fExact = TLaplaceExample1::ESing2D;
            break;
        case 8:
            config.exact.operator*().fExact = TLaplaceExample1::ESinMark;
            break;
        default:
            DebugStop();
            break;
    }

    config.k = preConfig.k;
    config.n = preConfig.n;
    config.problemname = preConfig.problem;
    config.exact->fmaxIter = preConfig.maxIter;
}
