//
// Created by victor on 16/03/2021.
//

#include "Output.h"
#include "InputTreatment.h"
#include "DataStructure.h"

void FlushTime(PreConfig &pConfig, clock_t start){
    float timer = float( clock () - start )/CLOCKS_PER_SEC;
    pConfig.timer << "Simulation time (" << pConfig.h << "x" << pConfig.h <<"): " << timer << "\n";
    pConfig.timer.flush();
}

void FlushTable(PreConfig &pConfig, char *argv[]){

    std::string plotname;
    plotname = pConfig.plotfile + "/" + pConfig.plotfile + ".csv";

    remove(plotname.c_str());
    ofstream table(plotname.c_str(), ios::app);

    ProblemConfig config;
    Configure(config, 0, pConfig, argv);

    pConfig.Erro.close();
    string file = "Erro.txt";
    CleanErrors(file);

    if(pConfig.mode == 1 || pConfig.mode == 2) InvertError(file);

    table << "Geometry" << "," << "Quadrilateral" << "\n";
    table << "Refinement" << "," << "Uniform" << "\n";
    table << "domain";

    if (pConfig.type != 2) table << "," << "[0 1]x[0 1]" << "\n";
    else table << "," << "[-1 1]x[-1 1]" << "\n";
    table << "Case" << "," << config.problemname << "\n";
    switch(pConfig.mode) {
        case 0:
            table << "Approximation" << "," << "H1" << "\n";
            table << "p order"  << "," << config.k << "\n";
            table << "---" << "," << "---" <<  "\n\n";
            table << "Norm" << "," << "H1" << "\n";
            break;
        case 1:
            table << "Approximation" << "," << "Hybrid" << "\n";
            table << "k order"  << "," << config.k << "\n";
            table << "Enrichment +n" << "," << config.n <<  "\n\n";
            table << "Norm" << "," << "Hybrid" << "\n";
            break;
        case 2:
            table << "Approximation" << "," << "Mixed" << "\n";
            table << "k order" << "," << config.k << "\n";
            table << "Enrichment +n" << "," << config.n<<  "\n\n";
            table << "Norm" << "," << "Mixed" << "\n";
            break;
    }
    FillErrors(table, file, pConfig.mode);
    table.close();
}

void InvertError(string file){
    std::vector<std::string> erro,rate;
    std::string sErro,sRate;
    int size = 0;

    std::ifstream iErro(file);
    std::ofstream temp("temp.txt");

    int it_count = -1, hash_count = 0;
    string Line;
    while(getline(iErro,Line)){
        it_count++;
        if (it_count == 0) {
            temp << Line << endl;
            continue;
        }

        hash_count++;
        if (Line.find("#") != string::npos) hash_count = 0;

        switch (hash_count) {
            case (1):
                sErro = Line;
                temp << Line << endl;
                continue;
            case (2):
                temp << sErro << endl;
                erro.push_back(Line);
                break;
            case (4):
                if (it_count == 5)
                    temp << Line << endl;
                else{
                    sRate = Line;
                    temp << Line << endl;
                }
                break;
            case (5):
                if (it_count == 6)
                    temp << Line << endl;
                else {
                    temp << sRate << endl;
                    rate.push_back(Line);
                }
                break;
            default:
                temp << Line << endl;
                break;
        }
    }
    temp.close();

    remove(file.c_str());
    std::ifstream itemp("temp.txt");
    std::ofstream Erro(file.c_str());

    it_count = -1; hash_count =0;
    int erro_counter =0, rate_counter =0;
    while(getline(itemp,Line)){
        it_count++;
        if (it_count == 0) {
            Erro << Line << endl;
            continue;
        }

        hash_count++;
        if (Line.find("#") != string::npos) hash_count = 0;

        switch (hash_count) {
            case (1):
                Erro << erro[erro_counter] << endl;
                erro_counter++;
                break;
            case (4):
                if (it_count == 5)
                    Erro << Line << endl;
                else{
                    Erro << rate[rate_counter] << endl;
                    rate_counter++;
                }
                break;
            default:
                Erro << Line << endl;
                continue;
        }
    }
    itemp.close();
    remove("temp.txt");
}

void FillErrors(ofstream &table,string f,int mode){
    std::ifstream iErro(f.c_str());

    int it_count = -1, hash_count = 0;
    string Line;
    while(getline(iErro,Line)) {
        it_count++;

        if (it_count == 0) {
            continue;
        }
        hash_count++;
        if (Line.find("#") != string::npos) hash_count = 0;

        FillLegend(table,hash_count,it_count);

        if(hash_count == 0) table <<"\n";
        else table << Line << "\n";
    }
}

void FillLegend(ofstream &table,int hash_count,int it_count){
    switch (hash_count) {
        case (0):
            table << "\n";
            break;
        case (1):
            table << "semiH1-error" << ",";
            break;
        case (2):
            table << "L2-error" << ",";
            break;
        case (3):
            table << "3rd-error" << ",";
            break;
        case (4):
            if (it_count == 5)
                table << "h" << ",";
            else
                table << "semiH1-rate" << ",";
            break;
        case (5):
            if (it_count == 6)
                table << "DOF" << ",";
            else
                table << "L2-rate" << ",";
            break;
        case (6):
            table << "3rd-rate" << ",";
            break;
        case (7):
            table << "h" << ",";
            break;
        case (8):
            table << "DOF" << ",";
            break;
        default:
            break;
    }
}

void CleanErrors(string file){
    std::ifstream iErro(file);
    std::ofstream temp("temp.txt");
    size_t last_index;
    string Line, residue, other = "other";


    int counter = -1;
    while (getline(iErro, Line)) {
        size_t found = Line.find(other);
        if (found != string::npos) continue;

        last_index = Line.find("=");

        if (last_index == string::npos) {
            last_index = Line.find(":");
            if(last_index == string::npos)
                temp << Line << endl;
            else{
                residue = Line.substr(last_index+2);
                temp << residue <<endl;
            }
        }
        else {
            residue = Line.substr(last_index+2);
            temp << residue <<endl;
        }
    }
    iErro.close(); temp.close();
    remove(file.c_str());
    rename("temp.txt", file.c_str());
}