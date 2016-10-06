#ifndef SIMU_H
#define SIMU_H
#include "settings.h"
//#include "model.h"
#include <string>
#include <map>

extern const std::string configFile;            // used in simu.cpp [simu-specific prefix added], initialized in main.cpp



class Model;

class Simu {

    public:
    Simu(Model&, std::string);
    //Simu(std::string);
    void createConfigFile() ;
    int launchAdmixem() const;
    void writeSummaryStat() ;
    void printSimu(std::ostream &) const;
    ~Simu();
    
    private:
    Model& simu_model;
    std::map<std::string,std::string> simu_parameters;
    int nHybPop;
    std::string simIterN;          // sth like 1_10 // => the 10th simu for the run on CPU 1
//    std::vector<double> list_distance;
    // then all the simu specific parameters
    



};

std::ostream &operator<<( std::ostream &, Simu const &);

#endif
