#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
#include <vector>
#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h> /* draw from distribution */
/* "Global" variables used in many scripts during simulation but initialized in main.cpp */
extern std::string model;
extern std::string geneF;
extern std::string markerF;
extern std::string recF;
extern std::string chromoF;
extern int nDist;
extern int nSelGenes;

extern int seed;

extern gsl_rng * rd;                                // set in main.cpp

extern std::vector<double> distance_list;            // set in main.cpp (re-sorted in simu.cpp)
extern std::vector<int> pop_sampling_list;           // set in main.cpp
extern std::vector<int> nbr_pop_list;                // set in main.cpp


extern  std::string controlFile;
 
#endif


