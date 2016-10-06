#ifndef _CREATECONFIG
#define _CREATECONFIG

#include <vector>

//variables defined in the main:

extern double randomSeed;
extern int numbersOfMarkers;
//extern int chromoSize;
extern int nChromo;
extern std::vector <int> allChromoSize;

extern float pop1AlleleFreqAvg;
extern float pop1AlleleFreqStdDev;
extern float pop2AlleleFreqAvg;
extern float pop2AlleleFreqStdDev;
extern int nRecombMeiosis;                  // recombination file
extern int factSamplePoint;


// function prototypes
void makeAdmixemFiles();

void makeCustomFiles(std::string);


// RETRIEVE VARIABLES DEFINED ELSEWHERE FOR MATH FUNCTIONS
//extern Normal NormalGen; 
//extern Uniform UniformGen;

#endif
