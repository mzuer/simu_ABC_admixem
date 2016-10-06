#ifndef SUMMARYSTAT_H
#define SUMMARYSTAT_H

#include <string>

/* "extern" variables, used in model.cpp, initialized in main.cpp */
extern const std::string geneFile;


// gene output of admixem, marker output of admixem for the last simulation, simulationID needed for output file name
void writeSummaryTables(std::string admixGeneFile, std::string admixMarkerFile, std::string simIterN);

void writeChromosomeBlocks(std::string simIterN);



#endif
