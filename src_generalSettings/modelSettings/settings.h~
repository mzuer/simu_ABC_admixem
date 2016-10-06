
#ifndef _SETTINGS_
#define _SETTNGS_

/* THIS IS THE TEMPLATE FILE USED BY OTHER SCRIPTS FOR CREATING 
 * THE "settings.h" THAT IS USED TO CREATE THE FILES
 * NEEDED FOR RUNNING THE SIMULATION 
 * DO -> NOT <- WHAT FOLLOWED THE "//" AND IF YOU ADD
 * COMMENT, DO -> NOT <- USE THE "//" NOTATION 
 */

/* OUTPUT FILE NAMES */
//gene_file
std::string geneFile("genes.txt");      // will be changed into <model_name>_genes.txt
//rec_file
std::string recombFile("recomb.txt");
//marker_file
std::string markerFile("markers.txt");
//pheno_file
std::string phenotypeFile("pheno.txt");
//sexSel_file
std::string sexualSelFile("sexSel.txt");
//natSel_file
std::string naturalSelFile("natSel.txt");

// seed
double randomSeed(0.2);             // OK ????

/* pop names and labels */
//pop1_id
std::string pop1Name("rapa");
//pop1_lab
std::string pop1Label("r");
//pop2_id
std::string pop2Name("nigra");
//pop2_lab
std::string pop2Label("n");
//hyb_id
std::string hybName("hyb");

/* MARKER FILE */
//n_marker
int numbersOfMarkers(10000);                           /// CHANGE HERE : number of markers to generate
//n_chromo
int nChromo(4);                                      /// CHANGE HERE : number of haploid chromosomes
//chromo_size
int chromoSize(1000);                                /// CHANGE HERE : size of the chromosome
//pop1_freqAvg
float pop1AlleleFreqAvg(0.8);                        /// CHANGE HERE : average frequency most common allele pop1
//pop1_freqSd
float pop1AlleleFreqStdDev(0.1);                     /// CHANGE HERE : sd frequency most common allele pop1
//pop2_freqAvg
float pop2AlleleFreqAvg(0.8);                        /// CHANGE HERE : average frequency most common allele pop2
//pop2_freqSd
float pop2AlleleFreqStdDev(0.1);                     /// CHANGE HERE : sd frequency most common allele pop2

/* RECOMBINATION FILE */
//n_rec
int nRecombMeiosis(16);      // was suggested 4n            // ???
//factSamplePoint
// usage in code : factSamplePoint * nChromo 
int factSamplePoint(20);                                    // ???
    
#endif
