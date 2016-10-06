
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
std::string geneFile("genes_null_model.txt");
//rec_file
std::string recombFile("recomb_null_model.txt");
//marker_file
std::string markerFile("markers_null_model.txt");
//pheno_file
std::string phenotypeFile("pheno_null_model.txt");
//sexSel_file
std::string sexualSelFile("sexSel_null_model.txt");
//natSel_file
std::string naturalSelFile("natSel_null_model.txt");
//config_file
std::string configFile("null_model.cfg");

/* MARKER FILE */
//n_marker
int numbersOfMarkers(100);
//n_chromo
int nChromo(4);
//chromo_size
int chromoSize(1000);
//pop1_freqAvg
float pop1AlleleFreqAvg(0.8);
//pop1_freqSd
float pop1AlleleFreqStdDev(0.1);
//pop2_freqAvg
float pop2AlleleFreqAvg(0.8);
//pop2_freqSd
float pop2AlleleFreqStdDev(0.1);

/* RECOMBINATION FILE */
//n_rec
int nRecombMeiosis(16);      // was suggested 4n 
//factSamplePoint
int factSamplePoint(20);     // factSamplePoint * nChromo
    
/* CONFIGURATION FILE */
/* general settings */
// thread
int numThreads(20);
// seed
double randomSeed(0.2);             // OK ????
//YN_uniform
std::string uniformRec("yes");
//YN_ignMarkFreq
std::string igMarkerFreq("yes");
//OO_dumpNat
std::string dumpNatSel("Off");          // WHAT IS THIS ???
//out_fold
std::string outFolder("./admixOut");
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
/* pop size and sex ratio */
//pop1_initSize
int pop1InitSize(100);
//pop2_initSize
int pop2InitSize(100);
//pop1_limitSize
int pop1LimitSize(250);
//pop2_limitSize
int pop2LimitSize(250);
//male_ratio1
float maleRatio1(0.5);
//male_ratio2
float maleRatio2(0.5);
//hyb_limitSize
int hybLimitSize(250);
/* migration */
//YN_mig_gen1
std::string migrationOnlyGen1("yes");
//YN_migRate
std::string binomMigRate("yes");
//mig_p1H
float pop1ToHyb(0.1);
//mig_p2H
float pop2ToHyb(0.1);
//mig_Hp1
float hybToPop1(0);
//mig_Hp2
float hybToPop2(0);
//mig_p1p2
float pop1ToPop2(0);
//mig_gen1p1H
float gen1_pop1ToHyb(25);
//mig_gen1p2H
float gen1_pop2ToHyb(25);
/* female gamete and male choice */
//gam_avg
int avgFemGamete(10);                    // was 10 in the example       ???
//gam_sd
float stdFemGamete(3.16);               // was 3.16 in the example     ???
//kid_fct
std::string kidsFunc("Poisson");         // was Poisson in the example  ???
//male_samp
int maleSampling(50);                    // was 50 in the example       ???
/* generation, sampling and output */
//n_gen
int nGen(500);
//sample_freq
int sampFreq(50);
//OO_outMaker
string markerOut("On");
//gen_toSamp
string gensToSampleFile("");

#endif
