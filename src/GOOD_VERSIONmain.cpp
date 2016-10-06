#include "simu.h"       /* call simu fction */
#include "model.h"      /* call model fction */
#include <getopt.h>     /* parsing cmd line */
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <cstdlib>      /* atof */
//#include <gsl/gsl_randist.h> /* draw from distribution */

#ifdef _OPENMP
    # include <omp.h>
#endif


using namespace std;

char * getTime();
///// WARNING: COMPILE "-std=c++11", -fopenmp, and flags for lgsl !!!!!!!!!!!!!


//************************* => SHOULD DISTANCE HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// - if distance_list is not empty => number of demes = number of distance
// - if distance_list is empty and -d => -d demes uniformly distributed between 0 and 100
// - if distance_list is empty and no -d argument (or 0) => stop

// list of distances to coast (in the unit that you want, in the order that you want) 
vector<double> distance_list = {10, 30, 50, 20};


// !!!!!!!!!!!!!!!!!!!!!!! <= SET DISTANCE HERE ****************************************//



/* GLOBAL VARIABLES (DECLARED EXTERN IN THE OTHER FILES) */
// CAN BE EQUALLY CHANGED HERE OR VIA COMMAND LINE (OVERWRITTEN BY COMMAND LINE)
// EXCEPT "model" CAN ONLY BE SET VIA COMMAND LINE
string model;                               // -M => name of the model (see below for the available models)
string markerF("data/markers.txt");         // -m => path to marker file for Admixem
string recF("data/rec.txt");                // -r => path to recombination file for Admixem
string chromoF("data/chromo.txt");         // -l => path to chromosome length information for creating selection file for Admixem
string gffF;                                // -g => gff like file, needed if must general data files must be prepared
int nDemes(1);                              // -d => number of hybrid population
int nSimu(1000);                            //-n => number of simulation
int cpu(10);                               // -c => number of cores available for parallelization
int nSelGenes(50);                          // -s => number of genes contributing to the trait under selection


//TODO: rajouter minSel maxSel ???
const string natSelFile("natSel.txt"); // used in model.cpp [model-specific prefix added], declared extern in model.h 
const string sexSelFile("sexSel.txt");  // used in model.cpp [model-specific prefix added], declared extern in model.h
const string phenoFile("phenotype.txt");          // used in model.cpp [model-specific prefix added], declared extern in model.h    
const string geneFile("genes.txt");          // used in model.cpp [model-specific prefix added], declared extern in model.h    
const string configFile(".cfg");            // used in simu.cpp [simu-specific prefix added], declared extern in simu.h    

// => FOR DEBUG
const string controlFile("debug_control_test.txt");    
// <= FOR DEBUG

// ./runABC -C -M <model> -n <simulations> -c <cpu> -d <nDemes> -g <gff-like file> -m <marker_file> -r <rec_file> -l <chromo_file> -s <nSelGenes>]

// DEFINE GLOBALY THE NAMES OF THE MODELS THAT ARE AUTHORIZED < ============================= IF NEEDED UPDATE HERE POSSIBLE MODELSS
set<string> available_models = {"model0", "model1"};

int seed;
int main (int argc, char * const argv[]) 
{  

    srand(time(NULL));          
    seed = rand() % 1000000; // seed

    string command;
    int foo;
    /**/ // FOR DEBUG => 
    ofstream outControl;
    outControl.open(controlFile.c_str());
    /**/ // <= FOR DEBUG
    
    bool help(false);                   // if -h => print the help
    bool makeData(false);               // if -C => run createGeneralSettings
    int option_char;

    // RETRIEVE FROM COMMAND LINE    
    while ((option_char = getopt (argc , argv, "ChM:n:c:d:g:m:r:l:s:")) != -1)
    {
        switch (option_char)
        {  
            case 'C': makeData = true;break;
            case 'h': help = true;break;
            case 'M': model = optarg;break;
            case 'n': nSimu = atoi(optarg);break;
            case 'c': cpu = atoi(optarg);break;
            case 'd': nDemes = atoi(optarg);break;
            case 'g': gffF = optarg; break;
            case 'm': markerF = optarg; break;
            case 'r': recF = optarg; break;
            case 'l': chromoF = optarg;break;
            case 's': nSelGenes = atoi(optarg);break;
        }
    }    
    if(argc == 1 || help) //help requested or no argument
    {
        cout << "...... HELP  - runABC - Usage: " << endl;
        cout << "./runABC -C -M <model> -n <simulations> -c <cpu> -d <nDemes> -g <gff-like file> -m <marker_file> -r <rec_file> -l <chromo_file> -s <nSelGenes>]" << endl;
        cout << "... -C => if present, run createGeneralSettings to prepare markers.txt, rec.txt and chromo.txt before simulations" << endl;
        cout << "... -M => valid model name (e.g. model0) [no default value: stop if not present and -C absent]" << endl;
        cout << "... -n => number of simulations [default: 1000]" << endl;
        cout << "... -c => number of cores available [default: 10]" << endl;
        cout << "... -d => the number of demes (number of hybrid populations to initialize) [default: 1]" << endl;
        cout << "... -g => the gff-like file needed for creating custom data with createGeneralSettings (if not pass: randomly simulated markers)" << endl;
        cout << "... -m => the marker file for Admixem [default: data/markers.txt]" << endl;
        cout << "... -r => the recombination file for Admixem [default: data/rec.txt]" << endl;
        cout << "... -l => the file with information about chromosome length (as outputed by createGeneralSettings)" << endl;
        cout << "... -s => the number of genes under selection" << endl;
        cout << "... -h => print the help and stop programm execution (if present, ignore all other arguments)" << endl;
        exit(EXIT_SUCCESS);
    }
    if(model.empty() and !makeData)
    {
        if(makeData) cout << "... WARNING: no model selected, will stop after creating data files." << endl;
        else 
        {
            cout << "Nothing to do and no model selected - bye bye :-) " << endl;
            exit(EXIT_SUCCESS);
        }
    }
    
    if(nDemes == 0 && distance_list.empty())
    {
        cout << "... STOP: Distance list is empty and number of demes not > 0" << endl;
        exit(EXIT_SUCCESS);
    }
    
    char* start = getTime();
    cout << "Start running at:\t" << start ;    



    /* SET THE RANDOM ONE AND FOR ALL */
/*
    gsl_rng * rng = gsl_rng_alloc (gsl_rng_default);        //use default one for the moment   
    gsl_rng_env_setup();
    long seed2 = 21;
    gsl_rng_set (rng, seed2);                  // set seed
*/

    /* STEP 1: if needed (-C), create data that are shared across all models and all simulations */
    // first create the general data if needed
    if(makeData)
    {
        int x;
        command = "mkdir data";
        x = system(command.c_str());
        
        if(gffF.empty())            // no gff like file => generate the files "randomly" as Admixem did
        {
             x = system("./createGeneralSettings -A");
            
        } else
        {
            command = "./createGeneralSettings -g " + gffF;
            x = system(command.c_str());
        }
        if(x != 0)
        {
            cout << "...... STOP! Problem during execution of \"createGeneralSettings\"" << endl;
            cout << "...... Error during general data preparation ! " << endl;
            cout << "End running at:\t" << getTime();
            exit(EXIT_FAILURE);
        }
    }

    cout << "...... Create general data for simulation: DONE." << endl;
    
    // some checks before running simulations
    if(makeData && model.empty()) // no model selected but "-C" -> just wanted to create data files
    {
        cout << "...... No model selected, stop here ! Bye bye :-) " << endl;
        cout << "End running at:\t" << getTime();
        exit(EXIT_SUCCESS);
    }
    if(available_models.find(model) == available_models.end())
    {
        cout << "...... \"" << model << "\" is not valid model name. Set of possible models may need to be updated in \"main.cpp\"." << endl;
        cout << "End running at:\t" << getTime();
        exit(EXIT_SUCCESS);
    }

    // create folders for the simulation
    command = "mkdir simu";
    foo = system(command.c_str());
    command = "mkdir simu/"+model;
    foo = system(command.c_str());
    command = "mkdir data/"+model;
    foo = system(command.c_str());
    
    /* STEP 2: create the model object and generate the data that are model specific */
    
    if(!distance_list.empty()) nDemes = distance_list.size();
    
    Model current_model(model, nDemes, nSelGenes);  //initialize model object
    /**/ // FOR DEBUG => 
    outControl << current_model;
    /**/ // <= FOR DEBUG 
    current_model.writeModelFiles();              // write out model specific files 


    //NOW RUN SIMULATIONS IN PARALLEL ("nSimu" to run in "cpu" cores)

    #ifdef _OPENMP
        omp_set_num_threads(cpu);
    #endif
    
    /* STEP 3: create the simulation object, generate simulation specific files and launch admixem */
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for(int simu_nbr = 0; simu_nbr < nSimu; simu_nbr++)
    {
        string simuID = to_string(simu_nbr);    
        #ifdef _OPENMP
            simuID = to_string(omp_get_thread_num()) +"_"+to_string(simu_nbr);
        #endif    
        command = "mkdir simu/"+model+"/"+simuID;
        foo = system(command.c_str());        
        
        Simu current_simu(current_model, nDemes, simuID);
        #ifdef _OPENMP
            #pragma omp critical
        #endif
        {
            /**/ // FOR DEBUG => 
            outControl << current_simu;
            /**/ // <= FOR DEBUG 
        }
        current_simu.createConfigFile(); 
        current_simu.launchAdmixem();  
        current_simu.writeSummaryStat();
    }
    /**/ // FOR DEBUG => 
    outControl.close(); 
    /**/ // <= FOR DEBUG 
    
 
    cout << "Program started at:\t" << start ;    
    cout << "Finished at:\t" << getTime();            
    cout << "************************************************************************************" << endl;
    cout << "*************************************** DONE ***************************************" << endl;
    cout << "************************************************************************************" << endl;
    
    outControl.close();
    
    return 0;
}
    
//    #include <cctype> 
    
char * getTime()
{
    time_t now = time(0);
    char* start = ctime(&now);
    return start;
}    
