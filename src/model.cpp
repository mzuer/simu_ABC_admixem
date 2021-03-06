#include "settings.h"       /* declare extern variables: model, geneF, markerF, recF, chromoF, nDemes, run, nSelGenes*/
#include "model.h"          /* class definition, definition of output file names*/
#include "my_math.h"        /* get from distribution*/

#include <cmath>            /* round()*/
#include <string>
#include <vector>           /* temporarily store length of chromosomes, */
#include <map>              /* map of the genes under selection (chr_nbr : vector of positions) */
#include <algorithm>        /* sort genes position for the gene file */
#include <gsl/gsl_randist.h> /* draw from distribution */
#include <iostream>
#include <fstream>
#include <sstream>      /* sstream */
#include <cstdlib>     


/* FOLLOWING FILES ARE CREATED VIA THIS SCRIPT:

_genes.txt      => file with list of genes under selection
_natSel.txt     => file with natural selection info like "hyb1	-1	exp(-pow(Pheno0+Pheno1+Pheno2+Pheno3+Pheno4-5, 2)/pow(0.125467, 2))"
_sexSel.txt     => file with sexual selection
_phenotype.txt  => file with phenotype info like "Pheno0	if(chr5_220==2,1,if(chr5_220==1,0.5,if(chr5_220==0,0,0)))"

*/


/* AVAILABLE FROM setttings.h AS INITIALIZED IN main.cpp
string model;
string geneF;
string markerF;
string recF;
string chromoF;
int nDemes;
int run;
int nSelGenes;
 */
/* AVAILABLE FROM model.h 
string natSelFile;
string sexSelFile;
string phenoFile;
 */


/*********************** NAME OF THE POPULATIONS
parental population names: pop1 and pop2
parental population  labels: x and y respectively
hybrid popuulation: depends on the number of deme
->  names:  hybFoo, hyb1, hyb2, ...                       
->  labels: a,b,c...

hybFoo serves as "parental" pop to other hybrid pop instead of rapa because we want to avoid
that rapa (which has optimal phenotype) totally invades the hybrid pop

*********************/
/**** DRAW FROM PRIORS 

** MODEL 0: constant gene flow, selection (omega) increases in direction of the coast
** MODEL 1: selection (omega) and gene flow increase in direction of the coast

hyb1 => is the closest to the sea (the smallest distance)
-> highest gene flow in model 1
-> highest selection pressure in models 0 and 1

###### for the gene flow (MODEL 1 ONLY) : => set in simu.cpp (config file)

hybFoo --- m*d1 ---> hyb1
hybFoo --- m*d2 ---> hyb2
hybFoo --- m*d3 ---> hyb3

=> draw m in uniform distribution 

how to find d1, d2, d3 ?

equation from Austerlitz -> distribution defined by a single parameter 'beta'

then function of distance : d1 = fct_beta(hyb1_distance) etc.

=> draw beta in uniform distribution [0-10'000]


###### for the selection (MODEL 0 and MODEL 1):  => set in natSel.txt

omega * d1  -> omega_hyb1     => omega should be very small - sel +++     
omega * d2  -> omega_hyb2     => increasing omega  - sel ++
omega * d3  -> omega_hyb3     => omega should be high - sel +

where d comes from a line that goes through 0;1 et dist_max;0  // y=a*x where x is distance d is y

high omega (e.g. 100) => neutral selection
small omega (e.g. 0.04) => strong selection

[Rhoné: 1/omega is the intensity of selection, 0.01 neutral, 25 very strong]

****/
/*
const double shape_omega

*/
using namespace std;


string allParametersLine("");

/// ??????????????????????????????????
inline double getOmegaModif(int hyb_pop)  // hyb_pop = 1 to i_pop < nDemes+1
{
    //sort(distance_list.begin(), distance_list.end());  already sorted ??!!
    // the curve y = x, should go through 0;0 and be 1 at max dist
    double slope = 1.0/distance_list.back();
    double dist = (distance_list[hyb_pop-1]/distance_list.back());  //relative distance
    return (dist*slope);
}

//nDemes was set to: nDemes = distance_list.size();

/* CONSTRUCTORS */
Model::Model() 
{
    Model(model, nDist, nSelGenes);    // by default retrieve "globally" defined parameters (from settings.h, initialized in main.cpp)
}
Model::Model(string model) 
{
    Model(model, nDist, nSelGenes);    // by default retrieve "globally" defined parameters (from settings.h, initialized in main.cpp)
}

Model::Model(string model, int dist, int selGenes) : model_name(model), n_dist(dist), n_sel_genes(selGenes)
{
    /**/ // => FOR DEBUG
    ofstream outControl;
    outControl.open(controlFile.c_str(), fstream::app|fstream::out);  // extern controlFile from settings.h, set in main.cpp
    outControl << "Start initializing model with following parameters:" << endl;
    outControl << "Model: " << model << endl << "Number of dist: " << dist << endl;
    outControl << "Number of genes under selection: " << selGenes << endl;
    /**/ // <= FOR DEBUG
    
    //DEFINE vector<string> list_pop_name, //DEFINE vector<char> list_pop_lab // no need for lab ?
    
    // for each distance (n_dist), I have many pop !! 1a 1b, etc. 
    
    int k(1);
    while(k < distance_list.size() + 1)  //hyb0 is now 'hybFoo' with 'z' lab  // +2 or +1 ???
    {
        //cout << "nbr_pop_list[k-1] : " << nbr_pop_list[k-1] << endl;
        for(int j = 0; j < nbr_pop_list[k-1]; j++)
        {
            char ith_pop = (j < 26) ? 97 + j: 65 + j - 26;        
            
            list_pop_name.push_back("hyb"+to_string(k)+ith_pop);   //hyb1a, hyb1b,..., hyb1A etc., hyb2a, ...
      //  list_pop_lab.push_back('a'+(k-1));                 //a, b, c, ...
            outControl << "In Class Model: added to list_pop_name : " << list_pop_name.back() << endl;   
        }        
      // outControl << "In Class Model: added to list_pop_lab : " << list_pop_lab.back() << endl;        
        ++k;
    }             
    // AIM: FILL vector<string> list_genes_sel;
    // chromoF looks like: should be ordered (chr1, chr2...) and length in second fields, tab-separated
    //chr1	41381587
    //chr2	53008256
    // retrieve the length of the chromosome in chomoF -> this will allow to find correct "random" positions of the genes
    vector<int> chr_len;
    ifstream inFile;
    string line;
    int nField(0);
	inFile.open(chromoF.c_str());     // chromoF => extern, initialized in main.cpp
    while (getline(inFile,line))
    {        
        stringstream ss(line);        
        string field;
        nField = 0;
        while(getline(ss, field, '\t'))         // split the line in tab-separated fields
        {
            if (nField == 1)   // chromosome length at second position
            {
                chr_len.push_back(stoi(field));
                outControl << "In Class Model: added to chr_len: " << to_string(chr_len.back()) << endl;
                break;
             }
            nField++;
        } 
    } 
    inFile.close();
    // now I have to put in list_genes_sel string like chr1_30 chrNumber_position
    // vectore must be ordered !!!
    map<int, vector<int>>temp_map_genes;
    int randChr, randPos;
    int foundGenes(0);
    while(foundGenes < selGenes +1) { // +1 because 1 mandatory for sex !! (and +1 because < not <=)
        //find randomly the chromosome on which the gene is located
        seed += temp_map_genes.size();
//        randChr = (int) round(getFromUnif(1, chr_len.size()));     
        randChr = (int) round(gsl_ran_flat(rd, 1, chr_len.size())); //  rd extern in settings.h
        outControl << "In Class Model: search position on chromosome nbr: " << to_string(randChr) << endl;  
        // find random position on this chromosome that has not already a gene
        int t = 0;
        do
        { 
            outControl << "In Class Model: ... searching position ... " << endl;
            //randPos = (int) round(getFromUnif(2, chr_len[randChr]));                     //!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            randPos = (int) round(gsl_ran_flat(rd, 2, chr_len[randChr])); //  rd extern in settings.h
            outControl << "FIND POS: " << randPos << endl;
            t++;
            seed += t;
        } while(find(begin(temp_map_genes[randChr]), end(temp_map_genes[randChr]), randPos) != end(temp_map_genes[randChr]));
        //now we are sure the position is unique, add this to the map
        temp_map_genes[randChr].push_back(randPos);
        foundGenes++;
        outControl << "In Class Model: add to temp_map_genes[" << to_string(randChr) << "]: " << to_string(temp_map_genes[randChr].back()) << endl;  
    }
    outControl << "In Class Model: finish to add to temp_map_genes " << endl;
    // now we have to put the list of genes ordered by chromosome and length into list_genes_sel     
    for(int n_chr = 0; !temp_map_genes.empty(); n_chr++)
    {
        sort(temp_map_genes[n_chr].begin(), temp_map_genes[n_chr].end());   // positions must be sorted on the chromo        
        for (auto gene_pos: temp_map_genes[n_chr])
        {
            string geneID("chr"+to_string(n_chr)+"_"+to_string(gene_pos));
            list_genes_sel.push_back(geneID);   
            outControl << "In Class Model: add to list_genes_sel: " << list_genes_sel.back() << endl;        
        }    
        temp_map_genes.erase(n_chr);
    }   
    /**/ // => FOR DEBUG
    outControl.close();
    /**/ // <= FOR DEBUG
}
/* GETTERS */
// Function to get model name
string Model::getModelName() const
{
    return model_name;
}
/* FUNCTIONS TO WRITE CONFIGURATION FILES */
// Function to write the model specific files 
void Model::writeModelFiles() const
{
    this -> writeNaturalSelectionFile();
    this -> writePhenotypeFile();
    this -> writeSexualSelectionFile();
    this -> writeGeneFile();
}
/* Sexual selection file */ // => same for all models ??????????
void Model::writeSexualSelectionFile() const
{
    /* FOR ALL MODELS: no sexual selection */
    string sexOutFile("data/"+ model_name +"/"+model_name +"_"+sexSelFile);
    ofstream fOutFile;
    fOutFile.open(sexOutFile.c_str());
    fOutFile << "Population\tGen\tSelection" << endl;  //~~ print the header
    for(auto pop: list_pop_name)
    {
        fOutFile << pop << "\t-1\t1" << endl;
    }
	cout << "Sexual selection file written in: " << sexOutFile << endl;
    return;
}
/* Gene file */    
/* WARNING: must be sorted  by chromosome and then position! */
void Model::writeGeneFile() const
{
    // iterate over the genes to write the gene file
    string geneOutFile("data/"+model_name+"/"+model_name+"_"+geneFile);
	ofstream fOutFile;
	fOutFile.open(geneOutFile.c_str());
	fOutFile << "name\tChromosome\tposition\tDominant Allele\tDominant Allele Value\tRecessive Allele\tRecessive Allele Value\tMode\tDominant Freq Pop1\tDominant Freq Pop2" << endl; // header line
    // all with allele values 1 and 0 -> it is ok for the moment
    // all additive except sex chromosome
    // by deifinition we say that the first gene in list_genes_sel is the one for sex -> Hemizygous
    int count_char(0);
    int shift_ascii;
    char allele1, allele2;
    string mode, freq;
    for(auto gene: list_genes_sel)
    {
        string geneID(gene), chromo, pos, part;
        stringstream ss(geneID);
        int chr_field(0);
        while(getline(ss, part, '_'))         // split the geneID chr1_10 -> chr1 and 10
        {
            if (chr_field == 0) chromo = part.replace(0, 3, "");    // retrieve 1 from chr1
            else if (chr_field == 1) pos = part;
            chr_field++;
        } 
        if(count_char == 32) count_char += 32 ;
        allele1 = 65 + count_char;      //ASCII code A and a, B and b, etc. // PROBLEM IF TOO MUCH GENES !!!!!!!
        allele2 = allele1 + 32;
        mode = count_char == 0 ? "Hemizygous": "Additive"; 
        freq = count_char == 0? "\t0.5\t0.5": "\t0\t1";
// DOMINANT ALLELE VALUE HAS FREQ OF 0 IN POP1 AND A FREQ OF 1 IN POP2
// BECAUSE OF THE WAY WE DEFINED SELECTION, WE PUT SELECTION TOWARDS BEING HOMOZYGOTE TOWARDS DOMINANT ALLELE
// AND SO WE WANT TO HAVE SELECTION IN THE DIRECTION OF RAPA (POP2)
// CAPITAL LETTER -> 1 AND FREQ OF CAPITAL LETTER -> 1 IN POP2 ("RAPA")
        fOutFile << gene << "\t" << chromo << "\t" << pos << "\t" << allele1 << "\t1\t" << allele2 << "\t0\t" << mode << freq << endl;
        count_char++;
    }	
	cout << "Gene file written in: " << geneOutFile << endl;
    return;
}
/* Natural and phenotype selection file */
void Model::writePhenotypeFile() const
{
    // BY DEFINITION, FOR ALL MODELS: WE SAY THAT THE SEX CHROMOSOME IS THE FIRST OF THE list_genes_sel VECTOR
    string phenotypeOutFile("data/"+model_name +"/"+ model_name +"_"+phenoFile);
    ofstream phenoOutFile;
    phenoOutFile.open(phenotypeOutFile.c_str());

    // print the header for phenotype file
    // and print also the 1st line for sex, by convention, we say that it is the first gene of list_genes_sel
    phenoOutFile << "Phenotypes\tFormula" << endl;   // ~~ print the header
    phenoOutFile << "Sex\t" << list_genes_sel[0] << endl;        // needed to start simulations !!!     

    // IN BOTH MODEL, JUST DEFINE THAT AA -> 1, Aa -> 0.5, aa -> 0 FOR ALL GENES UNDER SELECTION
    // "A" WAS SET TO HAVE A FREQ OF 0 IN POP1 (NIGRA) AND 1 IN POP2 (RAPA)
    // SO THAT OPTIMAL PHENOTYPE WILL JUST BE THE TOTAL NUMBER OF GENES UNDER SELECTION - EX: 
    // Pheno0   if(chr1_10==2, 1, if(chr1_10==1, 0.5, if(chr1_10==0, 0, 0)))
    // SO IN NATURAL SELECTION WE SUM THE PhenoX AND COMPARE THE VALUE TO 1*THE NUMBER OF SELECTED GENES (= OPTIMAL PHENOTYPE)
    for(int i_g = 1; i_g < list_genes_sel.size(); i_g++)      
    {
        phenoOutFile << "Pheno" << to_string(i_g-1) << "\t";  // start with Pheno0
        phenoOutFile << "if(" + list_genes_sel[i_g] + "==2,1,";
        phenoOutFile << "if(" + list_genes_sel[i_g] + "==1,0.5,";
        phenoOutFile << "if(" + list_genes_sel[i_g] + "==0,0,0)))" << endl;         
    }
    phenoOutFile.close();
	cout << "Phenotype file written in: " << phenotypeOutFile << endl;    
    return;
}
    
void Model::writeNaturalSelectionFile() const
{    
    // NATURAL SELECTION FILE:       
    string naturalOutFile("data/"+model_name+"/"+model_name+"_"+natSelFile);  //natSelFile -> "natSel.txt" main.cpp, used in simu.cpp cfg
    
    ofstream natOutFile;
    natOutFile.open(naturalOutFile.c_str());

    // print the header natural selection file
    natOutFile << "Population\tGen\tSelection" << endl;  //~~ print the header    
    
    // print the parental pop1 and pop2 lines and hybFoo -> no selection
    for(int i = 0; i < 3; i++)
    {
        natOutFile << list_pop_name[i] << "\t-1\t1" << endl;
    }       
    // for the hybrid pops, use the formula (P-P_opt)²/w² => where we draw omega (flat curve far from the sea, tight near see)
    //AS I USED 1 for AA -> OPTIMAL PHENOTYPE IS JUSTE THE NUMBER OF GENES UNDER SELECTION
    string sumOfPheno("Pheno0");
    for(int i_g = 2; i_g < list_genes_sel.size(); i_g++)
    {
        sumOfPheno += ("+Pheno"+to_string(i_g-1));          // P = Pheno0+Pheno1+Pheno2+..., if all homozygotes rapa -> sum to 1*nbr of genes
    }      

//    double omega0 = getFromUnif(0.04, 100.0);             // Rhoné: 1/w = 0.01 neutral, 1/w = 25 very strong

    double omega0 = gsl_ran_flat(rd, 0.04, 100.0);   // Rhoné: 1/w = 0.01 neutral, 1/w = 25 very strong

    // WRITE NATURAL SELECTION -> ADVANTAGE FOR THE RAPA GENE -> THIS DOES NOT CHANGE ACROSS SIMU
    // BUT FOR ALL MODELS, THE VALUE OF OMEGA CHANGES ACROSS THE POPULATIONS
    // THE SAME VALUE OF OMEGA FOR ALL POPULATIONS THAT ARE AT THE SAME DISTANCE


    for(int i_pop = 1; i_pop < distance_list.size() + 1; i_pop++)  //hyb pop start with hyb1, hyb2, etc.
    {
        string pop = "hyb"+to_string(i_pop);
        string omega = to_string(getOmegaModif(i_pop)*omega0);  // omega is the same for all pop that are at the same distance
        
        for(int j = 0; j < nbr_pop_list[i_pop-1]; j++)  // retrieve the number of pop that are at that distance
        {
            char ith_pop = (j < 26) ? 97 + j: 65 + j - 26;        // hyb1a, hyb1b, etc.
            natOutFile << pop+ith_pop << "\t-1\t";
            natOutFile << "exp(-pow(" << sumOfPheno << "-" << list_genes_sel.size()-1 << ", 2)/"; //-1 because of sex gene
            natOutFile << "pow(" << omega << ", 2))" << endl;
            
        }
        if(i_pop==1) allParametersLine += "omega_"+pop+"="+omega;   // this is to retrieve parameter values for the ABC
        else allParametersLine += ",omega_"+pop+"="+omega;
    }
    
    natOutFile.close();    
    
	cout << "Natural selection file written in: " << naturalOutFile << endl;

    return;
}




/* COMMODITY FCTS, VARIA */

// Print out the model info
void Model::printModel(ostream &stream) const
{
    stream << "Current model:\t" << model_name << endl;
    stream << "Nbr of classes of distance:\t" << n_dist << endl;
//    stream << "Following populations are in the model (Name\tLabel):" << endl;
    stream << "Following populations are in the model (Name):" << endl;          // without lab
    for(int i = 0; i < list_pop_name.size() ; i++)
    {
        //stream << list_pop_name[i] << "\t" << list_pop_lab[i] << endl;  // no need lab ?
        stream << list_pop_name[i] << endl;
    }
    stream << "Nbr of genes under selection:\t" << n_sel_genes << endl;
    stream << "Genes under selection are at the following positions:" << endl;
    for(auto gene: list_genes_sel)
    {
        stream << gene << endl;
    }
}

int Model::getTotalHybPop() const
{
    return (list_pop_name.size() - 3);
}


//destructor: nothing to do for the moment, no memory allocation
Model::~Model() 
{}





//********************************************************* NOT CLASS MEMBER FUNCTIONS (stream operator overload)

ostream &operator<<( ostream &stream, Model const &model)
{
    model.printModel(stream) ;
    return stream;
}
