#include "settings.h"
#include "summarystat.h"
#include "my_math.h"        /* draw from unif for sampling */
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <vector>
#include <algorithm>  /* sort*/
#include <gsl/gsl_randist.h> /* draw from distribution */
#include <sstream>      /* sstream */
#include <fstream>
#include <string>
#include <stdlib.h>     
#include <stdio.h>


using namespace std;

//TODO : subsampling !!!!!!!

///////////////////////////////////// ASSUMPTIONS
// genotype file: gene alleles written A/a, B/b, C/c etc. [65+i/97+i, etc.]
// parental pop labels: x and y; outputed is the frequency of "x" ancestry
// in genotype file, genes ordered by chromosome then position (required by admixem)
// marker data file, ordered by chromosome then position


///////////////////////////////////// ONLY USED IN THIS FILE 
inline string popNbr2Str(int nbr) {   
	return (nbr == 1 || nbr == 2) ? "pop"+to_string(nbr): (nbr == 3) ? "hybFoo" : "hyb"+to_string(nbr - 3);
}
struct Locus
{
    string id;
    string chromo;
    string pos;
};

struct IndivPopId
{
    string pop;
    string ind;
};


void printLocusCount(map <int, map<string, int>> &, queue<Locus> &, ostream &);
void printLocusValues(map <string, queue <string>> & locusM, queue<Locus> & locusQ, ostream & out);
queue<IndivPopId> subsampleByPop(string admixGeneFile);
bool cmp_loci( const Locus &, const Locus &);  
bool cmp_id( const IndivPopId & , const IndivPopId & );  //(l1, l2) return true if l1 smaller than l2

/////////////////////////////////////////////////////////////////////////////////////////// CALLED FROM simu.cpp
/* FUNCTION TO WRITE POPULATION THE CHROMOSOME BLOCKS BY CHROMOSOME AND INDIVIDUALS - TODO WHAT IS THE SUMMARY STAT ???? */
// WILL NEED CHANGES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// AS IMPLEMENTED NOW, ANCESTRY xx or yy TREATED IN THE SAME MANNER
void writeChromosomeBlocks(string simIterN) 
{
    ofstream outControl;
    outControl.open(controlFile.c_str(), fstream::app|fstream::out);  // extern controlFile from settings.h, set in main.cpp
    
    string outBlockFile("simu/"+model+"/"+simIterN+"/"+simIterN+"_SS_chromosomeBlockFile.txt");   
    ofstream outFile;
    // input is the allele value table output of from writeSummaryTables
    string valueTableFile("simu/"+model+"/"+simIterN+"/"+simIterN+"_SS_indValueTable.txt");   
    bool newChromo;
    ifstream inFile;
    string line, field, currChromo, chromoID;
    string geneID;
    int nField(0), nLine(0);
	inFile.open(valueTableFile.c_str());     
	
	if(!inFile.is_open())
	{
	    cout << "Unable to open: " << valueTableFile << endl;
	    exit(EXIT_FAILURE);
	}
	
    vector<string> individuals, current_values;     // store list of individuals ID retrieve from line 0
    // all_blocks[chromoID][individualID] = vector_size_of_homozygous_blocks
    map < string, map <string, vector<int>>> all_blocks;
    //string * current_values;     // temporary table will hold only the value of 1 line (sizie not known at compilation -> size of <individuals>
    while (getline(inFile,line))
    {        
        outControl << "In chromosome blocks - starting with line: " << to_string(nLine) << endl;
        stringstream ss(line);        
        nField = 0;
        if(nLine == 1) current_values = individuals; // so I have a vector of the correct size...
        while(getline(ss, field, '\t'))        
        {
            if(nField == 0) outControl << field << endl;
                
            if(nLine == 0 && nField ==0) //"Loci" -> do nothing
            {
                nField++;
                continue; //go to next field
            }
            if(nLine == 0 && nField > 0) // in the header line, store names of the individuals
            {   
                individuals.push_back(field); 
                outControl << "Individual: " << field << endl;

            } else
            {
                if(nField==0)  //locus ID field -> retrieve chromosome ID // ASSUME locus ID like "chr1_100" -> "chr1"
                {
                    currChromo = field.replace(field.find("_"), field.size()-field.find("_"), ""); //"chr#_#" -> "chr#"
                    newChromo = !(currChromo == chromoID);
                    chromoID = currChromo;
                    if(newChromo) outControl << "Find chromo: " << chromoID << endl;
                } else // fields with allele values 0, 0.5 or 1
                {
                    if(newChromo) //if new chromo, necessarily new blocks -> initialize a new block only if homozygous
                    {

outControl << "new chromo - individual: " << individuals[nField-1] << " see: " << field << endl;
                        if(field == "1")  all_blocks[chromoID][individuals[nField-1]] = {1};
                        else if(field == "0")  all_blocks[chromoID][individuals[nField-1]] = {1};
                        else all_blocks[chromoID][individuals[nField-1]] = {0};
                    } else
                    {
outControl << "2nd chromo - individual: " << individuals[nField-1] << " previous was: " << current_values[nField-1]<< " see: " << field << endl;
                        if(field != "0.5") // do sth only if homozygous
                        {                    
                            // same genotypes as previously -> +1 size of the block (no empty() test, because cannot be empty at this point)
                            if(field == current_values[nField-1]) all_blocks[chromoID][individuals[nField-1]].back()++;
                            // we enter a new homozygous block -> add a block of 1 size    
                            else all_blocks[chromoID][individuals[nField-1]].push_back(1);      
                        } 
                    } 
                    current_values[nField-1] = field; //whatever homo/heterozygous, new chromosome or not -> store it for next iteration
                } 
            } 
            nField++;
        } 
        nLine++;
    } 
    inFile.close();
//    delete[] current_values;
    // WRITE IN A FILE -> THIS IS NOT A TABLE !!!! //
	outFile.open(outBlockFile.c_str());     
    // new chromosome -> 1 line with the chromosome ID, then 1 line per individual, tab-separated the size of the blocks
    // no difference between homozygous xx or yy (1 or 0)
    for(auto chromo : all_blocks)
    {
        outFile << chromo.first << endl;
        for(auto ind : chromo.second)
        {
            outFile << ind.first;
            for(auto block : ind.second)
            {
                if(block > 0) outFile << "\t" << block;
            }
            outFile << endl;
        }
    }
    outFile.close();
    outControl.close();
    return;
}

/////////////////////////////////////////////////////////////////////////////////////////// CALLED FROM simu.cpp
/* FUNCTION TO WRITE POPULATION ALLELE FREQUENCY FROM 'x' FOR ALL LOCI AND WRITE ALLELE TABLE FOR ALL INDIVIDUALS */
// print allele freq by pop and the allele table for all individuals in the same function
// because they use the same loci queues 
void writeSummaryTables(string admixGeneFile, string admixMarkerFile, string simIterN)
{
    ofstream outControl;
    outControl.open(controlFile.c_str(), fstream::app|fstream::out);  // extern controlFile from settings.h, set in main.cpp
    /********************************************************************/
    /******* READ THE GENE DATA FILE TO STORE GENE LOCI IN QUEUE ********/
    /********************************************************************/
    // gene file written in model.cpp 
    queue<Locus> all_genesQ;
    string line, field, nChr, pos;
    int nField(0), nLine(0);
    string geneDataFile("data/"+model+"/"+model+"_"+geneFile);
    ifstream geneDataINfile;        
	geneDataINfile.open(geneDataFile.c_str());     
	if(!geneDataINfile.is_open())
	{
	    cout << "Unable to open: " << geneDataFile << endl;
	    exit(EXIT_FAILURE);
	}
    while (getline(geneDataINfile,line))
    {        
        if(nLine == 0){
            nLine++;
            continue;    // skip the header
        }
        stringstream ss(line);        
        nField = 0;
        while(getline(ss, field, '\t'))         
        {
            if (nField == 1)   nChr = field;    // chromosome number
            if (nField == 2)   pos = field;     // position on the chromosome
            if(nField == 3)
            {
                all_genesQ.push(Locus {field, nChr, pos});
                break;                          // no need to see next fields of the line
             }
            nField++;
        } 
        nLine++;
    } 
    geneDataINfile.close();    
    
    outControl << "PRINT THE QUEUE OF GENES INFO" << endl;
    queue<Locus> genesQprint = all_genesQ;
    while(!genesQprint.empty())
    {
        Locus i = genesQprint.front();
        outControl << "chromo\t" << i.chromo << "\tpos\t" << i.pos << endl;
        genesQprint.pop();
    }
    /************************************************************************/
    /******* READ THE MARKER DATA FILE TO STORE MARKER LOCI IN QUEUE ********/    
    /************************************************************************/
    // not model specific, created in createGeneralSettings    ("markerF" from settings.h, set in main.cpp)
    queue<Locus> all_markersQ;
    nField = nLine = 0;
	int chromo_nbr(0), nMarker(0);           // for Locus: need chromo and ID (no marker ID, so just increment nMarker)
    bool chromoline(false);
    string chromo_id;
    ifstream markerDataINfile;
	markerDataINfile.open(markerF.c_str());     
	
	if(!markerDataINfile.is_open())
	{
	    cout << "Unable to open: " << markerF << endl;
	    exit(EXIT_FAILURE);
	}	
	
    while (getline(markerDataINfile,line))
    {        
        if(nLine < 6 || line.length() < 1)      //skip header lines of the file or empty lines
        {
            nLine++;
            continue;
        }
        if(line[0] == ':')
        {
            chromo_nbr++;
            chromoline = true;                      // we will have to retrieve in which chromosome we are from the 1st field
        }
        nField = 0;
        stringstream ss(line);
        string field;        
        while(getline(ss, field, '\t'))         
        {
            if(field == "Num") break;                // header line of next chromosome, skip the line
            if(chromoline)
            {
                chromo_id = field.replace(0, 5, "");  //:chr 1 -> 1, store in which chromo we are and skip the line
                chromoline = false;
                break;                          
            }
            if(nField == 2)
            {
                nMarker++;      
                all_markersQ.push(Locus{to_string(nMarker), chromo_id, field}); 
                break;            
            }        
            nField++;
        } 
        nLine++;    
    } 
    markerDataINfile.close();

    outControl << "PRINT THE QUEUE OF MARKERS INFO" << endl;
    queue<Locus> markersQprint = all_markersQ;
    while(!markersQprint.empty())
    {
        Locus i = markersQprint.front();
        outControl << "chromo\t" << i.chromo << "\tpos\t" << i.pos << endl;
        markersQprint.pop();
    }


/*
    /********************************************************************/
    /*********************** WORK FOR FREQUENCY TABLE *******************/
    /********************************************************************/  
     
    /******* READ THE GENES DATA FROM ADMIXEM, STORE BY POP IN A MAP ********/
    //all_genes_count[popID][geneID] = count_x_ancestry_of_geneID_in_popID
    // admixGeneFile received in parameter
    string individualID, geneID, ancestry, indID;
    //``` => HERE: for the sampling effort:  <=
    int sample_idx(0);  
    const queue<IndivPopId> indivToSample = subsampleByPop(admixGeneFile);  // already ordered by 
    
    queue<IndivPopId> indiv_for_freqG(indivToSample);

    queue<IndivPopId> indiv_for_freq_temp(indivToSample);
    
    string outSamplingFile("simu/"+model+"/"+simIterN+"/"+simIterN+"_SS_sampling.txt");   
    ofstream outSamp;
    outSamp.open(outSamplingFile.c_str());   
    
    outSamp << "pop\tind" << endl;
    while(!indiv_for_freq_temp.empty())

    {
        outSamp << indiv_for_freq_temp.front().pop << "\t" << indiv_for_freq_temp.front().ind << endl;    
        indiv_for_freq_temp.pop();
    }
   
    int popID;
    nField = nLine = 0;
    ifstream admixFileG;
	admixFileG.open(admixGeneFile.c_str());         //<pop, <gene, count>>
	
	if(!admixFileG.is_open())
	{
	    cout << "Unable to open: " << admixGeneFile << endl;
	    exit(EXIT_FAILURE);
	}	
	
	bool newInd;
	
    map <int, map<string, int>> all_genes_count; //popID as integer -> easier to iterate over for writing in file
//100	1	1	447319	447437	-1	Ax	0.000000	-1	-1	-1	-1	Bx	0.000000	-1	-1	-1	-1	-1	-1	-1	-1	Cy    
    while (getline(admixFileG,line))
    {        
        outControl << "Admixem gene file - reading line: " << nLine << endl;
        stringstream ss(line);        
        nField = 0;
        while(getline(ss, field, '\t'))         
        {             
            if(nField == 0)    // look if the current individual is in subsampling
            {       
                newInd = !(field == indID); // ? false: true;            
                indID = field;
                nField++; 
                continue;


            }

            if(nField == 1)                     // second column is the popID
            {
                popID = stoi(field);             //better to integer for writing final file 
                if(popID < 4) break;              // we don't want to sample pop1 (1), pop2 (2), hybFoo(3) !!!!!!!!!!!!!                
                
                if(indiv_for_freqG.empty()) break;
                
                outControl << ">LINE " << nLine << " see_freqG " << field << "_" << indID << endl;
                outControl << "CURRENTLY_HOLD " << indiv_for_freqG.front().pop << "_" << indiv_for_freqG.front().ind  << endl;
                
                if(field == indiv_for_freqG.front().pop && indID == indiv_for_freqG.front().ind)
                {
                    outControl << "=>to_sample_freqG " << endl;
                    nField++;
                    if(!newInd)   // if we see the 2nd line
                    { 
                        indiv_for_freqG.pop();
                    }
                    continue;
                } else {
                    break; // go to next line
                }
                
            }            
            
            //skip the fields before 1st gene, the ones with allelic value, and -1 (indicates new chromo)        
            if(nField < 6 || field.length() > 2 || field == "-1" )
            {
                nField++;
                continue;
            }
            // im interested in ancestry, Cx or cx -> does not matter, so convert c -> C before storing (in the queue, C is stored)
            geneID = field[0] > 96 ? string(1, field[0]-32): string(1, field[0]); 

            ancestry = field.replace(0,1,"");       //retrieve ancestry: Ax -> x
            if(ancestry == "x") all_genes_count[popID][geneID]++;    //store the count, freq retrieved later
            else all_genes_count[popID][geneID] +=0 ;           
            outControl << "Pushed in map: pop " << popID << "-gene " << geneID << "++" << endl;
            
            nField++;
        } 
        nLine++;
    } 
    admixFileG.close();      
   
    
/// WARNING: IF NO ALLELE 'X' AT ALL FOR A GIVEN GENE IN A POP, THIS GENE WILL NOT BE IN THE all_genes_count MAP !!!    
    
    outControl << "*** FQCY BY POP Admixem genes - Iterate over all_genes_count ***" << endl;
    for(auto i : all_genes_count)
    {
        outControl << "For pop : " << i.first << endl;
        for(auto j: i.second)
        {
            outControl << "for gene: " << j.first << "\t number of alleles found: " << j.second <<  endl;
        }
    }

outControl << "***********************************************************************************" << endl;
  
    /******* READ THE MARKERS DATA FROM ADMIXEM, STORE BY POP IN A MAP ********/
    //all_markers_count[popID][markerNbr] = count_x_ancestry_of_markerNbr_in_popID
    // admixMarkerFile received in parameter

    queue<IndivPopId> indiv_for_freqM(indivToSample);
    
    ifstream admixFileM;
	admixFileM.open(admixMarkerFile.c_str());
	
	if(!admixFileM.is_open())
	{
	    cout << "Unable to open: " << admixMarkerFile << endl;
	    exit(EXIT_FAILURE);
	}		
	
	outControl << "size of ind to sample: " << to_string(indivToSample.size()) << endl;
	
    map<int, map<string, int>> all_markers_count;   // <pop, <marker, count>>
	sample_idx = nField = nLine = 0;
//100	1	1	49088	49363	-1	Ax	-1	Ax	-1	-1	-1	Ay	-1	-1	-1	-1	-1	Ax	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	Ay	-1	
    while (getline(admixFileM,line))
    {        
        outControl << "Admixem marker file - read line : " << nLine << endl;
        stringstream ss(line);        
        nMarker = nField = 0;
        while(getline(ss, field, '\t'))
        {
        
            if(nField == 0)    // look if the current individual is in subsampling
            {                             
                newInd = !(field == indID); // ? false: true;
                indID = field;
                nField++; 
                continue;
            }
            if(nField == 1)                     // second column is the popID
            {
                popID = stoi(field);             //better to integer for writing final file 
                if(popID < 4) break;              // we don't want to sample pop1 (1), pop2 (2), hybFoo(3) !!!!!!!!!!!!!                
                
                if(indiv_for_freqM.empty()) break;
                
                outControl << ">LINE " << nLine << " see_freqM " << field << "_" << indID << endl;
                outControl << "CURRENTLY_HOLD " << indiv_for_freqM.front().pop << "_" << indiv_for_freqM.front().ind  << endl;
                
                if(field == indiv_for_freqM.front().pop && indID == indiv_for_freqM.front().ind)
                {
                    outControl << "=>to_sample_freqM " << endl;
                    nField++;
                    if(!newInd)   // if we see the 2nd line
                    { 
                        indiv_for_freqM.pop();
                    }
                    continue;
                } else {
                    break; // go to next line
                }
                
            }        
                       
            if(nField < 6 || field == "-1")     // skip fields before the 1st marker and the -1 (indicates new chromosome)    
            {
                nField++;
                continue;
            }
            ancestry = field.replace(0,1,"");       //Ax -> only want ancestry x
            if(ancestry == "x") all_markers_count[popID][to_string(nMarker+1)]++;    //1-based numbering of the markers in the queue
            else all_markers_count[popID][to_string(nMarker+1)] +=0 ;
            outControl << "ADDED to marker map - pop " << popID << "marker number " << to_string(nMarker+1) << "++" << endl;
            
            nField++;
            nMarker++;
        } 
        nLine++;    
    } 
    admixFileM.close();    
    
    outControl << "*** FQCY BY POP - Admixem Iterate over all_markers_count ***" << endl;
    for(auto i : all_markers_count)
    {
        outControl << "For pop : " << i.first << endl;
        for(auto j: i.second)
        {
            outControl << "for marker: " << j.first << "\t number of alleles found: " << j.second <<  endl;
        }
    }    
    

    outControl << "Write frequency table : " << endl;  
    /******* NOW WRITE THE FREQUENCY TABLE ********/   
    // we want to have them ordred by position
    // 2 queues of Loci (genes and markers) are available   
    // check if we have to place a locus or a marker, pop the queue
    //simIterN received in parameter, store in simu/model0/simu1001/simu1001_SS_popFreqTable.txt
    string outputfile("simu/"+model+"/"+simIterN+"/"+simIterN+"_SS_popFreqTable.txt");   
    ofstream outTable;
    outTable.open(outputfile.c_str());   
   // print a header with population name
    outTable << "Loci";
    outControl << "Loci"; 
    for(int i = 0; i < all_genes_count.size(); i++) // we don't want to sample pop1 (1), pop2 (2), hybFoo(3) !!!!!!!!!!!!!
    {
        outTable << "\t" << popNbr2Str(i+4);
        outControl << "\t" << popNbr2Str(i+4);
    }   
    outTable << endl;
    outControl << endl;
    

    // work on copy of the queues because also needed later
    queue <Locus> fqcy_markersQ = all_markersQ;
    queue <Locus> fqcy_genesQ = all_genesQ;
    // iterate over the queues until one is empty
    while(!(fqcy_markersQ.empty() || fqcy_genesQ.empty()))
    {
        if(cmp_loci(fqcy_markersQ.front(), fqcy_genesQ.front()))   //if the next to put is a marker
        {     
            printLocusCount(all_markers_count, fqcy_markersQ, outTable);
        } else
        {
            printLocusCount(all_genes_count, fqcy_genesQ, outTable);
        }
    }
    while(!fqcy_markersQ.empty())            // if there are remaining markers
    {
        printLocusCount(all_markers_count, fqcy_markersQ, outTable);
    } 
    while(!fqcy_genesQ.empty())               //if there are remaining genes
    {
        printLocusCount(all_genes_count, fqcy_genesQ, outTable);
    }
    
    outTable.close();
    
    cout << "Allele frequency table for all loci per population written in: " << outputfile << endl;


    outControl << "********************************* STARTING ALLELE VALUE TABLE **********************" << endl;



    /********************************************************************/
    /*********************** WORK FOR ALLELE VALUE TABLE ****************/
    /********************************************************************/
     
    /******* READ THE GENES DATA FROM ADMIXEM, STORE IN A MAP ALL INDIVIDUAL VALUES ********/
    //all_genes_values[geneID] = queue_of_allele_values
    // admixGeneFile received in parameter
    
    string myGene("mygenefile.txt");
    ofstream myG;
    myG.open(myGene.c_str());
    

    queue<IndivPopId> indiv_for_valueG(indivToSample);


    queue<IndivPopId> sampling_temp(indivToSample);
    
    while(!sampling_temp.empty())
    {
        cout << "pop: " << sampling_temp.front().pop << "\t";
        cout << "ind: " << sampling_temp.front().ind << endl;
        sampling_temp.pop();
    
    }

    
    ifstream admixFileG2;
    
    // nField, nLine, sample_idx;   ////////////////////////////////////////////////////////////////////////////////////////////////
//    string line, field;////////////////////////////////////////////////////////////////////////////////////////////////
    nField = nLine  = 0;
	admixFileG2.open(admixGeneFile.c_str());     
	
//	int popID; 
//	string indID, ancestry, geneID;   ///////////////////////////////////////////////////////////////////////////////////////////////
	
	if(!admixFileG2.is_open())
	{
	    cout << "Unable to open: " << admixGeneFile << endl;
	    exit(EXIT_FAILURE);
	}		
	
	//string indID;
    //bool newInd;
    map<string, queue<string>> all_genes_values;
    queue <string> temp_ancestryG;
    
    outControl << "admixFileG2: " << admixGeneFile << endl;
    
    while (getline(admixFileG2,line))
    {        
        myG << line << endl;
        stringstream ss(line);        
        nField = 0;
                      	           
        while(getline(ss, field, '\t'))         // split the line in tab-separated fields
        {
            if(nField == 0)    // look if the current individual is in subsampling
            {                             
                newInd = !(field == indID); // ? false: true;
                indID = field;
                nField++; 
                continue;
            }
            if(nField == 1)                     // second column is the popID
            {
                popID = stoi(field);             //better to integer for writing final file 
                if(popID < 4) break;              // we don't want to sample pop1 (1), pop2 (2), hybFoo(3) !!!!!!!!!!!!!                
                
                if(indiv_for_valueG.empty()) break;
                
                outControl << ">LINE " << nLine << " see_valueG " << field << "_" << indID << endl;
                outControl << "CURRENTLY_HOLD " << indiv_for_valueG.front().pop << "_" << indiv_for_valueG.front().ind  << endl;
                
                if(field == indiv_for_valueG.front().pop && indID == indiv_for_valueG.front().ind)
                {
                    outControl << "=>to_sample_valueG " << endl;
                    nField++;
                    if(!newInd)   // if we see the 2nd line
                    { 
                        indiv_for_valueG.pop();
                    }
                    continue;
                } else {
                    break; // go to next line
                }
                
            }
        // skip the 6 fields, the allelic value (>2) and the chromosome without genes (-1)
            if(nField < 6 || field.length() > 2 || field == "-1" )//|| nField > 6+2*lastChromo )  
            {
                nField++;
                continue;
            }
            geneID = field[0] > 96? string(1, field[0]-32): string(1, field[0]); // I only store 1 allele, the 1st one in geneFile (65-96)
            ancestry = field.replace(0,1,"");       //string  if Am => only want ancestry info "m"            
            if(newInd)
            {
                temp_ancestryG.push(ancestry);
            } else{
                if(ancestry != temp_ancestryG.front())
                {
                    all_genes_values[geneID].push("0.5");
                } else 
                {
                    if(ancestry == "x")
                    {
                        all_genes_values[geneID].push("1");
                    } else
                    {
                        all_genes_values[geneID].push("0");
                    }
                }
                temp_ancestryG.pop();
            }
            nField++;
        } 
        nLine++;
    } 
    admixFileG2.close();      
    
outControl << "***** VALUE TABLE - Admixem genes ******" << endl;
for(auto i: all_genes_values)
{
    outControl << "FOR GENE " << i.first << " VALUES: " << endl;
    queue<string> tempQ = i.second;
    while(!tempQ.empty())
    {
        outControl << "\t" << tempQ.front();
        tempQ.pop();
    }
    outControl << endl;
}


    
    /******* READ THE MARKERS DATA FROM ADMIXEM, STORE IN A MAP WITH VALUES FOR EACH GENE ********/
    //all_markers_values[markerNbr] = queue_of_values
    // admixMarkerFile received in parameter
    
    string myMarker("mymarkerfile.txt");
    ofstream myM;
    myM.open(myMarker.c_str());
    
    queue<IndivPopId> indiv_for_valueM(indivToSample);    
    
    ifstream admixFileM2;
	admixFileM2.open(admixMarkerFile.c_str());
	
	    outControl << "admixFileM2: " << admixMarkerFile << endl;
	
	if(!admixFileM2.is_open())
	{
	    cout << "Unable to open: " << admixMarkerFile << endl;
	    exit(EXIT_FAILURE);
	}	
	
    map<string, queue<string>> all_markers_values;
	indID = "";
    nField = nLine = 0;
    queue <string> temp_ancestryM;            
    vector<string> indiv_header;
    while (getline(admixFileM2,line))
    {
        myM << line << endl;        
        nField = 0;
        stringstream ss(line);        
        nMarker = 0;
        while(getline(ss, field, '\t'))         // split the line in tab-separated fields
        {

            if(nField == 0)    // look if the current individual is in subsampling
            {                             
                newInd = !(field == indID); // ? false: true;
                indID = field;
                nField++; 
                continue;
            }
            if(nField == 1)                     // second column is the popID
            {
                popID = stoi(field);             //better to integer for writing final file 
                if(popID < 4) break;              // we don't want to sample pop1 (1), pop2 (2), hybFoo(3) !!!!!!!!!!!!!                
                                
                if(indiv_for_valueM.empty()) break;
                
                outControl << ">LINE " << nLine << " see_valueM " << field << "_" << indID << endl;
                outControl << "CURRENTLY_HOLD " << indiv_for_valueM.front().pop << "_" << indiv_for_valueM.front().ind  << endl;
                
                if(field == indiv_for_valueM.front().pop && indID == indiv_for_valueM.front().ind)
                {
//if(newInd) 
                    outControl << "=>to_sample_valueG " << endl;
                    nField++;
                    if(!newInd)   // if we see the 2nd line
                    { 
                        indiv_for_valueM.pop();
                    } else 
                    {
                        indiv_header.push_back(field+"_"+indID);
                    }
                    continue;
                } else {
                    break; // go to next line
                }
                
            }             
        
            if(nField < 6 || field == "-1")
            {
                nField++;
                continue;
            }
            ancestry = field.replace(0,1,"");       //string  if Am => only want ancestry info "m"
            if(newInd)
            {
                temp_ancestryM.push(ancestry);
            
            } else{
                if(ancestry != temp_ancestryM.front())
                {
                    all_markers_values[to_string(nMarker+1)].push("0.5");
                } else 
                {
                    if(ancestry == "x")
                    {
                        all_markers_values[to_string(nMarker+1)].push("1");
                    } else
                    {
                        all_markers_values[to_string(nMarker+1)].push("0");
                    }
                }
                temp_ancestryM.pop();
            }
            nField++;
            nMarker++;
        } 
        nLine++;    
    } 
    admixFileM2.close();    

myM.close();
myG.close();


outControl << "***** VALUE TABLE - Admixem markers ******" << endl;
for(auto i: all_markers_values)
{
    outControl << "FOR MARKER " << i.first << "VALUES: " << endl;
    queue<string> tempQ = i.second;
    while(!tempQ.empty())
    {
        outControl << "\t" << tempQ.front();
        tempQ.pop();
    }
    outControl << endl;
}

  
    
    
    
    /******* NOW WRITE THE ALLELE VALUE PER INDIVIDUAL TABLE ********/   
    string outputfileValues("simu/"+model+"/"+simIterN+"/"+simIterN+"_SS_indValueTable.txt");   
    ofstream outTableM;
    outTableM.open(outputfileValues.c_str());
    
    string indivPopFile("simu/"+model+"/"+simIterN+"/"+simIterN+"_indiv_pop_info.txt");   
    ofstream outIndivPop;
    outIndivPop.open(indivPopFile.c_str());
    
    
    
   // print a header
   outTableM << "Loci";
   for(auto i : indiv_header)   //4_1234 hyb1 120
   {
        outTableM << "\t" << i;
        int popNbr = stoi(i.substr(0, i.find("_")));  //inline string popNbr2Str(int nbr) {   
        outIndivPop << i << "\t" << popNbr2Str(popNbr) << "\t" << distance_list[popNbr-4] << endl;
        
    }   
    outTableM << endl;
    
    outIndivPop.close();
    
    // work on copy of the queues because also needed later and iterate over the queues until one is empty 
    queue <Locus> val_markersQ = all_markersQ;
    queue <Locus> val_genesQ = all_genesQ;
   
    while(!(val_markersQ.empty() || val_genesQ.empty()))
    {
        if(cmp_loci(val_markersQ.front(), val_genesQ.front()))   //if the next to put is a marker
        {     
            printLocusValues(all_markers_values, val_markersQ, outTableM);
        } else
        {
            printLocusValues(all_genes_values, val_genesQ, outTableM);
        }
    }  
    while(!val_markersQ.empty())            // if there are remaining markers
    {
        printLocusValues(all_markers_values, val_markersQ, outTableM);
    } 
    while(!val_genesQ.empty())               //if there are remaining genes
    {
        printLocusValues(all_genes_values, val_genesQ, outTableM);
    }
    cout << "Allele frequency table for all loci per population written in: " << outputfileValues << endl;   
    
    outControl.close();
    outTableM.close();
    return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////// COMMODITY FUNCTIONS
/* COMPARE TWO LOCI */

bool cmp_loci( const Locus & l1, const Locus & l2)  //(l1, l2) return true if l1 smaller than l2
{

//cout << "l1 chromo: " << l1.chromo << " l1 pos: " << l1.pos << endl;
//cout << "l2 chromo: " << l2.chromo << " l2 pos: " << l2.pos << endl;
  return l1.chromo != l2.chromo ? stoi(l1.chromo) < stoi(l2.chromo): stoi(l1.pos) < stoi(l2.pos);
}

//* COMPARE 2 INDIVIUDAL ID */

// sort: as in marker and gene Admixem files: 1st pop then individual ID // TO CONTROL
bool cmp_id( const IndivPopId & l1, const IndivPopId & l2)  //(l1, l2) return true if l1 smaller than l2
{
  return l1.pop != l2.pop ? stoi(l1.pop) < stoi(l2.pop): stoi(l1.ind) < stoi(l2.ind);
}



/* PRINT LOCI BY POP*/
//            printLocusCount(all_markers_count, fqcy_markersQ, outTable);
//all_markers_count[popID][to_string(nMarker+1)]
//    map<int, map<string, int>> all_markers_count;   // <pop, <marker, count>>
void printLocusCount(map <int, map<string, int>> & locusM, queue<Locus> & locusQ, ostream & out)
{
    ofstream outControl;
    outControl.open(controlFile.c_str(), fstream::app|fstream::out);  // extern controlFile from settings.h, set in main.cpp

    out << "chr" << locusQ.front().chromo << "_" << locusQ.front().pos ;
  //  cout << "chr" << locusQ.front().chromo << "_" << locusQ.front().pos ;
    for(int i = 0; i < locusM.size(); i++)  //i=0 -> pop 4
    {
        int sample_size = pop_sampling_list[i];
        outControl << "nominateur : " << (double) locusM[i+4][locusQ.front().id] << endl;
        outControl << "denominateur: " << to_string(sample_size*2) << endl;

        outControl << ((double) locusM[i+4][locusQ.front().id] /  (sample_size*2)) << endl;

        out << "\t" << ((double) locusM[i+4][locusQ.front().id] / (sample_size*2));
      //  cout << "\t" << ((double) locusM[i+1][locusQ.front().id] / (sample_size*2));
    }
    locusQ.pop();   
    out << endl;
   // cout << endl;
    outControl.close();
    return;
}

/* PRINT LOCI VALUES BY INDIVIDUALS */
//            printLocusValues(all_markers_values, val_markersQ, outTableM);
//    map<string, queue<string>> all_markers_values;
//                    all_markers_values[to_string(nMarker+1)].push("0.5");
void printLocusValues(map <string, queue <string>> & locusM, queue<Locus> & locusQ, ostream & out)
{
    out << "chr" << locusQ.front().chromo << "_" << locusQ.front().pos ;
    while(!locusM[locusQ.front().id].empty())
    {
        out << "\t" << locusM[locusQ.front().id].front() ;
        locusM[locusQ.front().id].pop();
    }
    locusQ.pop();
    out << endl;
    return;
}

///////////////////// DRAW sample_size INDIVIUDAL *PER POP* 
// RETURN A VECTOR OF IndivPopId, SORTED BY POP THEN INDIVIDUAL ID


queue<IndivPopId> subsampleByPop(string admixGeneFile)
{
cout << ">>>> YOU ARE IN SUBSAMPLE <<<< " << endl;
ofstream outControl;
outControl.open(controlFile.c_str(), fstream::app|fstream::out);  // extern controlFile from settings.h, set in main.cpp

    string indID, popID, field, line;
    int nField(0), nLine(0);
    ifstream admixFileG;
	admixFileG.open(admixGeneFile.c_str());         
    map <string, vector<string>> ind_by_pop; //ind_by_pop[popNbr] = {ind1, ind2, ...}

outControl << "START SUBSAMPLING " << endl;    

    while (getline(admixFileG,line))
    {        
        if(nLine % 2 !=0)       //read 1 of 2 lines because an individual = 2 lines
        {
            nLine++;
            continue;
        }
        stringstream ss(line);        
        nField = 0;
        while(getline(ss, field, '\t'))         
        {
            if(nField == 0) indID = field;
            if(nField == 1) 
            {
                popID = field;
                if(stoi(popID) < 4) break;              // we don't want to sample pop1 (1), pop2 (2), hybFoo(3) !!!!!!!!!!!!!
                ind_by_pop[popID].push_back(indID);
                break;
            }
            nField++;
        } 
        nLine++;
    } 
    admixFileG.close(); 
    //  nsample by pop ??????????????????????
    
    // the number of individuals to sample in each pop is stored in pop_sampling_list (sorted by pop: less distant to more distant)
    
    
    vector<IndivPopId> sampling;
    
    // pop number start at 4
    for(int i = 0; i < ind_by_pop.size(); i++)
    {
        string pop_id = to_string(i+4);
        int size_sample_pop = pop_sampling_list[i];
        outControl << "For pop: " << pop_id << ", number of individuals to sample: " << to_string(size_sample_pop) << endl;
        set<int> rand_nbr_set;

        while(rand_nbr_set.size() < size_sample_pop)
        {
//            rand_nbr_set.insert((int) getFromUnif(0, ind_by_pop[pop_id].size()-1));
            rand_nbr_set.insert((int) gsl_ran_flat(rd,0, ind_by_pop[pop_id].size()-1));
//            outControl << "Inserted individual in rand_nbr_set " << endl;
        }
        for (set<int>::iterator it=rand_nbr_set.begin(); it!=rand_nbr_set.end(); ++it)
        {
            sampling.push_back(IndivPopId{pop_id, ind_by_pop[pop_id][*it]});
//            outControl << "Push back to sampling" << endl;
        }
    }
    // return sorted vector
    //sort(sampling.begin(), sampling.end(), cmp_id);
    
    queue<IndivPopId> my_sampling;
    
    outControl << "***INDIVIDUALS_SAMPLED: " << endl;
    for(auto i: sampling)
    {
        outControl << i.pop << "\t" << i.ind << endl;
        my_sampling.push(i);
    }
    outControl.close();    
    return (my_sampling);
    
}

