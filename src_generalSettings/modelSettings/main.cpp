#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <string>
#include "maths.h"
#include "settings.h"

#pragma once

/*
 * GENERATE RECOMBINATION AND MARKER FILES => CODE RETRIEVED FROM ADMIXEM (Cui et al. 2015)
 */

using namespace std;

// RETRIEVE VARIABLES DEFINED ELSEWHERE FOR MATH FUNCTIONS
extern Normal NormalGen; 
extern Uniform UniformGen;


// RETRIEVE VARIABLES FROM "settings.h"                     => THIS IS WHAT WILL CHANGED FROM ONE MODEL TO ANOTHER !!!
extern std::string geneFile;                // file names
extern std::string recombFile;
extern std::string markerFile;
extern std::string phenotypeFile;
extern std::string sexualSelFile;
extern std::string naturalSelFile;

extern int numbersOfMarkers;                // markers file
extern int chromoSize;
extern int nChromo;
extern float pop1AlleleFreqAvg;
extern float pop1AlleleFreqStdDev;
extern float pop2AlleleFreqAvg;
extern float pop2AlleleFreqStdDev;
extern int nRecombMeiosis;                  // recombination file
extern int factSamplePoint;


// FUNCTION PROTOTYPES
void makeNewMarkerFile(string model);
void generateRecombinationFreqMap(int nChromosomes, double nLargestChromSize, double * pChromosomeSizes, double * pCentromerePos, string model);
void generateGeneFile(string model);
void generateNaturalSelFile(string model);
void generateSexualSelFile(string model);	
void generatePhenoFile(string model);

/*******************************
****** MAIN FUNCTION HERE ******
********************************/
int main (int argc, char * const argv[]) {

    string model = argv[1];

    // ~~ Generate the marker file and the recombination files (called within the function)
    makeNewMarkerFile(model); 

    // ~~ Generate the gene selection file
	generateGeneFile(model);
	
	// ~~ Generate the sexual selection file
	generateSexualSelFile(model);

	// ~~ Generate the natural selection file
	generateNaturalSelFile(model);
	
	// ~~ Generate the phenotype file
	generatePhenoFile(model);
	

    return 0;
}

/*******************************
***** FUNCTION DEFINITIONS *****
********************************/
/*
 * Fonction to generate the marker file
 */
void makeNewMarkerFile(string model) 
{
    //~~ Choice of some parameters for our simulations
    double nRandSeed = randomSeed;                        // number between 0 and 1 passed by user        
    int nHaploidChromosomeNum (nChromo);                  // number of haploid chromosome -> say 4 since B nigra 8n
    double nChromosomeSize (chromoSize);                   // size of the chromosome
 

    int nTotalMarkers (numbersOfMarkers);                      // total number of markers
    double nPop1AvgAlleleFreq (pop1AlleleFreqAvg);                //average freq. of the more abundant allele in pop 1 (0.5-1.0)
    double nPop1AvgAlleleFreqStdDev (pop1AlleleFreqStdDev);          //std. dev. of the freq of the more abundant allele in pop 1 (0-0.5)
    double nPop2AvgAlleleFreq (pop2AlleleFreqAvg);                //average freq. of the more abundant allele in pop 2 (0.5-1.0)
    double nPop2AvgAlleleFreqStdDev (pop2AlleleFreqStdDev);          //std. dev. of the freq of the more abundant allele in pop 2 (0-0.5)
    string sOutFile ("simu/"+model+"/"+model+"_"+markerFile);   // name of output file
    
    //double nCentromerePos;
    int nCurrChromosome (0);                        // index for filling chromosome size array

    // ~~ Set random seed     // random seed to repeat the simulation (-> no assignment is ok)
    Random::Set(nRandSeed);
   
    double * pChromosomeSizes = new double[nHaploidChromosomeNum];           // ~~ declaration was in def.h
	double * pCentromerePos = new double[nHaploidChromosomeNum];	
	
	// ~~ Initialize the size of the chromosomes -> all the same size
	while(nCurrChromosome < nHaploidChromosomeNum)
	{
		pChromosomeSizes[nCurrChromosome] = nChromosomeSize;	
		/*
		Now do not allow users to move centromere any more since centromere can be simply a marker, and it's forced 
		to the end of the chromosome.
        */
		pCentromerePos[nCurrChromosome] = nChromosomeSize;    // ~~ as implemented, not changed 
		nCurrChromosome++;
    }	
		
	// now prepare file:
	//Figure out how to distribute the markers
	
	ofstream fOutFile;
	fOutFile.open(sOutFile.c_str());
	fOutFile << "HaploidChromosomeNum = " << nHaploidChromosomeNum << endl; //Write haploid n.
	fOutFile << "RandomSeed = " << nRandSeed << endl; //Write random seed n.
	fOutFile << "Pop1AvgAlleleFreq = " << nPop1AvgAlleleFreq << endl; //
	fOutFile << "Pop1AvgAlleleFreqStdev = " << nPop1AvgAlleleFreqStdDev << endl; //
	fOutFile << "Pop2AvgAlleleFreq = " << nPop2AvgAlleleFreq << endl; //
	fOutFile << "Pop2AvgAlleleFreqStdev = " << nPop2AvgAlleleFreqStdDev << endl; //
		
		
	
	double nSumLens = fnSum(pChromosomeSizes , nHaploidChromosomeNum); // get sum of chromosome lengths (~~ defined in maths.cpp)
	double nLargestChromSize = 0;
	
	for (int i=0;  i < nHaploidChromosomeNum; i++) { // generate markers chromosome by chromosome
		
		double nChrLen = pChromosomeSizes[i];
		int	   nNumMarkers = (int) ceil( ((double) nTotalMarkers) * (nChrLen / nSumLens));
		
		nLargestChromSize = (nLargestChromSize < nChrLen)? nChrLen : nLargestChromSize;
		
		fOutFile << ":chr " << (i+1) << "\tlen = " << nChrLen << "\tcentromere = " << pCentromerePos[i] << endl; // chromosome number
		fOutFile << "\nNum\tPosition Percentage\tPosition Abs.\tFreq. Pop 1\tFreq. Pop 2\n" ;
		
		double * pPos = new double[nNumMarkers];
		for (int j=0;j<nNumMarkers;j++) {
			pPos[j] = UniformGen.Next(); // create random set of uniformly distributed markers
		}
		
		//Sort markers according to position on chromosome:
		
		qsort (pPos , nNumMarkers , sizeof(double), fnCompare);
		
		for (int j=0;j<nNumMarkers;j++) { 
			fOutFile << (j + 1);
			fOutFile << '\t';
			fOutFile << pPos[j] ;
			fOutFile << '\t';
			fOutFile << pPos[j] * nChrLen ;
			fOutFile << '\t';
			fOutFile << NormalExt(nPop1AvgAlleleFreq ,nPop1AvgAlleleFreqStdDev , 0.5, 1.0);         // ~~ defined in maths.cpp
			fOutFile << '\t';
			fOutFile <<  NormalExt( 1.0 - nPop2AvgAlleleFreq ,nPop2AvgAlleleFreqStdDev , 0, 0.5);
			fOutFile << endl; 
		}
		
		delete[] pPos;
	}
	
	fOutFile.close();
	
    // ~~ Generate the recombination probabilities
	generateRecombinationFreqMap(nHaploidChromosomeNum , nLargestChromSize, pChromosomeSizes, pCentromerePos, model); 
}


/*
 * Fonction to generate recombination map
 */
void generateRecombinationFreqMap(int nChromosomes, double nLargestChromSize, double * pChromosomeSizes, double * pCentromerePos, string model) 
{ 
	double nSumLens = fnSum(pChromosomeSizes , nChromosomes); // get sum of chromosome lengths
    // ~~ Use uniform probabilities rate
	char nRate(1); //1 - uniform, 2 - parabola  // ~~ for the moment we never change this !!!
    // ~~ Set the total expected recombination events in the whole genome per meiosis "4n should be fine" -> 16
    double nTotalExpectedRecombPerMeiosis(nRecombMeiosis);
    // ~~ Number of rec points
	double nTotalSamplePoint (factSamplePoint * nChromosomes);  //How often to put a point. // ~~ -> was 20 in the original
	//~~ Fold difference between female:male => set to 1, no diff
    double nFemFold (1); // ~~ for the moment we never change this !!!
	// ~~ Set output file name
	string sOutFile("simu/"+model+"/"+model+"_"+recombFile);

	// now prepare file:	~~ => THIS HAS MOSTLY NOT BEEN CHANGED
	ofstream fOutFile;
	fOutFile.open(sOutFile.c_str());
	
	for (int i=0;i< nChromosomes;i++) {		
		int nChrIndex =i;
		double nChrLen = pChromosomeSizes[nChrIndex];
		double nCentromerePos =  pCentromerePos[nChrIndex];	
		
		int	   nNumSamplePoints = (int) ceil( ((double) nTotalSamplePoint) * (nChrLen / nSumLens));

		double nArm1Len =  nCentromerePos;
		double nArm2Len =  nChrLen - nCentromerePos;
			
		fOutFile << ":chr " << (nChrIndex+1) << endl;
		
		double nExpectedMaleRecPerMeiosisArm1 = nTotalExpectedRecombPerMeiosis * ( nArm1Len/nSumLens) / 2;
		double nExpectedMaleRecPerMeiosisArm2 = nTotalExpectedRecombPerMeiosis * ( nArm2Len/nSumLens) / 2;
		double nExpectedFemaleRecPerMeiosisArm1 = nTotalExpectedRecombPerMeiosis * nFemFold * ( nArm1Len/nSumLens) / 2;
		double nExpectedFemaleRecPerMeiosisArm2 = nTotalExpectedRecombPerMeiosis * nFemFold * ( nArm2Len/nSumLens) / 2;
		double nExpectedRecPerMeiosisArm1 = ( nExpectedMaleRecPerMeiosisArm1 + nExpectedFemaleRecPerMeiosisArm1 ) / 2;
		double nExpectedRecPerMeiosisArm2 = ( nExpectedMaleRecPerMeiosisArm2 + nExpectedFemaleRecPerMeiosisArm2 ) / 2;
		
		fOutFile << ":ExpectedMaleRecPerMeiosisArm1 = " << nExpectedMaleRecPerMeiosisArm1 << endl;
		fOutFile << ":ExpectedMaleRecPerMeiosisArm2 = " << nExpectedMaleRecPerMeiosisArm2 << endl;
		fOutFile << ":ExpectedFemaleRecPerMeiosisArm1 = " << nExpectedFemaleRecPerMeiosisArm1 << endl;
		fOutFile << ":ExpectedFemaleRecPerMeiosisArm2 = " << nExpectedFemaleRecPerMeiosisArm2 << endl;
		fOutFile << "\tmale"<<"\tmale_accumulative" << "\tfemale\tfemale_accumulative\tMale_recombination_fraction\tFemale_recombination_fraction\tAvg_recombination_fraction\tKosambi_male_interval\tKosambi_male_pos\tKosambi_female_interval\tKosambi_female_pos\tKosambi_both_interval\tKosambi_both_pos" << endl;
	
		std::map<double, double *>  mSamplePoints;

		for (int j=0; j < nNumSamplePoints; j++) {
			double nPos = nChrLen * UniformGen.Next(); // random position
			double nDistanceFromCentromere = ( nPos - nCentromerePos ) / nLargestChromSize; //  distance from centromere standarized by the longest chromosome
			nDistanceFromCentromere = (nDistanceFromCentromere >= 0 )? nDistanceFromCentromere: -nDistanceFromCentromere;
			double nArmLen = (nPos > nCentromerePos)? nChrLen - nCentromerePos : nCentromerePos;
			double nDistanceFromCentromereOnArm = (nPos > nCentromerePos)? (nPos - nCentromerePos)/nArmLen : (nCentromerePos - nPos)/nArmLen ;
			double nProbAtPosMale(0.1); //~~ 
			//~~ double nProbAtPosMale =  (nRate=='1')? 0.1 : pow(nParabolaA * nDistanceFromCentromere, 2); 
			//double nProbAtPosFemale = nProbAtPosMale * nFemCorrA  + nFemCorrB;
			double nProbAtPosFemale(0.1*nFemFold); //~~ 
			//~~ double nProbAtPosFemale = (nRate=='1')? 0.1 * nFemFold : pow(nParabolaAFemale * nDistanceFromCentromere, 2) * nFemFold ; 
			double * pProbs = new double[4];
			pProbs[0] = NormalExt(nProbAtPosMale , nProbAtPosMale / 2, 0, 1);
			//pProbs[1] = nMaxAvgRecombRate * pProbs[0] ;
			pProbs[2] = NormalExt(nProbAtPosFemale , nProbAtPosFemale / 2, 0, 1);
			//pProbs[3] = nMaxAvgRecombRate * pProbs[2] ;
			if (mSamplePoints.find(nPos) == mSamplePoints.end()) {
				mSamplePoints[nPos] = pProbs ;
			}
			else {
				j--;
				continue; // won't work because same position already exists.
			}
		}
		double nMaleArm1Accum = 0;
		double nFemaleArm1Accum = 0;
		double nMaleArm2Accum = 0;
		double nFemaleArm2Accum = 0;
		
		for (std::map<double,double *>::iterator it = mSamplePoints.begin(); it != mSamplePoints.end(); ++it) {
			double nMaleAccum, nFemaleAccum;
			if (it->first < nCentromerePos) { //on arm 1
				nMaleArm1Accum += (it->second)[0];
				nFemaleArm1Accum += (it->second)[2];
			}
			else {
				 nMaleArm2Accum += (it->second)[0];
				 nFemaleArm2Accum += (it->second)[2];
			}
		}
		double nScaledMaleArm1Accum = 0;
		double nScaledFemaleArm1Accum = 0;
		double nScaledMaleArm2Accum = 0;
		double nScaledFemaleArm2Accum = 0;

		double nScaledMaleProb = 0;
		double nScaledFemaleProb = 0;
		
		double nPrevScaledMaleArm1Accum =0;
		double nPrevScaledFemaleArm1Accum =0;
		double nPrevScaledMaleArm2Accum =0;
		double nPrevScaledFemaleArm2Accum =0;
		
		double nKosambiMapDistance_male = 0;
		double nKosambiMapDistance_female = 0;
		double nKosambiMapDistance_both = 0;			
		
		double nKosambiMapPos_male = 0;
		double nKosambiMapPos_female = 0;
		double nKosambiMapPos_both = 0;
		
		double nRecombFraction_male = 0;
		double nRecombFraction_female = 0;
		double nRecombFraction_both = 0;

		for (std::map<double,double *>::iterator it = mSamplePoints.begin(); it != mSamplePoints.end(); ++it) {
			double nMaleAccum, nFemaleAccum;
			if (it->first < nCentromerePos) { //on arm 1
				nScaledMaleProb = (it->second)[0] / nMaleArm1Accum;
				nScaledFemaleProb = (it->second)[2] / nFemaleArm1Accum;
				nMaleAccum = nScaledMaleArm1Accum += nScaledMaleProb;
				nFemaleAccum = nScaledFemaleArm1Accum += nScaledFemaleProb;
				
				 nRecombFraction_male = nExpectedMaleRecPerMeiosisArm1 * nScaledMaleProb  ;
				 nKosambiMapPos_male += nKosambiMapDistance_male = 0.25 * log( (1+2 * nRecombFraction_male)/(1-2 * nRecombFraction_male) ) ;
				 nRecombFraction_female = nExpectedFemaleRecPerMeiosisArm1 * nScaledFemaleProb ;
				 nKosambiMapPos_female += nKosambiMapDistance_female = 0.25 * log( (1+2 * nRecombFraction_female)/(1-2 * nRecombFraction_female) ) ;
			}
			else {
				nScaledMaleProb = (it->second)[0] / nMaleArm2Accum;
				nScaledFemaleProb = (it->second)[2] / nFemaleArm2Accum;
				nMaleAccum = nScaledMaleArm2Accum += nScaledMaleProb;
				nFemaleAccum = nScaledFemaleArm2Accum += nScaledFemaleProb;
				
				nRecombFraction_male = nExpectedMaleRecPerMeiosisArm2 * nScaledMaleProb  ;
				 nKosambiMapPos_male += nKosambiMapDistance_male = 0.25 * log( (1+2 * nRecombFraction_male)/(1-2 * nRecombFraction_male) ) ;
				 nRecombFraction_female = nExpectedFemaleRecPerMeiosisArm2 * nScaledFemaleProb  ;
				 nKosambiMapPos_female += nKosambiMapDistance_female = 0.25 * log( (1+2 * nRecombFraction_female)/(1-2 * nRecombFraction_female) ) ;
			}
			nRecombFraction_both = (nRecombFraction_male + nRecombFraction_female)/2;
			nKosambiMapPos_both += nKosambiMapDistance_both = 0.25 * log( (1+2 * nRecombFraction_both)/(1-2 * nRecombFraction_both) ) ;
			
			fOutFile << setprecision (10) << it->first << '\t' << nScaledMaleProb << '\t' << nMaleAccum << '\t' << nScaledFemaleProb << '\t' << nFemaleAccum;
			fOutFile << setprecision (10) << '\t' << nRecombFraction_male << '\t' << nRecombFraction_female << '\t' << nRecombFraction_both << '\t' << nKosambiMapDistance_male << '\t' << nKosambiMapPos_male << '\t' << nKosambiMapDistance_female << '\t' << nKosambiMapPos_female << '\t' << nKosambiMapDistance_both << '\t' << nKosambiMapPos_both << endl;
		
			double nPrevScaledMaleArm1Accum = nScaledMaleArm1Accum;
			double nPrevScaledFemaleArm1Accum = nScaledFemaleArm1Accum;
			double nPrevScaledMaleArm2Accum = nScaledMaleArm2Accum;
			double nPrevScaledFemaleArm2Accum = nScaledFemaleArm2Accum;
		}
	}
	fOutFile.close();
}

//Note that sex NEEDS to be specified in the phenotype file. The sex ratios in the configuration only acts as filters in generation 0. 
//If your genetic and phenotypic settings for sex (say, only females are produced) cannot meet these sex ratio filters (say 0.5), then the program will not stop but try to create more individuals. 
//For a typical XY organism, specify an imaginary male determining locus in the gene file:
//name	Chromosome	position	Dominant Allele	Dominant Allele Value	Recessive Allele	Recessive Allele Value	Mode	Dominant Freq Pop1	Dominant Freq Pop2
//chr10_3000000	10	3000000	A	1	a	0	Hemizygous	0.5	0.5

//And then specify a phenotype called "Sex" (case-sensitive) in the phenotype file:
//Phenotype	Formula
//Sex		chr10_3000000

/*
 * ~~ Function to generate the gene file
 */
void generateGeneFile(string model)
{
    string geneOutFile("simu/"+model+"/"+model+"_"+geneFile);
	ofstream fOutFile;
	fOutFile.open(geneOutFile.c_str());
	
	fOutFile << "name\tChromosome\tposition\tDominant Allele\tDominant Allele Value\tRecessive Allele\tRecessive Allele Value\tMode\tDominant Freq Pop1\tDominant Freq Pop2" << endl; // header line
	
	fOutFile << "chr1_10\t1\t10\tA\t1\ta\t0\tHemizygous\t0.5\t0.5" << endl; 
	// needed to start simulations !!! chr1_10 as foo (if change, change also in phenoFile)
	
	// ~~ within the script -> print the basis as needed for running simulation
	// ~~ then if needed in another script change this to fine-tune selection
	
    return;
}

/*
 * ~~ Function to generate the natural selection file	
 */
// Example:
    //Population	Gen	Selection
    //nigra	-1	1
    //rapa	-1	1
    //hyb	-1	1
    // -1 in Gen => apply to all generation
    // Selection => probability between 0 and 1
void generateNaturalSelFile(string model)
{
    string naturalOutFile("simu/"+model+"/"+model+"_"+naturalSelFile);
	ofstream fOutFile;
	fOutFile.open(naturalOutFile.c_str());
	fOutFile << "Population\tGen\tSelection" << endl;  //~~ print the header
	fOutFile << pop2Name << "\t-1\t1" << endl;                                
	fOutFile << pop1Name <<"\t-1\t1" << endl; 
	fOutFile << hybName << "\t-1\t1" << endl; 
    return;
}


/*
 * ~~ Function to generate the sexual selection file	
 */
// Example:
    //Population	Gen	Selection
    //nigra	-1	1
    //rapa	-1	1
    //hyb	-1	1
    // -1 in Gen => apply to all generation
    // Selection => probability between 0 and 1
void generateSexualSelFile(string model)
{
    string sexOutFile("simu/"+model+"/"+model+"_"+sexualSelFile);
	ofstream fOutFile;
	fOutFile.open(sexOutFile.c_str());
	
	fOutFile << "Population\tGen\tSelection" << endl;  //~~ print the header
	fOutFile << pop2Name << "\t-1\t1" << endl;                                
	fOutFile << pop1Name <<"\t-1\t1" << endl; 
	fOutFile << hybName << "\t-1\t1" << endl; 
	
    return;
}
	
/*	
 * ~~ Function to generate the phenotype file
 */
//And then specify a phenotype called "Sex" (case-sensitive) in the phenotype file:
//Phenotype	Formula
//Sex		chr10_3000000

void generatePhenoFile(string model)
{
    string phenoOutFile("simu/"+model+"/"+model+"_"+phenotypeFile);
	ofstream fOutFile;
	fOutFile.open(phenoOutFile.c_str());
	
	fOutFile << "Phenotypes\tFormula" << endl;   // ~~ print the header
	fOutFile << "Sex\tchr1_10" << endl;        // needed to start simulations !!! chr1_10 as foo (if change, change also in geneFile)
	
    return;
}   

