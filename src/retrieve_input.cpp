#include "settings.h"
#include <tuple>
#include <string>
#include <vector>
#include <tuple>
#include <utility>      /* pair */
#include <algorithm>    /* sort*/
#include <iostream>
#include <fstream>
#include <sstream>      /* strinstream */
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <cstdlib>      /* atof */

using namespace std;

void retrieveABCinput()
{

    //=> START RETRIEVING INFO FROM runABC.input
    cout << "... Retrieve information from input file ..." << endl;
    string inputF("runABC.input");
    ifstream input;
    int nField1(0);
    int nField2(0);
    string line, field1, param, field2;
    vector<double> distance_list_temp;
    vector<int> pop_sampling_list_temp;
    vector<int> nbr_pop_list_temp;
	input.open(inputF.c_str());     
	
	if(!input.is_open())
	{
	    cout << "Unable to open: " << inputF << endl;
	    exit(EXIT_FAILURE);
	}		
    while (getline(input,line))
    {            
        if(line.substr(0,1) == "#")
        {
            continue;
        }
        stringstream ss(line);        
        nField1 = 0;
                      	           
        while(getline(ss, field1, '='))         // split the line around the "="   distances=10,20,30
        {    
            if(nField1 == 0)
            {
                param = field1;
                nField1++;
                continue;
            }
            if(nField1 == 1)
            {
                stringstream ss2(field1);
                nField2 = 0;
                while(getline(ss2, field2, ','))         // split the line around the ","
                {                
                    if(param == "distances") distance_list_temp.push_back(stod(field2));
                    else if(param == "sampling_pop") pop_sampling_list_temp.push_back(stoi(field2));
                    else if(param == "nbr_pop") nbr_pop_list_temp.push_back(stoi(field2));
                }  
                break;     
            }
        }                
    }          
    input.close();
    
    cout << "Read from runABC.input :" << endl;
    cout << "Distances :" ;
    for(auto i: distance_list_temp) 
    {
        cout << i << "\t";
    }
    cout << "\nSampling size :" ;
    for(auto i: pop_sampling_list_temp) 
    {
        cout << i << "\t";
    }    
    cout << "\nNumber of pop :" ;
    for(auto i: nbr_pop_list_temp) 
    {
        cout << i << "\t";
    }   
    cout << endl;
    
    // set number of pop, number of individuals to sample by pop and number of populations at each distance (none is empty)
    if(!distance_list_temp.empty() && !nbr_pop_list_temp.empty() && !pop_sampling_list_temp.empty())
    {

        if(! (distance_list_temp.size() == nbr_pop_list_temp.size() && distance_list_temp.size() == pop_sampling_list_temp.size()) )
        {
            cout << "ERROR in runABC.input: not same size of distances, sampling and pop" << endl;
            exit(EXIT_FAILURE);
        }

        vector<tuple <double, int, int>> unsorted_all_vect;
        
        for(int i = 0 ; i < distance_list_temp.size(); i++)
        {
            unsorted_all_vect.push_back(make_tuple(distance_list_temp[i], pop_sampling_list_temp[i], nbr_pop_list_temp[i]));
        }
        sort(unsorted_all_vect.begin(), unsorted_all_vect.end(), [](tuple<double, int, int> const &p1, tuple<double, int, int> const &p2) {
            return get<0>(p1) < get<0>(p2);
            });
            
        //clear from the default values and populate with values retrieved from input file
        distance_list.clear(); // = {};
        pop_sampling_list.clear(); // = {};
        nbr_pop_list.clear(); 
        
        cout << "... Found following distances, subsampling size and number of pops:" << endl;
        for(auto i: unsorted_all_vect)
        {
            cout << get<0>(i) << "\t" << get<1>(i) << "\t" << get<2>(i)  << endl;
            distance_list.push_back(get<0>(i));
            pop_sampling_list.push_back(get<1>(i));
            nbr_pop_list.push_back(get<2>(i));
        } 
    // set the number of pop (in the same order as the distance) and the distance) => number of individuals to sample is empty
    } else if(!distance_list_temp.empty() && !nbr_pop_list_temp.empty() && pop_sampling_list_temp.empty())
    {
        if(!(distance_list_temp.size() == nbr_pop_list_temp.size()))
        {
            cout << "ERROR in runABC.input: not same size of distances and subsampling" << endl;
            exit(EXIT_FAILURE);
        }

        vector<pair<double, int>> unsorted_vectB;                    // sort by distance without losing corresponding sample size
        for(int i = 0 ; i < distance_list_temp.size(); i++)
        {
            unsorted_vectB.push_back(make_pair(distance_list_temp[i], nbr_pop_list_temp[i]));
        }
        sort(unsorted_vectB.begin(), unsorted_vectB.end(), [](pair<double, int> const &p1, pair<double, int> const &p2) {
            return p1.first < p2.first;
            });
        //clear from the default values and populate with values retrieved from input file
        distance_list.clear();  
        nbr_pop_list.clear();
        cout << "... Found following distances and number of pops:" << endl;
        for(auto i: unsorted_vectB)
        {
            cout << to_string(i.first) << "\t" << to_string(i.second) << endl;
            distance_list.push_back(i.first);
            nbr_pop_list.push_back(i.second);
        } 
    } 
    // set the sampling size and distance lists (number of at each distance is empty)
    else if(!distance_list_temp.empty() && !pop_sampling_list_temp.empty() && nbr_pop_list_temp.empty())
    {
         if(!(distance_list_temp.size() == pop_sampling_list_temp.size()))
        {
            cout << "ERROR in runABC.input: not same size of distances and subsampling size" << endl;
            exit(EXIT_FAILURE);
        }

        vector<pair<double, int>> unsorted_vect;                    // sort by distance without losing corresponding sample size
        for(int i = 0 ; i < distance_list_temp.size(); i++)
        {
            unsorted_vect.push_back(make_pair(distance_list_temp[i], pop_sampling_list_temp[i]));
        }
        sort(unsorted_vect.begin(), unsorted_vect.end(), [](pair<double, int> const &p1, pair<double, int> const &p2) {
            return p1.first < p2.first;
            });
        distance_list.clear();// = {};
        pop_sampling_list.clear(); // = {};
        cout << "... Found following distances and number of individuals to subsample:" << endl;
        for(auto i: unsorted_vect)
        {
            cout << to_string(i.first) << "\t" << to_string(i.second) << endl;
            distance_list.push_back(i.first);
            pop_sampling_list.push_back(i.second);
        }
    } 
    // only distance list is given -> need to be sorted
    else if(!distance_list_temp.empty() && pop_sampling_list_temp.empty() && nbr_pop_list_temp.empty())
    {
        cout << "... Found following distances:" << endl;
        sort(distance_list_temp.begin(), distance_list_temp.end());
        for(auto i: distance_list_temp)
        {
            cout << i << "\t" ;
        }
        cout << endl;
        distance_list = distance_list_temp;
    } 
    // only size of sampling is given 
    else if(!pop_sampling_list_temp.empty())
    {
        cout << "... Found following number of populations for each distance:" << endl;

        for(auto i: pop_sampling_list_temp)
        {
            cout << i << "\t" ;
        }
        cout << endl;
        pop_sampling_list = pop_sampling_list_temp;        
    } 
    // only 
    else if(!nbr_pop_list_temp.empty())
    {
        cout << "... Found following number of individuals to subsample:" << endl;

        for(auto i: nbr_pop_list_temp)
        {
            cout << i << "\t" ;
        }
        cout << endl;
        nbr_pop_list = nbr_pop_list_temp;    
    }        
    // <= END RETRIEVING INFO FROM runABC.input
    
    
    return;
}
