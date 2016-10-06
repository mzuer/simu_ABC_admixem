#include "my_math.h"
#include "settings.h"
#include <math.h>       /* exp */
#include <time.h>
#include <iostream>
#include <gsl/gsl_randist.h> /* draw from distribution */

 
// draw from distribution: use gsl-lib
//https://www.math.utah.edu/software/gsl/gsl-ref_281.html#SEC281
//http://ampl.github.io/amplgsl/ran-beta.html
// how to use gsl fct: https://github.com/jstac/hh_sampling/blob/master/hh_cftp.c
// how to compile: https://ubuntuforums.org/showthread.php?t=270924

using namespace std;


    
// NB: might be changed, but for the moment use extern seed initialized in the main

double getFromUnif(const double min, const double max)
{
   /*
    int static i = 0;
    srand(time(NULL));          
    int seed = rand() % 1000000+i; // seed
    i++;
    const gsl_rng_type * T;
    gsl_rng * rd;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    rd = gsl_rng_alloc(T);
    gsl_rng_set(rd, seed);*/
gsl_rng * rng;
    return(gsl_ran_flat(rng, min, max)); // returns a random variate from the flat (uniform) distribution from a to b. //was rd not rng
}

double getFromBeta(const double a, const double b)
{
/*
    int static i = 0;
    srand(time(NULL));          
    int seed = rand() % 1000000 + i; // seed
    i++;
    const gsl_rng_type * T;
    gsl_rng * rd;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    rd = gsl_rng_alloc(T);
    gsl_rng_set(rd, seed);*/
    gsl_rng * rng;
    return(gsl_ran_beta(rng, a, b)); // returns a random variate from the beta distribution. //was rd not rng
} 


double getFromGaussian(const double mean, const double sigma)
{
    /*const gsl_rng_type * T;
    gsl_rng * rd;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    rd = gsl_rng_alloc(T);
    gsl_rng_set(rd, seed);*/
    //was rd not rng
    // gsl_ran_gaussian: returns a Gaussian random variate, with mean zero and standard deviation sigma
    gsl_rng * rng;
    return (gsl_ran_gaussian(rng, sigma) + mean); // Use the transformation z = \mu + x to obtain a Gaussian distribution with mean \mu. 
    
}







