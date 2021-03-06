#ifndef MY_MATH_H
#define MY_MATH_H

#include <stdio.h>
#include <stdlib.h>

// draw from distribution: use gsl-lib
// https://www.math.utah.edu/software/gsl/gsl-ref_281.html#SEC281
// http://ampl.github.io/amplgsl/ran-beta.html
// how to use gsl fct: https://github.com/jstac/hh_sampling/blob/master/hh_cftp.c
// how to compile: https://ubuntuforums.org/showthread.php?t=270924

//double getFromUnif(const double, const double);
double getFromNormal(const double, const double);
//double getFromBeta(const double, const double);
//double getFromGaussian(const double, const double );

#endif
