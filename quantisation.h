#ifndef _QUANTH_
#define _QUANTH_

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_complex_math.h>
#include <math.h>
#include "createMat.h"
#include "consoleIN.h"

double cot(double);
void creatediscsys(int,int,double,int,gsl_matrix_complex**,gsl_vector_complex **
	,gsl_complex*,gsl_complex*,gsl_vector_complex*,gsl_vector*,gsl_vector_complex*);
void realtocomplexmat(int,int,gsl_matrix*,gsl_matrix_complex**,gsl_matrix*,
	gsl_matrix_complex**);
gsl_vector* getrealvector(int,gsl_vector_complex*);
void linearsystem(int,int,gsl_matrix_complex*,gsl_vector_complex*,gsl_matrix_complex*,
	gsl_matrix_complex*,gsl_vector_complex**);
void solveI(int,int,gsl_matrix**,double,int,gsl_matrix*,gsl_matrix*,gsl_vector*,
	gsl_vector_complex*,gsl_vector_complex*);
bool checkifzerocomp(gsl_complex);
//these are not BLAS calls, so may be slow

#endif