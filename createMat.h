#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#ifndef _MATH_
#define _MATH_
#include "dsGraph.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector_complex_double.h>

gsl_matrix * matfromgraph(graphNode,int,int);
gsl_matrix * removegnd(gsl_matrix*,int,int);
void printmat(gsl_matrix *, int, int);
void printvect(gsl_vector *,int);
void printcomplvect(gsl_vector_complex *,int);
void makeQBP(gsl_matrix *,int, int,gsl_matrix **,gsl_matrix **, gsl_permutation **);
void makeZK(graphNode,int, int,gsl_vector_complex **Z,gsl_vector **Zsymbolic,gsl_vector_complex **K, gsl_permutation *);
/*traverses graph and creates impedance matrix as a diagonal vector,
Zsymbolic tracks the variable s, ie if res,cap or ind,
permutation reorders the linkedNodeIndex to match pivoting of A,Q,B
K is the constant vector to accomodate voltage source
*/
extern int gnd;
#endif