#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "consoleIN.h"
#include "createMat.h"
#include "quantisation.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector_complex_double.h>
//curretn limitations: no gui, cannot place reactances in parallel
//(fix with small series resistance), no curretn sources
//(fix with include V_i as independant variabel and create new large matrix)
//initial conditions=0 (fix with adding constants into K vector, but then
//need to take input, else larger fix is to create new matrix with even
//more uknowns and if need 0 then set zero else solve by setting respective
//vairables of components to zero) 

int maxargs = 10;
bool running = true;
extern int graphNodecount;
extern int edgescount;

int main()
{
	char input[100];
	char din[maxargs][10];
	test(graphHead);
	gsl_matrix * A,*Q,*B,*ib;
	gsl_vector *Zsymbolic;
	gsl_vector_complex *Z,*K;
	gsl_permutation * p;
	while (running) 
	{
		printf(">");
		if (fgets(input, sizeof(input), stdin) != NULL) handleinput(input,din,ib);
		if(strcmp(input, "solve") == 0) 
		{
			A = matfromgraph(graphHead,graphNodecount,edgescount);
			makeQBP(A,graphNodecount,edgescount,&Q,&B,&p);
			makeZK(graphHead,graphNodecount,edgescount,&Z,&Zsymbolic,&K,p);
			solveI(graphNodecount,edgescount,&ib,0.3,1000,Q,B,Zsymbolic,K,Z);
		}
	}
	freeall();
	gsl_matrix_free(A);
	gsl_matrix_free(Q);
	gsl_matrix_free(B);
	gsl_vector_complex_free(Z);
	gsl_vector_free(Zsymbolic);
	gsl_vector_complex_free(K);
	gsl_permutation_free(p);
	return 0;
}