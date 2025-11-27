#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_complex_math.h>
#include "createMat.h"
#include "consoleIN.h"

gsl_matrix * matfromgraph(graphNode graphHead,int n, int b)
{
	gsl_matrix *A = gsl_matrix_alloc(n, b);
	gsl_matrix_set_zero(A);
	int edgeindex = 0;
	inigraph(graphHead);
	for(int i=0;i<n;i++)
	{
		graphNode tmp = searchgraphNode(i+1);
		linkedNode templink = tmp->nextlink;
		while(templink!=NULL)
		{
			if(templink->branchIndex == 0 && templink->isforwardlink==true) 
			{
				templink->branchIndex = ++edgeindex;
				gsl_matrix_set(A, i, edgeindex-1, 1);
				gsl_matrix_set(A, templink->destindex-1, edgeindex-1, -1);
			}
			templink = templink->nextlink;
		}
	}
	printmat(A,n,b);
	return A;
}

void printmat(gsl_matrix * A, int n, int b)
{
	for(int i =0;i<n;i++)
	{
		for(int j=0;j<b;j++)
		{
			printf("%g ", gsl_matrix_get(A, i, j));
		}
		printf("\n");
	}
}

void printvect(gsl_vector *v,int n)
{
	for(int i=0;i<n;i++)
	{
		printf("%g ", gsl_vector_get(v,i));
	}
	printf("\n");
}

void printcomplvect(gsl_vector_complex *v,int n)
{
    for (int i = 0; i < n; i++) 
    {
        gsl_complex z = gsl_vector_complex_get(v, i);
        printf("%.10g + i%.10g, ",GSL_REAL(z),GSL_IMAG(z));
    }
}

gsl_matrix * removegnd(gsl_matrix* A,int n,int b)
{
	gsl_matrix * Ar = gsl_matrix_alloc(n-1, b);
	int r = 0;
	for (int i = 0; i < n; i++) 
	{
		if (i == gnd-1) continue;
		for (int j = 0; j < b; j++) 
		{
			gsl_matrix_set(Ar, r, j, gsl_matrix_get(A, i, j));
		}
		r++;
	}
	printmat(Ar,n-1,b);
	return Ar;
}

void makeQBP(gsl_matrix *A, int n, int b, gsl_matrix **Qptr, gsl_matrix **Bptr, gsl_permutation **pptr) 
{
	gsl_matrix *Q = gsl_matrix_alloc(n-1, b);
	gsl_matrix *B = gsl_matrix_alloc(b-n+1, b);
	gsl_matrix *P = gsl_matrix_alloc(b, b);
	gsl_matrix *check1 = gsl_matrix_alloc(n-1, b-n+1);
	gsl_matrix *check2 = gsl_matrix_alloc(b-n+1, n-1);
	gsl_matrix *Bt = gsl_matrix_alloc(b, b-n+1);
	gsl_matrix *Qt = gsl_matrix_alloc(b, n-1);
	gsl_matrix *Aperm = gsl_matrix_alloc(n, b);
	gsl_matrix *Atrans = gsl_matrix_alloc(b, n);
	gsl_matrix *At = gsl_matrix_alloc(n-1, n-1);
	gsl_matrix *Al = gsl_matrix_alloc(n-1, b-n+1);
	gsl_matrix *Atinv = gsl_matrix_alloc(n-1, n-1);
	gsl_matrix *commM = gsl_matrix_alloc(n-1, b-n+1);
	gsl_matrix *commT = gsl_matrix_alloc(b-n+1, n-1);
	gsl_permutation *p = gsl_permutation_alloc(b);
	gsl_permutation *pt = gsl_permutation_alloc(n-1);
	int signum;

	gsl_matrix_set_zero(P);
	gsl_matrix_transpose_memcpy(Atrans, A);
	gsl_linalg_LU_decomp(Atrans, p, &signum);

	for (int i = 0; i < b; i++)	
	{
		gsl_matrix_set(P, i, gsl_permutation_get(p, i), 1.0);
	}
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, A, P, 0.0, Aperm);

	for (int i = 0; i < n-1; i++) 
	{
		for (int j = 0; j < n-1; j++)
		{
			gsl_matrix_set(At, i, j, gsl_matrix_get(Aperm, i, j));
		}
		for (int j = n-1; j < b; j++)
		{
			gsl_matrix_set(Al, i, j-n+1, gsl_matrix_get(Aperm, i, j));
		}
	}

	gsl_linalg_LU_decomp(At, pt, &signum);
	gsl_linalg_LU_invert(At, pt, Atinv);
	
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Atinv, Al, 0.0, commM);

	for (int i = 0; i < n-1; i++) 
	{
		for (int j = 0; j < b-n+1; j++)
		{
			gsl_matrix_set(Q, i, j, gsl_matrix_get(commM, i, j));
		}
		for (int j = b-n+1; j < b; j++)
		{
			gsl_matrix_set(Q, i, j, (i==j-(b-n+1)) ? 1.0 : 0.0);
		}
	}

	gsl_matrix_set_zero(B);
	for (int i = 0; i < b-n+1; i++)
	{
		gsl_matrix_set(B, i, i, 1.0);
	}

	gsl_matrix_transpose_memcpy(commT, commM);
	for (int i = 0; i < b-n+1; i++)
	{
		for (int j = b-n+1; j < b; j++)
		{
			gsl_matrix_set(B, i, j, -gsl_matrix_get(commT, i, j-b+n-1));
		}
	}

	*Qptr = Q;
	*Bptr = B;
	*pptr = p;
	gsl_matrix_transpose_memcpy(Bt, B);
	gsl_matrix_transpose_memcpy(Qt, Q);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Q, Bt, 0.0, check1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, B, Qt, 0.0, check2);

	/*printf("\nPermuted A =\n");
	printmat(Aperm,n,b);

	printf("\nQ =\n");
	printmat(Q,n-1,b);

	printf("\nB =\n");
	printmat(B,b-n+1,b);

	printf("\nQBt =\n");
	printmat(check1,n-1,b-n+1);

	printf("\nBQt =\n");
	printmat(check2,b-n+1,n-1);*/

	gsl_matrix_free(P);
	gsl_matrix_free(Aperm);
	gsl_matrix_free(At);
	gsl_matrix_free(Al);
	gsl_matrix_free(Atinv);
	gsl_matrix_free(commM);
	gsl_matrix_free(commT);
	gsl_matrix_free(Qt);
	gsl_matrix_free(Bt);
	gsl_matrix_free(check1);
	gsl_matrix_free(check2);
	gsl_permutation_free(pt);
}

void makeZK(graphNode graphHead,int n, int b,gsl_vector_complex **Zptr,gsl_vector **Zsymbolicptr,gsl_vector_complex **Kptr, gsl_permutation *p)
{
	gsl_vector_complex *Z = gsl_vector_complex_alloc(b);
	gsl_vector *Zsymbolic = gsl_vector_alloc(b);
	gsl_vector_complex *K = gsl_vector_complex_alloc(b);
	gsl_vector_complex_set_zero(Z);
	gsl_vector_set_zero(Zsymbolic);
	gsl_vector_complex_set_zero(K);

	graphNode tmp = graphHead;
	while(tmp!=NULL)
	{
		linkedNode tmplink = tmp->nextlink;
		while(tmplink!=NULL)
		{
			if(tmplink->isforwardlink)
			{
				int branchIndex = tmplink->branchIndex - 1;
				tmplink->branchIndex = 1+gsl_permutation_get(p,branchIndex);
				gsl_complex z = gsl_complex_rect(tmplink->value,0);
				gsl_complex zc = gsl_complex_rect(1.0/tmplink->value,0);
				switch(tmplink->comp)
				{
					case Resistor :
						gsl_vector_complex_set(Z,branchIndex,z);
						break;
					case VoltSource :
						gsl_vector_complex_set(K,branchIndex,z);
						gsl_vector_set(Zsymbolic,branchIndex,-1.0);
						break;
					case Capacitor :
						gsl_vector_complex_set(Z,branchIndex,zc);
						gsl_vector_set(Zsymbolic,branchIndex,-1.0);
						break;
					case Inductor :
						gsl_vector_complex_set(Z,branchIndex,z);
						gsl_vector_set(Zsymbolic,branchIndex,1.0);
						break;
				}
			}
			tmplink=tmplink->nextlink;
		}
		tmp=tmp->next;
	}
	*Zptr = Z;
	*Zsymbolicptr = Zsymbolic;
	*Kptr = K;
	//printf("\nZ = \n");
	//printcomplvect(Z,b);
	/*printf("\nZsym = \n");
	printvect(Zsymbolic,b);
	printf("\nK = \n");
	printcomplvect(K,b);*/
}