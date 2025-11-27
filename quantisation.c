#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <math.h>
#include "createMat.h"
#include "consoleIN.h"
#include "createMat.h"
#include "quantisation.h"

double cot(double theta)
{
	double cot_theta; //for small theta use taylors
	if (fabs(theta) < 1e-3)	cot_theta = 1.0/theta - theta/3.0 - (theta*theta*theta)/45.0;
	else cot_theta = cos(theta)/sin(theta);
	return cot_theta;
}

bool checkifzerocomp(gsl_complex z)
{
	if(gsl_complex_abs(z) < 1e-13) return true;
	return false;
}

void realtocomplexmat(int n, int b, gsl_matrix * Q, gsl_matrix_complex** Qcompptr,
	gsl_matrix * B, gsl_matrix_complex ** Bcompptr)
{
	gsl_matrix_complex * Qcomp = gsl_matrix_complex_alloc(n-1,b);
	gsl_matrix_complex * Bcomp = gsl_matrix_complex_alloc(b-n+1,b);
	gsl_matrix_complex_set_zero(Qcomp);
	gsl_matrix_complex_set_zero(Bcomp);
	for (int i = 0; i < n-1; i++)
	{
		for (int j = 0; j < b; j++)
		{
			double val = gsl_matrix_get(Q, i, j);
			gsl_matrix_complex_set(Qcomp, i, j, gsl_complex_rect(val, 0.0));
		}
	}
	for (int i = 0; i < b-n+1; i++)
	{
		for (int j = 0; j < b; j++)
		{
			double val = gsl_matrix_get(B, i, j);
			gsl_matrix_complex_set(Bcomp, i, j, gsl_complex_rect(val, 0.0));
		}
	}
	*Qcompptr = Qcomp;
	*Bcompptr = Bcomp;
}

gsl_vector * getrealvector (int b,gsl_vector_complex* vcomp)
{
	gsl_vector * vreal = gsl_vector_alloc(b);
	for(int i=0;i<b;i++)
	{
		gsl_vector_set(vreal,i,GSL_REAL(gsl_vector_complex_get(vcomp,i)));
	}
	return vreal;
}

/*void creatediscsys(int N, int b, double tj, int k, gsl_matrix_complex **Zdiscptr,
 gsl_vector_complex **Kdiscptr,gsl_complex * s_kptr,gsl_complex * weight_kptr, 
	gsl_vector_complex* Z, gsl_vector * Zsymbolic,gsl_vector_complex* K)
{
	gsl_matrix_complex * Zdisc = gsl_matrix_complex_alloc(b,b);
	gsl_vector_complex * Kdisc = gsl_vector_complex_alloc(b);
	gsl_matrix_complex_set_zero(Zdisc);
	gsl_vector_complex_set_zero(Kdisc);

	double theta = M_PI*k/N;
	double real1,compl1,real2,compl2,yfactor;
	gsl_complex s_k, s_kinv, weight_k;
	if(k==0)
	{
		real1 = 0.4*N;
		compl1 = 0.0;
		real2 = 0.5*exp(real1);
		compl2 = 0.0;
	}
	else
	{
		real1 = 0.4*k*M_PI*cot(theta);
		compl1 = 0.4*k*M_PI;
		if (fabs(theta) < 1e-3) yfactor = (4.0/45.0) * theta*theta*theta;
		else yfactor = (theta) * (1.0 + cos(theta)*cos(theta)/(sin(theta)*sin(theta))) - cos(theta)/sin(theta);
		real2 = exp(real1)*(cos(compl1)-sin(compl1)*(yfactor));
		compl2 = exp(real1)*(cos(compl1)+sin(compl1)*(yfactor));;
	}
	s_k = gsl_complex_rect(real1/tj,compl1/tj);
	weight_k = gsl_complex_rect(real2,compl2);
	double abss = real1*real1+compl1*compl1;
	s_kinv = gsl_complex_inverse(s_k);;

	for(int i=0;i<b;i++)
	{
		int x = (int) gsl_vector_get(Zsymbolic,i);
		gsl_complex lval = gsl_complex_mul(gsl_vector_complex_get(Z,i), s_k);
		gsl_complex cval = gsl_complex_mul(gsl_vector_complex_get(Z,i), s_kinv);
		gsl_complex vval = gsl_complex_mul(gsl_vector_complex_get(K,i), s_kinv);
		switch(x)
		{
		case 0:
			gsl_matrix_complex_set(Zdisc,i,i,gsl_vector_complex_get(Z,i));
			break;
		case 1:
			gsl_matrix_complex_set(Zdisc,i,i,lval);
			break;
		case -1:
			if(checkifzerocomp(gsl_vector_complex_get(K,i))) gsl_matrix_complex_set(Zdisc,i,i,cval);
			else gsl_vector_complex_set(Kdisc,i,vval);	
			break;
		}
	}
	
	*Zdiscptr = Zdisc;
	*Kdiscptr = Kdisc;
	*s_kptr = s_k;
	*weight_kptr = weight_k;
}*/

void linearsystem(int n, int b,	gsl_matrix_complex *Z,
	gsl_vector_complex *K, gsl_matrix_complex *Q, gsl_matrix_complex* B,
	gsl_vector_complex **Xptr)
{
	gsl_matrix_complex * M = gsl_matrix_complex_alloc(b,b);
	gsl_vector_complex * X = gsl_vector_complex_alloc(b);
	gsl_vector_complex * C = gsl_vector_complex_alloc(b);
	gsl_vector_complex_set_zero(C);
	gsl_matrix_complex * temp1 = gsl_matrix_complex_alloc(b-n+1,b);
	gsl_vector_complex * temp2 = gsl_vector_complex_alloc(b-n+1);
	gsl_vector_complex * temp3 = gsl_vector_complex_alloc(b);
	gsl_complex alpha = gsl_complex_rect(1.0,0.0);
	gsl_complex beta = gsl_complex_rect(0.0,0.0);
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,alpha,B,Z,beta,temp1);
	gsl_blas_zgemv(CblasNoTrans,alpha,B,K,beta,temp2);

	for(int i = 0;i<b-n+1;i++)
	{
		gsl_vector_complex_set(C,i+n-1,gsl_vector_complex_get(temp2,i));
	}
	for(int j=0;j<b;j++)
	{
		for(int i=0;i<b;i++)
		{
			if(i<n-1) gsl_matrix_complex_set(M,i,j,gsl_matrix_complex_get(Q,i,j));
			else gsl_matrix_complex_set(M,i,j,gsl_matrix_complex_get(temp1,i-n+1,j));
		}
	}
	gsl_permutation *p = gsl_permutation_alloc(b);
	int signum;
	gsl_linalg_complex_LU_decomp(M, p, &signum);
	gsl_linalg_complex_LU_solve(M, p, C, X);
	*Xptr = X;
	//gsl_blas_zgemv(CblasNoTrans,alpha,M,X,beta,temp3);
	//gsl_vector_complex_sub(temp3, C);
	//double norm = gsl_blas_dznrm2(temp3);
	//if(norm>0.4) printf("%lf\n",norm);
	gsl_matrix_complex_free(temp1);
	gsl_matrix_complex_free(M);
	gsl_vector_complex_free(temp2);
	gsl_vector_complex_free(temp3);
	gsl_vector_complex_free(C);
	gsl_permutation_free(p);
}

/*void solveI(int n,int b,gsl_matrix ** ibptr,double timewindow, int timecount,
	gsl_matrix * Q, gsl_matrix *B, gsl_vector * Zsymbolic, gsl_vector_complex * K
	,gsl_vector_complex * Z)
{
	double timestep = timewindow/timecount;
	int N=32;
	gsl_matrix* ib = gsl_matrix_alloc(b,timecount);
	gsl_matrix_complex *Zdisc,*Bcomp,*Qcomp;
	realtocomplexmat(n,b,Q, &Qcomp,B,&Bcomp);
	gsl_vector * ibtempsum = gsl_vector_alloc(b);


	gsl_vector_complex *X,*Kdisc;
	gsl_complex s_k,weight_k;
	gsl_vector * ireal;

	for(int i=0;i<timecount;i++)
	{
		gsl_vector_set_zero(ibtempsum);
		double tj = (i+1)*timestep;
		for(int k=0;k<N;k++)
		{
			creatediscsys(N,b,tj,k,&Zdisc,&Kdisc,&s_k,&weight_k,Z,Zsymbolic,K);
			linearsystem(n, b, Zdisc, Kdisc, Qcomp, Bcomp, &X);
			gsl_complex tempcomp = gsl_complex_mul_real(s_k,tj);
			gsl_complex expcomp = gsl_complex_exp(tempcomp);
			gsl_complex scalar = gsl_complex_mul(weight_k,expcomp);
			gsl_vector_complex_scale(X,scalar);
			ireal = getrealvector(b,X);
			gsl_blas_daxpy(1.0,ireal,ibtempsum);
		}
		gsl_vector_scale(ibtempsum,2.0/tj);
		printvect(ibtempsum,b);
		for(int j=0;j<b;j++)
		{
			gsl_matrix_set(ib,j,i,gsl_vector_get(ibtempsum,j));
		}
	}
	*ibptr = ib;
}*/
void solveI(int n,int b,gsl_matrix ** ibptr,double timewindow, int timecount,
	gsl_matrix * Q, gsl_matrix *B, gsl_vector * Zsymbolic, gsl_vector_complex * K,
	gsl_vector_complex * Z)
{
	double timestep = timewindow / (double)timecount;
	int N = 32;
	gsl_matrix* ib = gsl_matrix_alloc(b,timecount);
	gsl_matrix_complex *Qcomp,*Bcomp,*Zdisc;
	gsl_vector_complex *Kdisc,*X;
	gsl_complex s_k, weight_k;
	realtocomplexmat(n, b, Q, &Qcomp, B, &Bcomp);
	gsl_vector * ibtempsum_real = gsl_vector_alloc(b);
	for (int it = 0; it < timecount; ++it)
	{
		gsl_vector_set_zero(ibtempsum_real);
		double tj = (it + 1) * timestep;
		for (int k = 0; k < N; ++k)
		{
			creatediscsys(N, b, tj, k, &Zdisc, &Kdisc, &s_k, &weight_k, Z, Zsymbolic, K);
			linearsystem(n, b, Zdisc, Kdisc, Qcomp, Bcomp, &X);
			gsl_vector_complex *Xscaled = gsl_vector_complex_alloc(b);
			gsl_vector_complex_memcpy(Xscaled, X);
			gsl_vector_complex_scale(Xscaled, weight_k);

			// convert to real part vector and accumulate
			gsl_vector *Xreal = getrealvector(b, Xscaled);
			gsl_blas_daxpy(1.0, Xreal, ibtempsum_real);
		}
		gsl_vector_scale(ibtempsum_real, 2.0 / (5.0 * tj));
		//printvect(ibtempsum_real, b);
		for (int j = 0; j < b; ++j)
		{
			gsl_matrix_set(ib, j, it, gsl_vector_get(ibtempsum_real, j));
		}
	}
	*ibptr = ib;
}

void creatediscsys(int N, int b, double tj, int k,
				   gsl_matrix_complex **Zdiscptr,
				   gsl_vector_complex **Kdiscptr,
				   gsl_complex *s_kptr,
				   gsl_complex *weight_kptr,
				   gsl_vector_complex* Z,
				   gsl_vector * Zsymbolic,
				   gsl_vector_complex* K)
{
	// allocate result containers
	gsl_matrix_complex * Zdisc = gsl_matrix_complex_alloc(b,b);
	gsl_vector_complex * Kdisc = gsl_vector_complex_alloc(b);
	gsl_matrix_complex_set_zero(Zdisc);
	gsl_vector_complex_set_zero(Kdisc);

	// Talbot parameters: N is M in formulas
	// handle k==0 special case
	if (k == 0)
	{
		double beta0 = 2.0 * (double)N / 5.0;          // beta_0 (real)
		double gamma0 = 0.5 * exp(beta0);             // gamma_0 (real)
		// s_k = beta_k / tj  (here beta_k real)
		*s_kptr = gsl_complex_rect(beta0 / tj, 0.0);
		*weight_kptr = gsl_complex_rect(gamma0, 0.0);
	}
	else
	{
		double theta = M_PI * (double)k / (double)N;   // θ = kπ/M
		double cot_theta;
		double Delta; // = (kπ/M)*(1+cot^2 θ) - cot θ   (use series for small theta)

		// small-angle guard to avoid catastrophic cancellation / overflow
		if (fabs(theta) < 1e-3)
		{
			// use series:
			// cot θ ≈ 1/θ - θ/3 - θ^3/45
			cot_theta = 1.0/theta - theta/3.0 - (theta*theta*theta)/45.0;
			// Delta ≈ (4/45) * θ^3  (from series simplification)
			Delta = (4.0/45.0) * theta*theta*theta;
		}
		else
		{
			cot_theta = cos(theta)/sin(theta);
			double cot2 = cot_theta * cot_theta;
			Delta = (theta) * (1.0 + cot2) - cot_theta;
		}

		// build beta_k (node) = (2kπ/5) * (cotθ + i)
		double beta_real = (2.0 * (double)k * M_PI / 5.0) * cot_theta;
		double beta_imag = (2.0 * (double)k * M_PI / 5.0);

		// exp(beta_real) and trig of beta_imag
		double exp_beta_real = exp(beta_real);
		double cos_b = cos(beta_imag);
		double sin_b = sin(beta_imag);

		// gamma_k = [ 1 + i * Delta ] * e^{beta_k}
		// Real(gamma) = e^{beta_real} * ( cos_b - Delta * sin_b )
		// Imag(gamma) = e^{beta_real} * ( sin_b + Delta * cos_b )
		double gamma_re = exp_beta_real * ( cos_b - Delta * sin_b );
		double gamma_im = exp_beta_real * ( sin_b + Delta * cos_b );

		// s_k = beta_k / t_j
		*s_kptr = gsl_complex_rect(beta_real / tj, beta_imag / tj);
		*weight_kptr = gsl_complex_rect(gamma_re, gamma_im);
	}

	// Now build discrete Zdisc and Kdisc using s_k (and 1/s_k where required)
	// compute 1/s_k safely
	gsl_complex s_k = *s_kptr;
	gsl_complex one = gsl_complex_rect(1.0, 0.0);
	gsl_complex s_k_inv;
	s_k_inv = gsl_complex_inverse(s_k);
	for (int i = 0; i < b; ++i)
	{
		int x = (int) gsl_vector_get(Zsymbolic, i);
		gsl_complex z_i = gsl_vector_complex_get(Z, i);
		gsl_complex k_i = gsl_vector_complex_get(K, i);
		gsl_complex lval = gsl_complex_mul(z_i, s_k);
		gsl_complex cval = gsl_complex_mul(z_i, s_k_inv);
		gsl_complex vval = gsl_complex_mul(k_i, s_k_inv);
		if (x == 0)	gsl_matrix_complex_set(Zdisc, i, i, z_i);
		else if (x == 1) gsl_matrix_complex_set(Zdisc, i, i, lval);
		else if (x == -1)
		{
			if (checkifzerocomp(k_i)) gsl_matrix_complex_set(Zdisc, i, i, cval);
			else gsl_vector_complex_set(Kdisc, i, vval);
		}
	}
	*Zdiscptr = Zdisc;
	*Kdiscptr = Kdisc;
}
