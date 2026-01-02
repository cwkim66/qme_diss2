#ifndef D_MATHLIB
#define D_MATHLIB

#ifndef D_GENLIB
#define D_GENLIB
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#endif

#ifndef D_MKL
#define D_MKL
#include "mkl.h"
#include "mkl_blas.h"
#include "mkl_lapack.h"
#endif

#ifndef D_MKL_COMPLEX16
#define D_MKL_COMPLEX16
#define complex MKL_Complex16
#endif

#include "rand.h"

void gauss_gen(double *r1, double *r2);
int dround(double d);
int factorial(int n);
int comb(int n, int r);
int multicomb(int n, int r);
void intcopy(int *X, int *Y, int ndim);
void d_init(double *darr, int nelem);
void normalize(double *X, int dim);
void Axv(double *Y, double *A, double *X, double alpha, double beta, int ndim);
double v1TAv2(double *v1, double *A, double *v2, int ndim);
void AxB(double *C, double *A, double *B, int m, int n, int common, int mode);
void trans(double *U, double *M, int ndim);
void trans_inv(double *U, double *M, int ndim);
void dmat_comm_liou(double *dmat, double *dliou, int ndim);
void AdpB(double *AdpB, double *A, int ncolA, int nrowA, double *B, int ncolB, int nrowB);
void diag(double *D, double *H_eig, double *H, int ndim);
void gauss_elim(double *A, double *x, double *b, int ndim);
void znormalize(complex *X, int dim);
void zAxzv(complex *Y, complex *A, complex *X, complex alpha, complex beta, int ndim);
void zATxzv(complex *Y, complex *A, complex *X, complex alpha, complex beta, int ndim);
void zsyAxzv(complex *Y, complex *A, complex *X, complex alpha, complex beta, int ndim);
void zAxzB(complex *C, complex *A, complex *B, int m, int n, int common, int mode);
void zdir_prod(complex *C, complex *A, complex *B, int dimC, int dimA, int dimB);
void zdiag_wrap(complex *D, complex *H_eig, complex *H, int ndim);
void ztrans(complex *U, complex *M, int ndim);
void ztrans_inv(complex *U, complex *M, int ndim);
void zAdpzB(complex *zAdpzB, complex *zA, int ncolA, int nrowA, complex *zB, int ncolB, int nrowB);
void zpart_trace(complex *zA, complex *zB, int nA, int nB);
void zprod_trace(complex *zA, complex *zB, int ndim, complex *val);
void zmat_comm_liou(complex *zmat, complex *zliou, int ndim);
void zmat_acomm_liou(complex *zmat, complex *zliou, int ndim);
void z_init(complex *zarr, int nelem);
void z_zero(complex *a);
double z_abs(complex a);
void z_conj(complex *b, complex a);
void z_mult(complex *c, complex a, complex b);
void z_multc(complex *c, complex a, complex b);
void z_div(complex *c, complex a, complex b);
void fft_forw(complex *data, int nelem);
void fft_back(complex *data, int nelem);

#endif // D_MATHLIB
