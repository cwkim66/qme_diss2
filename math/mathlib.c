#include "mathlib.h"

void gauss_gen(double *r1, double *r2)
{
	unsigned long m1 = randomMT();
	unsigned long m2 = randomMT();

	m1 = m1 ? m1 : m1+1;
	m2 = m2 ? m2 : m2+1;

	double d1 = (double)(m1) / ((double)UINT_MAX);
	double d2 = (double)(m2) / ((double)UINT_MAX);

	double pref = sqrt(-2.0 * log(d1));
	double arg  = 2.0 * M_PI * d2;

	*r1 = pref * cos(arg);
	*r2 = pref * sin(arg);
}

int dround(double d)
{
	if(d > 0.0) return (int)(d + 0.5);
	else return (int)(d - 0.5);
}

int factorial(int n)
{
	int i;
	int fact = 1;

	for(i=0; i<n; i++) fact *= (i+1);

	return fact;
}

int comb(int n, int r)
{
	int i;
	int nCr = 1;

	for(i=0; i<r; i++){
		nCr *= (n-i);
		nCr /= (i+1);
	}

	return nCr;
}

int multicomb(int n, int r)
{
	// Combination with repetition
	// nHr = (n+r-1)Cr

	if(r > (n+r-1)/2) return comb(n+r-1, n-1);
	else return comb(n+r-1, r);
}

void intcopy(int *X, int *Y, int ndim)
{
	int i;
	for(i=0; i<ndim; i++) Y[i] = X[i];
}

void d_init(double *darr, int nelem)
{
	int i;
	for(i=0; i<nelem; i++) darr[i] = 0.0;
}

void normalize(double *X, int dim)
{
	int i;
	double norm = 0.0;

	for(i=0; i<dim; i++){
		norm += X[i] * X[i];
	}

	norm = sqrt(norm);

	for(i=0; i<dim; i++){
		X[i] /= norm;
	}
}

void Axv(double *Y, double *A, double *X, double alpha, double beta, int ndim)
{
	int one = 1;
	char N = 'N', T = 'T';

	dgemv(&N, &ndim, &ndim, &alpha, A, &ndim, X, &one, &beta, Y, &one);
}

double v1TAv2(double *v1, double *A, double *v2, int ndim)
{
	int one = 1;
	double val;
	double *buf = malloc(ndim * sizeof(double));

	Axv(buf, A, v2, 1.0, 0.0, ndim);
	val = ddot(&ndim, v1, &one, buf, &one);

	free(buf);
	return val;
}

void AxB(double *C, double *A, double *B, int m, int n, int common, int mode)
{
	double doubI  = 1.0, doub0 = 0.0, dnegI = -1.0;
	char N = 'N', T = 'T';

	if(mode == 1) {
		//A(m*common)*B(common*n) = C(m*n)
		dgemm(&N, &N, &m, &n, &common, &doubI, A, &m, B, &common, &doub0, C, &m);
	}
	else if(mode == 2) {
		//AT(m*common)*B(common*n) = C(m*n)
		dgemm(&T, &N, &m, &n, &common, &doubI, A, &m, B, &common, &doub0, C, &m);
	}
	else if(mode == 3) {
		//-A(m*common)*B(common*n) = C(m*n)
		dgemm(&N, &N, &m, &n, &common, &dnegI, A, &m, B, &common, &doub0, C, &m);
	}
	else if(mode == 4) {
		//-AT(m*common)*B(common*n) = C(m*n)
		dgemm(&T, &N, &m, &n, &common, &dnegI, A, &m, B, &common, &doub0, C, &m);
	}
	else if(mode == 5) {
		//A(m*common)*BT(common*n) = C(m*n)
		dgemm(&N, &T, &m, &n, &common, &doubI, A, &m, B, &common, &doub0, C, &m);
	}
	else if(mode == 6) {
		//AT(m*common)*BT(common*n) = C(m*n)
		dgemm(&T, &T, &m, &n, &common, &doubI, A, &m, B, &common, &doub0, C, &m);
	}
}

void trans(double *U, double *M, int ndim)
{
	int i, j;

	double *M_tmp = malloc(ndim * ndim * sizeof(double));

	AxB(M_tmp, M, U, ndim, ndim, ndim, 1);
	AxB(M, U, M_tmp, ndim, ndim, ndim, 2);

	free(M_tmp);
}

void trans_inv(double *U, double *M, int ndim)
{
	int i, j;

	double *M_tmp = malloc(ndim * ndim * sizeof(double));

	AxB(M_tmp, M, U, ndim, ndim, ndim, 5);
	AxB(M, U, M_tmp, ndim, ndim, ndim, 1);

	free(M_tmp);
}

void dmat_comm_liou(double *dmat, double *dliou, int ndim)
{
	// Construct a Liouville space operator for [dmat, *]

	int i, j, k;
	int ndimsq = ndim * ndim;
	int ndimp4 = ndimsq * ndimsq;

	for(i=0; i<ndimp4; i++) dliou[i] = 0.0;

	for(i=0; i<ndim; i++) {
		for(j=0; j<ndim; j++) {
			for(k=0; k<ndim; k++) {
				dliou[(i + ndim*j) + ndimsq*(k + ndim*j)] += dmat[i + ndim*k];
				dliou[(i + ndim*j) + ndimsq*(i + ndim*k)] -= dmat[k + ndim*j];
			}
		}
	}
}

void AdpB(double *AdpB, double *A, int ncolA, int nrowA, double *B, int ncolB, int nrowB)
{
	int i, j;
	int nrowAB = nrowA * nrowB;
	int ncolAB = ncolA * ncolB;

	for(i=0; i<nrowAB; i++){
		for(j=0; j<ncolAB; j++){
			AdpB[i + nrowAB*j] = A[(i/nrowB) + nrowA*(j/ncolB)] * B[(i%nrowB) + nrowB*(j%ncolB)];
		}
	}
}

void diag(double *D, double *H_eig, double *H, int ndim)
{
	char V = 'V', U = 'U';
	int i,j,info;
	int intI = 1;
	int ndim2 = ndim*ndim;
	int lwork = 2*ndim2+6*ndim+2;
	int liwork = 5*ndim+4;
	double *evec = (double *)malloc(sizeof(double)*ndim2);
	double *eval = (double *)malloc(sizeof(double)*ndim);
	double *work  =(double *) malloc(sizeof(double)*lwork);
	int *iwork = (int *)malloc(sizeof(int)*liwork);

	//copy H to evec
	dcopy(&ndim2, H, &intI, evec, &intI);

	//upper triangular matrix.
	for(i=0; i<ndim; i++) {
		for(j=i+1; j<ndim; j++) {
			evec[ndim*i+j] = 0.0;
		}
	}
 
	//lots of diagonalization function at MKL
	//MKL user guide says that syevd is efficient for a single call..
	dsyevd(&V, &U, &ndim, evec, &ndim, eval, work, &lwork, iwork, &liwork, &info);

	dcopy(&ndim2, evec, &intI, D, &intI);
	dcopy(&ndim, eval, &intI, H_eig, &intI);

	free(work);
	free(iwork);
	free(evec);
	free(eval);
}

void gauss_elim(double *A, double *x, double *b, int ndim)
{
	// Solve Ax = b where x is the unknown

	int i, j, k;
	int ndimsq = ndim * ndim;

	double *Ac = malloc(ndimsq * sizeof(double));

	for(i=0; i<ndimsq; i++) Ac[i] = A[i];
	for(i=0; i<ndim; i++) x[i] = b[i];

	for(i=0; i<ndim; i++){
		for(j=i+1; j<ndim; j++){
			double fac = Ac[j + ndim*i] / Ac[i + ndim*i];
			for(k=i; k<ndim; k++){
				Ac[j + ndim*k] -= fac * Ac[i + ndim*k];
			}
			x[j] -= fac * x[i];
		}
	}
	// Ac is now upper triangular

	for(i=0; i<ndim; i++){
		for(j=i+1; j<ndim; j++){
			int ii = (ndim-1 - i);
			int jj = (ndim-1 - j);
			double fac = Ac[jj + ndim*ii] / Ac[ii + ndim*ii];
			Ac[jj + ndim*ii] -= fac * Ac[ii + ndim*ii];
			x[jj] -= fac * x[ii];
		}
	}
	// Ac is now diagonal

	for(i=0; i<ndim; i++){
		x[i] /= Ac[i + ndim*i];
		Ac[i + ndim*i] /= Ac[i + ndim*i];
	}
	// Ac is now identity, and we have solved the equation

#if 0
	for(i=0; i<ndim; i++){
		double sum = 0.0;
		for(j=0; j<ndim; j++){
			sum += A[i + ndim*j] * x[j];
		}
		fprintf(stderr, "%13.6le %13.6le\n", sum, b[i]);
	}
#endif

	free(Ac);
}

void znormalize(complex *X, int dim)
{
	int i;
	double norm = 0.0;

	for(i=0; i<dim; i++){
		norm += X[i].real * X[i].real + X[i].imag * X[i].imag;
	}

	norm = sqrt(norm);

	for(i=0; i<dim; i++){
		X[i].real /= norm;
		X[i].imag /= norm;
	}
}

void zAxzv(complex *Y, complex *A, complex *X, complex alpha, complex beta, int ndim)
{
	int one = 1;
	char N = 'N';

	zgemv(&N, &ndim, &ndim, &alpha, A, &ndim, X, &one, &beta, Y, &one);
}

void zATxzv(complex *Y, complex *A, complex *X, complex alpha, complex beta, int ndim)
{
	int one = 1;
	char T = 'T';

	zgemv(&T, &ndim, &ndim, &alpha, A, &ndim, X, &one, &beta, Y, &one);
}

void zsyAxzv(complex *Y, complex *A, complex *X, complex alpha, complex beta, int ndim)
{
	int one = 1;
	char U = 'U';

	zsymv(&U, &ndim, &alpha, A, &ndim, X, &one, &beta, Y, &one);
}

void zAxzB(complex *C, complex *A, complex *B, int m, int n, int common, int mode)
{
	complex zposRe, znegRe, zzero;
	zposRe.real = 1.0, zposRe.imag = 0.0;
	znegRe.real = -1.0, znegRe.imag = 0.0;
	zzero.real = 0.0, zzero.imag = 0.0;
	char N = 'N', CJ = 'C';

	if(mode == 1) {
		//A(m*common)*B(common*n) = C(m*n)
		zgemm(&N, &N, &m, &n, &common, &zposRe, A, &m, B, &common, &zzero, C, &m);
	}
	else if(mode == 2) {
		//A*(m*common)*B(common*n) = C(m*n)
		zgemm(&CJ, &N, &m, &n, &common, &zposRe, A, &m, B, &common, &zzero, C, &m);
	}
	else if(mode == 3) {
		//-A(m*common)*B(common*n) = C(m*n)
		zgemm(&N, &N, &m, &n, &common, &znegRe, A, &m, B, &common, &zzero, C, &m);
	}
	else if(mode == 4) {
		//-A*(m*common)*B(common*n) = C(m*n)
		zgemm(&CJ, &N, &m, &n, &common, &znegRe, A, &m, B, &common, &zzero, C, &m);
	}
	else if(mode == 5) {
		//A(m*common)*B*(common*n) = C(m*n)
		zgemm(&N, &CJ, &m, &n, &common, &zposRe, A, &m, B, &common, &zzero, C, &m);
	}
	else if(mode == 6) {
		//A*(m*common)*B*(common*n) = C(m*n)
		zgemm(&CJ, &CJ, &m, &n, &common, &zposRe, A, &m, B, &common, &zzero, C, &m);
	}
}

void zdir_prod(complex *C, complex *A, complex *B, int dimC, int dimA, int dimB)
{
	// Matrix direct product between square matrices

	int i, j, k, l;
	int start_c, start_r;
	int cn, rn;

	if(dimC != dimA * dimB){
		fprintf(stderr, "Inconsistency in dir_prod: dimC = %d, dimA = %d, dimB = %d\n", dimC, dimA, dimB);
		exit(EXIT_FAILURE);
	}

	for(i=0; i<dimA; i++){
		for(j=0; j<dimA; j++){
			start_c = i * dimB;
			start_r = j * dimB;
			for(k=0; k<dimB; k++){
				for(l=0; l<dimB; l++){
					cn = start_c + k;
					rn = start_r + l;

					z_mult(&C[cn * dimC + rn], A[i * dimA + j], B[k * dimB + l]);
				}
			}
		}
	}
}

void zdiag_wrap(complex *D, complex *H_eig, complex *H, int ndim)
{
	int i;
	int ndim2 = ndim * ndim;
	double *Dbuf    = malloc(ndim2 * sizeof(double));
	double *Heigbuf = malloc(ndim2 * sizeof(double));
	double *Hbuf    = malloc(ndim2 * sizeof(double));

	for(i=0; i<ndim2; i++) Hbuf[i] = H[i].real;

	diag(Dbuf, Heigbuf, Hbuf, ndim);

	for(i=0; i<ndim; i++){
		H_eig[i].real = Heigbuf[i];
		H_eig[i].imag = 0.0;
	}

	for(i=0; i<ndim2; i++){
		D[i].real = Dbuf[i];
		D[i].imag = 0.0;
	}

	free(Dbuf);
	free(Heigbuf);
	free(Hbuf);
}

void ztrans(complex *U, complex *M, int ndim)
{
	complex *M_tmp = malloc(ndim * ndim * sizeof(complex));

	zAxzB(M_tmp, M, U, ndim, ndim, ndim, 1);
	zAxzB(M, U, M_tmp, ndim, ndim, ndim, 2);

	free(M_tmp);
}

void ztrans_inv(complex *U, complex *M, int ndim)
{
	complex *M_tmp = malloc(ndim * ndim * sizeof(complex));

	zAxzB(M_tmp, M, U, ndim, ndim, ndim, 5);
	zAxzB(M, U, M_tmp, ndim, ndim, ndim, 1);

	free(M_tmp);
}

void zAdpzB(complex *zAdpzB, complex *zA, int ncolA, int nrowA, complex *zB, int ncolB, int nrowB)
{
	int i, j;
	int nrowAB = nrowA * nrowB;
	int ncolAB = ncolA * ncolB;

	for(i=0; i<nrowAB; i++){
		for(j=0; j<ncolAB; j++){
			z_mult(&zAdpzB[i + nrowAB*j], zA[(i/nrowB) + nrowA*(j/ncolB)], zB[(i%nrowB) + nrowB*(j%ncolB)]);
		}
	}
}

void zpart_trace(complex *zA, complex *zB, int nA, int nB)
{
	int i, j, k;
	int count;

	if(nA%nB) {
		fprintf(stderr, "Error in %s: nA %d is not divisible by nB %d\n", __func__, nA, nB);
		exit(EXIT_FAILURE);
	}

	z_init(zB, nB*nB);
	int nT = nA/nB;

	for(i=0; i<nB; i++) {
		for(j=0; j<nB; j++) {
			count = (i*nT) + nA*(j*nT);

			for(k=0; k<nT; k++) {
				zB[i + nB*j].real += zA[count].real;
				zB[i + nB*j].imag += zA[count].imag;
				count += (nA + 1);
			}
		}
	}
}

void zprod_trace(complex *zA, complex *zB, int ndim, complex *val)
{
	int i, j;
	int count = 0;
	val->real = 0.0, val->imag = 0.0;

	complex *zprod = malloc(ndim * ndim * sizeof(complex));
	z_init(zprod, ndim * ndim);

	zAxzB(zprod, zA, zB, ndim, ndim, ndim, 1);

	for(i=0; i<ndim; i++) {
		val->real += zprod[count].real;
		val->imag += zprod[count].imag;
		count += (ndim + 1);
	}

	free(zprod);
}

void zmat_comm_liou(complex *zmat, complex *zliou, int ndim)
{
	// Construct a Liouville space operator for [zmat, *]

	int i, j, k;
	int ndimsq = ndim * ndim;
	int ndimp4 = ndimsq * ndimsq;

	z_init(zliou, ndimp4);

	for(i=0; i<ndim; i++) {
		for(j=0; j<ndim; j++) {
			for(k=0; k<ndim; k++) {
				zliou[(i + ndim*j) + ndimsq*(k + ndim*j)].real += zmat[i + ndim*k].real;
				zliou[(i + ndim*j) + ndimsq*(k + ndim*j)].imag += zmat[i + ndim*k].imag;

				zliou[(i + ndim*j) + ndimsq*(i + ndim*k)].real -= zmat[k + ndim*j].real;
				zliou[(i + ndim*j) + ndimsq*(i + ndim*k)].imag -= zmat[k + ndim*j].imag;
			}
		}
	}
}

void zmat_acomm_liou(complex *zmat, complex *zliou, int ndim)
{
	// Construct a Liouville space operator for {zmat, *}

	int i, j, k;
	int ndimsq = ndim * ndim;
	int ndimp4 = ndimsq * ndimsq;

	z_init(zliou, ndimp4);

	for(i=0; i<ndim; i++) {
		for(j=0; j<ndim; j++) {
			for(k=0; k<ndim; k++) {
				zliou[(i + ndim*j) + ndimsq*(k + ndim*j)].real += zmat[i + ndim*k].real;
				zliou[(i + ndim*j) + ndimsq*(k + ndim*j)].imag += zmat[i + ndim*k].imag;

				zliou[(i + ndim*j) + ndimsq*(i + ndim*k)].real += zmat[k + ndim*j].real;
				zliou[(i + ndim*j) + ndimsq*(i + ndim*k)].imag += zmat[k + ndim*j].imag;
			}
		}
	}
}

void z_init(complex *zarr, int nelem)
{
	int i;
	for(i=0; i<nelem; i++) z_zero(&zarr[i]);
}

void z_zero(complex *a)
{
	a->real = 0.0;
	a->imag = 0.0;
}

double z_abs(complex a)
{
	return sqrt(a.real * a.real + a.imag * a.imag);
}

void z_conj(complex *b, complex a)
{
	b->real = a.real;
	b->imag = -a.imag;
}

void z_mult(complex *c, complex a, complex b)
{
	// c = a * b
	c->real = a.real * b.real - a.imag * b.imag;
	c->imag = a.real * b.imag + a.imag * b.real;
}

void z_multc(complex *c, complex a, complex b)
{
	// c = (a*) * b
	c->real = a.real * b.real + a.imag * b.imag;
	c->imag = a.real * b.imag - a.imag * b.real;
}

void z_div(complex *c, complex a, complex b)
{
	// c = a / b
	c->real = b.real * b.real + b.imag * b.imag;
	c->imag = (a.imag * b.real - a.real * b.imag) / c->real;
	c->real = (a.real * b.real + a.imag * b.imag) / c->real;
}

void fft_forw(complex *data, int nelem)
{
	MKL_LONG errormessage;
	DFTI_DESCRIPTOR_HANDLE fft;

	errormessage = DftiCreateDescriptor(&fft, DFTI_DOUBLE, DFTI_COMPLEX, 1, nelem);
	errormessage = DftiCommitDescriptor(fft);
	errormessage = DftiComputeForward(fft, data);
	errormessage = DftiFreeDescriptor(&fft);
}

void fft_back(complex *data, int nelem)
{
	MKL_LONG errormessage;
	DFTI_DESCRIPTOR_HANDLE fft;

	errormessage = DftiCreateDescriptor(&fft, DFTI_DOUBLE, DFTI_COMPLEX, 1, nelem);
	errormessage = DftiCommitDescriptor(fft);
	errormessage = DftiComputeBackward(fft, data);
	errormessage = DftiFreeDescriptor(&fft);
}
