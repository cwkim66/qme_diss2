#ifndef D_QME
#define D_QME

#ifndef D_GENLIB
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#define D_GENLIB
#endif

//mkl
#ifndef D_MKL
#define D_MKL
#include "mkl.h"
#include "mkl_blas.h"
#include "mkl_lapack.h"
#define complex MKL_Complex16
#endif

static int one = 1;

typedef struct {
	double freq, spread;
	double temp, coth;
	double *gamma, *gexci;

	double *Idiss, *Kdiss, *Jdiss;
} t_osc;

typedef struct {
	double freq;
	double temp, coth;

	int spdn, oscn;
	int nsys, nsyssq;
	int nqst, nqstsq;
	int ndim, ndimsq;

	complex *aHpa;
	double *gamma, *gexci;
	double **FC;
	complex **FCU;
	complex **Req, **H, **U;
} t_qosc;

typedef struct {
	int unit;
	int nsys, nsyssq;

	double dt, dtout, T;
	int n, nout, nsto;

	double *Hsys, *Hsys0, *Ediag;

	int nspd;
	int *nosc;
	t_osc **osc;

	double temp_eq;
	double *reorg;
	complex **gt, **g1t, **g2t, **Gt;
	double *rate;

	double temp_uni;
	bool bTEMP_UNI;

	bool bSPLIT;
} t_fret;

typedef struct {
	int unit;
	int nsys, nsyssq;

	double dt, dtout, T;
	int n, nout, nsto;

	double *Hsys, *Usys, *Eexci, *Eexci0;

	int nspd;
	int *nosc;
	t_osc **osc;

	double *reaabb, *reabbb;
	complex **gtaabb, **gtabbb, **g1abbb, **g2abab;
	complex **Gt, **Zt, **Nt;
	double *rate;

	double temp_uni;
	bool bTEMP_UNI;

	bool bSPLIT;
} t_mrt;

typedef struct {
	int jtype, unit;
	int nsys, nsyssq;
	bool bEXCI, bSPLIT;

	double dt, dtout, T;
	int n, nout, nsto;

	double *rate_pop;
	double *pop0, *popt, *popt_sto;
	double *E0, *Hevec;

	int nspd;
	int *nosc, nosc_all;
	double *freq, **Idiss, **rate_diss;
	double *disst;
} t_dyn_inco;

typedef struct {
	int jtype, unit;
	int nsys, nsyssq;

	int nosc_tr;
	t_qosc *qosc;

	double dtout, T;
	int nout;

	complex ***tr1, ***tr2;
} t_tr;

typedef struct {
	int jtype, unit;
	int nsys, nsyssq, nstd;
	bool bDISS, bHSYS_ADJUST;

	int ntraj;
	int type_fcn;
	int rand_start;
	double freq_split, pref_split;

	int nout;
	double dtout;
	double *popt_avg, **disst_avg;

	double *Hsys_prist, *Hsys_offset, *Hsys_std, **std;

	int nspd;
	int *nosc, nosc_all;
	t_osc **osc0;

	t_fret *fret;
	t_mrt *mrt;
	t_dyn_inco *dyn;
} t_split;

typedef struct {
	/* General parameters */
	int jtype, unit;
	bool bINCO, bDISS, bEXCI, bDEBUG, bSPLIT;

	double dt1, dtout1, T1;
	double dt2, dtout2, T2;

	/* System Hamiltonian and density matrix */
	int nsys, nsyssq;
	int initexc, initcond;
	double *Hsys, *Hevec, *Heval;

	double *pop0, *popt;
	complex *rho0, *rhot;

	/* Bath oscillators */
	int nspd_ana, nspd_num, nspd;
	int *nosc;
	t_osc **osc;

	/* QME construction */
	t_fret *fret;
	t_mrt *mrt;

	/* Direct trace evaluation */
	t_tr *tr;
	int nosc_tr;
	int *spdn_tr, *oscn_tr;
	double *cut_tr;

	/* Split into low and high frequency + optimization */
	t_split *split;

	int mem;
} t_qme;

#include "constants.h"
#include "mathlib.h"
#include "set_param.h"
#include "fret.h"
#include "mrt.h"
#include "calc_tr.h"
#include "dyn_inco.h"
#include "util.h"

int main(int argc, char *argv[]);

#endif
