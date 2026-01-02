#include "dyn_inco.h"

void dyn_inco(t_qme *qme)
{
	t_dyn_inco *dyn = malloc(sizeof(t_dyn_inco));
	init_dyn_inco(qme, dyn);

	if(!(qme->bDISS)) {
		prop_pop_inco(dyn);
	} else {
		prop_diss_inco(dyn, -1);
	}
}

void conv_popt_inco(t_dyn_inco *dyn, double *popt_prt)
{
	int i;
	int nsys   = dyn->nsys;
	int nsyssq = dyn->nsyssq;

	double *popt = dyn->popt;
	double *rho  = calloc(nsyssq, sizeof(double));

	for(i=0; i<nsys; i++) {
		rho[i + nsys*i] = popt[i];
	}

	trans_inv(dyn->Hevec, rho, nsys);

	for(i=0; i<nsys; i++) {
		popt_prt[i] = rho[i + nsys*i];
	}

	free(rho);
}

void prop_pop_inco(t_dyn_inco *dyn)
{
	int i, j;
	int one = 1;

	int nsys = dyn->nsys;
	int n    = dyn->n;
	int nsto = dyn->nsto;

	double dt   = dyn->dt;
	double tfac = dyn->unit == 1 ? TUNIT_TO_FS : 1.0;

	double *E0        = dyn->E0;
	double *pop0      = dyn->pop0;
	double *popt      = dyn->popt;
	double *popt_prt  = malloc(nsys * sizeof(double));
	double **popt_sto = dyn->popt_sto;

	FILE *f_pop = fopen("pop.txt", "w");

	// Initial population
	dcopy(&nsys, pop0, &one, popt, &one);

	for(i=0; i<n; i++) {
		if(i % nsto == 0) {
			// Time
			fprintf(f_pop, "%13.6lf", i * dt * tfac);

			// Electronic energy
			double Eelec0 = 0.0;
			for(j=0; j<nsys; j++) {
				Eelec0 += E0[j] * popt[j];
			}

			fprintf(f_pop, "%20.10lf", Eelec0);

			// Population
			if(dyn->bEXCI) {
//				conv_popt_inco(dyn, popt_prt);
				dcopy(&nsys, popt, &one, popt_prt, &one);
			} else {
				dcopy(&nsys, popt, &one, popt_prt, &one);
			}

			for(j=0; j<nsys; j++) {
				fprintf(f_pop, "%13.6lf", popt_prt[j]);
			}
			fprintf(f_pop, "\n");

			if(dyn->bSPLIT) {
				dcopy(&nsys, popt_prt, &one, &popt_sto[(i/nsto) * nsys], &one);
			}
		}
		prop_pop_RK4(popt, dyn->rate_pop, dt, nsys);
	}

	fclose(f_pop);
	free(popt_prt);
}

void prop_diss_inco(t_dyn_inco *dyn, int flag)
{
	int i, j, k, l;
	int one = 1;
	char fname[256];

	int nsys = dyn->nsys;
	int n    = dyn->n;
	int nsto = dyn->nsto;

	double dt   = dyn->dt;
	double tfac = dyn->unit ? TUNIT_TO_FS : 1.0;

	int nspd     = dyn->nspd;
	int *nosc    = dyn->nosc;
	int nosc_all = dyn->nosc_all;

	double *pop0     = dyn->pop0;
	double *popt     = dyn->popt;
	double *popt_prt = malloc(nsys * sizeof(double));

	double *disst = dyn->disst;
	double *corr  = malloc(nosc_all * sizeof(double));

	if(flag >= 0) {
		sprintf(fname, "pop_%04d.txt", flag);
	} else {
		sprintf(fname, "pop.txt");
	}

	FILE *f_pop = fopen(fname, "w");

	// Initialize
	dcopy(&nsys, pop0, &one, popt, &one);
	d_init(disst, nosc_all);

	for(i=0; i<n; i++) {
		if(i % nsto == 0) {
			fprintf(f_pop, "%13.6lf", i * dt * tfac);

			if(dyn->bEXCI) {
//				conv_popt_inco(dyn, popt_prt);
				dcopy(&nsys, popt, &one, popt_prt, &one);
			} else {
				dcopy(&nsys, popt, &one, popt_prt, &one);
			}

			for(j=0; j<nsys; j++) {
				fprintf(f_pop, "%13.6lf", popt_prt[j]);
			}
			fprintf(f_pop, "\n");
		}
		prop_diss_RK4(popt, disst, dyn->rate_pop, dyn->rate_diss, dt, nsys, nosc_all);

		// For removing the linear drift
		if(i == n-2) {
			for(j=0; j<nosc_all; j++) {
				corr[j] = disst[j];
			}
		}
	}

	fclose(f_pop);

	// Print the dissipation info on files
	int count = 0;
	double *freq   = dyn->freq;
	double **Idiss = dyn->Idiss;
	double **Jdiss = dyn->rate_diss;

	double wfac = dyn->unit ? 1.0 / WAVENO_TO_EUNIT : 1.0;
	double Ifac = dyn->unit ? TUNIT_TO_FS : 1.0;
	double Jfac = dyn->unit ? 1.0 / TUNIT_TO_FS : 1.0;

	for(i=0; i<nspd; i++) {
		sprintf(fname, "diss_spd%02d.txt", i+1);
		FILE *f_diss = fopen(fname, "w");

		fprintf(f_diss, "   Freq(waveno)");

		for(j=0; j<nsys; j++) {
			for(k=0; k<j; k++) {
				fprintf(f_diss, "   Idiss %2d<-%2d   Idiss %2d<-%2d   Jdiss %2d<-%2d   Jdiss %2d<-%2d", j+1, k+1, k+1, j+1, j+1, k+1, k+1, j+1);
			}
		}
		fprintf(f_diss, "           diss\n");

		for(j=0; j<nosc[i]; j++) {
			disst[count] -= ((disst[count] - corr[count]) * (n-1));

			fprintf(f_diss, "%15.6le", freq[count] * wfac);

			for(k=0; k<nsys; k++) {
				for(l=0; l<k; l++) {
					fprintf(f_diss, "%15.6le%15.6le", Idiss[count][l + nsys*k] * Ifac, Idiss[count][k + nsys*l] * Ifac);
					fprintf(f_diss, "%15.6le%15.6le", Jdiss[count][l + nsys*k] * Jfac, Jdiss[count][k + nsys*l] * Jfac);
				}
			}
			fprintf(f_diss, "%15.6le\n", disst[count]);

			count++;
		}
		fclose(f_diss);
	}

	free(popt_prt);
}

void calc_dpop(double *in, double *out, double *rate, double dt, int n)
{
	int i, j;

	for(i=0; i<n; i++) {
		out[i] = 0.0;

		for(j=0; j<n; j++) {
			out[i] += rate[i + n*j] * in[j];
			out[i] -= rate[j + n*i] * in[i];
		}

		out[i] *= dt;
	}
}

void calc_ddiss(double *pop, double *out, double **rate, double dt, int nsys, int nosc)
{
	int i, j, k;

	for(i=0; i<nosc; i++) {
		out[i] = 0.0;

		for(j=0; j<nsys; j++) {
			for(k=0; k<j; k++) {
				out[i] += rate[i][j + nsys*k] * pop[k];
				out[i] += rate[i][k + nsys*j] * pop[j];
			}
		}

		out[i] *= dt;
	}
}

void prop_pop_RK4(double *pop, double *rate, double dt, int n)
{
	int one = 1;
	double d_one   = 1.0;
	double d_half  = 1.0 / 2.0;
	double d_third = 1.0 / 3.0;
	double d_sixth = 1.0 / 6.0;

	double *pop_buf = malloc(n * sizeof(double));
	double *dpop1   = malloc(n * sizeof(double));
	double *dpop2   = malloc(n * sizeof(double));
	double *dpop3   = malloc(n * sizeof(double));
	double *dpop4   = malloc(n * sizeof(double));

    // Propagates the auxiliary DM by 4th order Runge-Kutta method.
    // dpop(t) / dt = f(pop(t))

    // dpop1 = dt * f(pop(t))
    calc_dpop(pop, dpop1, rate, dt, n);

    // dpop2 = dt * f(pop(t) + 0.5 * dpop1)
    dcopy(&n, pop, &one, pop_buf, &one);
    daxpy(&n, &d_half, dpop1, &one, pop_buf, &one);
    calc_dpop(pop_buf, dpop2, rate, dt, n);

    // dpop3 = dt * f(pop(t) + 0.5 * dpop2)
    dcopy(&n, pop, &one, pop_buf, &one);
    daxpy(&n, &d_half, dpop2, &one, pop_buf, &one);
    calc_dpop(pop_buf, dpop3, rate, dt, n);

    // dpop4 = dt * f(pop(t) + dpop3)
    dcopy(&n, pop, &one, pop_buf, &one);
    daxpy(&n, &d_one, dpop3, &one, pop_buf, &one);
    calc_dpop(pop_buf, dpop4, rate, dt, n);

    // pop(t+dt) = pop(t) + (dpop1 + 2*dpop2 + 2*dpop3 + dpop4)/6;
    daxpy(&n, &d_sixth, dpop1, &one, pop, &one);
    daxpy(&n, &d_third, dpop2, &one, pop, &one);
    daxpy(&n, &d_third, dpop3, &one, pop, &one);
    daxpy(&n, &d_sixth, dpop4, &one, pop, &one);

	free(pop_buf);
	free(dpop1);
	free(dpop2);
	free(dpop3);
	free(dpop4);
}

void prop_diss_RK4(double *pop, double *diss, double *rate_pop, double **rate_diss, double dt, int nsys, int nosc)
{
	int one = 1;
	double d_one   = 1.0;
	double d_half  = 1.0 / 2.0;
	double d_third = 1.0 / 3.0;
	double d_sixth = 1.0 / 6.0;

	double *pop_buf = malloc(nsys * sizeof(double));
	double *dpop1   = malloc(nsys * sizeof(double));
	double *dpop2   = malloc(nsys * sizeof(double));
	double *dpop3   = malloc(nsys * sizeof(double));
	double *dpop4   = malloc(nsys * sizeof(double));

	double *ddiss1  = malloc(nosc * sizeof(double));
	double *ddiss2  = malloc(nosc * sizeof(double));
	double *ddiss3  = malloc(nosc * sizeof(double));
	double *ddiss4  = malloc(nosc * sizeof(double));

    // Propagates the DE by 4th order Runge-Kutta method.
	// Population and dissipation are simultaneously propagated.
    // dpop(t) / dt = f(pop(t))

    // dpop1 = dt * f(pop(t))
    calc_dpop(pop, dpop1, rate_pop, dt, nsys);
    calc_ddiss(pop, ddiss1, rate_diss, dt, nsys, nosc);

    // dpop2 = dt * f(pop(t) + 0.5 * dpop1)
    dcopy(&nsys, pop, &one, pop_buf, &one);
    daxpy(&nsys, &d_half, dpop1, &one, pop_buf, &one);
    calc_dpop(pop_buf, dpop2, rate_pop, dt, nsys);
    calc_ddiss(pop_buf, ddiss2, rate_diss, dt, nsys, nosc);

    // dpop3 = dt * f(pop(t) + 0.5 * dpop2)
    dcopy(&nsys, pop, &one, pop_buf, &one);
    daxpy(&nsys, &d_half, dpop2, &one, pop_buf, &one);
    calc_dpop(pop_buf, dpop3, rate_pop, dt, nsys);
    calc_ddiss(pop_buf, ddiss3, rate_diss, dt, nsys, nosc);

    // dpop4 = dt * f(pop(t) + dpop3)
    dcopy(&nsys, pop, &one, pop_buf, &one);
    daxpy(&nsys, &d_one, dpop3, &one, pop_buf, &one);
    calc_dpop(pop_buf, dpop4, rate_pop, dt, nsys);
    calc_ddiss(pop_buf, ddiss4, rate_diss, dt, nsys, nosc);

    // pop(t+dt) = pop(t) + (dpop1 + 2*dpop2 + 2*dpop3 + dpop4)/6;
    daxpy(&nsys, &d_sixth, dpop1, &one, pop, &one);
    daxpy(&nsys, &d_third, dpop2, &one, pop, &one);
    daxpy(&nsys, &d_third, dpop3, &one, pop, &one);
    daxpy(&nsys, &d_sixth, dpop4, &one, pop, &one);

    daxpy(&nosc, &d_sixth, ddiss1, &one, diss, &one);
    daxpy(&nosc, &d_third, ddiss2, &one, diss, &one);
    daxpy(&nosc, &d_third, ddiss3, &one, diss, &one);
    daxpy(&nosc, &d_sixth, ddiss4, &one, diss, &one);

	free(pop_buf);
	free(dpop1);
	free(dpop2);
	free(dpop3);
	free(dpop4);

	free(ddiss1);
	free(ddiss2);
	free(ddiss3);
	free(ddiss4);
}
