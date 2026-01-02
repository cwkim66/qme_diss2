#include "calc_tr.h"

void comp_tr(t_qme *qme)
{
	int jtype = qme->jtype;

	init_tr(qme);

	if(jtype == 20) {
		fprintf(stderr, "---- Direct vibrational trace calculation for FRET ----\n");
		comp_tr_fret(qme->tr, 0);
		comp_tr_fret(qme->tr, 1);
	} else if(jtype == 21) {
		fprintf(stderr, "---- Direct vibrational trace calculation for MRT ----\n");
//		comp_tr_mrt(qme->tr, 0);
//		comp_tr_mrt(qme->tr, 1); // deprecate, as it is based on wrong expressions
		comp_tr_mrt_new(qme->tr);
	}
}

void comp_tr_fret(t_tr *tr, int mode)
{
	int i, j, k, l;
	char fname[256];

	int nout  = tr->nout;
	int nosc  = tr->nosc_tr;
	double dt = tr->dtout;

	for(i=0; i<nosc; i++) {
		complex **tr1 = tr->tr1[i];
		complex **tr2 = tr->tr2[i];
		t_qosc *qosc  = &tr->qosc[i];

		if(mode == 0) {
			for(j=0; j<nout; j++) {
				double t = j * dt;
				calc_tr_num_fret_pop(qosc, t, tr1[j]);
				calc_tr_ana_fret_pop(qosc, t, tr2[j]);
			}
			sprintf(fname, "tr_fret_%03d_pop.txt", i+1);
			print_tr(tr, fname, tr1, tr2);			
		} else if(mode == 1){
			for(j=0; j<nout; j++) {
				double t = j * dt;
				calc_tr_num_fret_diss(qosc, t, tr1[j]);
				calc_tr_ana_fret_diss(qosc, t, tr2[j]);
			}
			sprintf(fname, "tr_fret_%03d_diss.txt", i+1);
			print_tr(tr, fname, tr1, tr2);			
		}
	}
}

void comp_tr_mrt(t_tr *tr, int mode)
{
	int i, j;
	char fname[256];

	int nout  = tr->nout;
	int nosc  = tr->nosc_tr;
	double dt = tr->dtout;

	for(i=0; i<nosc; i++) {
		complex **tr1 = tr->tr1[i];
		complex **tr2 = tr->tr2[i];
		t_qosc *qosc  = &tr->qosc[i];

		if(mode == 0) {
			for(j=0; j<nout; j++) {
				double t = (double)j * dt;
				calc_tr_num_mrt_pop(qosc, t, tr1[j]);
				calc_tr_ana_mrt_pop(qosc, t, tr2[j]);
			}
			sprintf(fname, "tr_mrt_%03d_pop.txt", i+1);
			print_tr(tr, fname, tr1, tr2);

			for(j=0; j<nout; j++) {
				double t = (double)j * dt;
				calc_tr_num_thermal_mrt_pop(qosc, t, tr2[j]);
			}
			sprintf(fname, "tr_mrt_%03d_pop_wick.txt", i+1);
			print_tr(tr, fname, tr1, tr2);
		} else if(mode == 1) {
			for(j=0; j<nout; j++) {
				double t = (double)j * dt;
//				calc_tr_num1_mrt_dpdt(qosc, t, tr1[j]);
				calc_tr_num2_mrt_dpdt(qosc, t, tr1[j]);
//				calc_tr_ana1_mrt_dpdt(qosc, t, tr2[j]);
				calc_tr_ana2_mrt_dpdt(qosc, t, tr2[j]);
			}
			sprintf(fname, "tr_mrt_%03d_dpdt.txt", i+1);
			print_tr(tr, fname, tr1, tr2);

			for(j=0; j<nout; j++) {
				double t = (double)j * dt;
				calc_tr_num2_thermal_mrt_dpdt(qosc, t, tr2[j]);
			}
			sprintf(fname, "tr_mrt_%03d_dpdt_wick.txt", i+1);
			print_tr(tr, fname, tr1, tr2);
		}
	}
}

void print_tr(t_tr *tr, char *fname, complex **tr1, complex **tr2)
{
	int i, j, k;
	double tfac = tr->unit ? TUNIT_TO_FS : 1.0;

	int nsys  = tr->nsys;
	int nout  = tr->nout;
	double dt = tr->dtout;

	FILE *f = fopen(fname, "w");

	fprintf(f, "           Time");
	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				fprintf(f, "  Tr1 %2d-%2d Re", i+1, j+1);
				fprintf(f, "  Tr1 %2d-%2d Im", i+1, j+1);
				fprintf(f, "  Tr2 %2d-%2d Re", i+1, j+1);
				fprintf(f, "  Tr2 %2d-%2d Im", i+1, j+1);
			}
		}
	}
	fprintf(f, "\n");

	for(i=0; i<nout; i++) {
		fprintf(f, "%15.6lf ", i * dt * tfac);

		for(j=0; j<nsys; j++) {
			for(k=0; k<nsys; k++) {
				if(j != k) {
					fprintf(f, "%13.6le %13.6le ", tr1[i][j + nsys*k].real, tr1[i][j + nsys*k].imag);
					fprintf(f, "%13.6le %13.6le ", tr2[i][j + nsys*k].real, tr2[i][j + nsys*k].imag);
				}
			}
		}
		fprintf(f, "\n");
	}

	fclose(f); 
}

void calc_tr_num_fret_pop(t_qosc *qosc, double t, complex *tr_num)
{
	int i, j;

	int nsys   = qosc->nsys;
	int nqst   = qosc->nqst;
	int nqstsq = qosc->nqstsq;

	complex *zbuf = malloc(nqstsq * sizeof(complex));

	set_U_diag(qosc, t);

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				// Calculate Tr[(U_a^H) (U_b) (R_a^eq)]
				// a <-> j, b <-> i
				complex *Ua = qosc->U[j + nsys*j];
				complex *Ub = qosc->U[i + nsys*i];
				complex *Ra = qosc->Req[j + nsys*j];

				z_init(zbuf, nqstsq);
				zAxzB(zbuf, Ua, Ub, nqst, nqst, nqst, 2);
				zprod_trace(zbuf, Ra, nqst, &tr_num[i + nsys*j]);
			}
		}
	}

	free(zbuf);
}

void calc_tr_ana_fret_pop(t_qosc *qosc, double t, complex *tr_ana)
{
	int i, j;
	complex uni, g;
	complex f0_wsq;

	int nsys = qosc->nsys;

	double coth   = qosc->coth;
	double w      = qosc->freq;
	double *gamma = qosc->gamma;

	f0_wsq.real = coth * (1.0 - cos(w*t));
	f0_wsq.imag = sin(w*t) - w*t;

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				double gamma_diff = gamma[i + nsys*i] - gamma[j + nsys*j];
				double hrfac = gamma_diff * gamma_diff;

				uni.real = 0.0;
				uni.imag = -w * hrfac * t;

				g.real = hrfac * f0_wsq.real;
				g.imag = hrfac * f0_wsq.imag;

				double exp = pow(M_E, uni.real - g.real);
				double arg = uni.imag - g.imag;

				tr_ana[i + nsys*j].real = exp * cos(arg);
				tr_ana[i + nsys*j].imag = exp * sin(arg);
			}
		}
	}
}

void calc_tr_num_fret_diss(t_qosc *qosc, double t, complex *tr_num)
{
	int i, j, k;

	int nsys   = qosc->nsys;
	int nqst   = qosc->nqst;
	int nqstsq = qosc->nqstsq;

	complex *Hdiff = malloc(nqstsq * sizeof(complex));
	complex *zbuf1 = malloc(nqstsq * sizeof(complex));
	complex *zbuf2 = malloc(nqstsq * sizeof(complex));

	set_U_diag(qosc, t);

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				// Calculate Tr[(U_a^H) (H_b - H_a) (U_b) (R_a^eq)]
				// a <-> j, b <-> i
				complex *Ua = qosc->U[j + nsys*j];
				complex *Ub = qosc->U[i + nsys*i];
				complex *Ha = qosc->H[j + nsys*j];
				complex *Hb = qosc->H[i + nsys*i];
				complex *Ra = qosc->Req[j + nsys*j];

				for(k=0; k<nqstsq; k++) {
					Hdiff[k].real = Hb[k].real - Ha[k].real;
					Hdiff[k].imag = Hb[k].imag - Ha[k].imag;
				}

				z_init(zbuf1, nqstsq);
				z_init(zbuf2, nqstsq);
				zAxzB(zbuf1, Hdiff, Ub, nqst, nqst, nqst, 1);
				zAxzB(zbuf2, Ua, zbuf1, nqst, nqst, nqst, 2);
				zprod_trace(zbuf2, Ra, nqst, &tr_num[i + nsys*j]);
			}
		}
	}

	free(zbuf1);
	free(zbuf2);
}

void calc_tr_ana_fret_diss(t_qosc *qosc, double t, complex *tr_ana)
{
	int i, j;
	complex f, z;

	int nsys   = qosc->nsys;
	int nsyssq = qosc->nsyssq;

	double coth   = qosc->coth;
	double w      = qosc->freq;
	double *gamma = qosc->gamma;

	complex *tr_pop = malloc(nsyssq * sizeof(complex));
	
	f.real = cos(w*t);
	f.imag = -coth * sin(w*t);

	calc_tr_ana_fret_pop(qosc, t, tr_pop);

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				double gamma_diff = gamma[i + nsys*i] - gamma[j + nsys*j];
				double hrfac = gamma_diff * gamma_diff;
				double reorg = w * hrfac;

				z.real = f.real * reorg;
				z.imag = f.imag * reorg;

				z_mult(&tr_ana[i + nsys*j], z, tr_pop[i + nsys*j]);
			}
		}
	}

	free(tr_pop);
}

void calc_tr_num_mrt_pop(t_qosc *qosc, double t, complex *tr_num)
{
	int i, j;

	int nsys   = qosc->nsys;
	int nqst   = qosc->nqst;
	int nqstsq = qosc->nqstsq;

	complex *zbuf1 = malloc(nqstsq * sizeof(complex));
	complex *zbuf2 = malloc(nqstsq * sizeof(complex));

	set_U_diag(qosc, t);

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				// Calculate Tr[(U_a^H) (H_ab) (U_b) (H_ba) (R_a^eq)]
				// a <-> j, b <-> i
				complex *Ua  = qosc->U[j + nsys*j];
				complex *Ub  = qosc->U[i + nsys*i];
				complex *Hab = qosc->H[j + nsys*i];
				complex *Hba = qosc->H[i + nsys*j];
				complex *Ra  = qosc->Req[j + nsys*j];

				z_init(zbuf1, nqstsq);
				zAxzB(zbuf1, Ua,    Hab, nqst, nqst, nqst, 2);
				zAxzB(zbuf2, zbuf1, Ub,  nqst, nqst, nqst, 1);
				zAxzB(zbuf1, zbuf2, Hba, nqst, nqst, nqst, 1);
				zprod_trace(zbuf1, Ra, nqst, &tr_num[i + nsys*j]);
			}
		}
	}

	free(zbuf1);
	free(zbuf2);
}

void calc_tr_num_thermal_mrt_pop(t_qosc *qosc, double t, complex *tr_num)
{
	int i, j;
	int one = 1;

	int nsys   = qosc->nsys;
	int nqst   = qosc->nqst;
	int nqstsq = qosc->nqstsq;

	complex *zbuf1 = malloc(nqstsq * sizeof(complex));
	complex *zbuf2 = malloc(nqstsq * sizeof(complex));

	complex **U  = malloc(nsys * sizeof(complex *));
	complex **UH = malloc(nsys * sizeof(complex *));

	for(i=0; i<nsys; i++) {
		U[i]  = malloc(nqstsq * sizeof(complex));
		UH[i] = malloc(nqstsq * sizeof(complex));
	}

	set_U_diag_wick(qosc, t, true);
	for(i=0; i<nsys; i++) zcopy(&nqstsq, qosc->U[i + nsys*i], &one, U[i], &one);

	set_U_diag_wick(qosc, t, false);
	for(i=0; i<nsys; i++) zcopy(&nqstsq, qosc->U[i + nsys*i], &one, UH[i], &one);

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				// Calculate Tr[(U_a^H) (H_ab) (U_b) (H_ba) (R_a^eq)] with wick rotation t -> t - i * hbar * beta
				// a <-> j, b <-> i
				complex *Ua  = UH[j];
				complex *Ub  = U[i];
				complex *Hab = qosc->H[j + nsys*i];
				complex *Hba = qosc->H[i + nsys*j];
				complex *Ra  = qosc->Req[j + nsys*j];

				z_init(zbuf1, nqstsq);
				zAxzB(zbuf1, Ua,    Hab, nqst, nqst, nqst, 1);
				zAxzB(zbuf2, zbuf1, Ub,  nqst, nqst, nqst, 1);
				zAxzB(zbuf1, zbuf2, Hba, nqst, nqst, nqst, 1);
				zprod_trace(zbuf1, Ra, nqst, &tr_num[i + nsys*j]);
			}
		}
	}

	free(zbuf1);
	free(zbuf2);

	for(i=0; i<nsys; i++) {
		free(U[i]);
		free(UH[i]);
	}

	free(U);
	free(UH);
}

void calc_tr_ana_mrt_pop(t_qosc *qosc, double t, complex *tr_ana)
{
	int i, j;
	complex zpref, zexp;
	complex f0_wsq, f1_wsq, f2_wsq;

	int nsys = qosc->nsys;

	double coth   = qosc->coth;
	double w      = qosc->freq;
	double *gexci = qosc->gexci;

	f0_wsq.real = coth * (1.0 - cos(w*t));
	f0_wsq.imag = sin(w*t) - w*t;
	f1_wsq.real = w * coth * sin(w*t);
	f1_wsq.imag = w * (cos(w*t) - 1.0);
	f2_wsq.real = w * w * coth * cos(w*t);
	f2_wsq.imag = -w * w * sin(w*t);

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				double gaa    = gexci[j + nsys*j];
				double gbb    = gexci[i + nsys*i];
				double gabgba = gexci[j + nsys*i] * gexci[i + nsys*j];
				double hrfac  = (gbb-gaa) * (gbb-gaa);

				double exp = pow(M_E, -hrfac * f0_wsq.real);
				double arg = -hrfac * (w*t + f0_wsq.imag);

				zexp.real = exp * cos(arg);
				zexp.imag = exp * sin(arg);

				zpref.real = (gbb-gaa) * f1_wsq.real;
				zpref.imag = (gbb-gaa) * f1_wsq.imag - 2.0 * w * gaa;
				z_mult(&zpref, zpref, zpref);

				zpref.real = gabgba * (f2_wsq.real - zpref.real);
				zpref.imag = gabgba * (f2_wsq.imag - zpref.imag);

				z_mult(&tr_ana[i + nsys*j], zpref, zexp);
			}
		}
	}
}

void calc_tr_num1_mrt_dpdt(t_qosc *qosc, double t, complex *tr_num)
{
	int i, j;
	double dt = 1.0e-06;

	int nsys   = qosc->nsys;
	int nsyssq = qosc->nsyssq;

	complex *trm = malloc(nsyssq * sizeof(complex));
	complex *trp = malloc(nsyssq * sizeof(complex));
	
	calc_tr_ana_mrt_pop(qosc, t-dt, trm);
	calc_tr_ana_mrt_pop(qosc, t+dt, trp);

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				tr_num[i + nsys*j].real = (trp[i + nsys*j].real - trm[i + nsys*j].real) / (2.0 * dt);
				tr_num[i + nsys*j].imag = (trp[i + nsys*j].imag - trm[i + nsys*j].imag) / (2.0 * dt);
			}
		}
	}

	free(trm);
	free(trp);
}

void calc_tr_num2_mrt_dpdt(t_qosc *qosc, double t, complex *tr_num)
{
	int i, j, k;
	complex z, zmone;
	zmone.real = 0.0, zmone.imag = 1.0;

	int nsys   = qosc->nsys;
	int nqst   = qosc->nqst;
	int nqstsq = qosc->nqstsq;

	complex *zbuf1 = malloc(nqstsq * sizeof(complex));
	complex *zbuf2 = malloc(nqstsq * sizeof(complex));

	set_U_diag(qosc, t);

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				// Calculate i * Tr[(U_a^H) (H_aa H_ab - H_ab H_bb) (U_b) (H_ba) (R_a^eq)]
				// a <-> j, b <-> i
				complex *Ua  = qosc->U[j + nsys*j];
				complex *Ub  = qosc->U[i + nsys*i];
				complex *Haa = qosc->H[j + nsys*j];
				complex *Hab = qosc->H[j + nsys*i];
				complex *Hba = qosc->H[i + nsys*j];
				complex *Hbb = qosc->H[i + nsys*i];
				complex *Ra  = qosc->Req[j + nsys*j];

				z_init(zbuf1, nqstsq);
				zAxzB(zbuf1, Haa,   Hab,   nqst, nqst, nqst, 1);
				z_init(zbuf2, nqstsq);
				zAxzB(zbuf2, Hab,   Hbb,   nqst, nqst, nqst, 1);

				for(k=0; k<nqstsq; k++) {
					zbuf1[k].real -= zbuf2[k].real;
					zbuf1[k].imag -= zbuf2[k].imag;
				}

				z_init(zbuf2, nqstsq);
				zAxzB(zbuf2, Ua,    zbuf1, nqst, nqst, nqst, 2);
				zAxzB(zbuf1, zbuf2, Ub,    nqst, nqst, nqst, 1);
				zAxzB(zbuf2, zbuf1, Hba,   nqst, nqst, nqst, 1);
				zprod_trace(zbuf2, Ra, nqst, &z);

				z_mult(&tr_num[i + nsys*j], zmone, z);
			}
		}
	}

	free(zbuf1);
	free(zbuf2);
}

void calc_tr_num2_thermal_mrt_dpdt(t_qosc *qosc, double t, complex *tr_num)
{
	int i, j, k;
	int one = 1;
	complex z, zmone;
	zmone.real = 0.0, zmone.imag = 1.0;

	int nsys   = qosc->nsys;
	int nqst   = qosc->nqst;
	int nqstsq = qosc->nqstsq;

	complex *zbuf1 = malloc(nqstsq * sizeof(complex));
	complex *zbuf2 = malloc(nqstsq * sizeof(complex));

	complex **U  = malloc(nsys * sizeof(complex *));
	complex **UH = malloc(nsys * sizeof(complex *));

	for(i=0; i<nsys; i++) {
		U[i]  = malloc(nqstsq * sizeof(complex));
		UH[i] = malloc(nqstsq * sizeof(complex));
	}

	set_U_diag_wick(qosc, t, true);
	for(i=0; i<nsys; i++) zcopy(&nqstsq, qosc->U[i + nsys*i], &one, U[i], &one);

	set_U_diag_wick(qosc, t, false);
	for(i=0; i<nsys; i++) zcopy(&nqstsq, qosc->U[i + nsys*i], &one, UH[i], &one);

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				// Calculate i * Tr[(U_a^H) (H_aa H_ab - H_ab H_bb) (U_b) (H_ba) (R_a^eq)] with wick rotation t -> t - i * hbar * beta
				// a <-> j, b <-> i
				complex *Ua  = UH[j];
				complex *Ub  = U[i];
				complex *Haa = qosc->H[j + nsys*j];
				complex *Hab = qosc->H[j + nsys*i];
				complex *Hba = qosc->H[i + nsys*j];
				complex *Hbb = qosc->H[i + nsys*i];
				complex *Ra  = qosc->Req[j + nsys*j];

				z_init(zbuf1, nqstsq);
				zAxzB(zbuf1, Haa,   Hab,   nqst, nqst, nqst, 1);
				z_init(zbuf2, nqstsq);
				zAxzB(zbuf2, Hab,   Hbb,   nqst, nqst, nqst, 1);

				for(k=0; k<nqstsq; k++) {
					zbuf1[k].real -= zbuf2[k].real;
					zbuf1[k].imag -= zbuf2[k].imag;
				}

				z_init(zbuf2, nqstsq);
				zAxzB(zbuf2, Ua,    zbuf1, nqst, nqst, nqst, 1);
				zAxzB(zbuf1, zbuf2, Ub,    nqst, nqst, nqst, 1);
				zAxzB(zbuf2, zbuf1, Hba,   nqst, nqst, nqst, 1);
				zprod_trace(zbuf2, Ra, nqst, &z);

				z_mult(&tr_num[i + nsys*j], zmone, z);
			}
		}
	}

	free(zbuf1);
	free(zbuf2);

	for(i=0; i<nsys; i++) {
		free(U[i]);
		free(UH[i]);
	}

	free(U);
	free(UH);
}

void calc_tr_ana1_mrt_dpdt(t_qosc *qosc, double t, complex *tr_ana)
{
	int i, j;
	complex zpref1, zpref2, zbuf, zexp;
	complex f0_wsq, f1_wsq, f2_wsq, f3_wsq;

	int nsys = qosc->nsys;

	double coth   = qosc->coth;
	double w      = qosc->freq;
	double *gexci = qosc->gexci;

	f0_wsq.real = coth * (1.0 - cos(w*t));
	f0_wsq.imag = sin(w*t) - w*t;
	f1_wsq.real = w * coth * sin(w*t);
	f1_wsq.imag = w * (cos(w*t) - 1.0);
	f2_wsq.real = w * w * coth * cos(w*t);
	f2_wsq.imag = -w * w * sin(w*t);
	f3_wsq.real = -w * w * w * coth * sin(w*t);
	f3_wsq.imag = -w * w * w * cos(w*t);

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				double gaa    = gexci[j + nsys*j];
				double gbb    = gexci[i + nsys*i];
				double gabgba = gexci[j + nsys*i] * gexci[i + nsys*j];
				double hrfac  = (gbb-gaa) * (gbb-gaa);

				double exp = pow(M_E, -hrfac * f0_wsq.real);
				double arg = -hrfac * (w*t + f0_wsq.imag);

				zexp.real = exp * cos(arg);
				zexp.imag = exp * sin(arg);

				// 1st term
				zpref1.real = (gbb-gaa) * f1_wsq.real;
				zpref1.imag = (gbb-gaa) * f1_wsq.imag - 2.0 * w * gaa;
				z_mult(&zpref1, zpref1, zpref1);

				zpref1.real = f2_wsq.real - zpref1.real;
				zpref1.imag = f2_wsq.imag - zpref1.imag;

				zbuf.real = -hrfac * f1_wsq.real;
				zbuf.imag = -hrfac * (w + f1_wsq.imag);

				z_mult(&zpref1, zbuf, zpref1);

				// 2nd term
				zpref2.real = (gbb-gaa) * f1_wsq.real;
				zpref2.imag = (gbb-gaa) * f1_wsq.imag - 2.0 * w * gaa;

				zbuf.real = 2.0 * (gbb-gaa) * f2_wsq.real;
				zbuf.imag = 2.0 * (gbb-gaa) * f2_wsq.imag;

				z_mult(&zpref2, zbuf, zpref2);

				zpref2.real = f3_wsq.real - zpref2.real;
				zpref2.imag = f3_wsq.imag - zpref2.imag;

				// Overall
				zpref1.real = gabgba * (zpref1.real + zpref2.real);
				zpref1.imag = gabgba * (zpref1.real + zpref2.imag);

				z_mult(&tr_ana[i + nsys*j], zpref1, zexp);
			}
		}
	}
}

void calc_tr_ana2_mrt_dpdt(t_qosc *qosc, double t, complex *tr_ana)
{
	int i, j;
	complex zterm1, zterm2, zbuf, zdeno, znume;
	complex f0_wsq, f1_wsq, f2_wsq, f3_wsq;

	int nsys   = qosc->nsys;
	int nsyssq = qosc->nsyssq;

	double coth   = qosc->coth;
	double w      = qosc->freq;
	double *gexci = qosc->gexci;

	f1_wsq.real = w * coth * sin(w*t);
	f1_wsq.imag = w * (cos(w*t) - 1.0);
	f2_wsq.real = w * w * coth * cos(w*t);
	f2_wsq.imag = -w * w * sin(w*t);
	f3_wsq.real = -w * w * w * coth * sin(w*t);
	f3_wsq.imag = -w * w * w * cos(w*t);

	complex *tr_pop = malloc(nsyssq * sizeof(complex));

	calc_tr_ana_mrt_pop(qosc, t, tr_pop);

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				double gaa    = gexci[j + nsys*j];
				double gbb    = gexci[i + nsys*i];
				double gabgba = gexci[j + nsys*i] * gexci[i + nsys*j];
				double hrfac  = (gbb-gaa) * (gbb-gaa);

				// term 1
				zterm1.real = -hrfac * f1_wsq.real;
				zterm1.imag = -hrfac * (w + f1_wsq.imag);

				// term 2 numerator
				znume.real = (gbb-gaa) * f1_wsq.real;
				znume.imag = (gbb-gaa) * f1_wsq.imag - 2.0 * w * gaa;

				zbuf.real = 2.0 * (gbb-gaa) * f2_wsq.real;
				zbuf.imag = 2.0 * (gbb-gaa) * f2_wsq.imag;

				z_mult(&znume, zbuf, znume);

				znume.real = f3_wsq.real - znume.real;
				znume.imag = f3_wsq.imag - znume.imag;

				// term 2 denominator
				zdeno.real = (gbb-gaa) * f1_wsq.real;
				zdeno.imag = (gbb-gaa) * f1_wsq.imag - 2.0 * w * gaa;

				z_mult(&zdeno, zdeno, zdeno);

				zdeno.real = f2_wsq.real - zdeno.real;
				zdeno.imag = f2_wsq.imag - zdeno.imag;

				// Finalize
				z_div(&zterm2, znume, zdeno);

				zterm1.real += zterm2.real;
				zterm1.imag += zterm2.imag;

				z_mult(&tr_ana[i + nsys*j], zterm1, tr_pop[i + nsys*j]);

// vvvvvvvv
#if 1
				if((t < (100.0 / TUNIT_TO_FS)) && (i == 0)) {
					if(t < 1.0e-10) fprintf(stderr, "gaa %13.6le   gbb %13.6le   gabgba %13.6le   hrfac %13.6le\n", gaa, gbb, gabgba, hrfac);
					fprintf(stderr, "Gt %13.6le %13.6le   zdeno %13.6le %13.6le   znume %13.6le %13.6le   t %13.6le %13.6le\n",
					tr_pop[i + nsys*j].real, tr_pop[i + nsys*j].imag, zdeno.real, zdeno.imag, znume.real, znume.imag, tr_ana[i + nsys*j].real, tr_ana[i + nsys*j].imag);
#endif
// ^^^^^^^^
				}
			}
		}
	}

	free(tr_pop);
}
