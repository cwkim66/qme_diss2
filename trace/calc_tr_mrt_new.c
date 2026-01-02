#include "calc_tr_mrt_new.h"

void comp_tr_mrt_new(t_tr *tr)
{
	int i, j, k, l;
	char fname[256];

	int nsyssq = tr->nsyssq;
	int nout   = tr->nout;
	int nosc   = tr->nosc_tr;
	double dt  = tr->dtout;

	// Allocate
	complex **zt1 = malloc(nout * sizeof(complex *));
	complex **zt2 = malloc(nout * sizeof(complex *));

	for(i=0; i<nout; i++) {
		zt1[i] = malloc(nsyssq * sizeof(complex));
		zt2[i] = malloc(nsyssq * sizeof(complex));
	}

	complex ***mtr = malloc(4 * sizeof(complex **));

	for(i=0; i<4; i++) {
		mtr[i] = malloc(nosc * sizeof(complex *));

		for(j=0; j<nosc; j++) {
			mtr[i][j] = malloc(nsyssq * sizeof(complex));
		}
	}

	complex ***Ct = malloc(4 * sizeof(complex **));

	for(i=0; i<4; i++) {
		Ct[i] = malloc(nout * sizeof(complex *));

		for(j=0; j<nout; j++) {
			Ct[i][j] = malloc(nsyssq * sizeof(complex));
		}
	}

	if(nosc < 3) {
		fprintf(stderr, "Error: to execute detailed analysis on mrt,\nthe number of oscillators must be larger than 3 (currently %d)\n", nosc);
		exit(EXIT_FAILURE);
	}

	// The last mode will be treated explicitly
	for(i=0; i<nout; i++) {
		double t = i * dt;
		calc_tr_ana_mrt_mult(tr, t, zt1[i]);
		calc_mrt_micro_tr_num(tr, t, mtr);
		calc_mrt_micro_tr_ana(tr, t, mtr);

// Expression 1
#if 0
		calc_tr_micro_mrt_mult(tr, t, mtr, zt2[i], 0);
#endif

// Expression 2-1
#if 1
		calc_Ct(tr, t, mtr, Ct, i);
		calc_tr_ana_from_Ct(tr, Ct, mtr, i, zt2[i], 0);
#endif

// Expression 2-2
#if 0
		calc_Ct(tr, t, mtr, Ct, i);
		calc_tr_ana_from_Ct(tr, Ct, mtr, i, zt2[i], 1);
#endif
	}

	sprintf(fname, "tr_mrt_pop.txt");
	print_tr(tr, fname, zt1, zt2);

	// Time derivatives for the last mode
	for(i=0; i<4; i++) {
		for(j=0; j<nout; j++) {
			double t = j * dt;
//			calc_mrt_micro_tr_deriv_num1(tr, t, mtr, zt1[j], i);
			calc_mrt_micro_tr_deriv_num2(tr, t, zt2[j], i);
			calc_mrt_micro_tr_deriv_ana(tr, t, zt1[j], i);
		}
		sprintf(fname, "tr_mrt_deriv%d.txt", i);
		print_tr(tr, fname, zt1, zt2);
	}

	// Energy kernel
	for(i=0; i<nout; i++) {
		double t = i * dt;
		calc_mrt_micro_tr_ana(tr, t, mtr);
		calc_Ct(tr, t, mtr, Ct, i);
		calc_eker_from_Ct(tr, Ct, i, zt1[i]);
		calc_eker_ana(tr, t, zt2[i]);
	}

	sprintf(fname, "tr_mrt_eker.txt");
	print_tr(tr, fname, zt1, zt2);

	for(i=0; i<nout; i++) {
		free(zt1[i]);
		free(zt2[i]);
	}

	free(zt1);
	free(zt2);

	for(i=0; i<4; i++) {
		for(j=0; j<nosc; j++) {
			free(mtr[i][j]);
		}
		free(mtr[i]);
	}
	free(mtr);

	for(i=0; i<4; i++) {
		for(j=0; j<nout; j++) {
			free(Ct[i][j]);
		}
		free(Ct[i]);
	}
	free(Ct);
}

void calc_tr_ana_mrt_mult(t_tr *tr, double t, complex *tr_ana)
{
	int i, j, k;
	complex zpref1, zpref2, zpref3;
	complex zexp, zfac;
	complex f0, f1, f2;

	int nsys = tr->nsys;
	int nosc = tr->nosc_tr;

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				zpref1.real = 0.0, zpref1.imag = 0.0;
				zpref2.real = 0.0, zpref2.imag = 0.0;
				zpref3.real = 0.0, zpref3.imag = 0.0;
				zexp.real   = 0.0, zexp.imag   = 0.0;

				for(k=0; k<nosc; k++) {
					t_qosc *qosc = &tr->qosc[k];
			
					double coth   = qosc->coth;
					double w      = qosc->freq;
					double wsq    = w*w;
					double *gexci = qosc->gexci;

					f0.real = coth * (1.0 - cos(w*t)) / wsq;
					f0.imag = (sin(w*t) - w*t) / wsq;
					f1.real = coth * sin(w*t) / w;
					f1.imag = (cos(w*t) - 1.0) / w;
					f2.real = coth * cos(w*t);
					f2.imag = -sin(w*t);

					double gaa   = w * gexci[j + nsys*j];
					double gbb   = w * gexci[i + nsys*i];
					double gab   = w * gexci[j + nsys*i];
					double gba   = w * gexci[i + nsys*j];
					double Gba   = gbb - gaa;
					double Gbasq = Gba * Gba;

					// zpref1
					zpref1.real += gab * gba * f2.real;
					zpref1.imag += gab * gba * f2.imag;

					// zpref2
					zpref2.real += gab * Gba * f1.real;
					zpref2.imag += gab * (Gba * f1.imag - 2.0 * gaa / w);

					// zpref3
					zpref3.real += gba * Gba * f1.real;
					zpref3.imag += gba * (Gba * f1.imag - 2.0 * gaa / w);

					// zexp
					zexp.real -= Gbasq * f0.real;
					zexp.imag -= Gbasq * (t/w + f0.imag);
				}

				z_mult(&zpref2, zpref2, zpref3);

				zpref1.real -= zpref2.real;
				zpref1.imag -= zpref2.imag;

				double exp = pow(M_E, zexp.real);
				double arg = zexp.imag;
				zfac.real = exp * cos(arg);
				zfac.imag = exp * sin(arg);

				z_mult(&tr_ana[i + nsys*j], zpref1, zfac);
			}
		}
	}
}

void calc_mrt_micro_tr_num(t_tr *tr, double t, complex ***mtr)
{
	int i, j, k;
	int one = 1;

	int nsys = tr->nsys;
	int nosc = tr->nosc_tr;

	complex **tr0  = mtr[0];
	complex **tr1  = mtr[1];
	complex **tr1p = mtr[2];
	complex **tr2  = mtr[3];

	for(i=0; i<nosc; i++) {
		t_qosc *qosc = &tr->qosc[i];

		int nqst   = qosc->nqst;
		int nqstsq = qosc->nqstsq;

		complex *aHpa  = malloc(nqstsq * sizeof(complex));
		complex *zbuf1 = malloc(nqstsq * sizeof(complex));
		complex *zbuf2 = malloc(nqstsq * sizeof(complex));

		set_U_diag(qosc, t);

		for(j=0; j<nsys; j++) {
			for(k=0; k<nsys; k++) {
				if(j != k) {
					complex *FCU = qosc->FCU[k + nsys*k];
					complex *Ua  = qosc->U[k + nsys*k];
					complex *Ub  = qosc->U[j + nsys*j];
					complex *Ra  = qosc->Req[k + nsys*k];

					// Transformation to ground state basis
					zcopy(&nqstsq, qosc->aHpa, &one, aHpa, &one);
					ztrans_inv(FCU, aHpa, nqst);

					// tr0: Tr[(U_a^H) (U_b) (R_a^eq)]
					z_init(zbuf1, nqstsq);
					zAxzB(zbuf1, Ua, Ub, nqst, nqst, nqst, 2);
					zprod_trace(zbuf1, Ra, nqst, &tr0[i][j + nsys*k]);

					// tr1: Tr[(U_a^H) (a^H + a) (U_b) (R_a^eq)]
					z_init(zbuf1, nqstsq);
					z_init(zbuf2, nqstsq);
					zAxzB(zbuf1, Ua, aHpa, nqst, nqst, nqst, 2);
					zAxzB(zbuf2, zbuf1, Ub, nqst, nqst, nqst, 1);
					zprod_trace(zbuf2, Ra, nqst, &tr1[i][j + nsys*k]);

					// tr1p: Tr[(U_a^H) (U_b) (a^H + a) (R_a^eq)]
					z_init(zbuf1, nqstsq);
					z_init(zbuf2, nqstsq);
					zAxzB(zbuf1, Ua, Ub, nqst, nqst, nqst, 2);
					zAxzB(zbuf2, zbuf1, aHpa, nqst, nqst, nqst, 1);
					zprod_trace(zbuf2, Ra, nqst, &tr1p[i][j + nsys*k]);

					// tr2: Tr[(U_a^H) (a^H + a) (U_b) (a^H + a) (R_a^eq)]
					z_init(zbuf1, nqstsq);
					z_init(zbuf2, nqstsq);
					zAxzB(zbuf1, Ua, aHpa, nqst, nqst, nqst, 2);
					zAxzB(zbuf2, zbuf1, Ub, nqst, nqst, nqst, 1);
					z_init(zbuf1, nqstsq);
					zAxzB(zbuf1, zbuf2, aHpa, nqst, nqst, nqst, 1);
					zprod_trace(zbuf1, Ra, nqst, &tr2[i][j + nsys*k]);
				}
			}
		}

		free(aHpa);
		free(zbuf1);
		free(zbuf2);
	}
}

void calc_mrt_micro_tr_ana(t_tr *tr, double t, complex ***mtr)
{
	int i, j, k;
	complex zbuf1, zbuf2, zexp;
	complex f0, f1, f2;

	int nsys = tr->nsys;
	int nosc = tr->nosc_tr;

	complex **tr0  = mtr[0];
	complex **tr1  = mtr[1];
	complex **tr1p = mtr[2];
	complex **tr2  = mtr[3];

	for(i=0; i<nosc; i++) {
		t_qosc *qosc = &tr->qosc[i];

		double coth   = qosc->coth;
		double w      = qosc->freq;
		double wsq    = w * w;
		double *gexci = qosc->gexci;

		f0.real = coth * (1.0 - cos(w*t)) / wsq;
		f0.imag = (sin(w*t) - w*t) / wsq;
		f1.real = coth * sin(w*t) / w;
		f1.imag = (cos(w*t) - 1.0) / w;
		f2.real = coth * cos(w*t);
		f2.imag = -sin(w*t);

		for(j=0; j<nsys; j++) {
			for(k=0; k<nsys; k++) {
				if(j != k) {
					double gaa   = w * gexci[k + nsys*k];
					double gbb   = w * gexci[j + nsys*j];
					double Gba   = gbb-gaa;
					double Gbasq = Gba * Gba;

					// tr0
					zexp.real = -Gbasq * f0.real;
					zexp.imag = -Gbasq * (t/w + f0.imag);
					double exp = pow(M_E, zexp.real);
					double arg = zexp.imag;
					tr0[i][j + nsys*k].real = exp * cos(arg);
					tr0[i][j + nsys*k].imag = exp * sin(arg);

					// tr1
					zbuf1.real = Gba * f1.real;
					zbuf1.imag = Gba * f1.imag;
					z_mult(&zbuf2, zbuf1, tr0[i][j + nsys*k]);

					zbuf1.real = 0.0;
					zbuf1.imag = -1.0;
					z_mult(&tr1[i][j + nsys*k], zbuf2, zbuf1);

					// tr1p
					tr1p[i][j + nsys*k].real = tr1[i][j + nsys*k].real;
					tr1p[i][j + nsys*k].imag = tr1[i][j + nsys*k].imag;

					// tr2
					z_mult(&zbuf1, f1, f1);
					zbuf2.real = f2.real - Gbasq * zbuf1.real;
					zbuf2.imag = f2.imag - Gbasq * zbuf1.imag;
					z_mult(&tr2[i][j + nsys*k], zbuf2, tr0[i][j + nsys*k]);
				}
			}
		}
	}
}

void calc_tr_micro_mrt_mult(t_tr *tr, double t, complex ***mtr, complex *tr_num, int mode)
{
	int i, j, k, l, m;
	int fin;
	complex zbuf1, zbuf2;

	int nsys   = tr->nsys;
	int nsyssq = tr->nsyssq;
	int nosc   = tr->nosc_tr;

	complex **tr0  = mtr[0];
	complex **tr1  = mtr[1];
	complex **tr1p = mtr[2];
	complex **tr2  = mtr[3];

	z_init(tr_num, nsyssq);

	switch(mode) {
		case 0: fin = nosc;
				break;
		case 1: fin = nosc-1;
				break;
		default: fprintf(stderr, "Error at %s\n", __func__);
				exit(EXIT_FAILURE);
				break;
	}

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				int naa = j + nsys*j;
				int nbb = i + nsys*i;
				int nab = j + nsys*i;
				int nba = i + nsys*j;

				// We add reorganization E for all modes
				double rebaaa = 0.0, reabaa = 0.0;

				for(k=0; k<nosc; k++) {
					t_qosc *qosc = &tr->qosc[k];

					double w      = qosc->freq;
					double *gexci = qosc->gexci;

					rebaaa += w * gexci[nba] * gexci[naa];
					reabaa += w * gexci[nab] * gexci[naa];
				}

				// Term 1
				for(k=0; k<fin; k++) {
					t_qosc *qosck = &tr->qosc[k];

					double wk  = qosck->freq;
					double *gk = qosck->gexci;

					double fac = wk * wk * gk[nab] * gk[nba];

					zbuf1.real = tr2[k][nba].real;
					zbuf1.imag = tr2[k][nba].imag;

					for(l=0; l<fin; l++) {
						if(l != k) {
							z_mult(&zbuf1, zbuf1, tr0[l][nba]);
						}
					}

					tr_num[nba].real += fac * zbuf1.real;
					tr_num[nba].imag += fac * zbuf1.imag;
				}

				// Term 2
				for(k=0; k<fin; k++) {
					for(l=0; l<fin; l++) {
						if(k != l) {
							t_qosc *qosck = &tr->qosc[k];
							t_qosc *qoscl = &tr->qosc[l];

							double wk = qosck->freq;
							double wl = qoscl->freq;

							double *gk = qosck->gexci;
							double *gl = qoscl->gexci;

							double fac = wk * wl * gk[nab] * gl[nba];

							zbuf1.real = tr1[k][nba].real;
							zbuf1.imag = tr1[k][nba].imag;

							z_mult(&zbuf1, zbuf1, tr1p[l][nba]);

							for(m=0; m<fin; m++) {
								if((m != k) && (m != l)) {
									z_mult(&zbuf1, zbuf1, tr0[m][nba]);
								}
							}

							tr_num[nba].real += fac * zbuf1.real;
							tr_num[nba].imag += fac * zbuf1.imag;
						}
					}
				}

				// Term 3
				for(k=0; k<fin; k++) {
					t_qosc *qosck = &tr->qosc[k];

					double wk   = qosck->freq;
					double *gk  = qosck->gexci;

					double fac = -2.0 * rebaaa * wk * gk[nab];

					zbuf1.real = tr1[k][nba].real;
					zbuf1.imag = tr1[k][nba].imag;
					
					for(l=0; l<fin; l++) {
						if(l != k) {
							z_mult(&zbuf1, zbuf1, tr0[l][nba]);
						}
					}

					tr_num[nba].real += fac * zbuf1.real;
					tr_num[nba].imag += fac * zbuf1.imag;
				}

				// Term 4
				for(k=0; k<fin; k++) {
					t_qosc *qosck = &tr->qosc[k];

					double wk   = qosck->freq;
					double *gk  = qosck->gexci;

					double fac = -2.0 * reabaa * wk * gk[nba];

					zbuf1.real = tr1p[k][nba].real;
					zbuf1.imag = tr1p[k][nba].imag;
					
					for(l=0; l<fin; l++) {
						if(l != k) {
							z_mult(&zbuf1, zbuf1, tr0[l][nba]);
						}
					}

					tr_num[nba].real += fac * zbuf1.real;
					tr_num[nba].imag += fac * zbuf1.imag;
				}

				// Term 5
				double fac = 4.0 * reabaa * rebaaa;
				zbuf1.real = 1.0;
				zbuf1.imag = 0.0;

				for(k=0; k<fin; k++) {
					z_mult(&zbuf1, zbuf1, tr0[k][nba]);
				}

				tr_num[nba].real += fac * zbuf1.real;
				tr_num[nba].imag += fac * zbuf1.imag;
			}
		}
	}
}

void calc_Ct(t_tr *tr, double t, complex ***mtr, complex ***Ct, int n)
{
	int i, j, k, l, m;
	double fac;
	complex zbuf1, zbuf2;

	int nsys   = tr->nsys;
	int nsyssq = tr->nsyssq;
	int nosc   = tr->nosc_tr;

	complex **tr0  = mtr[0];
	complex **tr1  = mtr[1];
	complex **tr1p = mtr[2];
	complex **tr2  = mtr[3];

	complex *Ct1 = Ct[0][n];
	complex *Ct2 = Ct[1][n];
	complex *Ct3 = Ct[2][n];
	complex *Ct4 = Ct[3][n];

	z_init(Ct1, nsyssq);
	z_init(Ct2, nsyssq);
	z_init(Ct3, nsyssq);
	z_init(Ct4, nsyssq);

	// Ct1
	calc_tr_micro_mrt_mult(tr, t, mtr, Ct1, 1);

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				int naa = j + nsys*j;
				int nbb = i + nsys*i;
				int nab = j + nsys*i;
				int nba = i + nsys*j;

				// We add reorganization E for all modes
				double rebaaa = 0.0, reabaa = 0.0;

				for(k=0; k<nosc; k++) {
					t_qosc *qosc = &tr->qosc[k];

					double w      = qosc->freq;
					double *gexci = qosc->gexci;

					rebaaa += w * gexci[nba] * gexci[naa];
					reabaa += w * gexci[nab] * gexci[naa];
				}

				// We denote the last mode as x
				t_qosc *qoscx = &tr->qosc[nosc-1];
				double wx  = qoscx->freq;
				double *gx = qoscx->gexci;

				// C2t term 1
				for(k=0; k<nosc-1; k++) {
					t_qosc *qosck = &tr->qosc[k];

					double wk  = qosck->freq;
					double *gk = qosck->gexci;

					fac = wk * gk[nab] * wx * gx[nba];

					zbuf1.real = tr1[k][nba].real;
					zbuf1.imag = tr1[k][nba].imag;
					
					for(l=0; l<nosc-1; l++) {
						if(l != k) {
							z_mult(&zbuf1, zbuf1, tr0[l][nba]);
						}
					}

					Ct2[nba].real += fac * zbuf1.real;
					Ct2[nba].imag += fac * zbuf1.imag;
				}

				// C2t term 2
				fac = -2.0 * reabaa * wx * gx[nba];
				zbuf1.real = 1.0;
				zbuf1.imag = 0.0;

				for(k=0; k<nosc-1; k++) {
					z_mult(&zbuf1, zbuf1, tr0[k][nba]);
				}

				Ct2[nba].real += fac * zbuf1.real;
				Ct2[nba].imag += fac * zbuf1.imag;

				// C3t term 1
				for(k=0; k<nosc-1; k++) {
					t_qosc *qosck = &tr->qosc[k];
					t_qosc *qoscx = &tr->qosc[nosc-1];

					double wk  = qosck->freq;
					double *gk = qosck->gexci;

					fac = wx * gx[nab] * wk * gk[nba];

					zbuf1.real = tr1p[k][nba].real;
					zbuf1.imag = tr1p[k][nba].imag;
					
					for(l=0; l<nosc-1; l++) {
						if(l != k) {
							z_mult(&zbuf1, zbuf1, tr0[l][nba]);
						}
					}

					Ct3[nba].real += fac * zbuf1.real;
					Ct3[nba].imag += fac * zbuf1.imag;
				}

				// C3t term 2
				fac = -2.0 * wx * gx[nab] * rebaaa;
				zbuf1.real = 1.0;
				zbuf1.imag = 0.0;

				for(k=0; k<nosc-1; k++) {
					z_mult(&zbuf1, zbuf1, tr0[k][nba]);
				}

				Ct3[nba].real += fac * zbuf1.real;
				Ct3[nba].imag += fac * zbuf1.imag;

				// C4t
				fac = wx * wx * gx[nab] * gx[nba];
				zbuf1.real = 1.0;
				zbuf1.imag = 0.0;

				for(k=0; k<nosc-1; k++) {
					z_mult(&zbuf1, zbuf1, tr0[k][nba]);
				}

				Ct4[nba].real += fac * zbuf1.real;
				Ct4[nba].imag += fac * zbuf1.imag;
			}
		}
	}
}

void calc_tr_ana_from_Ct(t_tr *tr, complex ***Ct, complex ***mtr, int n, complex *tr_ana, int mode)
{
	int i, j, k, l, m;
	double fac;
	complex zbuf;

	int nsys   = tr->nsys;
	int nsyssq = tr->nsyssq;
	int nosc   = tr->nosc_tr;

	z_init(tr_ana, nsyssq);

	if(mode == 0) {
		for(i=0; i<nsys; i++) {
			for(j=0; j<nsys; j++) {
				if(i != j) {
					for(k=0; k<4; k++) {
						z_mult(&zbuf, Ct[k][n][i + nsys*j], mtr[k][nosc-1][i + nsys*j]);
						tr_ana[i + nsys*j].real += zbuf.real;
						tr_ana[i + nsys*j].imag += zbuf.imag;

// vvvvv
						if((j == 0) && (n == 0)) {
							fprintf(stderr, "Ct  %13.6le %13.6le  ", Ct[k][0][2].real, Ct[k][0][2].imag);
							fprintf(stderr, "mtr %13.6le %13.6le  \n", mtr[k][nosc-1][2].real, mtr[k][nosc-1][2].imag);
						}
// ^^^^^
					}
				}
			}
		}
	} else if(mode == 1) {
		for(i=0; i<nsys; i++) {
			for(j=0; j<nsys; j++) {
				if(i != j) {
					for(k=0; k<4; k++) {
						z_conj(&zbuf, mtr[k][nosc-1][i + nsys*j]);
						z_multc(&zbuf, Ct[k][n][i + nsys*j], zbuf);
						tr_ana[i + nsys*j].real += zbuf.real;
						tr_ana[i + nsys*j].imag += zbuf.imag;
					}
					z_conj(&tr_ana[i + nsys*j], tr_ana[i + nsys*j]);
				}
			}
		}
	}
}

void calc_mrt_micro_tr_deriv_num1(t_tr *tr, double t, complex ***mtr, complex *tr_num, int mode)
{
	if((mode < 0) || (mode >= 4)) exit(EXIT_FAILURE);

	int i;
	int one = 1;
	double dt = 1.0e-06;
	complex zmone;
	zmone.real = 0.0, zmone.imag = -1.0;

	int nsyssq = tr->nsyssq;
	int nosc   = tr->nosc_tr;

	complex *mtrp = malloc(nsyssq * sizeof(complex));
	complex *mtrm = malloc(nsyssq * sizeof(complex));

	calc_mrt_micro_tr_ana(tr, t+dt, mtr);
	zcopy(&nsyssq, mtr[mode][nosc-1], &one, mtrp, &one);
	
	calc_mrt_micro_tr_ana(tr, t-dt, mtr);
	zcopy(&nsyssq, mtr[mode][nosc-1], &one, mtrm, &one);

	for(i=0; i<nsyssq; i++) {
		tr_num[i].real = (mtrp[i].real - mtrm[i].real) / (2.0 * dt);
		tr_num[i].imag = (mtrp[i].imag - mtrm[i].imag) / (2.0 * dt);
		z_mult(&tr_num[i], zmone, tr_num[i]);
	}

	free(mtrp);
	free(mtrm);
}

void calc_mrt_micro_tr_deriv_num2(t_tr *tr, double t, complex *tr_num, int mode)
{
	if((mode < 0) || (mode >= 4)) exit(EXIT_FAILURE);

	int i, j, k;
	int one = 1;

	int nsys = tr->nsys;
	int nosc = tr->nosc_tr;

	t_qosc *qosc = &tr->qosc[nosc-1];

	int nqst   = qosc->nqst;
	int nqstsq = qosc->nqstsq;

	complex *aHpa   = malloc(nqstsq * sizeof(complex));
	complex *Hdiff  = malloc(nqstsq * sizeof(complex));
	complex *Hdiffp = malloc(nqstsq * sizeof(complex));
	complex *zbuf1  = malloc(nqstsq * sizeof(complex));
	complex *zbuf2  = malloc(nqstsq * sizeof(complex));

	set_U_diag(qosc, t);

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				complex *FCU = qosc->FCU[j + nsys*j];
				complex *Ua  = qosc->U[j + nsys*j];
				complex *Ub  = qosc->U[i + nsys*i];
				complex *Ha  = qosc->H[j + nsys*j];
				complex *Hb  = qosc->H[i + nsys*i];
				complex *Ra  = qosc->Req[j + nsys*j];

				// Transformation to ground state basis
				zcopy(&nqstsq, qosc->aHpa, &one, aHpa, &one);
				ztrans_inv(FCU, aHpa, nqst);

				z_init(zbuf1, nqstsq);
				z_init(zbuf2, nqstsq);
				zAxzB(zbuf1, Ha, aHpa, nqst, nqst, nqst, 1);
				zAxzB(zbuf2, aHpa, Hb, nqst, nqst, nqst, 1);

				for(k=0; k<nqstsq; k++) {
					Hdiff[k].real = Ha[k].real - Hb[k].real;
					Hdiff[k].imag = Ha[k].imag - Hb[k].imag;

					Hdiffp[k].real = zbuf1[k].real - zbuf2[k].real;
					Hdiffp[k].imag = zbuf1[k].imag - zbuf2[k].imag;
				}

				switch(mode) {
					case 0:
					// dtr0/dt: Tr[(U_a^H) (H_a - H_b) (U_b) (R_a^eq)]
					z_init(zbuf1, nqstsq);
					z_init(zbuf2, nqstsq);
					zAxzB(zbuf1, Ua, Hdiff, nqst, nqst, nqst, 2);
					zAxzB(zbuf2, zbuf1, Ub, nqst, nqst, nqst, 1);
					zprod_trace(zbuf2, Ra, nqst, &tr_num[i + nsys*j]);
					break;

					case 1:
					// tr1: Tr[(U_a^H) [H_a(a^H + a) - (a^H + a)H_b] (U_b) (R_a^eq)]
					z_init(zbuf1, nqstsq);
					z_init(zbuf2, nqstsq);
					zAxzB(zbuf1, Ua, Hdiffp, nqst, nqst, nqst, 2);
					zAxzB(zbuf2, zbuf1, Ub, nqst, nqst, nqst, 1);
					zprod_trace(zbuf2, Ra, nqst, &tr_num[i + nsys*j]);
					break;

					case 2:
					// tr1p: Tr[(U_a^H) (H_a - H_b) (U_b) (a^H + a) (R_a^eq)]
					z_init(zbuf1, nqstsq);
					z_init(zbuf2, nqstsq);
					zAxzB(zbuf1, Ua, Hdiff, nqst, nqst, nqst, 2);
					zAxzB(zbuf2, zbuf1, Ub, nqst, nqst, nqst, 1);
					zAxzB(zbuf1, zbuf2, aHpa, nqst, nqst, nqst, 1);
					zprod_trace(zbuf1, Ra, nqst, &tr_num[i + nsys*j]);
					break;

					case 3:
					// tr2: Tr[(U_a^H) [H_a(a^H + a) - (a^H + a)H_b] (U_b) (a^H + a) (R_a^eq)]
					z_init(zbuf1, nqstsq);
					z_init(zbuf2, nqstsq);
					zAxzB(zbuf1, Ua, Hdiffp, nqst, nqst, nqst, 2);
					zAxzB(zbuf2, zbuf1, Ub, nqst, nqst, nqst, 1);
					z_init(zbuf1, nqstsq);
					zAxzB(zbuf1, zbuf2, aHpa, nqst, nqst, nqst, 1);
					zprod_trace(zbuf1, Ra, nqst, &tr_num[i + nsys*j]);
					break;

					default: break;
				}
			}
		}
	}

	free(aHpa);
	free(Hdiff);
	free(Hdiffp);
	free(zbuf1);
	free(zbuf2);
}

void calc_mrt_micro_tr_deriv_ana(t_tr *tr, double t, complex *tr_ana, int mode)
{
	int i, j, k;
	complex zbuf1, zbuf2, zbuf3, zexp;
	complex f0, f1, f2, f3;
	complex zimag, zmimag;
	zimag.real  = 0.0, zimag.imag  =  1.0;
	zmimag.real = 0.0, zmimag.imag = -1.0;

	int nsys = tr->nsys;
	int nosc = tr->nosc_tr;

	t_qosc *qosc = &tr->qosc[nosc-1];

	double coth   = qosc->coth;
	double w      = qosc->freq;
	double wsq    = w * w;
	double *gexci = qosc->gexci;

	f0.real = coth * (1.0 - cos(w*t)) / wsq;
	f0.imag = (sin(w*t) - w*t) / wsq;
	f1.real = coth * sin(w*t) / w;
	f1.imag = (cos(w*t) - 1.0) / w;
	f2.real = coth * cos(w*t);
	f2.imag = -sin(w*t);
	f3.real = -w * coth * sin(w*t);
	f3.imag = -w * cos(w*t);

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				double gaa   = w * gexci[j + nsys*j];
				double gbb   = w * gexci[i + nsys*i];
				double Gba   = gbb-gaa;
				double Gbasq = Gba * Gba;
				double Gba3  = Gbasq * Gba;
				double Gba4  = Gba3 * Gba;

				zbuf1.real = -Gbasq * f0.real;
				zbuf1.imag = -Gbasq * (t/w + f0.imag);
				double exp = pow(M_E, zbuf1.real);
				double arg = zbuf1.imag;
				zexp.real  = exp * cos(arg);
				zexp.imag  = exp * sin(arg);
				
				switch(mode) {
					case 0:
					zbuf1.real = Gbasq * f1.real;
					zbuf1.imag = Gbasq * f1.imag;
					z_mult(&zbuf1, zimag, zbuf1);
					zbuf1.real -= Gbasq / w;

					z_mult(&tr_ana[i + nsys*j], zbuf1, zexp);
					break;

					case 1:
					z_mult(&zbuf1, f1, f1);
					zbuf1.real *= Gba3;
					zbuf1.imag *= Gba3;

					zbuf2.real = Gba3 / w * f1.real;
					zbuf2.imag = Gba3 / w * f1.imag;
					z_mult(&zbuf2, zimag, zbuf2);
					zbuf1.real += zbuf2.real - Gba * f2.real;
					zbuf1.imag += zbuf2.imag - Gba * f2.imag;

					z_mult(&tr_ana[i + nsys*j], zbuf1, zexp);
					break;

					case 2:
					z_mult(&zbuf1, f1, f1);
					zbuf1.real *= Gba3;
					zbuf1.imag *= Gba3;

					zbuf2.real = Gba3 / w * f1.real;
					zbuf2.imag = Gba3 / w * f1.imag;
					z_mult(&zbuf2, zimag, zbuf2);
					zbuf1.real += zbuf2.real - Gba * f2.real;
					zbuf1.imag += zbuf2.imag - Gba * f2.imag;

					z_mult(&tr_ana[i + nsys*j], zbuf1, zexp);
					break;
					
					case 3:
					z_mult(&zbuf1, f1, f1);
					z_mult(&zbuf1, zbuf1, f1);
					zbuf1.real *= Gba4;
					zbuf1.imag *= Gba4;

					z_mult(&zbuf3, f1, f2);
					zbuf1.real += f3.real - 3.0 * Gbasq * zbuf3.real;
					zbuf1.imag += f3.imag - 3.0 * Gbasq * zbuf3.imag;

					z_mult(&zbuf1, zmimag, zbuf1);

					z_mult(&zbuf2, f1, f1);
					zbuf2.real *= -Gbasq;
					zbuf2.imag *= -Gbasq;
					zbuf2.real += f2.real;
					zbuf2.imag += f2.imag;

					zbuf1.real -= Gbasq / w * zbuf2.real;
					zbuf1.imag -= Gbasq / w * zbuf2.imag;

					z_mult(&tr_ana[i + nsys*j], zbuf1, zexp);
					break;

					default: break;
				}
			}
		}
	}
}

void calc_eker_from_Ct(t_tr *tr, complex ***Ct, int n, complex *tr_num)
{
	int i, j, k;
	double fac;
	complex zbuf;

	int nsys   = tr->nsys;
	int nsyssq = tr->nsyssq;
	int nosc   = tr->nosc_tr;
	double t   = tr->dtout * n;

	complex *mtr = malloc(nsyssq * sizeof(complex));

	z_init(tr_num, nsyssq);

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				for(k=0; k<4; k++) {
					calc_mrt_micro_tr_deriv_ana(tr, t, mtr, k);
					z_mult(&zbuf, Ct[k][n][i + nsys*j], mtr[i + nsys*j]);
					tr_num[i + nsys*j].real += zbuf.real;
					tr_num[i + nsys*j].imag += zbuf.imag;

// vvvvv
					if((j == 0) && (n == 0)) {
						fprintf(stderr, "Ct  %13.6le %13.6le  ", Ct[k][0][2].real, Ct[k][0][2].imag);
						fprintf(stderr, "mtr %13.6le %13.6le  \n", mtr[2].real, mtr[2].imag);
					}
// ^^^^^
				}
			}
		}
	}

	free(mtr);
}

void calc_eker_ana(t_tr *tr, double t, complex *tr_ana)
{
	int i, j, k;
	complex zbuf1, zbuf2, zt1, zt2;
	complex zexp, zfac;
	complex f0, f1, f2, f3;
	complex zimag, zmimag;
	zimag.real  = 0.0, zimag.imag  =  1.0;
	zmimag.real = 0.0, zmimag.imag = -1.0;

	int nsys   = tr->nsys;
	int nsyssq = tr->nsyssq;
	int nosc   = tr->nosc_tr;

	z_init(tr_ana, nsyssq);

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				// Quantities involving all modes
				zbuf1.real = 0.0, zbuf1.imag = 0.0;
				zt1.real   = 0.0, zt1.imag   = 0.0;
				zt2.real   = 0.0, zt2.imag   = 0.0;
				zexp.real  = 0.0, zexp.imag  = 0.0;

				for(k=0; k<nosc; k++) {
					t_qosc *qosc = &tr->qosc[k];
			
					double coth   = qosc->coth;
					double w      = qosc->freq;
					double wsq    = w*w;
					double *gexci = qosc->gexci;

					f0.real = coth * (1.0 - cos(w*t)) / wsq;
					f0.imag = (sin(w*t) - w*t) / wsq;
					f1.real = coth * sin(w*t) / w;
					f1.imag = (cos(w*t) - 1.0) / w;
					f2.real = coth * cos(w*t);
					f2.imag = -sin(w*t);
					f3.real = -w * coth * sin(w*t);
					f3.imag = -w * cos(w*t);

					double gaa   = w * gexci[j + nsys*j];
					double gbb   = w * gexci[i + nsys*i];
					double gab   = w * gexci[j + nsys*i];
					double gba   = w * gexci[i + nsys*j];
					double Gba   = gbb - gaa;
					double Gbasq = Gba * Gba;

					zbuf1.real += gab * gba * f2.real;
					zbuf1.imag += gab * gba * f2.imag;

					zt1.real += gab * Gba * f1.real;
					zt1.imag += gab * (Gba * f1.imag - 2.0 * gaa / w);

					zt2.real += gba * Gba * f1.real;
					zt2.imag += gba * (Gba * f1.imag - 2.0 * gaa / w);

					zexp.real -= Gbasq * f0.real;
					zexp.imag -= Gbasq * (t/w + f0.imag);
				}

				double exp = pow(M_E, zexp.real);
				double arg = zexp.imag;
				zfac.real = exp * cos(arg);
				zfac.imag = exp * sin(arg);

				// Explicitly treated last mode
				t_qosc *qosc = &tr->qosc[nosc-1];

				double coth   = qosc->coth;
				double w      = qosc->freq;
				double wsq    = w*w;
				double *gexci = qosc->gexci;

				f0.real = coth * (1.0 - cos(w*t)) / wsq;
				f0.imag = (sin(w*t) - w*t) / wsq;
				f1.real = coth * sin(w*t) / w;
				f1.imag = (cos(w*t) - 1.0) / w;
				f2.real = coth * cos(w*t);
				f2.imag = -sin(w*t);
				f3.real = -w * coth * sin(w*t);
				f3.imag = -w * cos(w*t);

				double gaa   = w * gexci[j + nsys*j];
				double gbb   = w * gexci[i + nsys*i];
				double gab   = w * gexci[j + nsys*i];
				double gba   = w * gexci[i + nsys*j];
				double Gba   = gbb - gaa;
				double Gbasq = Gba * Gba;

				// The prefactor will be stored in zbuf1
				// Term 1
				z_mult(&zbuf2, zt1, zt2);

				zbuf1.real -= zbuf2.real;
				zbuf1.imag -= zbuf2.imag;

				z_mult(&zbuf2, zimag, f1);
				zbuf2.real *= Gbasq;
				zbuf2.imag *= Gbasq;
				zbuf2.real -= Gbasq / w;

				z_mult(&zbuf1, zbuf1, zbuf2);

				// Term 2
				zbuf2.real = Gba * (gba * zt1.real + gab * zt2.real);
				zbuf2.imag = Gba * (gba * zt1.imag + gab * zt2.imag);
				z_mult(&zbuf2, f2, zbuf2);
				z_mult(&zbuf2, zimag, zbuf2);

				zbuf1.real += zbuf2.real;
				zbuf1.imag += zbuf2.imag;

				// Term 3
				zbuf2.real = gab * gba * f3.real;
				zbuf2.imag = gab * gba * f3.imag;
				z_mult(&zbuf2, zmimag, zbuf2);

				zbuf1.real += zbuf2.real;
				zbuf1.imag += zbuf2.imag;

				// Final evaluation
				z_mult(&tr_ana[i + nsys*j], zbuf1, zfac);
			}
		}
	}
}
