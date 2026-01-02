#include "fret.h"

void fret_prep(t_qme *qme)
{
	int jtype = qme->jtype;

	qme->fret = malloc(sizeof(t_fret));
	t_fret *fret = qme->fret;
	init_fret(qme, fret);

	fret_calc_Hsys0_reorg(fret);

	if(jtype == 0 || jtype == 10) {
		fret_calc_g(fret);
		fret_calc_rates(fret);
	}

	if(jtype == 10) {
		fprintf(stderr, "---- Dissipation calculation requested ----\n");
		fret_calc_diss(fret);
	}
}

void fret_noisy_Hsys(t_fret *fret, t_split *split, double *std)
{
	int i, j;
	int nsys = split->nsys;
	int count = 0;

	double *Hsys_prist = split->Hsys_prist;
	double *Hsys_std   = split->Hsys_std;
	double *Hsys_noise = fret->Hsys;

	for(i=0; i<nsys; i++) {
		for(j=0; j<i+1; j++) {
			Hsys_noise[i + nsys*j] =  Hsys_prist[i + nsys*j];
			Hsys_noise[i + nsys*j] += std[count] * Hsys_std[i + nsys*j];
			count++;
		}
	}

	for(i=0; i<nsys; i++) {
		for(j=0; j<i; j++) {
			Hsys_noise[j + nsys*i] = Hsys_noise[i + nsys*j];
		}
	}
}

void fret_calc_Hsys0_reorg(t_fret *fret)
{
	int i, j, k, l;
	int one = 1;
	double w, *gamma;

	double *Hsys  = fret->Hsys;
	double *Hsys0 = fret->Hsys0;
	double *reorg = fret->reorg;
	double *Ediag = fret->Ediag;

	int nsys   = fret->nsys;
	int nsyssq = fret->nsyssq;
	int nspd   = fret->nspd;
	int *nosc  = fret->nosc;

	d_init(reorg, nsyssq);

	for(i=0; i<nspd; i++) {
		t_osc *osc = fret->osc[i];

		for(j=0; j<nosc[i]; j++) {
			w = osc[j].freq;
			gamma = osc[j].gamma;

			for(k=0; k<nsys; k++) {
				for(l=0; l<k+1; l++) {
					reorg[k + nsys*l] += w * gamma[k + nsys*k] * gamma[l + nsys*l];
				}
			}
		}
	}

	// Symmetrize
	for(i=0; i<nsys; i++) {
		for(j=0; j<i; j++) {
			reorg[j + nsys*i] = reorg[i + nsys*j];
		}
	}

	// Construct Hsys0 and Ediag
	dcopy(&nsyssq, Hsys, &one, Hsys0, &one);

	for(i=0; i<nsys; i++) {
		Hsys0[i + nsys*i] = Hsys[i + nsys*i] - reorg[i + nsys*i];
		Ediag[i] = Hsys0[i + nsys*i];
	}
}

void fret_calc_gderiv(t_fret *fret)
{
	int i, j, k, l, m;
	int one = 1;
	double arg, hrgen;
	complex f1wsq, f2wsq;

	int nsys     = fret->nsys;
	int nsyssq   = fret->nsyssq;
	int nspd     = fret->nspd;
	int *nosc    = fret->nosc;
	t_osc **osc  = fret->osc;

	int n     = fret->n;
	double dt = fret->dt;

	for(i=0; i<nsyssq; i++) {
		z_init(fret->g1t[i], n);
		z_init(fret->g2t[i], n);
	}

	for(i=0; i<nspd; i++) {
		for(j=0; j<nosc[i]; j++) {
			double w = osc[i][j].freq;
			double coth = osc[i][j].coth;
			double *gamma = osc[i][j].gamma;

			for(k=0; k<nsys; k++) {
				for(l=0; l<k+1; l++) {
					hrgen = gamma[k + nsys*k] * gamma[l + nsys*l];
					complex *g1t = fret->g1t[k + nsys*l];
					complex *g2t = fret->g2t[k + nsys*l];

					for(m=0; m<n; m++) {
						arg = w * m * dt;

						f1wsq.real = w * coth * sin(arg);
						f1wsq.imag = w * (cos(arg) - 1.0);

						f2wsq.real = w * w * coth * cos(arg);
						f2wsq.imag = -w * w * sin(arg);

						g1t[m].real += hrgen * f1wsq.real;
						g1t[m].imag += hrgen * f1wsq.imag;

						g2t[m].real += hrgen * f2wsq.real;
						g2t[m].imag += hrgen * f2wsq.imag;
					}
				}
			}
		}
	}

	for(i=0; i<nsys; i++) {
		for(j=0; j<i; j++) {
			zcopy(&n, fret->gt[i + nsys*j], &one, fret->gt[j + nsys*i], &one);
		}
	}
}

void fret_calc_g(t_fret *fret)
{
	int i, j, k, l, m;
	int one = 1;
	double arg, hrgen;
	complex f0wsq;

	int nsys     = fret->nsys;
	int nsyssq   = fret->nsyssq;
	int nspd     = fret->nspd;
	int *nosc    = fret->nosc;
	t_osc **osc  = fret->osc;

	int n     = fret->n;
	double dt = fret->dt;

	for(i=0; i<nsyssq; i++) {
		z_init(fret->gt[i], n);
	}

	for(i=0; i<nspd; i++) {
		for(j=0; j<nosc[i]; j++) {
			double w = osc[i][j].freq;
			double coth = osc[i][j].coth;
			double *gamma = osc[i][j].gamma;

			for(k=0; k<nsys; k++) {
				for(l=0; l<k+1; l++) {
					hrgen = gamma[k + nsys*k] * gamma[l + nsys*l];
					complex *gt = fret->gt[k + nsys*l];

					for(m=0; m<n; m++) {
						arg = w * m * dt;

						f0wsq.real = coth * (1.0 - cos(arg));
						f0wsq.imag = sin(arg) - arg;

						gt[m].real += hrgen * f0wsq.real;
						gt[m].imag += hrgen * f0wsq.imag;
					}
				}
			}

			if((!fret->bSPLIT) && (j % 100 == 0)) fprintf(stderr, "SPD %2d  OSC %6d / %6d\n", i+1, j, nosc[i]);
		}
		if(!fret->bSPLIT) fprintf(stderr, "\n");
	}
	if(!fret->bSPLIT) fprintf(stderr, "\n");

	for(i=0; i<nsys; i++) {
		for(j=0; j<i; j++) {
			zcopy(&n, fret->gt[i + nsys*j], &one, fret->gt[j + nsys*i], &one);
		}
	}
}

void fret_calc_rates(t_fret *fret)
{
	int i, j, k;
	double exp_re;
	complex mzone, ztwo;
	mzone.real = -1.0, mzone.imag = 0.0;
	ztwo.real = 2.0,   ztwo.imag = 0.0;

	int nsys   = fret->nsys;
	int nsyssq = fret->nsyssq;

	int n     = fret->n;
	double dt = fret->dt;

	double *Hsys0 = fret->Hsys0;
	double *reorg = fret->reorg;
	complex **gt  = fret->gt;
	double *rate  = fret->rate;

	complex *t_ji = malloc(n * sizeof(complex));
	complex *t_ij = malloc(n * sizeof(complex));

	d_init(rate, nsyssq);

	for(i=0; i<nsys; i++) {
		for(j=0; j<i; j++) {
			z_init(t_ij, n);
			z_init(t_ji, n);

			// Identical components
			zaxpy(&n, &mzone, gt[i + nsys*i], &one, t_ij, &one);
			zaxpy(&n, &ztwo,  gt[i + nsys*j], &one, t_ij, &one);
			zaxpy(&n, &mzone, gt[j + nsys*j], &one, t_ij, &one);

			double reorg_sum = reorg[i + nsys*i] - 2.0 * reorg[i + nsys*j] + reorg[j + nsys*j];
//			fprintf(stderr, "reorg_sum %d %d %13.6lf %13.6lf %13.6lf  %13.6lf\n", i, j, reorg[i + nsys*i], reorg[i + nsys*j], reorg[j + nsys*j], reorg_sum);
//			fprintf(stderr, "Hsys0 %13.6lf %13.6lf\n", Hsys0[j + nsys*j], Hsys0[i + nsys*i]);

			for(k=0; k<n; k++) {
				t_ij[k].imag -= reorg_sum * k * dt;
			}

			zcopy(&n, t_ij, &one, t_ji, &one);

			// Different components
			double sitee_diff = Hsys0[j + nsys*j] - Hsys0[i + nsys*i];

			for(k=0; k<n; k++) {
				t_ij[k].imag += sitee_diff * k * dt;
				t_ji[k].imag -= sitee_diff * k * dt;
			}

			complex *Gt_ij = fret->Gt[i + nsys*j];
			complex *Gt_ji = fret->Gt[j + nsys*i];

			for(k=0; k<n; k++) {
				exp_re = pow(M_E, t_ij[k].real);
				Gt_ij[k].real = exp_re * cos(t_ij[k].imag);
				Gt_ij[k].imag = exp_re * sin(t_ij[k].imag);

				exp_re = pow(M_E, t_ji[k].real);
				Gt_ji[k].real = exp_re * cos(t_ji[k].imag);
				Gt_ji[k].imag = exp_re * sin(t_ji[k].imag);
			}

			// Integration by trapezoidal method
			rate[i + nsys*j] += Gt_ij[0].real / 2.0;
			rate[j + nsys*i] += Gt_ji[0].real / 2.0;

			for(k=1; k<n-1; k++) {
				rate[i + nsys*j] += Gt_ij[k].real;
				rate[j + nsys*i] += Gt_ji[k].real;
			}

			rate[i + nsys*j] += Gt_ij[n-1].real / 2.0;
			rate[j + nsys*i] += Gt_ji[n-1].real / 2.0;

			// Final multiplication
			double twoVsq = 2.0 * Hsys0[i + nsys*j] * Hsys0[i + nsys*j];

			rate[i + nsys*j] *= twoVsq * dt;
			rate[j + nsys*i] *= twoVsq * dt;
		}
	}

// vvvvv
#if 0
	double time;
	complex **g1t = fret->g1t;
	complex **g2t = fret->g2t;
	complex **Gt  = fret->Gt;
	FILE *f_func = fopen("func.txt", "w");
	for(i=0; i<n; i++) {
		time = i * dt * TUNIT_TO_FS;	
		if(time < 10000.001) fprintf(f_func, "%13.6lf %13.6le %13.6le %13.6le %13.6le\n", time, g2t[0][i].real, g1t[0][i].real, gt[0][i].real, Gt[1][i].real);
	}
	fclose(f_func);
#endif
// ^^^^^

	if(fret->unit) {
		fprintf(stderr, "---- Population Transfer Rate Constants (1/fs) ----\n");
		for(i=0; i<nsys; i++) {
			for(j=0; j<nsys; j++) {
				fprintf(stderr, "%15.6le", rate[i + nsys*j] / TUNIT_TO_FS);
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");
	} else {
		fprintf(stderr, "---- Population Transfer Rate Constants (a.u.) ----\n");
		for(i=0; i<nsys; i++) {
			for(j=0; j<nsys; j++) {
				fprintf(stderr, "%15.6le", rate[i + nsys*j]);
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");
	}

	if(fret->bTEMP_UNI) {
		fprintf(stderr, "---- Detailed Balance ----\n");
		for(i=0; i<nsys; i++) {
			for(j=0; j<i; j++) {
				fprintf(stderr, "%3d <-> %3d  %13.6le %13.6le\n", j+1, i+1, rate[j + nsys*i] / rate[i + nsys*j], pow(M_E, -(Hsys0[j + nsys*j] - Hsys0[i + nsys*i]) / fret->temp_uni));
			}
		}
		fprintf(stderr, "\n");
	}

	free(t_ji);
	free(t_ij);
}

void fret_calc_diss(t_fret *fret)
{
	int i, j, k, l, m;
	double arg;
	complex zfac;

	int nsys = fret->nsys;

	double *Hsys0 = fret->Hsys0;

	int n     = fret->n;
	double dt = fret->dt;

	int nspd    = fret->nspd;
	int *nosc   = fret->nosc;
	t_osc **osc = fret->osc;

	complex *t_lk = malloc(n * sizeof(complex));
	complex *t_kl = malloc(n * sizeof(complex));

	fprintf(stderr, "---- Calculating FRET dissipation parameters ----\n");

	for(i=0; i<nspd; i++) {
		for(j=0; j<nosc[i]; j++) {
			double w      = osc[i][j].freq;
			double coth   = osc[i][j].coth;
			double spread = osc[i][j].spread;
			double *gamma = osc[i][j].gamma;
			double *Idiss = osc[i][j].Idiss;
			double *Kdiss = osc[i][j].Kdiss;
			double *Jdiss = osc[i][j].Jdiss;

			for(k=0; k<nsys; k++) {
				for(l=0; l<k; l++) {
					complex *Gt_kl = fret->Gt[k + nsys*l];
					complex *Gt_lk = fret->Gt[l + nsys*k];

					z_init(t_kl, n);
					z_init(t_lk, n);

					for(m=0; m<n; m++) {
						arg = w * m * dt;
						zfac.real = cos(arg);
						zfac.imag = -coth * sin(arg);
						z_mult(&t_kl[m], Gt_kl[m], zfac);
						z_mult(&t_lk[m], Gt_lk[m], zfac);
					}

					// Integration by trapezoidal method
					Idiss[k + nsys*l] = t_kl[0].real / 2.0;
					Idiss[l + nsys*k] = t_lk[0].real / 2.0;

					for(m=1; m<n-1; m++) {
						Idiss[k + nsys*l] += t_kl[m].real;
						Idiss[l + nsys*k] += t_lk[m].real;
					}

					Idiss[k + nsys*l] += t_kl[n-1].real / 2.0;
					Idiss[l + nsys*k] += t_lk[n-1].real / 2.0;

					Idiss[k + nsys*l] *= dt;
					Idiss[l + nsys*k] *= dt;

					// Calculate Jdiss
					double twoVsq = 2.0 * Hsys0[k + nsys*l] * Hsys0[k + nsys*l];
					double gamma_diff = gamma[k + nsys*k] - gamma[l + nsys*l];
					double reorg_gen = w * gamma_diff * gamma_diff;

					Kdiss[k + nsys*l] = twoVsq * reorg_gen * Idiss[k + nsys*l];
					Kdiss[l + nsys*k] = twoVsq * reorg_gen * Idiss[l + nsys*k];

					Jdiss[k + nsys*l] = Kdiss[k + nsys*l] / spread;
					Jdiss[l + nsys*k] = Kdiss[l + nsys*k] / spread;
				}
			}

			if((!fret->bSPLIT) && (j % 100 == 0)) fprintf(stderr, "SPD %2d  OSC %6d / %6d\n", i+1, j, nosc[i]);
		}
		if(!fret->bSPLIT) fprintf(stderr, "\n");
	}

#if 1
	check_fret_physics(fret);
#endif
	
	free(t_lk);
	free(t_kl);
}

void check_fret_physics(t_fret *fret)
{
	// Assumes that the temperature of all modes are the same
	int i, j, k, l;

	int nsys = fret->nsys;
	int nspd = fret->nspd;
	int *nosc   = fret->nosc;
	t_osc **osc = fret->osc;

	double *rate  = fret->rate;
	double *Ediag = fret->Ediag;
	double T      = osc[0][0].temp;

	// Detailed Balance
	fprintf(stderr, "Detailed Balance for Dissipation\n");

	for(i=0; i<nsys; i++) {
		for(j=0; j<i; j++) {
			for(k=0; k<nspd; k++) {
				for(l=0; l<nosc[k]; l++) {
					double w      = osc[k][l].freq;
					double *Jdiss = osc[k][l].Jdiss;
					double diss_ratio = Jdiss[i + nsys*j] == 0.0 ? 0.0 : -Jdiss[j + nsys*i] / Jdiss[i + nsys*j];
					if(l % 50 == 0) {
						fprintf(stderr, "SPD %3d   %2d <-> %2d  OSC %5d   freq %13.6le   diss_ratio %13.6le   Boltzmann_ratio %13.6le\n",
								k+1, i+1, j+1, l, w, diss_ratio, pow(M_E, -(Ediag[j] - Ediag[i]) / T));
					}
				}
				fprintf(stderr, "\n");
			}
		}
	}

	// Energy Conservation
	fprintf(stderr, "Energy Conservation\n");

	for(i=0; i<nsys; i++) {
		for(j=0; j<i; j++) {
			double Edec_elec_ij = -(Ediag[i] - Ediag[j]) * rate[i + nsys*j];
			double Edec_elec_ji = -(Ediag[j] - Ediag[i]) * rate[j + nsys*i];
			double Einc_diss_ij = 0.0;
			double Einc_diss_ji = 0.0;

			for(k=0; k<nspd; k++) {
				for(l=0; l<nosc[k]; l++) {
					Einc_diss_ji += osc[k][l].Kdiss[j + nsys*i];
					Einc_diss_ij += osc[k][l].Kdiss[i + nsys*j];
				}
			}

			fprintf(stderr, "%2d -> %2d  Elec %13.6le  Diss %13.6le\n", j+1, i+1, Edec_elec_ij, Einc_diss_ij);
			fprintf(stderr, "%2d -> %2d  Elec %13.6le  Diss %13.6le\n", i+1, j+1, Edec_elec_ji, Einc_diss_ji);
		}
	}
}
