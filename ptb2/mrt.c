#include "mrt.h"

void mrt_prep(t_qme *qme)
{
	int jtype = qme->jtype;

	qme->mrt = malloc(sizeof(t_mrt));
	t_mrt *mrt = qme->mrt;
	init_mrt(qme, mrt);

	mrt_calc_H_reorg(mrt);

	if(jtype == 1 || jtype == 11) {
		mrt_calc_g(mrt);
		mrt_calc_rates(mrt);
	}

	if(jtype == 11) {
		fprintf(stderr, "---- Dissipation calculation requested ----\n");
		mrt_calc_diss(mrt);
	}
}

void mrt_calc_H_reorg(t_mrt *mrt)
{
	int i, j, k, l;
	int one = 1;
	double w, *gamma, *gexci;

	double *Hsys   = mrt->Hsys;
	double *Usys   = mrt->Usys;
	double *Eexci  = mrt->Eexci;
	double *Eexci0 = mrt->Eexci0;

	double *reaabb = mrt->reaabb;
	double *reabbb = mrt->reabbb;

	int nsys   = mrt->nsys;
	int nsyssq = mrt->nsyssq;
	int nspd   = mrt->nspd;
	int *nosc  = mrt->nosc;

	diag(Usys, Eexci, Hsys, nsys);

	for(i=0; i<nspd; i++) {
		t_osc *osc = mrt->osc[i];

		for(j=0; j<nosc[i]; j++) {
			w = osc[j].freq;
			gamma = osc[j].gamma;
			gexci = osc[j].gexci;
			dcopy(&nsyssq, gamma, &one, gexci, &one);
			trans_inv(Usys, gexci, nsys);

			// Symmetric component
			for(k=0; k<nsys; k++) {
				for(l=0; l<k+1; l++) {
					reaabb[k + nsys*l] += w * gexci[l + nsys*l] * gexci[k + nsys*k];
				}
			}

			// Asymmetric component
			for(k=0; k<nsys; k++) {
				for(l=0; l<nsys; l++) {
					reabbb[k + nsys*l] += w * gexci[k + nsys*l] * gexci[l + nsys*l];
				}
			}
		}
	}

	// Symmetrize
	for(i=0; i<nsys; i++) {
		for(j=0; j<i; j++) {
			reaabb[j + nsys*i] = reaabb[i + nsys*j];
		}
	}

	// Construct Eexci0
	for(i=0; i<nsys; i++) {
		Eexci0[i] = Eexci[i] - reaabb[i + nsys*i];
	}

// vvvvv
#if 0
	for(i=0; i<nsys; i++) {
		fprintf(stderr, "Eexci %13.6le  Eexci0 %13.6le\n", Eexci[i], Eexci0[i]);
	}
#endif
// ^^^^^
}

void mrt_calc_g(t_mrt *mrt)
{
	int i, j, k, l, m;
	int one = 1;
	double arg;

	int nsys    = mrt->nsys;
	int nsyssq  = mrt->nsyssq;
	int nspd    = mrt->nspd;
	int *nosc   = mrt->nosc;
	t_osc **osc = mrt->osc;

	int n     = mrt->n;
	double dt = mrt->dt;

	for(i=0; i<nsyssq; i++) {
		z_init(mrt->gtaabb[i], n);
		z_init(mrt->gtabbb[i], n);

		z_init(mrt->g1abbb[i], n);
		z_init(mrt->g2abab[i], n);
	}

	for(i=0; i<nspd; i++) {
		for(j=0; j<nosc[i]; j++) {
			double w = osc[i][j].freq;
			double coth = osc[i][j].coth;
			double *gexci = osc[i][j].gexci;

			// Symmetric components
			for(k=0; k<nsys; k++) {
				for(l=0; l<k+1; l++) {
					double hraabb    = gexci[k + nsys*k] * gexci[l + nsys*l];
					double hrabab_w2 = gexci[k + nsys*l] * gexci[k + nsys*l] * w * w;

					complex *gtaabb = mrt->gtaabb[k + nsys*l];
					complex *g2abab = mrt->g2abab[k + nsys*l];

					for(m=0; m<n; m++) {
						arg = w * m * dt;

						gtaabb[m].real += hraabb * coth * (1.0 - cos(arg));
						gtaabb[m].imag += hraabb * (sin(arg) - arg);

						g2abab[m].real += hrabab_w2 * coth * cos(arg);
						g2abab[m].imag -= hrabab_w2 * sin(arg);
					}
				}
			}

			// Asymmetric components
			for(k=0; k<nsys; k++) {
				for(l=0; l<nsys; l++) {
					double hrabbb    = gexci[k + nsys*l] * gexci[l + nsys*l];
					double hrabbb_w  = hrabbb * w;

					complex *gtabbb = mrt->gtabbb[k + nsys*l];
					complex *g1abbb = mrt->g1abbb[k + nsys*l];

					for(m=0; m<n; m++) {
						arg = w * m * dt;

						gtabbb[m].real += hrabbb * coth * (1.0 - cos(arg));
						gtabbb[m].imag += hrabbb * (sin(arg) - arg);

						g1abbb[m].real += hrabbb_w * coth * sin(arg);
						g1abbb[m].imag += hrabbb_w * (cos(arg) - 1.0);
					}
				}
			}
			if(j % 100 == 0) fprintf(stderr, "SPD %2d  OSC %6d / %6d\n", i+1, j, nosc[i]);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "\n");

	for(i=0; i<nsys; i++) {
		for(j=0; j<i; j++) {
			zcopy(&n, mrt->gtaabb[i + nsys*j], &one, mrt->gtaabb[j + nsys*i], &one);
			zcopy(&n, mrt->g2abab[i + nsys*j], &one, mrt->g2abab[j + nsys*i], &one);
		}
	}
}

void mrt_calc_rates(t_mrt *mrt)
{
	int i, j, k;
	double exp_re;
	complex mzone, ztwo, zbuf;
	mzone.real = -1.0, mzone.imag = 0.0;
	ztwo.real = 2.0,   ztwo.imag = 0.0;

	int nsys = mrt->nsys;

	int n     = mrt->n;
	double dt = mrt->dt;

	double *Eexci0   = mrt->Eexci0;
	double *reaabb   = mrt->reaabb;
	double *reabbb   = mrt->reabbb;

	complex **gtaabb = mrt->gtaabb;
	complex **gtabbb = mrt->gtabbb;
	complex **g1abbb = mrt->g1abbb;
	complex **g2abab = mrt->g2abab;

	double *rate     = mrt->rate;

	complex *zexp = malloc(n * sizeof(complex));
	complex *zker = malloc(n * sizeof(complex));

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				z_init(zexp, n);

				// i <-> b, j <-> a

				// Identical components
				zaxpy(&n, &mzone, gtaabb[j + nsys*j], &one, zexp, &one);
				zaxpy(&n, &ztwo,  gtaabb[j + nsys*i], &one, zexp, &one);
				zaxpy(&n, &mzone, gtaabb[i + nsys*i], &one, zexp, &one);

				double reorg_sum = reaabb[j + nsys*j] - 2.0 * reaabb[j + nsys*i] + reaabb[i + nsys*i];

				for(k=0; k<n; k++) {
					zexp[k].imag -= reorg_sum * k * dt;
				}

				// Different components
				double Eexci_diff = Eexci0[i] - Eexci0[j];

				for(k=0; k<n; k++) {
					zexp[k].imag -= Eexci_diff * k * dt;
				}

				complex *Gt = mrt->Gt[i + nsys*j];

				for(k=0; k<n; k++) {
					exp_re = pow(M_E, zexp[k].real);
					Gt[k].real = exp_re * cos(zexp[k].imag);
					Gt[k].imag = exp_re * sin(zexp[k].imag);
				}

				complex *Nt = mrt->Nt[i + nsys*j];
				complex *Zt = mrt->Zt[i + nsys*j];

				zcopy(&n, g2abab[i + nsys*j], &one, Nt, &one);

				for(k=0; k<n; k++) {
					Zt[k].real = g1abbb[j + nsys*i][k].real - g1abbb[i + nsys*j][k].real;
					Zt[k].imag = g1abbb[j + nsys*i][k].imag - g1abbb[i + nsys*j][k].imag - 2.0 * reabbb[i + nsys*j];
					z_mult(&zbuf, Zt[k], Zt[k]);

					Nt[k].real -= zbuf.real;
					Nt[k].imag -= zbuf.imag;
					z_mult(&zker[k], Gt[k], Nt[k]);
				}

				// Integration by trapezoidal method
				rate[i + nsys*j] = zker[0].real / 2.0;

				for(k=1; k<n-1; k++) {
					rate[i + nsys*j] += zker[k].real;
				}

				rate[i + nsys*j] += zker[n-1].real / 2.0;

				// Final multiplication
				double twodt = 2.0 * dt;
				rate[i + nsys*j] *= twodt;
			}
		}
	}

	fprintf(stderr, "---- Population Transfer Rate Constants (a.u.) ----\n");
	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
            fprintf(stderr, "%15.6le", rate[i + nsys*j]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");

	fprintf(stderr, "---- Population Transfer Rate Constants (1/fs) ----\n");
	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
            fprintf(stderr, "%15.6le", rate[i + nsys*j] / TUNIT_TO_FS);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");

	if(nsys == 2) {
		fprintf(stderr, "Total transfer rate: %15.6le a.u. (%15.6le fs^-1)\n\n", rate[1] + rate[2], (rate[1] + rate[2]) / TUNIT_TO_FS);
	}

// vvvvv
#if 0
	fprintf(stderr, "---- Detailed Balance ----\n");
	for(i=0; i<nsys; i++) {
		for(j=0; j<i; j++) {
			fprintf(stderr, "%13.6le %13.6le\n", rate[j + nsys*i] / rate[i + nsys*j], pow(M_E, -(Eexci0[j] - Eexci0[i]) / 1.0));
		}
	}
#endif
// ^^^^^

	free(zexp);
	free(zker);
}

void mrt_calc_diss(t_mrt *mrt)
{
	int i, j, k, l, m;
	complex zbuf1, zbuf2;
	complex zimag, zmimag;
	zimag.real  = 0.0, zimag.imag  =  1.0;
	zmimag.real = 0.0, zmimag.imag = -1.0;

	int nsys = mrt->nsys;

	int n     = mrt->n;
	double dt = mrt->dt;
	double t;

	int nspd  = mrt->nspd;
	int *nosc = mrt->nosc;

	double *Eexci0 = mrt->Eexci0;

	complex *zker = malloc(n * sizeof(complex));
	complex *f1   = malloc(n * sizeof(complex));
	complex *f2   = malloc(n * sizeof(complex));
	complex *f3   = malloc(n * sizeof(complex));

	fprintf(stderr, "---- Calculating MRT dissipation parameters ----\n");

	// i <-> b, j <-> a
	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				complex *g1abbb = mrt->g1abbb[j + nsys*i];
				complex *g1abaa = mrt->g1abbb[i + nsys*j];
				complex *g2abab = mrt->g2abab[i + nsys*j];
				complex *Gt     = mrt->Gt[i + nsys*j];
				complex *Zt     = mrt->Zt[i + nsys*j];
				complex *Nt     = mrt->Nt[i + nsys*j];

				// Explicit treatment of a particular mode
				for(k=0; k<nspd; k++) {
					for(l=0; l<nosc[k]; l++) {
						t_osc *osc    = &(mrt->osc[k][l]);
						double coth   = osc->coth;
						double w      = osc->freq;
						double wsq    = w*w;
						double *gexci = osc->gexci;

						double gaa   = w * gexci[j + nsys*j];
						double gbb   = w * gexci[i + nsys*i];
						double gab   = w * gexci[j + nsys*i];
						double gba   = w * gexci[i + nsys*j];
						double Gba   = gbb - gaa;
						double Gbasq = Gba * Gba;

						for(m=0; m<n; m++) {
							t = dt * (double)m;

							f1[m].real = coth * sin(w*t) / w;
							f1[m].imag = (cos(w*t) - 1.0) / w;

							f2[m].real = coth * cos(w*t);
							f2[m].imag = -sin(w*t);

							f3[m].real = -w * coth * sin(w*t);
							f3[m].imag = -w * cos(w*t);

							// The prefactor will be stored in zbuf1
							// Term 1
							zbuf1.real = Nt[m].real;
							zbuf1.imag = Nt[m].imag;

							z_mult(&zbuf2, zimag, f1[m]);
							zbuf2.real *= Gbasq;
							zbuf2.imag *= Gbasq;
							zbuf2.real -= Gbasq / w;

							z_mult(&zbuf1, zbuf1, zbuf2);

							// Term 2
							zbuf2.real = Gba * (gba * Zt[m].real + gab * Zt[m].real);
							zbuf2.imag = Gba * (gba * Zt[m].imag + gab * Zt[m].imag);
							z_mult(&zbuf2, f2[m], zbuf2);
							z_mult(&zbuf2, zimag, zbuf2);

							zbuf1.real += zbuf2.real;
							zbuf1.imag += zbuf2.imag;

							// Term 3
							zbuf2.real = gab * gba * f3[m].real;
							zbuf2.imag = gab * gba * f3[m].imag;
							z_mult(&zbuf2, zmimag, zbuf2);

							zbuf1.real += zbuf2.real;
							zbuf1.imag += zbuf2.imag;

							// Final evaluation
							z_mult(&zker[m], zbuf1, Gt[m]);
						}

						double reorg_exci = mrt->reaabb[j + nsys*j] - 2.0 * mrt->reaabb[j + nsys*i] + mrt->reaabb[i + nsys*i];
						double *Idiss = osc->Idiss; // we will make it zero in MRT
						double *Kdiss = osc->Kdiss;
						double *Jdiss = osc->Jdiss;
						double spread = osc->spread;

						// Integration by trapezoidal method
						Kdiss[i + nsys*j] = zker[0].real / 2.0;

						for(m=1; m<n-1; m++) {
							Kdiss[i + nsys*j] += zker[m].real;
						}

						Kdiss[i + nsys*j] += zker[n-1].real / 2.0;

						// Final evaluation and convert to continuous density
						Kdiss[i + nsys*j] *= -(2.0 * dt);
						Idiss[i + nsys*j] = Kdiss[i + nsys*j] / reorg_exci;
						Jdiss[i + nsys*j] = Kdiss[i + nsys*j] / spread;

						if(l % 100 == 0) fprintf(stderr, "%3d -> %3d    SPD %2d    OSC %6d / %6d\n", j+1, i+1, k+1, l, nosc[k]);
					}
				}
			}
		}
	}

#if 1
	check_mrt_physics(mrt);
#endif

	free(zker);

	free(f1);
	free(f2);
	free(f3);
}

void check_mrt_physics(t_mrt *mrt)
{
	// Assumes that the temperature of all modes are the same
	int i, j, k, l;

	int nsys    = mrt->nsys;
	int nspd    = mrt->nspd;
	int *nosc   = mrt->nosc;
	t_osc **osc = mrt->osc;

	double *rate   = mrt->rate;
	double *Eexci0 = mrt->Eexci0;
	double T       = osc[0][0].temp;

	// Detailed Balance
	fprintf(stderr, "Detailed Balance for Dissipation\n");

	for(i=0; i<nsys; i++) {
		for(j=0; j<i; j++) {
			for(k=0; k<nspd; k++) {
				for(l=0; l<nosc[k]; l++) {
					double w      = osc[k][l].freq;
					double *Jdiss = osc[k][l].Jdiss;
					if(l % 50 == 0) {
						fprintf(stderr, "SPD %3d   OSC %5d   freq %13.6le   diss_ratio %13.6le   Boltzmann_ratio %13.6le\n",
								k+1, l, w, -Jdiss[j + nsys*i] / Jdiss[i + nsys*j], pow(M_E, -(Eexci0[j] - Eexci0[i]) / T));
					}
				}
				fprintf(stderr, "\n");
			}
		}
	}

#if 1
	// Energy Conservation
	fprintf(stderr, "Energy Conservation\n");

	for(i=0; i<nsys; i++) {
		for(j=0; j<i; j++) {
			double Edec_elec_ij = -(Eexci0[i] - Eexci0[j]) * rate[i + nsys*j];
			double Edec_elec_ji = -(Eexci0[j] - Eexci0[i]) * rate[j + nsys*i];
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
#endif
}
