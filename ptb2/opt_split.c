#include "opt_split.h"

void avg_split(t_qme *qme)
{
	t_split *split = qme->split;

	init_split(qme, split);

	if(split->jtype == 30) { // FRET
		split->fret = malloc(sizeof(t_fret));
		init_fret(qme, split->fret);
	} else { // MRT
		split->mrt = malloc(sizeof(t_mrt));
		init_mrt(qme, split->mrt);
	}

	calc_Hsys_offset(split);

	split->dyn = malloc(sizeof(t_dyn_inco));
	init_dyn_inco(qme, split->dyn);

	avg_split_lowlvl(split);
}

void calc_Hsys_offset(t_split *split)
{
#if 0
	int nsys   = split->nsys;
	int nsyssq = nsys * nsys;

	bool b1 = split->bHSYS_ADJUST;
	bool b2 = split->fret->bTEMP_UNI;

	if(b1 && b2) {
		
	} else {
		d_init(
	}
#endif
}

void avg_split_lowlvl(t_split *split)
{
	int i, j, k, l;
	int one = 1;
	char fname[256];

	double tfac = split->unit ? TUNIT_TO_FS : 1.0;
	double wfac = split->unit ? 1.0 / WAVENO_TO_EUNIT : 1.0;
	double Jfac = split->unit ? 1.0 / TUNIT_TO_FS : 1.0;

	int ntraj  = split->ntraj;
	int nsys   = split->nsys;
	int nsyssq = nsys * nsys;
	int nout   = split->nout;
	int nelem  = nsys*nout;
	double dt  = split->dtout;

	t_fret *fret    = split->fret;
	t_dyn_inco *dyn = split->dyn;
	double **std    = split->std;

	int nspd    = split->nspd;
	int *nosc   = split->nosc;
	t_osc **osc = split->osc0;

	double *popt_avg = split->popt_avg;
	double *popt_one = dyn->popt_sto;

	double **disst_avg = split->disst_avg;
	double *disst_one = dyn->disst;

	fprintf(stderr, "Starting averaging the dynamics by spectral density partition\n");
	fprintf(stderr, "No. of trajectories: %d\n", ntraj);
	fprintf(stderr, "\n");

	split_spd_fret(split);
	d_init(split->popt_avg, nelem);

	for(i=0; i<ntraj; i++) {
		if(split->jtype == 30) { // FRET
			fprintf(stderr, "Noise realization %5d / %5d\n", i+1, ntraj);
			fret_noisy_Hsys(fret, split, std[i]);
			fret_calc_Hsys0_reorg(fret);
//			fret_calc_gderiv(fret);
			fret_calc_g(fret);
			fret_calc_rates(fret);

			if(split->bDISS) {
				fret_calc_diss(fret);
			}
		} else { // MRT
		}

		if(!split->bDISS) {
			prop_pop_inco(dyn);

			for(j=0; j<nelem; j++) {
				popt_avg[j] += popt_one[j];
			}
		} else {
			prop_diss_inco(dyn, i);
			
			int count = 0;
			for(j=0; j<nspd; j++) {
				for(k=0; k<nosc[j]; k++) {
					disst_avg[j][k] += disst_one[count];
					count++;
				}
			}
		}
	}

	// Sum to average
	if(!split->bDISS) {
		for(i=0; i<nelem; i++) {
			popt_avg[i] /= ntraj;
		}
	} else {
		for(i=0; i<nspd; i++) {
			for(j=0; j<nosc[i]; j++) {
				disst_avg[i][j] /= ntraj;
			}
		}
	}

	// Print
	if(!split->bDISS) {
		int count = 0;
		FILE *favg = fopen("popt_avg.txt", "w");
		for(i=0; i<nout; i++) {
			fprintf(favg, "%13.6lf", i * dt * tfac);
			for(j=0; j<nsys; j++) {
				fprintf(favg, "%13.6lf", popt_avg[count]);
				count++;
			}
			fprintf(favg, "\n");
		}
		fclose(favg);
	} else {
		for(i=0; i<nspd; i++) {
			sprintf(fname, "diss_spd%02d_avg.txt", i+1);
			FILE *f_diss = fopen(fname, "w");

			fprintf(f_diss, "           Freq           diss\n");

			for(j=0; j<nosc[i]; j++) {
				fprintf(f_diss, "%15.6le%15.6le\n", osc[i][j].freq * wfac, disst_avg[i][j]);
			}

			fclose(f_diss);
		}
	}

#if 0
	fprintf(stderr, "---- Last sampling for dissipation ----\n");
	fret_calc_diss(fret);
#endif
}

void split_spd_fret(t_split *split)
{
	int i, j, k;
	int one = 1;

	int nsys   = split->nsys;
	int nsyssq = split->nsyssq;
	int nspd   = split->nspd;
	int *nosc  = split->nosc;

	int type_fcn = split->type_fcn;

	double wsp       = split->freq_split;
	double pref      = split->pref_split;
	double *Hsys_std = split->Hsys_std;

	for(i=0; i<nspd; i++) {
		t_osc *osc0 = split->osc0[i];
		t_osc *osc  = split->fret->osc[i];

		for(j=0; j<nosc[i]; j++) {
			double w       = osc0[j].freq;
			double wsq     = w * w;
			double coth    = osc0[j].coth;
			double *gamma0 = osc0[j].gamma;
			double *gamma  = osc[j].gamma;

			if(type_fcn == 0) {
				if(w < wsp) {
					double scalefac = calc_scalefactor_2015castillo_gen(w, wsp, pref);
					double remain = sqrt(1.0 - scalefac);

					for(k=0; k<nsyssq; k++) {
						gamma[k] = remain * gamma0[k];
						Hsys_std[k] += scalefac * wsq * gamma0[k] * gamma0[k] * coth;
					}
				} else {
					dcopy(&nsyssq, gamma0, &one, gamma, &one);
				}
			} else if(type_fcn == 1) {
				double scalefac = calc_scalefactor_exponential(w, wsp, pref);
				double remain = sqrt(1.0 - scalefac);

				for(k=0; k<nsyssq; k++) {
					gamma[k] = remain * gamma0[k];
					Hsys_std[k] += scalefac * wsq * gamma0[k] * gamma0[k] * coth;
				}
			}
		}
	}

	for(i=0; i<nsyssq; i++) {
		Hsys_std[i] = sqrt(Hsys_std[i]);
	}

// vvvvv
#if 0
	double fac = (split->unit) ? 1.0 / WAVENO_TO_EUNIT : 1;

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			fprintf(stderr, "%15.6le ", split->Hsys0[i + nsys*j] * fac);
		}
		fprintf(stderr, "\n");
	}

	fprintf(stderr, "\n");

	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			fprintf(stderr, "%15.6le ", Hsys_std[i + nsys*j] * fac);
		}
		fprintf(stderr, "\n");
	}
#endif
// ^^^^^
}

double calc_scalefactor_2015castillo_gen(double w, double wsp, double pref)
{
	double fac = 1.0 - (w/wsp) * (w/wsp);
	return pref*fac*fac;
}

double calc_scalefactor_exponential(double w, double wsp, double pref)
{
	double fac = pow(M_E, -w/wsp);
	return pref*fac;
}
