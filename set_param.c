#include "set_param.h"

#define LARGE 10000
#define BYTE_TO_MB 1.0e-06

void get_option(t_qme *qme, FILE *f_opt)
{
	char *line = malloc(sizeof(char) * LARGE);

	load_defaults(qme);

	// read key parameters
	while(fscanf(f_opt, "%s", line) > 0) {
		upper2low(line);
		if(!strcmp(line, "$jobtype")){;
			fscanf(f_opt, "%s", line);
			if(!strcmp(line, "fret_pop")) {
				qme->jtype = 0;
				qme->bEXCI = false;
			} else if(!strcmp(line, "mrt_pop")) {
				qme->jtype = 1;
				qme->bEXCI = true;
			} else if(!strcmp(line, "fret_diss")) {
				qme->jtype = 10;
				qme->bDISS = true;
				qme->bEXCI = false;
			} else if(!strcmp(line, "mrt_diss")) {
				qme->jtype = 11;
				qme->bDISS = true;
				qme->bEXCI = true;
			} else if(!strcmp(line, "fret_trace")) {
				qme->jtype = 20;
				qme->bDEBUG = true;
				qme->bEXCI = false;
			} else if(!strcmp(line, "mrt_trace")) {
				qme->jtype = 21;
				qme->bDEBUG = true;
				qme->bEXCI = true;
			} else if(!strcmp(line, "fret_split")) {
				qme->jtype = 30;
				qme->bEXCI = false;
				qme->bDISS = false;
				qme->bSPLIT = true;
			} else if(!strcmp(line, "mrt_split")) {
				qme->jtype = 31;
				qme->bEXCI = true;
				qme->bDISS = false;
				qme->bSPLIT = true;
			} else {
				fprintf(stderr, "INPUT ERROR: unrecognized jobtype %s\n", line);
				fprintf(stderr, "Allowed jobtypes are: fret_pop, mrt_pop, fret_diss, mrt_diss, fret_trace, mrt_trace, fret_split, mrt_split\n");
			}
		}

		if(!strcmp(line, "$unit")) {
			fscanf(f_opt, "%s", line);
			if(!strcmp(line, "atomic")) qme->unit = 0;
			else if(!strcmp(line, "waveno_fs")) qme->unit = 1;
			else {
				fprintf(stderr, "INPUT ERROR: unrecognized unit %s\n", line);
				fprintf(stderr, "Allowed units are: atomic, waveno_fs\n");
				exit(EXIT_FAILURE);
			}
		}

		if(!strcmp(line, "$dt_lns")) {
			fscanf(f_opt, "%lf", &(qme->dt1));
		}
		if(!strcmp(line, "$dtout_lns")) {
			fscanf(f_opt, "%lf", &(qme->dtout1));
		}
		if(!strcmp(line, "$t_lns")) {
			fscanf(f_opt, "%lf", &(qme->T1));
		}
		if(!strcmp(line, "$dt_prop")) {
			fscanf(f_opt, "%lf", &(qme->dt2));
		}
		if(!strcmp(line, "$dtout_prop")) {
			fscanf(f_opt, "%lf", &(qme->dtout2));
		}
		if(!strcmp(line, "$t_prop")) {
			fscanf(f_opt, "%lf", &(qme->T2));
		}
	}
	rewind(f_opt);

	if(qme->unit){
		qme->dt1     /= TUNIT_TO_FS;
		qme->dtout1  /= TUNIT_TO_FS;
		qme->T1      /= TUNIT_TO_FS;
		qme->dt2     /= TUNIT_TO_FS;
		qme->dtout2  /= TUNIT_TO_FS;
		qme->T2      /= TUNIT_TO_FS;
	}

	get_option_sys(qme, f_opt);
	get_option_spd(qme, f_opt);
	get_option_init(qme, f_opt);

	if(qme->bDEBUG) {
		get_option_trace(qme, f_opt);
	}

	if(qme->bSPLIT) {
		get_option_split(qme, f_opt);
	}

	fprintf(stderr, "allocated memory:   %.1lf MB\n", (double)qme->mem * BYTE_TO_MB);

	free(line);
}

void get_option_sys(t_qme *qme, FILE *f_opt)
{
	int i, j;

	char *line = malloc(sizeof(char) * LARGE);

	load_defaults_sys(qme);

	while(fscanf(f_opt, "%s", line) > 0) {
		upper2low(line);
		if(!strcmp(line, "$nsys")) {
			fscanf(f_opt, "%d", &qme->nsys);
			qme->nsyssq = qme->nsys * qme->nsys;
		}
	}
	rewind(f_opt);

	int nsys   = qme->nsys;
	int nsyssq = qme->nsyssq;

	/* Electronic Hamiltonian and density matrices */
	/* The diagonal elements are vertical transition energies at ground state minimum */
	scalloc(qme->Hsys,  nsyssq, qme->mem);
	smalloc(qme->Heval, nsys,   qme->mem);
	smalloc(qme->Hevec, nsyssq, qme->mem);

	double fac = WAVENO_TO_EUNIT;

	while(fscanf(f_opt, "%s", line) > 0) {
		upper2low(line);
		if(!strcmp(line, "$hsys")) {
			for(i=0; i<nsys; i++) {
				for(j=0; j<i+1; j++) {
					// Read lower triangular part and symmetrize
					fscanf(f_opt, "%lf", &qme->Hsys[i + nsys*j]);
					qme->Hsys[j + nsys*i] = qme->Hsys[i + nsys*j];
				}
				fgets(line, 255, f_opt);
			}
			if(qme->unit) {
				dscal(&nsyssq, &fac, qme->Hsys, &one);
			}
		}
	}
	rewind(f_opt);

	diag(qme->Hevec, qme->Heval, qme->Hsys, nsys);

// vvvvv
#if 0
	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			fprintf(stderr, "%10.3le", qme->Hsys[i + nsys*j]);
		}
		fprintf(stderr, "\n");
	}
#endif

#if 0
	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			fprintf(stderr, "%10.3le", qme->Hevec[i + nsys*j]);
		}
		fprintf(stderr, "\n");
	}
#endif
// ^^^^^

	free(line);
}

void get_option_spd(t_qme *qme, FILE *f_opt)
{
	int i, j;

	char *line  = malloc(sizeof(char) * LARGE);
	char *fname = malloc(sizeof(char) * LARGE);

	while(fscanf(f_opt, "%s", line) > 0) {
		upper2low(line);
		if(!strcmp(line, "$nspd_ana")) {
			fscanf(f_opt, "%d", &qme->nspd_ana);
		}
		if(!strcmp(line, "$nspd_num")) {
			fscanf(f_opt, "%d", &qme->nspd_num);
		}
	}
	rewind(f_opt);

	int nspd_ana = qme->nspd_ana;
	int nspd_num = qme->nspd_num;
	qme->nspd    = nspd_ana + nspd_num;
	int nspd     = qme->nspd;
	int nsyssq   = qme->nsyssq;

	smalloc(qme->nosc, nspd, qme->mem);
	smalloc(qme->osc,  nspd, qme->mem);

	while(fscanf(f_opt, "%s", line) > 0) {
		upper2low(line);

		if(!strcmp(line, "$spd_ana_files")) {
			for(i=0; i<nspd_ana; i++) {
				int crit = 0;

				fscanf(f_opt, "%s", fname);
				FILE *f_spd_ana = fopen(fname, "r");

				while(fscanf(f_spd_ana, "%s", line) > 0) {
					if(!strcmp(line, "-spdtype")) {
						fscanf(f_spd_ana, "%s", line);
						if(!strcmp(line, "drude")) {
							qme->nosc[i] = const_drude_spd(f_spd_ana, fname, qme, &qme->osc[i]);
						} else if(!strcmp(line, "log-normal")) {
							qme->nosc[i] = const_lognormal_spd(f_spd_ana, fname, qme, &qme->osc[i]);
						} else if(!strcmp(line, "brownian")) {
							qme->nosc[i] = const_brownian_spd(f_spd_ana, fname, qme, &qme->osc[i]);
						} else if(!strcmp(line, "single-osc")) {
							qme->nosc[i] = const_singleosc_spd(f_spd_ana, fname, qme, &qme->osc[i]);
						} else {
							fprintf(stderr, "INPUT ERROR: unrecognized spdtype %s\n", line);
							fprintf(stderr, "Allowed spdtypes are: drude, log-normal, brownian and single-osc\n");
						}
						crit++;
					}
				}

				fclose(f_spd_ana);

				if(!crit) {
					fprintf(stderr, "No spdtype detected in %s\n", fname);
					exit(EXIT_FAILURE);
				}
			}
		}

		if(!strcmp(line, "$spd_num_files")) {
			for(i=nspd_ana; i<nspd; i++) {
				fscanf(f_opt, "%s", fname);
				FILE *f_spd_num = fopen(fname, "r");

				while(fscanf(f_spd_num, "%s", line) > 0) {
					qme->nosc[i] = const_num_spd(f_spd_num, fname, &qme->osc[i], qme->nsys);
				}

				fclose(f_spd_num);
			}
		}
	}
	rewind(f_opt);

	free(line);
	free(fname);
}

void get_option_init(t_qme *qme, FILE *f_opt)
{
	int jtype = qme->jtype;
	int nsys  = qme->nsys;

	char *line  = malloc(sizeof(char) * LARGE);

	load_defaults_init(qme);

	while(fscanf(f_opt, "%s", line) > 0) {
		upper2low(line);

		if(!strcmp(line, "$initexc")) {
			fscanf(f_opt, "%d", &qme->initexc);
			if(qme->initexc <= 0 || qme->initexc > nsys) {
				fprintf(stderr, "INPUT ERROR: initexc must be a positive number less than or equal to nsys\n");
				exit(EXIT_FAILURE);
			}
			qme->initexc--;
		}

		if(!strcmp(line, "$initcond")) {
			fscanf(f_opt, "%s", line);
			if(!strcmp(line, "site")) {
				qme->initcond = 0;
			} else if(!strcmp(line, "exci")) {
				qme->initcond = 1;
			} else if(!strcmp(line, "custom")) {
				qme->initcond = 2;
			} else {
				fprintf(stderr, "INPUT ERROR: unrecognized initcond %s\n", line);
				fprintf(stderr, "Allowed initcond: site, exci, custom\n");
				exit(EXIT_FAILURE);
			}
		}
	}
	rewind(f_opt);

	if(qme->bINCO) get_option_init_pop(qme, f_opt);
	else get_option_init_rho(qme, f_opt);

	free(line);
}

void get_option_init_pop(t_qme *qme, FILE *f_opt)
{
	int i, j;

	char *line  = malloc(sizeof(char) * LARGE);

	int nsys   = qme->nsys;
	int nsyssq = qme->nsyssq;
	int ie     = qme->initexc;
	int ic     = qme->initcond;

	scalloc(qme->pop0, nsys, qme->mem);

	double *rho = calloc(nsyssq, sizeof(double));
	rho[ie + nsys*ie] = 1.0;

	if(ic == 0) {
		if(qme->bEXCI) {
			trans(qme->Hevec, rho, nsys);
			for(i=0; i<nsys; i++) {
				qme->pop0[i] = rho[i + nsys*i];
			}
		} else {
			qme->pop0[ie] = 1.0;
		}
	} else if(ic == 1) {
		if(qme->bEXCI) {
			qme->pop0[ie] = 1.0;
		} else {
			trans_inv(qme->Hevec, rho, nsys);
			for(i=0; i<nsys; i++) {
				qme->pop0[i] = rho[i + nsys*i];
			}
		}
	} else if(ic == 2) {
		while(fscanf(f_opt, "%s", line) > 0) {
			upper2low(line);

			if(!strcmp(line, "$pop0_custom")) {
				for(i=0; i<nsys; i++) {
					fscanf(f_opt, "%lf", &qme->pop0[i]);
				}
			}
		}

		rewind(f_opt);
	}

	free(line);
	free(rho);
}

void get_option_init_rho(t_qme *qme, FILE *f_opt)
{
	int i, j;

	char *line = malloc(sizeof(char) * LARGE);

	int nsys   = qme->nsys;
	int nsyssq = qme->nsyssq;
	int ie     = qme->initexc;
	int ic     = qme->initcond;

	smalloc(qme->rho0, nsyssq, qme->mem);
	z_init(qme->rho0, nsyssq);

	if(ic == 0) {
		qme->rho0[ie + nsys * ie].real = 1.0;
	} else if(ic == 1) {
		// Unfinished
	} else if(ic == 2) {
		while(fscanf(f_opt, "%s", line) > 0) {
			upper2low(line);

			if(!strcmp(line, "$rho0_custom_real")) {
				for(i=0; i<nsys; i++) {
					for(j=0; j<i+1; j++) {
						fscanf(f_opt, "%lf%lf", &qme->rho0[i + nsys*j].real);
					}
				}
			}

			if(!strcmp(line, "$rho0_custom_imag")) {
				for(i=0; i<nsys; i++) {
					for(j=0; j<i+1; j++) {
						fscanf(f_opt, "%lf%lf", &qme->rho0[i + nsys*j].imag);
					}
				}
			}
		}

		for(i=0; i<nsys; i++) {
			for(j=0; j<i; j++) {
				qme->rho0[j + nsys*i].real = qme->rho0[i + nsys*j].real;
				qme->rho0[j + nsys*i].imag = qme->rho0[i + nsys*j].imag;
			}
		}

		rewind(f_opt);
	}

	free(line);
}

void get_option_trace(t_qme *qme, FILE *f_opt)
{
	int i, j;

	char *line = malloc(sizeof(char) * LARGE);

	load_defaults_trace(qme, false);

	while(fscanf(f_opt, "%s", line) > 0) {
		upper2low(line);

		if(!strcmp(line, "$nosc_tr")) {
			fscanf(f_opt, "%d", &qme->nosc_tr);
		}
	}
	rewind(f_opt);

	int nosc = qme->nosc_tr;

	smalloc(qme->spdn_tr, nosc, qme->mem);
	smalloc(qme->oscn_tr, nosc, qme->mem);
	smalloc(qme->cut_tr,  nosc, qme->mem);

	load_defaults_trace(qme, true);

	while(fscanf(f_opt, "%s", line) > 0) {
		upper2low(line);

		if(!strcmp(line, "$osc_tr")) {
			for(i=0; i<nosc; i++) {
				fscanf(f_opt, "%d %d %lf\n", &qme->spdn_tr[i], &qme->oscn_tr[i], &qme->cut_tr[i]);
				qme->spdn_tr[i]--;
				qme->oscn_tr[i]--;
			}
		}
	}
	rewind(f_opt);

	free(line);
}

void get_option_split(t_qme *qme, FILE *f_opt)
{
	int i, j;
	int buf;
	int rand_start;

	char *line = malloc(sizeof(char) * LARGE);

	qme->split = malloc(sizeof(t_split));
	t_split *split = qme->split;

	load_defaults_split(qme);

	while(fscanf(f_opt, "%s", line) > 0) {
		upper2low(line);

		if(!strcmp(line, "$ntraj")) {
			fscanf(f_opt, "%d", &(split->ntraj));
		}

		if(!strcmp(line, "$freq_split")) {
			fscanf(f_opt, "%lf", &(split->freq_split));
		}

		if(!strcmp(line, "$pref_split")) {
			fscanf(f_opt, "%lf", &(split->pref_split));
			if(split->pref_split > 1.0) {
				fprintf(stderr, "Error: pref_split must be smaller than 1.0\n");
				exit(EXIT_FAILURE);
			}
		}

		if(!strcmp(line, "$fcn_split")) {
			fscanf(f_opt, "%d", &(split->type_fcn));
			if((split->type_fcn < 0) || (split->type_fcn > 1)) {
				fprintf(stderr, "Error: type_fcn must be 0 or 1\n");
				exit(EXIT_FAILURE);
			}
		}

		if(!strcmp(line, "$rand_start")) {
			fscanf(f_opt, "%d", &(split->rand_start));
		}

		if(!strcmp(line, "$calc_diss")) {
			fscanf(f_opt, "%d", &buf);
			qme->bDISS = (bool)buf;
		}

		if(!strcmp(line, "$hsys_adj")) {
			fscanf(f_opt, "%d", &buf);
			split->bHSYS_ADJUST = (bool)buf;
		}
	}
	rewind(f_opt);

	if(qme->unit) {
		split->freq_split *= WAVENO_TO_EUNIT;
	}

	free(line);
}

void load_defaults(t_qme *qme)
{
	qme->jtype  = 0;
	qme->bINCO  = true;
	qme->bDISS  = false;
	qme->bDEBUG = false;
	qme->bSPLIT = false;
	qme->unit   = 1;

	qme->dt1    = 0.5;
	qme->dtout1 = 10.0;
	qme->T1     = 50000.0;

	qme->dt2    = 0.5;
	qme->dtout2 = 10.0;
	qme->T2     = 50000.0;

	qme->mem   = 0.0;
}

void load_defaults_sys(t_qme *qme)
{
	qme->nsys     = 2;
	qme->initexc  = 1;
	qme->initcond = 0;
}

void load_defaults_init(t_qme *qme)
{
	qme->initexc = 0;
	qme->initcond = 0;
}

void load_defaults_trace(t_qme *qme, bool bFLAG)
{
	if(!bFLAG) {
		qme->nosc_tr = 1;
	} else {
		qme->spdn_tr[0] = 0;
		qme->oscn_tr[0] = qme->nosc[0] - 1;
		qme->cut_tr[0]  = 1.0e-06;
	}
}

void load_defaults_split(t_qme *qme)
{
	qme->split->bDISS = false;
	qme->split->bHSYS_ADJUST = false;

	qme->split->ntraj       = 100;
	qme->split->freq_split  = 1.0e-15;
	qme->split->pref_split  = 0.99;
	qme->split->type_fcn    = 0;
	qme->split->rand_start  = 0;
}
