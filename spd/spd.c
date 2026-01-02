#include "spd.h"

t_osc *osc_init(t_qme *qme, int nosc)
{
	int i;
	int nsyssq = qme->nsyssq;

	t_osc *osc = malloc(nosc * sizeof(t_osc));

	for(i=0; i<nosc; i++) {
		osc[i].gamma = malloc(nsyssq * sizeof(double));
	}

	if(qme->bDISS) {
		for(i=0; i<nosc; i++) {
			osc[i].Idiss = calloc(nsyssq, sizeof(double));
			osc[i].Kdiss = calloc(nsyssq, sizeof(double));
			osc[i].Jdiss = calloc(nsyssq, sizeof(double));
		}
	}

	return osc;
}

void osc_copy(t_osc *osc_src, t_osc *osc_dest, int nsyssq)
{
	int i, j;
	int one = 1;

	osc_dest->freq   = osc_src->freq;
	osc_dest->spread = osc_src->spread;
	osc_dest->temp   = osc_src->temp;
	osc_dest->coth   = osc_src->coth;
}

int const_drude_spd(FILE *f_spd, char *f_name, t_qme *qme, t_osc **osc)
{
	int i, j, k;
	int nosc;
	int crit = 0;
	double reorg, wc, wmax, w;
	double rho, J;
	double gamma_raw, T;
	double reorg_sum = 0.0;

	int nsys   = qme->nsys;
	int nsyssq = qme->nsyssq;

	double *coupl = malloc(sizeof(double) * nsyssq);
	char *line    = malloc(sizeof(char) * LARGE);

	while(fscanf(f_spd, "%s", line) > 0) {
		upper2low(line);
		if(!strcmp(line, "-nosc")) {
			fscanf(f_spd, "%d", &nosc);
			crit++;
		}
		if(!strcmp(line, "-reorg")) {
			fscanf(f_spd, "%lf", &reorg);
			crit++;
		}
		if(!strcmp(line, "-w_cutoff")) {
			fscanf(f_spd, "%lf", &wc);
			crit++;
		}
		if(!strcmp(line, "-w_max")) {
			fscanf(f_spd, "%lf", &wmax);
			crit++;
		}
		if(!strcmp(line, "-temp")) {
			fscanf(f_spd, "%lf", &T);
			crit++; 
		}
		if(!strcmp(line, "-coupl")) {
			for(i=0; i<nsys; i++) {
				for(j=0; j<i+1; j++) {
					fscanf(f_spd, "%lf", &coupl[i + nsys*j]);
					coupl[j + nsys*i] = coupl[i + nsys*j];
				}
			}
			crit++;
		}
	}

	if(crit != 6) {
		fprintf(stderr, "Missing information in %s\n", f_name);
		fprintf(stderr, "Required keywords are: nosc, reorg, w_cutoff, w_max, temp, coupl\n");
		exit(EXIT_FAILURE);
	}

	*osc = osc_init(qme, nosc);

	if(qme->unit) {
		reorg *= WAVENO_TO_EUNIT;
		wc    *= WAVENO_TO_EUNIT;
		wmax  *= WAVENO_TO_EUNIT;
	}

	for(i=0; i<nosc; i++) {
		w = wmax * pow((double)(i+1) / (double)nosc, 2.0);
		rho = (double)nosc / (2.0 * sqrt(w * wmax));
		J   = (2.0 * reorg / M_PI) * w * wc / (wc * wc + w * w);
		gamma_raw = sqrt(J / rho) / w;
		reorg_sum += w * gamma_raw * gamma_raw;

		(*osc)[i].freq   = w;
		(*osc)[i].spread = 2.0 * w / (double)(i+1);
		(*osc)[i].temp   = T;
		(*osc)[i].coth   = 1.0 / tanh(w / (2.0 * T));

		for(j=0; j<nsyssq; j++) {
			(*osc)[i].gamma[j] = gamma_raw * coupl[j];
		}
	}

// vvvvv
#if 0
	for(i=0; i<20; i++) {
		for(j=0; j<nsyssq; j++) {
			fprintf(stderr, "%13.6le ", (*osc)[i].gamma[j]);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "\n");
#endif
// ^^^^^

	fprintf(stderr, "--------------------------------------------------\n"); 
	fprintf(stderr, "Spectral density file %s - Drude-Lorentz\n", f_name);
	fprintf(stderr, "Reorg E pristine  %15.6le a.u.   %15.6le waveno.\n", reorg, reorg / WAVENO_TO_EUNIT);
	fprintf(stderr, "Reorg E osc sum   %15.6le a.u.   %15.6le waveno.  (%6.2lf \%)\n", reorg_sum, reorg_sum / WAVENO_TO_EUNIT, reorg_sum/reorg*100.0);
	fprintf(stderr, "w_cutoff          %15.6le a.u.   %15.6le waveno.\n", wc, wc / WAVENO_TO_EUNIT);
	fprintf(stderr, "w_max             %15.6le a.u.   %15.6le waveno.\n", wmax, wmax / WAVENO_TO_EUNIT);
	fprintf(stderr, "Coupling\n");
	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			fprintf(stderr, "%8.3lf", coupl[i + nsys*j]);
		}
		fprintf(stderr, "\n");
	}

	free(coupl);
	free(line);

	return nosc;
}

int const_lognormal_spd(FILE *f_spd, char *f_name, t_qme *qme, t_osc **osc)
{
	// Buggy
	int i, j;
	int nosc;
	int crit = 0;
	double reorg, std, wc, wmax, w;
	double temp1, temp2;
	double gamma_raw, T;
	double reorg_sum = 0.0;

	int nsys   = qme->nsys;
	int nsyssq = qme->nsyssq;

	double *coupl = malloc(sizeof(double) * nsyssq);
	char *line    = malloc(sizeof(char) * LARGE);

	while(fscanf(f_spd, "%s", line) > 0) {
		upper2low(line);
		if(!strcmp(line, "-nosc")) {
			fscanf(f_spd, "%d", &nosc);
			crit++;
		}
		if(!strcmp(line, "-reorg")) {
			fscanf(f_spd, "%lf", &reorg);
			crit++;
		}
		if(!strcmp(line, "-w_cutoff")) {
			fscanf(f_spd, "%lf", &wc);
			crit++;
		}
		if(!strcmp(line, "-w_max")) {
			fscanf(f_spd, "%lf", &wmax);
			crit++;
		}
		if(!strcmp(line, "-std")) {
			fscanf(f_spd, "%lf", &std);
			crit++;
		}
		if(!strcmp(line, "-temp")) {
			fscanf(f_spd, "%lf", &T);
			crit++; 
		}
		if(!strcmp(line, "-coupl")) {
			for(i=0; i<nsys; i++) {
				for(j=0; j<i+1; j++) {
					fscanf(f_spd, "%lf", &coupl[i + nsys*j]);
					coupl[j + nsys*i] = coupl[i + nsys*j];
				}
			}
			crit++;
		}
	}

	if(crit != 7) {
		fprintf(stderr, "Missing information in %s\n", f_name);
		fprintf(stderr, "Required keywords are: nosc, reorg, std, w_cutoff, w_max, temp, coupl\n");
		exit(EXIT_FAILURE);
	}

	*osc = osc_init(qme, nosc);

	if(qme->unit) {
		reorg *= WAVENO_TO_EUNIT;
		wc    *= WAVENO_TO_EUNIT;
		wmax  *= WAVENO_TO_EUNIT;
	}

	double K = reorg / wc / pow(M_E, std * std / 2.0) / std / sqrt(2.0 * M_PI);

	for(i=0; i<nosc; i++) {
		w = wmax * pow((double)(i+1) / (double)nosc, 2.0);
		temp1 = 2.0 * K / (double)nosc * sqrt(wmax * w);
		temp2 = -pow((log(w / wc) / std), 2.0) / 2.0;

		(*osc)[i].freq   = w;
		(*osc)[i].spread = 2.0 * w / (double)(i+1);
		(*osc)[i].temp   = T;
		(*osc)[i].coth   = 1.0 / tanh(w / (2.0 * T));

		gamma_raw = sqrt(temp1 * pow(M_E, temp2) / w);
		reorg_sum += w * gamma_raw * gamma_raw;

		for(j=0; j<nsyssq; j++) {
			(*osc)[i].gamma[j] = gamma_raw * coupl[j];
		}
	}

	fprintf(stderr, "--------------------------------------------------\n"); 
	fprintf(stderr, "Spectral density file %s - log-normal\n", f_name);
	fprintf(stderr, "Reorg E pristine  %15.6le a.u.   %15.6le waveno.\n", reorg, reorg / WAVENO_TO_EUNIT);
	fprintf(stderr, "Reorg E osc sum   %15.6le a.u.   %15.6le waveno.  (%6.2lf \%)\n", reorg_sum, reorg_sum / WAVENO_TO_EUNIT, reorg_sum/reorg*100.0);
	fprintf(stderr, "w_cutoff          %15.6le a.u.   %15.6le waveno.\n", wc, wc / WAVENO_TO_EUNIT);
	fprintf(stderr, "w_max             %15.6le a.u.   %15.6le waveno.\n", wmax, wmax / WAVENO_TO_EUNIT);
	fprintf(stderr, "standard dev.     %15.6le\n", std);
	fprintf(stderr, "Coupling\n");
	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			fprintf(stderr, "%8.3lf", coupl[i + nsys*j]);
		}
		fprintf(stderr, "\n");
	}

	free(coupl);
	free(line);

	return nosc;
}

int const_brownian_spd(FILE *f_spd, char *f_name, t_qme *qme, t_osc **osc)
{
	int i, j;
	int nosc, nosc_half, index;
	int crit = 0;
	double reorg, w0, w1, wc, wmax, w;
	double wcsq, w0sq, wsq;
	double rho, J;
	double gamma_raw, T;
	double reorg_sum = 0.0;

	int nsys   = qme->nsys;
	int nsyssq = qme->nsyssq;

	double *coupl = malloc(sizeof(double) * nsyssq);
	char *line    = malloc(sizeof(char) * LARGE);

	while(fscanf(f_spd, "%s", line) > 0) {
		upper2low(line);
		if(!strcmp(line, "-nosc")) {
			fscanf(f_spd, "%d", &nosc);
			nosc -= (nosc % 2);
			nosc_half = nosc / 2;
			crit++;
		}
		if(!strcmp(line, "-reorg")) {
			fscanf(f_spd, "%lf", &reorg);
			crit++;
		}
		if(!strcmp(line, "-w0")) {
			fscanf(f_spd, "%lf", &w0);
			crit++;
		}
		if(!strcmp(line, "-w_cutoff")) {
			fscanf(f_spd, "%lf", &wc);
			crit++;
		}
		if(!strcmp(line, "-w_max")) {
			fscanf(f_spd, "%lf", &wmax);
			crit++;
		}
		if(!strcmp(line, "-temp")) {
			fscanf(f_spd, "%lf", &T);
			crit++; 
		}
		if(!strcmp(line, "-coupl")) {
			for(i=0; i<nsys; i++) {
				for(j=0; j<i+1; j++) {
					fscanf(f_spd, "%lf", &coupl[i + nsys*j]);
					coupl[j + nsys*i] = coupl[i + nsys*j];
				}
			}
			crit++;
		}
	}

	if(crit != 7) {
		fprintf(stderr, "Missing information in %s\n", f_name);
		fprintf(stderr, "Required keywords are: nosc, reorg, w0, w_cutoff, w_max, temp, coupl\n");
		exit(EXIT_FAILURE);
	}

	if(w0 > wmax) {
		fprintf(stderr, "For Brownian oscillator, w0 must be smaller than wmax in %s\n", f_name);
		exit(EXIT_FAILURE);
	}

	*osc = osc_init(qme, nosc);

	if(qme->unit) {
		reorg *= WAVENO_TO_EUNIT;
		w0    *= WAVENO_TO_EUNIT;
		wc    *= WAVENO_TO_EUNIT;
		wmax  *= WAVENO_TO_EUNIT;
	}

	w0sq = w0 * w0;
	wcsq = wc * wc;

	// Determine whether the range is whole or split
	if(w0sq < 2.0 * wcsq) { // whole
		for(i=0; i<nosc; i++) {
			w   = wmax * pow((double)(i+1) / (double)nosc, 2.0);
			wsq = w * w;
			rho = (double)nosc / (2.0 * sqrt(w * wmax));
			J   = (2.0 * reorg * wc / M_PI) * (2.0 * w0sq * w / ((wsq - w0sq) * (wsq - w0sq) + 4.0 * wcsq * wsq));

			gamma_raw = sqrt(J / rho) / w;
			reorg_sum += w * gamma_raw * gamma_raw;

			(*osc)[i].freq   = w;
			(*osc)[i].spread = 2.0 * w / (double)(i+1);
			(*osc)[i].temp   = T;
			(*osc)[i].coth   = 1.0 / tanh(w / (2.0 * T));

			for(j=0; j<nsyssq; j++) {
				(*osc)[i].gamma[j] = gamma_raw * coupl[j];
			}
		}
	} else { // split
		w1 = sqrt(w0sq - 2.0 * wcsq);

		// Lower half
		for(i=0; i<nosc_half-1; i++) {
			w   = (1.0 - pow(1.0 - (double)(i+1) / (double)nosc_half, 2.0)) * w1;
			wsq = w * w;
			rho = (double)nosc_half / (2.0 * sqrt((w1 - w) * w1));
			J   = (2.0 * reorg * wc / M_PI) * (2.0 * w0sq * w / ((wsq - w0sq) * (wsq - w0sq) + 4.0 * wcsq * wsq));

			gamma_raw = sqrt(J / rho) / w;
			reorg_sum += w * gamma_raw * gamma_raw;

			(*osc)[i].freq   = w;
			(*osc)[i].spread = (4.0 * (double)nosc - 8.0 * (double)(i+1)) / ((double)nosc * (double)nosc) * w1;
			(*osc)[i].temp   = T;
			(*osc)[i].coth   = 1.0 / tanh(w / (2.0 * T));

//			fprintf(stderr, "%13.6le %13.6le\n", w * w * gamma_raw * gamma_raw / (*osc)[i].spread, J);

			for(j=0; j<nsyssq; j++) {
				(*osc)[i].gamma[j] = gamma_raw * coupl[j];
			}
		}

		// w = w1
		w   = w1;
		wsq = w * w;
		J   = (2.0 * reorg * wc / M_PI) * (2.0 * w0sq * w / ((wsq - w0sq) * (wsq - w0sq) + 4.0 * wcsq * wsq));

		double reorg_w1 = reorg * w0sq / (M_PI * wc * (w0sq - wcsq)); // reorganization energy density at w1
		double spread   = 2.0 * wmax / ((double)nosc * (double)nosc);
		reorg_w1 *= spread; // reorganization energy allocated for the oscillator at w = w1
		reorg_sum += reorg_w1;
		gamma_raw = sqrt(reorg_w1 / w1);

		(*osc)[nosc_half-1].freq   = w1;
		(*osc)[nosc_half-1].spread = spread;
		(*osc)[nosc_half-1].temp   = T;
		(*osc)[nosc_half-1].coth   = 1.0 / tanh(w / (2.0 * T));

//		fprintf(stderr, "%13.6le %13.6le\n", w * w * gamma_raw * gamma_raw / (*osc)[nosc_half-1].spread, J);

		for(j=0; j<nsyssq; j++) {
			(*osc)[nosc_half-1].gamma[j] = gamma_raw * coupl[j];
		}

		// Upper half
		for(i=nosc_half; i<nosc; i++) {
			int I = i - nosc_half;
			w = w1 + (wmax - w1) * pow((double)(I+1) / (double)nosc_half, 2.0);
			wsq = w * w;
			rho = (double)nosc_half / (2.0 * sqrt((w - w1) * (wmax - w1)));
			J   = (2.0 * reorg * wc / M_PI) * (2.0 * w0sq * w / ((wsq - w0sq) * (wsq - w0sq) + 4.0 * wcsq * wsq));

			gamma_raw = sqrt(J / rho) / w;
			reorg_sum += w * gamma_raw * gamma_raw;

			(*osc)[i].freq   = w;
			(*osc)[i].spread = (8.0 * (double)(I+1)) / ((double)nosc * (double)nosc) * (wmax - w1);
			(*osc)[i].temp   = T;
			(*osc)[i].coth   = 1.0 / tanh(w / (2.0 * T));

//			fprintf(stderr, "%13.6le %13.6le\n", w * w * gamma_raw * gamma_raw / (*osc)[i].spread, J);

			for(j=0; j<nsyssq; j++) {
				(*osc)[i].gamma[j] = gamma_raw * coupl[j];
			}
		} 
	}

	fprintf(stderr, "--------------------------------------------------\n"); 
	fprintf(stderr, "Spectral density file %s - Brownian oscillator\n", f_name);
	fprintf(stderr, "Reorg E pristine  %15.6le a.u.   %15.6le waveno.\n", reorg, reorg / WAVENO_TO_EUNIT);
	fprintf(stderr, "Reorg E osc sum   %15.6le a.u.   %15.6le waveno.  (%6.2lf \%)\n", reorg_sum, reorg_sum / WAVENO_TO_EUNIT, reorg_sum/reorg*100.0);
	fprintf(stderr, "w0                %15.6le a.u.   %15.6le waveno.\n", w0, w0 / WAVENO_TO_EUNIT);
	fprintf(stderr, "w_cutoff          %15.6le a.u.   %15.6le waveno.\n", wc, wc / WAVENO_TO_EUNIT);
	fprintf(stderr, "w_max             %15.6le a.u.   %15.6le waveno.\n", wmax, wmax / WAVENO_TO_EUNIT);
	fprintf(stderr, "Coupling\n");
	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			fprintf(stderr, "%8.3lf", coupl[i + nsys*j]);
		}
		fprintf(stderr, "\n");
	}

	free(coupl);
	free(line);

	return nosc;
}

int const_singleosc_spd(FILE *f_spd, char *f_name, t_qme *qme, t_osc **osc)
{

}

int const_num_spd(FILE *f_spd, char *f_name, t_qme *qme, t_osc **osc)
{

}
