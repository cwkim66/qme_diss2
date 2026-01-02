#include "init_tr.h"

void init_tr(t_qme *qme)
{
	int i, j;

	qme->tr = malloc(sizeof(t_tr));

	t_tr *tr = qme->tr;

	tr->jtype = qme->jtype;
	tr->unit  = qme->unit;

	tr->dtout = qme->dtout1;
	tr->T     = qme->T1;

	tr->nout = dround(tr->T / tr->dtout) + 1;

	tr->nsys    = qme->nsys;
	tr->nsyssq  = qme->nsyssq;
	tr->nosc_tr = qme->nosc_tr;

	init_osc_tr(qme);

	int nosc   = tr->nosc_tr;
	int nout   = tr->nout;
	int nsyssq = tr->nsyssq;

	smalloc(tr->tr1, nosc, qme->mem);
	smalloc(tr->tr2, nosc, qme->mem);

	for(i=0; i<nosc; i++) {
		smalloc(tr->tr1[i], nout, qme->mem);
		smalloc(tr->tr2[i], nout, qme->mem);

		for(j=0; j<nout; j++) {
			smalloc(tr->tr1[i][j], nsyssq, qme->mem);
			smalloc(tr->tr2[i][j], nsyssq, qme->mem);
		}
	}
}

void init_osc_tr(t_qme *qme)
{
	int i, j;

	int jtype   = qme->jtype;
	int nosc_tr = qme->nosc_tr;
	t_tr *tr    = qme->tr;

	tr->qosc = malloc(nosc_tr * sizeof(t_qosc));

	t_qosc *qosc = tr->qosc;

	fprintf(stderr, "---- Direct trace evaluation for bath oscillators ----\n"); 

	for(i=0; i<nosc_tr; i++) {
		int spdn = qme->spdn_tr[i];
		int oscn = qme->oscn_tr[i];
		double w = qme->osc[spdn][oscn].freq;
		double T = qme->osc[spdn][oscn].temp;

		if(oscn < 0 || qme->nosc[spdn] <= oscn) {
			fprintf(stderr, "Error: oscn %d of entry no.%d is out of the range for the spd no.%d (1 - %d)\n", oscn+1, i+1, spdn+1, qme->nosc[spdn]);
			exit(EXIT_FAILURE);
		}

		qosc[i].freq = w;
		qosc[i].temp = T;
		qosc[i].coth = qme->osc[spdn][oscn].coth;

		qosc[i].spdn = qme->spdn_tr[i];
		qosc[i].oscn = qme->oscn_tr[i];

		qosc[i].nsys   = qme->nsys;
		qosc[i].nsyssq = qme->nsyssq;
		qosc[i].nqst   = set_nqst(qosc[i], qme->cut_tr[i]);
		qosc[i].nqstsq = qosc[i].nqst * qosc[i].nqst;
		qosc[i].ndim   = qme->nsys * qosc[i].nqst;
		qosc[i].ndimsq = qosc[i].ndim * qosc[i].ndim;

		qosc[i].gamma  = qme->osc[spdn][oscn].gamma;
		qosc[i].gexci  = qme->osc[spdn][oscn].gexci;

		int nsys   = qosc[i].nsys;
		int nsyssq = qosc[i].nsyssq;
		int nqstsq = qosc[i].nqstsq;

		smalloc(qosc[i].aHpa, nqstsq, qme->mem);

		smalloc(qosc[i].FC,  nsyssq, qme->mem);
		smalloc(qosc[i].FCU, nsyssq, qme->mem);
		smalloc(qosc[i].Req, nsyssq, qme->mem);
		smalloc(qosc[i].H,   nsyssq, qme->mem);
		smalloc(qosc[i].U,   nsyssq, qme->mem);
		for(j=0; j<nsyssq; j++) {
			smalloc(qosc[i].FC[j],  nqstsq, qme->mem);
			smalloc(qosc[i].FCU[j], nqstsq, qme->mem);
			smalloc(qosc[i].Req[j], nqstsq, qme->mem);
			smalloc(qosc[i].H[j],   nqstsq, qme->mem);
			smalloc(qosc[i].U[j],   nqstsq, qme->mem);
		}

		fprintf(stderr, "qosc entry %d   spd %3d   oscn %5d   freq %.3le a.u. (%.3le waveno.)   temp %.3le   coth %.3le   nstates %d\n",
				i+1, spdn+1, oscn+1, w, w / WAVENO_TO_EUNIT, T, qosc[i].coth, qosc[i].nqst);

		set_FC(&qosc[i], jtype);
		ortho_FC(&qosc[i]);
		set_Req(&qosc[i]);
		set_H_diag(&qosc[i]);

		if(jtype == 21) {
			set_H_offdiag(&qosc[i]);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "\n");
}

double set_nqst(t_qosc qosc, double cut)
{
	// Determine number of harmonic oscillator eigenstates
	int i, j;
	int nqst = 0;

	double w = qosc.freq;
	double T = qosc.temp;

	double pop_sum = 0.0;
	double partf = 1.0 / (2.0 * sinh(w / (2.0 * T)));

	while(pop_sum < 1.0 - cut) {
		double vibE = w * ((double)nqst + 0.5);
		pop_sum += pow(M_E, -vibE / T) / partf;
		nqst++;
	}
	nqst++;

	if(nqst < NMAX_DEF) nqst = NMAX_DEF;

	return nqst;
}

void set_FC(t_qosc *qosc, int jtype)
{
	int i, j, k, l;
	int one = 1;
	int sign;
	double S;

	int nsys   = qosc->nsys;
	int nqst   = qosc->nqst;
	int nqstsq = qosc->nqstsq;

	// Currently only diagonal parts are constructed
	// Skipped the variable j on purpose

	for(i=0; i<nsys; i++) {
		int n      = i + nsys*i;
		double *FC = qosc->FC[n];

		if(jtype == 20) {
			// FRET
			S = qosc->gamma[n] * qosc->gamma[n];
			sign = (qosc->gamma[n] > 0) ? 1 : -1;
		} else if(jtype == 21) {
			// MRT
			S = qosc->gexci[n] * qosc->gexci[n];
			sign = (qosc->gexci[n] > 0) ? 1 : -1;
		}
		fprintf(stderr, "Element %03d / %03d  S %13.6le\n", i+1, i+1, S);

		// <k_g|l_e> = FC[k + nqst*l]: overlap between k-th ground vibrational state and l-th excited vibrational state
		// The center of excited state PES is at -d where the Huang-Rhys factor S = wd^2/hbar

		for(k=0; k<nqst; k++) {
			FC[0 + nqst*k] = pow(M_E, -S / 2.0);

			for(l=1; l<k+1; l++) {
				FC[0 + nqst*k] *= (pow(S, 0.5) / sqrt(l));
				FC[0 + nqst*k] *= ((sign == -1) && (k % 2)) ? -1.0 : 1.0;
			}

			FC[k + nqst*0] = (k % 2) ? -FC[0 + nqst*k] : FC[0 + nqst*k];

			for(l=1; l<k+1; l++) {
				FC[k + nqst*l] = sqrt((double)k / (double)l) * FC[(k-1) + nqst*(l-1)] + sign * sqrt(S / (double)l) * FC[k + nqst*(l-1)];
				FC[l + nqst*k] = ((k+l) % 2) ? -FC[k + nqst*l] : FC[k + nqst*l];
			}
		}
	}

// vvvvv
#if 0
	for(i=0; i<nsys; i++) {
		double *FC = qosc->FC[i + nsys*i];

		for(k=0; k<nqst; k++) {
			for(l=0; l<nqst; l++) {
				fprintf(stderr, "%8.1le ", FC[k + nqst*l]);
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");
	}
#endif
// ^^^^^

}

void ortho_FC(t_qosc *qosc)
{
	// orthonormalize FC coefficients by symmetric orthonormalization
	// FCU = FC * X where X = S^(-1/2) and S is overlap matrix

	int i, j, k, l;
	int one = 1;

	int nsys   = qosc->nsys;
	int nqst   = qosc->nqst;
	int nqstsq = qosc->nqstsq;

	double *G    = malloc(nqstsq * sizeof(double));
	double *S    = malloc(nqstsq * sizeof(double));
	double *Seig = malloc(nqstsq * sizeof(double));
	double *U    = malloc(nqstsq * sizeof(double));
	double *X    = malloc(nqstsq * sizeof(double));

	for(i=0; i<nsys; i++) {
		int n      = i + nsys*i;
		double *FC = qosc->FC[n];

		d_init(X, nqstsq);

		for(k=0; k<nqst; k++) {
			for(l=0; l<nqst; l++) {
				S[k + nqst*l] = ddot(&nqst, &FC[nqst*k], &one, &FC[nqst*l], &one);
			}
		}

		diag(U, Seig, S, nqst);
		for(k=0; k<nqst; k++) X[k + nqst*k] = sqrt(1.0 / Seig[k]);
		trans_inv(U, X, nqst);
		AxB(G, FC, X, nqst, nqst, nqst, 1);

// vvvvv
#if 0
		for(k=0; k<nqst; k++) {
			for(l=0; l<nqst; l++) {
				fprintf(stderr, "%15.8le ", FC[k + nqst*l]);
			}
			fprintf(stderr, "\n");
		}

		fprintf(stderr, "\n");

		for(k=0; k<nqst; k++) {
			for(l=0; l<nqst; l++) {
				fprintf(stderr, "%15.8le ", G[k + nqst*l]);
			}
			fprintf(stderr, "\n");
		}

		fprintf(stderr, "\n");
#endif
// ^^^^^

		complex *FCU = qosc->FCU[n];
		z_init(FCU, nqstsq);
		for(k=0; k<nqstsq; k++) {
			FCU[k].real = G[k];
		}
	}

	free(G);
	free(S);
	free(Seig);
	free(U);
	free(X);
}

void set_Req(t_qosc *qosc)
{
	// Construct equilibrium density matrix for the oscillator

	int i, j, k;

	int nsys   = qosc->nsys;
	int nqst   = qosc->nqst;
	int nqstsq = qosc->nqstsq;

	double w = qosc->freq;
	double T = qosc->temp;

	for(i=0; i<nsys; i++) {
		double pop, pop_sum = 0.0;
		int n = i + nsys*i;

		complex *FCU = qosc->FCU[n];
		complex *Req = qosc->Req[n];
		z_init(Req, nqstsq);

		for(j=0; j<nqst; j++) {
			double vibE = w * ((double)j + 0.5);
			pop = pow(M_E, -vibE / T);
			Req[j + nqst*j].real = pop; 
			pop_sum += pop;
		}

		for(j=0; j<nqst; j++) {
			Req[j + nqst*j].real /= pop_sum;
		}

		// We have created the equilibrium density at the excited state
		// so, to represent it in the ground state basis, we use ztrans_inv
		ztrans_inv(FCU, Req, nqst);
	}

// vvvvv
#if 0
	for(i=0; i<nsys; i++) {
		complex *Req = qosc->Req[i + nsys*i];

		for(j=0; j<nqst; j++) {
			for(k=0; k<nqst; k++) {
				fprintf(stderr, "%13.6le ", Req[j + nqst*k].real);
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");
	}
#endif
// ^^^^^
}

void set_H_diag(t_qosc *qosc)
{
	// Construct vibrational Hamiltonian for diagonal electronic block

	int i, j, k;

	int nsys   = qosc->nsys;
	int nqst   = qosc->nqst;
	int nqstsq = qosc->nqstsq;

	double w = qosc->freq;

	for(i=0; i<nsys; i++) {
		int n = i + nsys*i;

		complex *FCU = qosc->FCU[n];
		complex *H   = qosc->H[n];
		z_init(H, nqstsq);

		for(j=0; j<nqst; j++) {
			H[j + nqst*j].real = w * ((double)j + 0.5);
		}

		// The Hamiltonian is diagonal at the excited state basis
		// so, to represent it in the ground state basis, we use ztrans_inv
		ztrans_inv(FCU, H, nqst);
	}

// vvvvv
#if 0
	for(i=0; i<nqst; i++) {
		for(j=0; j<nqst; j++) {
			fprintf(stderr, "%13.6le ", qosc->H[0][i + nqst*j].real - qosc->H[3][i + nqst*j].real);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "\n");
#endif
// ^^^^^


// vvvvv
#if 0
	for(i=0; i<nsys; i++) {
		complex *H = qosc->H[i + nsys*i];

		for(j=0; j<nqst; j++) {
			for(k=0; k<nqst; k++) {
				fprintf(stderr, "%13.6le ", H[j + nqst*k].real);
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");
	}
#endif
// ^^^^^
}

void set_H_offdiag(t_qosc *qosc)
{
	// Construct vibrational Hamiltonian for off-diagonal electronic block
	// Currently only works for MRT

	int i, j, k, l;
	int one = 1;

	int nsys   = qosc->nsys;
	int nqst   = qosc->nqst;
	int nqstsq = qosc->nqstsq;

	double w = qosc->freq;

	complex *aHpa = qosc->aHpa;
	z_init(aHpa, nqstsq);

	for(i=0; i<nqst-1; i++) {
		aHpa[i + nqst*(i+1)].real = sqrt(i+1);
		aHpa[(i+1) + nqst*i].real = sqrt(i+1);
	}

	for(i=0; i<nsys; i++) {
		for(j=0; j<i; j++) {
			int n       = i + nsys*j;
			double coup = w * qosc->gexci[n];
			complex *H  = qosc->H[n];
			z_init(H, nqstsq);

			for(k=0; k<nqst-1; k++) {
				H[k + nqst*(k+1)].real = coup * sqrt(k+1);
				H[(k+1) + nqst*k].real = coup * sqrt(k+1);
			}

			zcopy(&nqstsq, H, &one, qosc->H[j + nsys*i], &one);
		}
	}

// vvvvv
#if 0
	for(i=0; i<nsys; i++) {
		for(j=0; j<nsys; j++) {
			if(i != j) {
				complex *H = qosc->H[i + nsys*j];

				for(k=0; k<nqst; k++) {
					for(l=0; l<nqst; l++) {
						fprintf(stderr, "%13.6le ", H[k + nqst*l].real);
					}
					fprintf(stderr, "\n");
				}
				fprintf(stderr, "\n");
			}
		}
	}
#endif
// ^^^^^
}

void set_U_diag(t_qosc *qosc, double t)
{
	// Construct vibrational propagator for diagonal electronic block

	int i, j, k;
	double arg;

	int nsys   = qosc->nsys;
	int nqst   = qosc->nqst;
	int nqstsq = qosc->nqstsq;

	double w = qosc->freq;

	for(i=0; i<nsys; i++) {
		int n = i + nsys*i;

		complex *FCU = qosc->FCU[n];
		complex *U   = qosc->U[n];
		z_init(U, nqstsq);

		for(j=0; j<nqst; j++) {
			arg = -w * ((double)j + 0.5) * t;
			U[j + nqst*j].real = cos(arg);
			U[j + nqst*j].imag = sin(arg);
		}

		// The propagator is diagonal at the excited state basis
		// so, to represent it in the ground state basis, we use ztrans_inv
		ztrans_inv(FCU, U, nqst);
	}

// vvvvv
#if 0
	for(i=0; i<nsys; i++) {
		complex *U = qosc->U[i + nsys*i];

		for(j=0; j<nqst; j++) {
			for(k=0; k<nqst; k++) {
				fprintf(stderr, "%13.6le ", U[j + nqst*k].real);
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");
	}
#endif
// ^^^^^
}

void set_U_diag_wick(t_qosc *qosc, double t, bool mode)
{
	// Construct vibrational propagator for diagonal electronic block
	// subject to the wick rotation

	// mode = true:  forward rotation
	// mode = false: backward rotation

	int i, j, k;
	double exp, arg, ener;
	double sign = mode ? -1.0 : 1.0;

	int nsys = qosc->nsys;
	int nqst = qosc->nqst;
	int nqstsq = qosc->nqstsq;

	double w = qosc->freq;
	double T = qosc->temp;

	for(i=0; i<nsys; i++) {
		int n = i + nsys*i;

		complex *FCU = qosc->FCU[n];
		complex *U   = qosc->U[n];
		z_init(U, nqstsq);

		for(j=0; j<nqst; j++) {
			ener = w * ((double)j + 0.5);
			exp = pow(M_E, sign * ener / T);
			arg = ener * t;
			U[j + nqst*j].real = exp * cos(arg);
			U[j + nqst*j].imag = sign * exp * sin(arg);
		}

		// The propagator is diagonal at the excited state basis
		// so, to represent it in the ground state basis, we use ztrans_inv
		ztrans_inv(FCU, U, nqst);
	}
}
