#include "init_ptb2.h"

void init_fret(t_qme *qme, t_fret *fret)
{
	int i, j;
	int one = 1;

	fret->unit   = qme->unit;

	fret->nsys   = qme->nsys;
	fret->nsyssq = qme->nsyssq;

	fret->dt    = qme->dt1;
	fret->dtout = qme->dtout1;
	fret->T     = qme->T1;

	fret->n    = dround(fret->T / fret->dt) + 1;
	fret->nout = dround(fret->T / fret->dtout) + 1;
	fret->nsto = dround(fret->dtout / fret->dt);

	fret->nspd = qme->nspd;
	fret->nosc = qme->nosc;

	fret->bSPLIT = qme->bSPLIT;

	int nsys   = fret->nsys;
	int nsyssq = fret->nsyssq;
	int nspd   = fret->nspd;

	smalloc(fret->Hsys0, nsyssq, qme->mem);
	scalloc(fret->reorg, nsyssq, qme->mem);
	smalloc(fret->Ediag, nsys,   qme->mem);

	if(fret->bSPLIT) {
		smalloc(fret->Hsys, nsyssq, qme->mem);
		smalloc(fret->osc,  nspd,   qme->mem);
		init_osc_split(qme, fret->osc);
	} else {
		fret->Hsys = qme->Hsys;
		fret->osc  = qme->osc;
	}

	check_uniform_temp(fret->osc, fret->nspd, fret->nosc, &fret->temp_uni, &fret->bTEMP_UNI);

	smalloc(fret->gt,  nsyssq, qme->mem);
	smalloc(fret->g1t, nsyssq, qme->mem);
	smalloc(fret->g2t, nsyssq, qme->mem);
	smalloc(fret->Gt,  nsyssq, qme->mem);

	for(i=0; i<nsyssq; i++) {
		smalloc(fret->gt[i],  fret->n, qme->mem);
		smalloc(fret->g1t[i], fret->n, qme->mem);
		smalloc(fret->g2t[i], fret->n, qme->mem);
		smalloc(fret->Gt[i],  fret->n, qme->mem);
	}

	scalloc(fret->rate, nsyssq, qme->mem);

	fprintf(stderr, "---- FRET simulation requested ----\n");
	fprintf(stderr, "Off-diagonal system-bath coupling will be neglected\n");
	fprintf(stderr, "\n");
}

void init_mrt(t_qme *qme, t_mrt *mrt)
{
	int i, j;

	mrt->unit   = qme->unit;

	mrt->nsys   = qme->nsys;
	mrt->nsyssq = qme->nsyssq;

	mrt->dt    = qme->dt1;
	mrt->dtout = qme->dtout1;
	mrt->T     = qme->T1;

	mrt->n    = dround(mrt->T / mrt->dt) + 1;
	mrt->nout = dround(mrt->T / mrt->dtout) + 1;
	mrt->nsto = dround(mrt->dtout / mrt->dt);

	mrt->nspd = qme->nspd;
	mrt->nosc = qme->nosc;
	mrt->osc  = qme->osc;

	mrt->bSPLIT = qme->bSPLIT;

	int nsys    = mrt->nsys;
	int nsyssq  = mrt->nsyssq;
	int nspd    = mrt->nspd;
	int *nosc   = mrt->nosc;

	smalloc(mrt->Usys,   nsyssq, qme->mem);
	smalloc(mrt->Eexci,  nsys,   qme->mem);
	smalloc(mrt->Eexci0, nsys,   qme->mem);

	if(mrt->bSPLIT) {
		smalloc(mrt->Hsys, nsyssq, qme->mem);
		smalloc(mrt->osc,  nspd,   qme->mem);
		init_osc_split(qme, mrt->osc);
	} else {
		mrt->Hsys = qme->Hsys;
		mrt->osc  = qme->osc;
	}

	t_osc **osc = mrt->osc;

	for(i=0; i<nspd; i++) {
		for(j=0; j<nosc[i]; j++) {
			osc[i][j].gexci = malloc(nsyssq * sizeof(double));
		}
	}

	smalloc(mrt->reaabb, nsyssq, qme->mem);
	smalloc(mrt->reabbb, nsyssq, qme->mem);

	smalloc(mrt->gtabbb, nsyssq, qme->mem);
	smalloc(mrt->gtaabb, nsyssq, qme->mem);
	smalloc(mrt->g1abbb, nsyssq, qme->mem);
	smalloc(mrt->g2abab, nsyssq, qme->mem);
	smalloc(mrt->Gt,     nsyssq, qme->mem);
	smalloc(mrt->Zt,     nsyssq, qme->mem);
	smalloc(mrt->Nt,     nsyssq, qme->mem);

	for(i=0; i<nsyssq; i++) {
		smalloc(mrt->gtabbb[i], mrt->n, qme->mem);
		smalloc(mrt->gtaabb[i], mrt->n, qme->mem);
		smalloc(mrt->g1abbb[i], mrt->n, qme->mem);
		smalloc(mrt->g2abab[i], mrt->n, qme->mem);
		smalloc(mrt->Gt[i],     mrt->n, qme->mem);
		smalloc(mrt->Zt[i],     mrt->n, qme->mem);
		smalloc(mrt->Nt[i],     mrt->n, qme->mem);
	}

	scalloc(mrt->rate, nsyssq, qme->mem);

	fprintf(stderr, "---- MRT simulation requested ----\n");
	fprintf(stderr, "\n");
}

void init_split(t_qme *qme, t_split *split)
{
	int i, j;

	split->jtype = qme->jtype;
	split->unit  = qme->unit;
	split->bDISS = qme->bDISS;

	split->nsys   = qme->nsys;
	split->nsyssq = qme->nsyssq;
	split->nstd = qme->nsys * (qme->nsys+1) / 2;
	if(split->nstd % 2) split->nstd++;

	split->nspd  = qme->nspd;
	split->nosc  = qme->nosc;
	split->osc0  = qme->osc;

	split->dtout = qme->dtout2;
	split->nout  = dround(qme->T2 / qme->dtout2) + 1;

	int nout = split->nout;
	int nsys = split->nsys;

	smalloc(split->popt_avg, nout * nsys, qme->mem);

	split->Hsys_prist = qme->Hsys;
	scalloc(split->Hsys_offset, split->nsyssq, qme->mem);
	scalloc(split->Hsys_std,    split->nsyssq, qme->mem);

	int nskip = split->rand_start;
	int nstd  = split->nstd;
	int ntraj = split->ntraj;
	smalloc(split->std, ntraj, qme->mem);

	// Random numbers discarded
	double r[2];
	for(i=0; i<nskip; i++) {
		for(j=0; j<nstd/2; j++) {
			gauss_gen(&r[0], &r[1]);
		}
	}

	// Random numbers actually used
	for(i=0; i<ntraj; i++) {
		smalloc(split->std[i], nstd, qme->mem);

		for(j=0; j<nstd/2; j++) {
			gauss_gen(&split->std[i][j*2], &split->std[i][j*2+1]);
		}
	}
	
	if(split->bDISS) {
		int nspd   = split->nspd;
		int *nosc  = split->nosc;
		int nsyssq = nsys * nsys;

		smalloc(split->disst_avg, nspd, qme->mem);
		for(i=0; i<nspd; i++) {
			scalloc(split->disst_avg[i], nosc[i], qme->mem);
		}
	}
}

void init_osc_split(t_qme *qme, t_osc **osc)
{
	int i, j;

	int nsyssq = qme->nsyssq;

	int nspd     = qme->nspd;
	int *nosc    = qme->nosc;
	t_osc **osc0 = qme->osc;

	for(i=0; i<nspd; i++) {
		osc[i] = osc_init(qme, nosc[i]);
		for(j=0; j<nosc[i]; j++) {
			osc_copy(&osc0[i][j], &osc[i][j], nsyssq);
		}
	}
}

void check_uniform_temp(t_osc **osc, int nspd, int *nosc, double *temp_uni, bool *bTEMP_UNI)
{
	int i, j;
	(*bTEMP_UNI) = true;

	double temp_ref = osc[0][0].temp;

	for(i=0; i<nspd; i++) {
		for(j=0; j<nosc[i]; j++) {
			if(temp_ref - osc[i][j].temp > SMALL) (*bTEMP_UNI) = false;
		}
	}

	if(*bTEMP_UNI) (*temp_uni) = temp_ref;
}
