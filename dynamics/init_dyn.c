#include "init_dyn.h"

void init_dyn_inco(t_qme *qme, t_dyn_inco *dyn)
{
	int i, j;
	int one = 1;
	t_osc **osc;

	dyn->jtype  = qme->jtype;
	dyn->unit   = qme->unit;
	dyn->nsys   = qme->nsys;
	dyn->nsyssq = qme->nsyssq;
	dyn->bEXCI  = qme->bEXCI;
	dyn->bSPLIT = qme->bSPLIT;

	int jtype  = dyn->jtype;
	int nsys   = dyn->nsys;
	int nsyssq = dyn->nsyssq;

	dyn->dt    = qme->dt2;
	dyn->dtout = qme->dtout2;
	dyn->T     = qme->T2;

	dyn->n    = dround(dyn->T / dyn->dt) + 1;
	dyn->nout = dround(dyn->T / dyn->dtout) + 1;
	dyn->nsto = dround(dyn->dtout / dyn->dt);

	dyn->Hevec = qme->Hevec;
	dyn->pop0  = qme->pop0;
	scalloc(dyn->popt, nsys, qme->mem);

	int nout = dyn->nout;
	if(qme->bSPLIT) {
		smalloc(dyn->popt_sto, nout * nsys, qme->mem);
	}

	if(jtype % 10 == 0) { // fret
		t_fret *fret = qme->bSPLIT ? qme->split->fret : qme->fret;
		osc = qme->bSPLIT ? qme->split->fret->osc : qme->osc;
		dyn->rate_pop = fret->rate;
		dyn->E0 = fret->Ediag;
	} else if(jtype % 10 == 1) { // mrt
		t_mrt *mrt = qme->bSPLIT ? qme->split->mrt : qme->mrt;
		osc = qme->bSPLIT ? qme->split->mrt->osc : qme->osc;
		dyn->rate_pop = mrt->rate;
		dyn->E0 = mrt->Eexci0;
	}
		
	if(qme->bINCO && qme->bDISS) { // The required quantities are not affected by time scale separation
		dyn->nspd = qme->nspd;
		dyn->nosc = qme->nosc;
		dyn->nosc_all = 0;

		int nspd = dyn->nspd;
		int *nosc = dyn->nosc;

		for(i=0; i<nspd; i++) {
			dyn->nosc_all += dyn->nosc[i];
		}

		smalloc(dyn->freq,      dyn->nosc_all, qme->mem);
		smalloc(dyn->Idiss,     dyn->nosc_all, qme->mem);
		smalloc(dyn->rate_diss, dyn->nosc_all, qme->mem);
		smalloc(dyn->disst,     dyn->nosc_all, qme->mem);

		int count = 0;
		for(i=0; i<nspd; i++) {
			for(j=0; j<nosc[i]; j++) {
				dyn->freq[count]      = osc[i][j].freq;
				dyn->Idiss[count]     = osc[i][j].Idiss;
				dyn->rate_diss[count] = osc[i][j].Jdiss;
				count++;
			}
		}
	}
}
