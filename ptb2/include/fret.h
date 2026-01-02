#ifndef D_FRET
#define D_FRET

#include "qme.h"
#include "init_ptb2.h"

void fret_prep(t_qme *qme);
void fret_noisy_Hsys(t_fret *fret, t_split *split, double *std);
void fret_calc_Hsys0_reorg(t_fret *fret);
void fret_calc_gderiv(t_fret *fret);
void fret_calc_g(t_fret *fret);
void fret_calc_rates(t_fret *fret);
void fret_calc_diss(t_fret *fret);
void check_fret_physics(t_fret *fret);

#endif
