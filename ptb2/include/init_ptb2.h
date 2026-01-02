#ifndef D_INIT_PTB2
#define D_INIT_PTB2

#include "qme.h"

void init_fret(t_qme *qme, t_fret *fret);
void init_mrt(t_qme *qme, t_mrt *mrt);
void init_split(t_qme *qme, t_split *split);
void init_osc_split(t_qme *qme, t_osc **osc);
void check_uniform_temp(t_osc **osc, int nspd, int *nosc, double *temp_uni, bool *bTEMP_UNI);

#endif
