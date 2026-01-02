#ifndef D_OPT_SPLIT
#define D_OPT_SPLIT

#include "qme.h"
#include "init_ptb2.h"

void avg_split(t_qme *qme);
void avg_split_lowlvl(t_split *split);
void split_spd_fret(t_split *split);
double calc_scalefactor_2015castillo_gen(double w, double wsp, double pref);
double calc_scalefactor_exponential(double w, double wsp, double pref);

#endif
