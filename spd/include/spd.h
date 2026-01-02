#ifndef D_SPD
#define D_SPD

#include "qme.h"

t_osc *osc_init(t_qme *qme, int nosc);
void osc_copy(t_osc *osc_src, t_osc *osc_dest, int nsyssq);
int const_drude_spd(FILE *f_spd, char *f_name, t_qme *qme, t_osc **osc);
int const_lognormal_spd(FILE *f_spd, char *f_name, t_qme *qme, t_osc **osc);
int const_brownian_spd(FILE *f_spd, char *f_name, t_qme *qme, t_osc **osc);
int const_singleosc_spd(FILE *f_spd, char *f_name, t_qme *qme, t_osc **osc);
int const_num_spd(FILE *f_spd, char *f_name, t_qme *qme, t_osc **osc);
int calc_reorg(t_qme *qme);

#endif
