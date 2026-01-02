#ifndef D_MRT_POP
#define D_MRT_POP

#include "qme.h"
#include "init_ptb2.h"

void mrt_prep(t_qme *qme);
void mrt_calc_H_reorg(t_mrt *mrt);
void mrt_calc_g(t_mrt *mrt);
void mrt_calc_rates(t_mrt *mrt);
void mrt_calc_diss(t_mrt *mrt);
void check_mrt_physics(t_mrt *mrt);

#endif
