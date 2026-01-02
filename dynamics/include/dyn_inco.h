#ifndef D_DYN_INCO
#define D_DYN_INCO

#include "qme.h"
#include "init_dyn.h"

void dyn_inco(t_qme *qme);
void conv_popt_inco(t_dyn_inco *dyn, double *popt_prt);
void prop_pop_inco(t_dyn_inco *dyn);
void prop_diss_inco(t_dyn_inco *dyn, int flag);
void calc_dpop(double *input, double *output, double *k, double dt_pop, int n);
void calc_ddiss(double *pop, double *out, double **rate, double dt, int nsys, int nosc);
void prop_pop_RK4(double *pop, double *k, double dt_pop, int n);
void prop_diss_RK4(double *pop, double *diss, double *rate_pop, double **rate_diss, double dt, int nsys, int nosc);

#endif
