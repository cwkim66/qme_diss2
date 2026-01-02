#ifndef D_CALC_TR_MRT_NEW
#define D_CALC_TR_MRT_NEW

#include "qme.h"
#include "calc_tr.h"

void comp_tr_mrt_new(t_tr *tr);
void calc_tr_ana_mrt_mult(t_tr *tr, double t, complex *tr_ana);
void calc_mrt_micro_tr_num(t_tr *tr, double t, complex ***mtr);
void calc_mrt_micro_tr_ana(t_tr *tr, double t, complex ***mtr);
void calc_tr_micro_mrt_mult(t_tr *tr, double t, complex ***mtr, complex *tr_num, int mode);
void calc_Ct(t_tr *tr, double t, complex ***mtr, complex ***Ct, int n);
void calc_tr_ana_from_Ct(t_tr *tr, complex ***Ct, complex ***mtr, int n, complex *tr_ana, int mode);
void calc_mrt_micro_tr_deriv_num1(t_tr *tr, double t, complex ***mtr, complex *tr_num, int mode);
void calc_mrt_micro_tr_deriv_num2(t_tr *tr, double t, complex *tr_num, int mode);
void calc_mrt_micro_tr_deriv_ana(t_tr *tr, double t, complex *tr_ana, int mode);
void calc_eker_from_Ct(t_tr *tr, complex ***Ct, int n, complex *tr_num);
void calc_eker_ana(t_tr *tr, double t, complex *tr_ana);

#endif
