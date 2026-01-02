#ifndef D_CALC_TR
#define D_CALC_TR

#include "qme.h"
#include "init_tr.h"

void calc_tr(t_qme *qme);
void comp_tr_fret(t_tr *tr, int mode);
void comp_tr_mrt(t_tr *tr, int mode);
void print_tr(t_tr *tr, char *fname, complex **tr1, complex **tr2);
void calc_tr_num_fret_pop(t_qosc *qosc, double t, complex *tr_num);
void calc_tr_ana_fret_pop(t_qosc *qosc, double t, complex *tr_ana);
void calc_tr_num_fret_diss(t_qosc *qosc, double t, complex *tr_num);
void calc_tr_ana_fret_diss(t_qosc *qosc, double t, complex *tr_ana);
void calc_tr_num_mrt_pop(t_qosc *qosc, double t, complex *tr_num);
void calc_tr_ana_mrt_pop(t_qosc *qosc, double t, complex *tr_ana);
void calc_tr_num1_mrt_dpdt(t_qosc *qosc, double t, complex *tr_num);
void calc_tr_num2_mrt_dpdt(t_qosc *qosc, double t, complex *tr_num);
void calc_tr_num2_thermal_mrt_dpdt(t_qosc *qosc, double t, complex *tr_num);
void calc_tr_ana1_mrt_dpdt(t_qosc *qosc, double t, complex *tr_ana);
void calc_tr_ana2_mrt_dpdt(t_qosc *qosc, double t, complex *tr_ana);

#endif
