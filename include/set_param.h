#ifndef D_SET_PARAM
#define D_SET_PARAM

#include "qme.h"
#include "spd.h"

void get_option(t_qme *qme, FILE *f_opt);
void get_option_sys(t_qme *qme, FILE *f_opt);
void get_option_spd(t_qme *qme, FILE *f_opt);
void get_option_init(t_qme *qme, FILE *f_opt);
void get_option_init_pop(t_qme *qme, FILE *f_opt);
void get_option_init_rho(t_qme *qme, FILE *f_opt);
void get_option_trace(t_qme *qme, FILE *f_opt);
void load_defaults(t_qme *qme);
void load_defaults_sys(t_qme *qme);
void load_defaults_init(t_qme *qme);
void load_defaults_trace(t_qme *qme, bool bFLAG);
void load_defaults_split(t_qme *qme);

#endif
