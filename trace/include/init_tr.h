#ifndef D_INIT_TR
#define D_INIT_TR

#include "qme.h"

#define NMAX_DEF 5

void init_tr(t_qme *qme);
void init_osc_tr(t_qme *qme);
double set_nqst(t_qosc qosc, double cut);
void set_FC(t_qosc *qosc, int jtype);
void ortho_FC(t_qosc *qosc);
void set_Req(t_qosc *qosc);
void set_H_diag(t_qosc *qosc);
void set_H_offdiag(t_qosc *qosc);
void set_U_diag(t_qosc *qosc, double t);
void set_U_diag_wick(t_qosc *qosc, double t, bool mode);

#endif
