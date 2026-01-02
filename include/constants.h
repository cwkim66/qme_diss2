#ifndef D_CONSTANTS

#define PI 

// hbar and kb is set to be unity
#define HBAR_MKS 1.0545718e-34
#define kb_MKS 1.38064852e-23
#define C_CGS 2.997992458e+10
#define PS_TO_S 1.0e-12
#define FS_TO_S 1.0e-15
#define TWOPI (2.0 * M_PI)

#define TUNIT_TO_S (HBAR_MKS / kb_MKS)
#define TUNIT_TO_FS (TUNIT_TO_S / FS_TO_S)
#define TUNIT_TO_PS (TUNIT_TO_S / PS_TO_S)
#define WAVENO_TO_EUNIT (C_CGS * TWOPI * TUNIT_TO_S)
#define WAVENO_TO_OMEGA_FS (C_CGS * TWOPI * FS_TO_S)

#define SMALL 1.0e-08
#define SMALL_POP 1.0e-03
#define SMALL_NW 1.0e-15
#define LARGE 10000

#define XX 0
#define YY 1
#define ZZ 2
#define ID 3

#endif // D_CONSTANTS
