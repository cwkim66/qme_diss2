#ifndef D_UTIL
#define D_UTIL

#define BYTE_TO_MB 1.0e-06
#define smalloc(ptr,nelem,mem) (ptr)=save_malloc(#ptr,nelem,&mem,sizeof(*(ptr)))
#define scalloc(ptr,nelem,mem) (ptr)=save_calloc(#ptr,nelem,&mem,sizeof(*(ptr)))

#include "qme.h"

FILE* sfopen(char *filename, int mode);
void *save_calloc(char *name, int nelem, double *mem, size_t size);
void *save_malloc(char *name, int nelem, double *mem, size_t size);
int get_lines(char *filename);
void upper2low(char *line);
void prt_vib_dmat(int ndim, int nvib, double *vib_mat, double fac);
void prt_vib_dmat_sci(int ndim, int nvib, double *vib_mat, double fac);
void prt_vib_zmat(int ndim, int nvib, complex *vib_mat, double fac);
void prt_vib_zmat_sci(complex *vib_mat, int ndim, int nvib, double fac);
void prt_vib_zmat_imag(int ndim, int nvib, complex *vib_mat, double fac);

#endif // D_UTIL
