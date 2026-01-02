#ifndef D_GAUSS
#define D_GAUSS

#ifndef D_GENLIB
#define D_GENLIB
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#endif

typedef unsigned long uint32;

#define NN             (624)                 // length of state vector
#define MM             (397)                 // a period parameter
#define KK             (0x9908B0DFU)         // a magic constant
#define hiBit(u)       ((u) & 0x80000000U)   // mask all but highest   bit of u
#define loBit(u)       ((u) & 0x00000001U)   // mask all but lowest    bit of u
#define loBits(u)      ((u) & 0x7FFFFFFFU)   // mask     the highest   bit of u
#define mixBits(u, v)  (hiBit(u)|loBits(v))  // move hi bit of u to hi bit of v

void seedMT(uint32 seed);
uint32 reloadMT(void);
uint32 randomMT(void);

#endif
