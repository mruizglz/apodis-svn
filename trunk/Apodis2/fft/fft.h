#include <stdlib.h>
#include <math.h>

/* macros */
#define TWO_PI (6.2831853071795864769252867665590057683943L)

/* function prototypes */
void fft(int N, double (*x)[2], double (*X)[2]);
void fft_rec(int N, int offset, int delta,
             double (*x)[2], double (*X)[2], double (*XX)[2]);
void ifft(int N, double (*x)[2], double (*X)[2]);
