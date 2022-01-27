#include <math.h>
#include <stdlib.h>

float max_abs(int npts, float *vec);
float* merge_hor(int npts, float *N, float *E);
void shock_response(int npts, float *a, float dt, float wn, float xi,
    float **resp);
void calc_psa(int nfreq, float *f, int npts, float *a, float dt, float xi,
    float **psa);
void calc_psa_hor(int nfreq, float *f, int npts, float *N, float *E, 
    float dt, float xi, float **psa);
