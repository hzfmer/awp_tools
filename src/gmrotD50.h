#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "irvine2.h"

float gmrotD50(float *x, float*y, int nt, float dt, float xi, 
   float f, int nang);

float* logspace(float min, float max, int n);

float* gmrotIpp(float *xv, float *yv, int nt, float dt, float *freq, int nf, 
   float xi, float nang, float *minang);
