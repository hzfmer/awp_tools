#include "irvine2.h"

/* gets the maximum absolute value of a vector */
float max_abs(int npts, float *vec){
   int l;
   float ma=0;
   for (l=0;l<npts;l++){
      if (fabs(vec[l]) > ma) ma=fabs(vec[l]);
   }
   return(ma);
}

/* computes the absolute value of two vectors */
float* merge_hor(int npts, float *N, float *E){
     float *H;
     int k;
    
     H=(float*) calloc(npts, sizeof(float));
     for (k=0;k<npts;k++){
        H[k]=sqrt(pow(N[k],2) + pow(E[k],2));
     }
     return(H);
}

/* computes the shock response of a acceleration time series a with
   npts points, time step dt for the natural circular period omega
   with damping factor xi 
   the method is described in Irvine (2002), available at
   http://www.vibrationdata.com/tutorials2/srs_intr.pdf 
*/
void shock_response(int npts, float *a, float dt, float wn, float xi,
    float **resp){
    float wd = wn * sqrt(1-pow(xi,2));
    
    int n;   
    float *ra;
    ra=(float*) calloc (npts, sizeof(float) );
    for (n=2;n<npts;n++){
       ra[n]= 2 * exp (-xi*wn*dt) * cos(wd*dt) * ra[n-1]
          -exp(-2*xi*wn*dt)*ra[n-2] 
          +2*xi*wn*dt*a[n] 
          +wn*dt*exp(-xi*wn*dt) * ( (wn/wd * (1-2*pow(xi,2) ) ) *
             sin(wd*dt) -2*xi*cos(wd*dt) ) * a[n-1];
    }
    *resp = ra;    
}

/* computes peak of acceleration response for a vector of frequencies f,
   from time series a */
void calc_psa(int nfreq, float *f, int npts, float *a, float dt, float xi,
    float **psa){
    int k;

    float *psa0, *ra;
    float wn;
    psa0=(float*) calloc(nfreq, sizeof(float) );

    for (k=0; k<nfreq; k++){
       wn=f[k]* 2*M_PI;
       shock_response(npts, a, dt, wn, xi, &ra);
       psa0[k]=max_abs(npts, ra);
       free(ra);
    }

    *psa=psa0;
}

/* computes peak acceleration response for two horizontal vectors 
   N and E for and frequencies f */
void calc_psa_hor(int nfreq, float *f, int npts, float *N, float *E, 
    float dt, float xi, float **psa){

    int k;
    float wn;
    float *psaH, *respN, *respE, *respH;

    psaH=(float*) calloc(nfreq, sizeof(float) );

    for (k=0; k<nfreq; k++){
       wn=f[k]* 2*M_PI;
       shock_response(npts, N, dt, wn, xi, &respN);
       shock_response(npts, E, dt, wn, xi, &respE);
       respH=merge_hor(npts, respN, respE);
       free(respE);
       free(respN);
       psaH[k]=max_abs(npts, respH);
       free(respH);
    }
    *psa=psaH;
}
