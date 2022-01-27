/* Orientation-independent horizontal response spectrum,
   gmrotI50 (Boore et all., 2006, BSSA 96(4) 1502-1511) 
   Daniel Roten, Apr 2010 <droten@sciences.sdsu.edu> */

#include "gmrotD50.h"

/* a structure that contains the rotation angle and the geometric mean */
struct gm_ang {
   float gm;
   float theta;
};

float* logspace(float min, float max, int n){
    float *f;
    f=(float*) calloc(n, sizeof(float));
    int k;
    for (k=0; k<n; k++){
       f[k] = pow(10., (log10(max)-log10(min)) * (float)k/(n-1) + log10(min));
    }
    return(f);
}

typedef int (*compfn)(const void*, const void*);

/* invoke qsort on structure with geometric mean and angle */
int gm_cmp(struct gm_ang *a, struct gm_ang *b) {
    /* comparison of gm in struct: returns negative if a.gm > b.gm
    and positive if a.gm > b.gm */
    if (a->gm < b->gm ) return(-1);
    else if (a->gm > b->gm ) return(1);
    else return(0);
}

void rot_ts(float *x, float *y, float *xr, float *yr, int nt, float ang){
   int n;
   float dang;
   dang = ang / 180. * M_PI;
   for (n=0; n<nt; n++){
       xr[n] =  x[n] * cosf(dang) + y[n] * sinf(dang);
       yr[n] = -x[n] * sinf(dang) + y[n] * cosf(dang);
   }
}

float maxabs(float *h, int nt){
   int n;
   float mv=0;
   for (n=0; n<nt; n++){
      if (fabs(h[n]) > mv) mv=fabs(h[n]);
   }
   return(mv);
}

struct gm_ang* gmrotDpp(float *x, float *y, int nt, float dt, float xi, 
   float f, int nang){
   float *ox, *oy, *oxr, *oyr; /* oscillator repsonse time series */
   float wn;
   struct gm_ang *GMA;
   float xm, ym;
   int r; 

   wn=f*2*M_PI;

   oxr=(float*) calloc(nt, sizeof(float));
   oyr=(float*) calloc(nt, sizeof(float));

   GMA=(struct gm_ang*) calloc(nang, sizeof(struct gm_ang));

   shock_response(nt, x, dt, wn, xi, &ox);
   shock_response(nt, y, dt, wn, xi, &oy);

   for (r=0; r<nang; r++){
      GMA[r].theta = r * 1.;
      rot_ts(ox, oy, oxr, oyr, nt, GMA[r].theta);
      xm=maxabs(oxr, nt);
      ym=maxabs(oyr, nt);
      GMA[r].gm=sqrtf(xm * ym); 
   }

   free(ox);
   free(oy);
   free(oxr);
   free(oyr);

   return(GMA);
}

void get_median(struct gm_ang *GMA, int nang, float *median, float *angle){
   struct gm_ang *tmpgm;
   int k, medpos;
 
   /* make a local copy, to avoid altering input array */
   tmpgm = (struct gm_ang*) calloc(nang, sizeof(struct gm_ang));
   for (k=0; k<nang; k++){
      tmpgm[k].gm = GMA[k].gm;
      tmpgm[k].theta = GMA[k].theta;
   }
   medpos=(int) roundf((float) nang/2.);
   qsort(GMA, nang, sizeof(struct gm_ang), (compfn) gm_cmp);
   *median=GMA[medpos].gm;
   *angle=GMA[medpos].theta;
   free(tmpgm);
}

float gmrotD50(float *x, float*y, int nt, float dt, float xi, 
   float f, int nang) {
   struct gm_ang *gma;
   float median; 
   float med_angle;
   gma=gmrotDpp(x, y, nt, dt, xi, f, nang);
   get_median(gma, nang, &median, &med_angle);
   return(median);
}

float* penalty_function(struct gm_ang **GMA, int nf, int nang){
   int g, r;
   float *penalty;
   float median, med_angle;

   median=(int) roundf(nang / 2.);
   penalty=(float*) calloc(nang, sizeof(float));

   for (r=0; r<nang; r++){
      for (g=0; g<nf; g++){
         get_median(GMA[g], nang, &median, &med_angle);
         penalty[r] += powf(GMA[g][r].gm / median - 1., 2.);
      }
      penalty[r] /= nf;
   }
   return(penalty);
}

float* gmrotIpp(float *xv, float *yv, int nt, float dt, float *freq, int nf, 
   float xi, float nang, float *minang){
   struct gm_ang **GMA;
   float *iresp, *penalty;
   float minpen=1e9;
   int minidx=-1, g, r;

   GMA=(struct gm_ang**) calloc(nf, sizeof(struct gm_ang*));

   for (g=0; g<nf; g++){
      GMA[g]=gmrotDpp(xv, yv, nt, dt, xi, freq[g], nang);
   }

   /* locate minimum in penalty function */
   penalty=penalty_function(GMA, nf, nang);
   for (r=0; r<nang; r++){
      if (penalty[r] < minpen) {
         minpen = penalty[r];
         minidx = r; 
      }
   }
   *minang=GMA[0][minidx].theta;

   iresp=(float*) calloc(nf, sizeof(float));
   for (g=0; g<nf; g++){
      iresp[g] = GMA[g][minidx].gm;
   }
   for (g=0; g<nf; g++) free(GMA[g]);
   free(GMA);
   free(penalty);
   return(iresp);
}

float* load (char *fname, int nt){
   FILE *fid;
   int n;
   float *dat;
   dat=(float*) calloc(nt, sizeof(float));
   fid=fopen(fname, "r");
   if (fid==NULL){
      perror("could not load file.");
      exit(-1);
   }
   for (n=0; n<nt; n++){
      fscanf(fid, "%e\n", dat+n);
   }
   fclose(fid);
   return(dat);
}

int test () {
   float *xv, *yv;
   int nt=2300;
   float freq=0.50;
   int nang=91;
   float xi=0.05;
   float dt=0.025;
   float resp;

   xv=load("/Users/daniel/work/SLV/M7a_3c/L_5a.M_tp/output_sfc/SX96PS_0666_0889.dat", nt);
   yv=load("/Users/daniel/work/SLV/M7a_3c/L_5a.M_tp/output_sfc/SY96PS_0666_0889.dat", nt);

      /*fprintf(stdout, "%e %e %e %e %e %e %e\n", freq[g], 
         GMA[g][0].gm, GMA[g][46].gm, GMA[g][90].gm, 
         GMA[g][0].theta, GMA[g][46].theta, GMA[g][90].theta);*/
   /*freq=logspace(0.1, 1.0, nf);
   resp=gmrotIpp(xv, yv, nt, dt, freq, nf, xi, pp, nang, &minang);
   fprintf(stdout, "#angle of minimum: %e\n", minang);
   for (g=0; g<nf; g++) fprintf(stdout, "%e %e\n", freq[g], resp[g]);*/

   resp=gmrotD50(xv, yv, nt, dt, xi, freq, nang);
   fprintf(stdout, "resp=%e\n", resp);
   return(0);
}
