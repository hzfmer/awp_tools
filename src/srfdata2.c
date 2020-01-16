#include "srfdata2.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

struct srfdata srfdata_init(int nx, int ny, int nt, float dt, int itype, int wstep, 
        int ntiskp, char *fname){
   struct srfdata srf;
   int l;

   srf.nx=nx;
   srf.ny=ny;
   srf.nt=nt;
   srf.dt=dt;
   srf.itype=itype;
   srf.wstep=wstep;
   srf.ntiskp=ntiskp;
   srf.fname=fname;
  
   srf.V=(float**) calloc(ny, sizeof(float*));
   if (srf.V==NULL) perror("error using calloc()");
   for (l=0; l<ny; l++){
      srf.V[l]=(float*) calloc(nx, sizeof(float));
      if (srf.V[l]==NULL) perror("error using calloc()");
   }
  
   return(srf);
}

void srfdata_loadts(struct srfdata *srf, int t){
   FILE *fid;
   off_t nskip;
   char errmsg[200];
   int sv;
   size_t npts, nread=0;
   int  n, l; 
   int askip, fnum, eskip;
   char fname2[200];

   if ((t % srf->ntiskp) != 0) {
      fprintf(stderr, "Error. Timestep %d requested, but ntiskp is %d.\n", 
          t, srf->ntiskp);
      exit(1);
   }

   npts=(size_t) srf->nx * srf->ny;

   if (srf->itype==0){
      n=t / srf->ntiskp;
      fid=fopen(srf->fname, "r");
      if (fid==NULL){
         sprintf(errmsg, "Error opening %s", srf->fname);
         perror(errmsg);
         exit(1);
      }
      nskip=(off_t) srf->nx * srf->ny * n * sizeof(float);
      /*fprintf(stdout, "seeking to n=%d, seekpos=%ld\n", n, nskip);*/
      sv=fseeko(fid, nskip, SEEK_SET);
      if (sv!=0) {
         perror("Error using fseek()");
         fprintf(stderr, "fseek was trying to seek to pos %lld\n", nskip);
         exit(1);
      }
      
      for (l=0; l < srf->ny; l++){
         nread+=fread(srf->V[l], sizeof(float), srf->nx, fid);
      }

      if (nread != npts) {
         sprintf(errmsg, "%ld items expected in fread, but %ld read.", npts, nread);
         perror(errmsg);
         exit(1);
      }
      fclose(fid);
   }
   else {
      askip=srf->wstep * srf->ntiskp;
      fnum=(int) ceilf(t/ ((float) askip)) * askip;

      sprintf(fname2, "%s%07d", srf->fname, fnum);
      fid=fopen(fname2, "r");
      if (fid==NULL){
         sprintf(errmsg, "Error opening %s", fname2);
         perror(errmsg);
         exit(1);
      }

      eskip= (int) ((float) (t - fnum + askip -1 )) / srf->ntiskp;
      /*fprintf(stdout, " >> loading entry %d in %s\n", eskip, fname2);*/

      nskip=(size_t) eskip * npts * sizeof(float);
      sv=fseeko(fid, nskip, SEEK_SET);
      if (sv!=0) {
         perror("Error using fseek()");
         fprintf(stderr, "fseek was trying to seek to pos %lld\n", nskip);
         exit(1);
      }

      for (l=0; l < srf->ny; l++){
         nread+=fread(srf->V[l], sizeof(float), srf->nx, fid);
      }

      if (nread != npts) {
         sprintf(errmsg, "%ld items expected in fread, but %ld read.", npts, nread);
         perror(errmsg);
         exit(1);
      }
      fclose(fid);
  }


}

void srfdata_getts(struct srfdata *srf, float *D, int t){
   FILE *fid;
   char fname2[200];
   off_t nskip;
   char errmsg[200];
   int n;
   int sv;
   size_t npts, nread=0;
   int askip, fnum, eskip;

   if ((t % srf->ntiskp) != 0) {
      fprintf(stderr, "Error. Timestep %d requested, but write_step is %d.\n", 
          t, srf->ntiskp);
      exit(1);
   }
   npts=(size_t) srf->nx * srf->ny;

   if (srf->itype==0){
      n=t / srf->ntiskp;
      fid=fopen(srf->fname, "r");
      if (fid==NULL){
         sprintf(errmsg, "Error opening %s", srf->fname);
         perror(errmsg);
         exit(1);
      }
      nskip=(off_t) srf->nx * srf->ny * n * sizeof(float);
      /*fprintf(stdout, "seeking to n=%d, seekpos=%ld\n", n, nskip);*/
      sv=fseeko(fid, nskip, SEEK_SET);
      if (sv!=0) {
         perror("Error using fseek()");
         fprintf(stderr, "fseek was trying to seek to pos %lld\n", nskip);
         exit(1);
      }
      
      nread=fread(D, sizeof(float), npts, fid);

      if (nread != npts) {
         sprintf(errmsg, "%ld items expected in fread, but %ld read.", npts, nread);
         perror(errmsg);
         exit(1);
      }
      fclose(fid);
   }
   else { 
      askip=srf->wstep * srf->ntiskp;
      fnum=(int) ceilf(t/((float) askip)) * askip;

      sprintf(fname2, "%s%07d", srf->fname, fnum);
      fid=fopen(fname2, "r");
      if (fid==NULL){
         sprintf(errmsg, "Error opening %s", fname2);
         perror(errmsg);
         exit(1);
      }

      eskip= (int) ((float) (t - fnum + askip -1 )) / srf->ntiskp;
      /*fprintf(stdout, " >> loading entry %d in %s\n", eskip, fname2);*/

      nskip=(size_t) eskip * npts * sizeof(float);
      sv=fseeko(fid, nskip, SEEK_SET);
      if (sv!=0) {
         perror("Error using fseek()");
         fprintf(stderr, "fseek was trying to seek to pos %lld\n", nskip);
         exit(1);
      }

      nread=fread(D, sizeof(float), npts, fid);

      if (nread != npts) {
         sprintf(errmsg, "%ld items expected in fread, but %ld read.", npts, nread);
         perror(errmsg);
         exit(1);
      }
      fclose(fid);
  }

}

float* srfdata_gettrace(struct srfdata *srf, int k, int l){
   size_t size_per_ts, nskip;
   size_t seekpos;
   int n;
   float *ts;
   int nt2;
   FILE *fid;
   char errmsg[200];

   nt2=srf->nt / srf-> ntiskp;

   if (srf->itype != 0){
     fprintf(stderr, "only IVELOCITY=0 is implemented as for now\n");
     return(NULL);
   }

   size_per_ts=srf->ny * srf->nx * sizeof(float);

   nskip= ( srf->nx * l + k ) * sizeof(float);

   ts=(float*) calloc(nt2, sizeof(float));

   fid=fopen(srf->fname, "r");
   if (fid==NULL){
      sprintf(errmsg, "could not open %s", srf->fname);
      perror(errmsg);
      return(NULL);
   }

   for (n=0; n<srf->nt; n++){
      seekpos = n * size_per_ts + nskip;
      fseek(fid, seekpos, SEEK_SET);
      fread(ts+n, sizeof(float), 1, fid);
   }
   fclose(fid);
   return(ts);
}

void srfdata_getrow(struct srfdata *srf, int l, float **D){
   size_t size_per_ts, nskip;
   size_t seekpos;
   int n, ix;
   int nt2;
   FILE *fid;
   char errmsg[200];
   int seekstat;
   float *tmp;

   nt2=srf->nt / srf-> ntiskp;
   
   tmp=(float*) calloc(srf->nx, sizeof(float));

   if (srf->itype != 0){
     fprintf(stderr, "only IVELOCITY=0 is implemented for now\n");
     return;
   }

   size_per_ts=srf->ny * srf->nx * sizeof(float);

   nskip= ( srf->nx * l ) * sizeof(float);

   fid=fopen(srf->fname, "r");
   if (fid==NULL){
      sprintf(errmsg, "could not open %s", srf->fname);
      perror(errmsg);
      return;
   }

   for (n=0; n<nt2; n++){
      /*fprintf(stdout, "\r read %d out of %d timesteps", n, nt2);
      fflush(stdout);*/
      seekpos = (size_t) n * size_per_ts + nskip;
      seekstat=fseeko(fid, seekpos, SEEK_SET);
      if (seekstat != 0){
        sprintf(errmsg, "could not seek to position %ld", seekpos);
        perror(errmsg);
        exit(-1);
      }
      fread(tmp, sizeof(float), srf->nx, fid);
      for (ix=0; ix<srf->nx; ix++){
        D[ix][n]=tmp[ix];
      }
   }
   /*fprintf(stdout, "\n");*/
   free(tmp);
   fclose(fid);
}

void srfdata_getrow_wstep(struct srfdata *srf, int l, float **D){
   size_t skipto_next_ts, nskip;
   int n=0, ix, si;
   FILE *fid;
   char errmsg[200];
   int seekstat;
   float *tmp;
   int askip, ti, ti2, nentries;
   char tfname[200];
   
   tmp=(float*) calloc(srf->nx, sizeof(float));

   skipto_next_ts=(srf->ny-1) * srf->nx * sizeof(float);
   nskip= ( srf->nx * l ) * sizeof(float);
   askip=srf->wstep * srf->ntiskp;

   for (ti=askip; ti < (srf->nt+askip); ti+=askip){
      if (ti > srf->nt){
         ti2=srf->nt-1;
         nentries= srf->nt % askip / srf->ntiskp;
      }
      else{
         ti2=ti;
         nentries=srf->wstep;
      }
      sprintf(tfname, "%s%07d", srf->fname, ti2);
      /*fprintf(stdout, "opening %s\n", tfname);*/
      fid=fopen(tfname, "r");
      if (fid==NULL){
         sprintf(errmsg, "could not open %s", tfname);
         perror(errmsg);
         exit(-1);
      }

      seekstat=fseeko(fid, nskip, SEEK_SET);
      if (seekstat != 0){
        sprintf(errmsg, "could not seek to position %ld", nskip);
        perror(errmsg);
        exit(-1);
      }
      for (si=0; si<nentries; si++){
         /*fprintf(stdout, "reading entry %d, saving to ts %d\n", si, n);*/
         fread(tmp, sizeof(float), srf->nx, fid);
         for (ix=0; ix<srf->nx; ix++){
            D[ix][n]=tmp[ix];
         }
         seekstat=fseeko(fid, skipto_next_ts, SEEK_CUR);
         if (seekstat != 0){
            sprintf(errmsg, "could not seek to position %ld", skipto_next_ts);
            perror(errmsg);
            exit(-1);
         }
         n++;
      }
      fclose(fid);
   }
   free(tmp);
}

void srfdata_getrows_wstep(struct srfdata *srf, int l, int nrows, float ***D){
   size_t skipto_next_ts, nskip;
   int n=0, ix, si;
   FILE *fid;
   char errmsg[200];
   int seekstat;
   float *tmp;
   int askip, ti, ti2, nentries;
   char tfname[200];
   int rn;
   
   tmp=(float*) calloc(srf->nx, sizeof(float));

   skipto_next_ts=(srf->ny-nrows) * srf->nx * sizeof(float);
   nskip= ( srf->nx * l ) * sizeof(float);
   askip=srf->wstep * srf->ntiskp;

   for (ti=askip; ti < (srf->nt+askip); ti+=askip){
      if (ti > srf->nt){
         ti2=srf->nt-1;
         nentries= srf->nt % askip / srf->ntiskp;
      }
      else{
         ti2=ti;
         nentries=srf->wstep;
      }
      sprintf(tfname, "%s%07d", srf->fname, ti2);
      fprintf(stdout, "opening %s\n", tfname);
      fid=fopen(tfname, "r");
      if (fid==NULL){
         sprintf(errmsg, "could not open %s", tfname);
         perror(errmsg);
         exit(-1);
      }

      seekstat=fseeko(fid, nskip, SEEK_SET);
      if (seekstat != 0){
        sprintf(errmsg, "could not seek to position %ld", nskip);
        perror(errmsg);
        exit(-1);
      }
      for (si=0; si<nentries; si++){
         /*fprintf(stdout, "reading entry %d, saving to ts %d\n", si, n);*/
         for (rn=0; rn<nrows; rn++){
            fread(tmp, sizeof(float), srf->nx, fid);
            for (ix=0; ix<srf->nx; ix++){
               D[rn][ix][n]=tmp[ix];
            }
         }
         seekstat=fseeko(fid, skipto_next_ts, SEEK_CUR);
         if (seekstat != 0){
            sprintf(errmsg, "could not seek to position %ld", skipto_next_ts);
            perror(errmsg);
            exit(-1);
         }
         n++;
      }
      fclose(fid);
   }
   free(tmp);
}


void srfdata_getchunk(struct srfdata *srf, int l1, int l2, float ***D){
   size_t size_per_ts, nskip;
   size_t seekpos;
   int n, ix;
   int nt2;
   FILE *fid;
   char errmsg[200];
   int seekstat;
   float *tmp;
   int nrows, iy;

   nt2=srf->nt / srf-> ntiskp;
   
   nrows=l2 - l1;
   tmp=(float*) calloc(srf->nx * nrows, sizeof(float));


   if (srf->itype != 0){
     fprintf(stderr, "only IVELOCITY=0 is implemented for now\n");
     return;
   }

   size_per_ts=srf->ny * srf->nx * sizeof(float);

   nskip= ( srf->nx * l1 ) * sizeof(float);

   fid=fopen(srf->fname, "r");
   if (fid==NULL){
      sprintf(errmsg, "could not open %s", srf->fname);
      perror(errmsg);
      return;
   }

   for (n=0; n<nt2; n++){
      fprintf(stdout, "\r read %d out of %d timesteps", n+1, nt2);
      fflush(stdout);
      seekpos = (size_t) n * size_per_ts + nskip;
      seekstat=fseeko(fid, seekpos, SEEK_SET);
      if (seekstat != 0){
        sprintf(errmsg, "could not seek to position %ld", seekpos);
        perror(errmsg);
        exit(-1);
      }
      fread(tmp, sizeof(float), srf->nx*nrows, fid);
      for (iy=0; iy<nrows; iy++){
         for (ix=0; ix<srf->nx; ix++){
           D[iy][ix][n]=tmp[iy * srf->nx + ix];
         }
      }
   }
   fprintf(stdout, "\n");
   free(tmp);
   fclose(fid);
}

/* saves timestep in pre-allocated, continuous (i.e. 1-D) array */
/*void srfdata_getts(struct srfdata *srf, float *D, int t){
   srfdata_loadts(srf, t);
   int k, l, m=0;*/

   /*float *D; 
   size_t ne;
   ne = (size_t) srf->nx * srf->ny;

   D=(float*) calloc(ne, sizeof(float));*/

   /*
   srfdata_loadts(srf, t);
   for (l=0; l< srf->ny; l++){
      for (k=0; k< srf->nx; k++){
         D[m] = srf->V[l][k];
         m++;
      }
   }
}*/
