/* Computes peak ground velocity from AWP-ODC output.
   Improved version using MPI-IO.
   Daniel Roten, March 2015 <droten@sdsc.edu> */

/* This version also computes gmrotD50 and duration, in addition to 
 * PGVs. */

#include "awp_data.h"
#include "fd3dparam.h"
#include "xapiir.h"
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdlib.h>
#include "gmrotD50.h"

#define PI 3.1415927
#define G 9.8


int index_above(float *in, int nt, float thres) {
    int l = 0, r = nt - 1;
    while (l < r) {
        int mid = l + (r - l) / 2;
        if (in[mid] < thres) {
            l = mid + 1;
        } else {
            r = mid;
        }
    }
    return l;
}

void comp_cum(float *accX, float *accY, float *accZ, int nt, float dt, float low, float high, float *dur, float *ener){
    float *cumeng;
    float low_bound, high_bound;
    int l;

    cumeng = (float*) calloc(nt, sizeof(float));
    cumeng[0] = pow(accX[0], 2.) + pow(accY[0], 2.) + pow(accZ[0], 2.);
 
    for (l=1; l<nt; l++){
        cumeng[l] = cumeng[l - 1] + pow(accX[l], 2.) + pow(accY[l], 2.) + pow(accZ[l], 2.);
    }
    low_bound = cumeng[nt - 1] * low;
    high_bound = cumeng[nt - 1] * high;

    int low_idx = index_above(cumeng, nt, low_bound);
    int high_idx = index_above(cumeng, nt, high_bound);

    free(cumeng);
    *dur = (high_idx - low_idx) * dt;
    *ener = cumeng[nt - 1] / 3;

}

void vel2acc(float *vel, int nt, float dt, float *acc){
   int n;

   for (n=0; n<(nt-1); n++){
      acc[n] = (vel[n+1] - vel[n]) / dt;
   }
   acc[nt-1] = 0.;
}

int main (int argc, char *argv[] ){
   FD3D_param par;
   read_settings(&par, "IN3D.out");
   char *xfile, *yfile, *zfile;
   int nx=1, ny=1, nt=1, nz=1;
   float **x, **y, **z;
   float *accx, *accy, *accz, acc;
   float *bufx, *bufy, *bufz;
   float *PH, *PZ, *PGV, vh, vz, vel;
   float *AI, *PGA, *ENER, *DUR;
   MPI_File hfid, zfid, pgvfid, pgafid;
   MPI_File aifid, durfid, enerfid;
   int nchunks = 18, csize; //nchunks=15
   int rank, nprocs;
   int k, l, n;
   MPI_Offset off;
   long int s0;
   int xi, yi, zi;
   float dt2;

   /*parameters for xapiir*/
   int dofilt = 0;
   int comp_sa = 0;
   int iord=3, npas=1;
   float trbndw=0., a=0.; /* chebyshev parameters */
   char *aproto="BU";
   char ftype[3];
   float lp=0.15, hp=5.00; 
   float low=0.05, high=0.75;

   /* parameters for getopt */
   int c;
   extern char *optarg;

   char phfile[200], pzfile[200], pgvfile[200], pgafile[200], aifile[200], durfile[200], enerfile[200], gmfile[200];
   int err;

   float freq[]={0.1, 0.333333, 0.5, 0.6666666, 1.0, 1.5, 2.0, 3.0, 3.5};
   int nfreq=(int) sizeof(freq) / sizeof(freq[0]);

   float **GM;  //gmrotD50 and duration
   int f;
   float xid = 0.05; // damping constant for gmrotD50
   int nang=91;  // number of angles for gmrotD50

   MPI_File *gmfid;
 
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  while ((c = getopt (argc, argv, "fgl:h:n:s:e:")) != -1) {
     switch(c) {
        case 'f':
           dofilt = 1;
           break;
        case 'g':
           comp_sa = 1;
           break;
        case 'l':
           sscanf(optarg, "%f", &lp);
           break;
        case 'h':
           sscanf(optarg, "%f", &hp);
           break;
        case 'n':
           sscanf(optarg, "%d", &nchunks);
           break;
        case 's':
           sscanf(optarg, "%f", &low);
           break;
        case 'e':
           sscanf(optarg, "%f", &high);
           break;
        case '?':
           fprintf(stdout, "usage: %s [-f] [-l lp] [-h hp] [-n nchunks] [-s low] [-e high]\n", argv[0]);
           exit(0);
        default:
           abort();
     }
  }

  if (hp > 0.) strncpy(ftype, "BP\0", 3);
  else strncpy(ftype, "LP\0", 3);
  if (rank==0) fprintf(stdout, "dofilt=%d, lp=%f, hp=%f, ftype=%s, nchunks=%d, low=%f, high=%f\n", dofilt, lp, hp, ftype, nchunks, low, high);

   /* determine dimensions from these parameters */
  nx=(int) floorf( (par.nedx-par.nbgx)/par.nskpx + 1.);
  ny=(int) floorf( (par.nedy-par.nbgy)/par.nskpy + 1.);
  nz=(int) floorf( (par.nedz-par.nbgz)/par.nskpz + 1.);
  nt=(int) floorf( (par.tmax/par.dt/par.ntiskp));
  dt2=par.dt*par.ntiskp;
   
   /* output parameters for debug */
   if (rank==0) fprintf(stdout, "  >> nx=%d, ny=%d, nz=%d, nt=%d," 
                   " ivelocity=%d, ntiskp=%d, writestep=%d\n",
		    nx, ny, nz, nt, par.itype, par.ntiskp, par.writestep);

   xfile=(char*) calloc(100, sizeof(char) );
   yfile=(char*) calloc(100, sizeof(char) );
   zfile=(char*) calloc(100, sizeof(char) );
   sscanf(par.sxrgo, "%s", xfile);
   sscanf(par.syrgo, "%s", yfile);
   sscanf(par.szrgo, "%s", zfile);
   /* this removes the remaining ' at the end of the string, if any */
   xfile=strsep(&xfile, "\'");
   yfile=strsep(&yfile, "\'");
   zfile=strsep(&zfile, "\'");

   /* number of points to process per read operation */
   if (((nx*ny*nz) % (nprocs * nchunks)) != 0) {
      if (rank==0) fprintf(stdout, "number of points not divisible by number of cpus * buffer size\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
   }
 
   csize = nx*ny*nz / nprocs / nchunks;
   if (rank==0) fprintf(stdout, "%d: csize=%ld\n", rank, csize);

   /* arrays for peak velocities, displacements */
   PH=(float*) calloc(csize, sizeof(float));
   PZ=(float*) calloc(csize, sizeof(float));
   PGV=(float*) calloc(csize, sizeof(float));
   PGA=(float*) calloc(csize, sizeof(float));
   AI=(float*) calloc(csize, sizeof(float));
   ENER=(float*) calloc(csize, sizeof(float));
   DUR=(float*) calloc(csize, sizeof(float));

   GM=(float**) calloc(nfreq, sizeof(float*)); 
   for (f=0; f < nfreq, comp_sa > 0; f++) GM[f] = (float*) calloc(csize, sizeof(float));


   if (rank==0) fprintf(stdout, "%d: Allocating buffers...\n", rank);
   /* buffers for unsorted velocities */
   bufx = (float*) calloc(csize*nt, sizeof(float));
   bufy = (float*) calloc(csize*nt, sizeof(float));
   bufz = (float*) calloc(csize*nt, sizeof(float));

   /* 2D arrays for velocities */
   x = (float**) calloc (csize, sizeof(float*));
   y = (float**) calloc (csize, sizeof(float*));
   z = (float**) calloc (csize, sizeof(float*));

   if (rank==0) fprintf(stdout, "%d: Allocating 2D arrays ...\n", rank);
   for (l=0; l<csize; l++){
      x[l] = (float*) calloc(nt, sizeof(float));
      y[l] = (float*) calloc(nt, sizeof(float));
      z[l] = (float*) calloc(nt, sizeof(float));
   }
   if (z[csize-1]==NULL) perror("error using calloc()");

   if (rank==0) fprintf(stdout, "%d: defining output files ...\n", rank);
   if (dofilt == 0){
      sprintf(phfile, "peak_velocity_H.bin");
      sprintf(pzfile, "peak_velocity_Z.bin");
      sprintf(pgvfile, "pgv.bin");
      sprintf(pgafile, "pga.bin");
      sprintf(aifile, "arias.bin");
      sprintf(enerfile, "ener.bin");
      sprintf(durfile, "dur.bin");
   }
   else {
       sprintf(phfile, "peak_velocity_H_%05.2f_%05.2fHz.bin", lp, hp);
       sprintf(pzfile, "peak_velocity_Z_%05.2f_%05.2fHz.bin", lp, hp);
       sprintf(pgvfile, "pgv_%05.2f_%05.2fHz.bin", lp, hp);
       sprintf(pgafile, "pga_%05.2f_%05.2fHz.bin", lp, hp);
       sprintf(aifile, "arias_%05.2f_%05.2fHz.bin", lp, hp);
       sprintf(durfile, "dur_%05.2f_%05.2fHz.bin", lp, hp);
       sprintf(enerfile, "ener_%05.2f_%05.2fHz.bin", lp, hp);
   }
   MPI_Barrier(MPI_COMM_WORLD);

   if (rank==0) fprintf(stdout, "%d: opening output files ...\n", rank);

   err=MPI_File_open(MPI_COMM_SELF, phfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &hfid);
   error_check(err, "MPI_File_open phfile");

   err=MPI_File_open(MPI_COMM_SELF, pzfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &zfid);
   error_check(err, "MPI_File_open pzfile");

   err=MPI_File_open(MPI_COMM_SELF, pgvfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &pgvfid);
   error_check(err, "MPI_File_open pgvfile");

   err=MPI_File_open(MPI_COMM_SELF,pgafile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &pgafid);
   error_check(err, "MPI_File_open pgafile");

   err=MPI_File_open(MPI_COMM_SELF, aifile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &aifid);
   error_check(err, "MPI_File_open aifile");

   err=MPI_File_open(MPI_COMM_SELF, durfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &durfid);
   error_check(err, "MPI_File_open durfile");
   err=MPI_File_open(MPI_COMM_SELF, enerfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &enerfid);
   error_check(err, "MPI_File_open enerfile");

   gmfid=(MPI_File*) calloc(nfreq, sizeof(MPI_File));

   for (f=0; f<nfreq, comp_sa > 0; f++){
      sprintf(gmfile, "gmrotD50_%05.2fHz.bin", freq[f]);
      err=MPI_File_open(MPI_COMM_SELF, gmfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
	  MPI_INFO_NULL, gmfid+f);
      error_check(err, "MPI_File_open gmfile");
   }

   MPI_Barrier(MPI_COMM_WORLD);

   accx=(float*) calloc(nt, sizeof(float));
   accy=(float*) calloc(nt, sizeof(float));
   accz=(float*) calloc(nt, sizeof(float));

   for (k=0; k<nchunks; k++){
      if (rank==0) fprintf(stdout, "processing part %02d of %02d\n", k+1, nchunks);
      s0=(long int) rank* (long int) nchunks* (long int) csize + (long int) k*csize;
      zi=s0 / (nx*ny);
      yi=(s0 % (nx*ny)) / nx;
      xi=s0 % nx;
      if (k == 0 && rank % 10 == 0) {
        fprintf(stdout, "(%d): s0=%d, xi=%d, yi=%d, zi=%d\n", rank, s0, xi, yi, zi);
        fflush(stdout);
      }
      if (par.itype == 0) {
         if (rank==0) fprintf(stdout, "reading %s ...\n", xfile);
         read_awp_timeseries(xfile, nx, ny, nz, xi, yi, zi, nt, 
             csize, bufx);
         if (rank==0) fprintf(stdout, "reading %s ...\n", yfile);
         read_awp_timeseries(yfile, nx, ny, nz, xi, yi, zi, nt, 
             csize, bufy);
         if (rank==0) fprintf(stdout, "reading %s ...\n", zfile);
         read_awp_timeseries(zfile, nx, ny, nz, xi, yi, zi, nt, 
             csize, bufz);
         if (rank==0) fprintf(stdout, "done reading.\n");
      }
      else {
         read_awp_timeseries_multi(xfile, nx, ny, nz, xi, yi, zi, nt, 
            par.writestep, par.ntiskp, csize, bufx);
         read_awp_timeseries_multi(yfile, nx, ny, nz, xi, yi, zi, nt, 
            par.writestep, par.ntiskp, csize, bufy);
         read_awp_timeseries_multi(zfile, nx, ny, nz, xi, yi, zi, nt, 
            par.writestep, par.ntiskp, csize, bufz);
      }

      MPI_Barrier(MPI_COMM_WORLD);

      if (rank==0) fprintf(stdout, "reshaping data ...\n");
      for (n=0; n<nt; n++){
         for (l=0; l<csize; l++){
            x[l][n]=bufx[n*csize+l];
            y[l][n]=bufy[n*csize+l];
            z[l][n]=bufz[n*csize+l];
         }
         //if (k==0) fprintf(stdout, "%e\n", x[3851][n]);
         //if (k==1) fprintf(stdout, "%e\n", x[3871][n]);
         //if (k==87) fprintf(stdout, "%e\n", x[87][n]);
      }

      if (rank==0) fprintf(stdout, "filtering time series ...\n");

      for (l=0; l<csize; l++){
         PH[l] = PZ[l] = PGV[l] = PGA[l] = AI[l] = DUR[l] = ENER[l] = 0.f;

         vel2acc(x[l], nt, dt2, accx);
         vel2acc(y[l], nt, dt2, accy);
         vel2acc(z[l], nt, dt2, accz);

         #ifdef DEBUG
         if (l==16){
            FILE *debugfid;
            debugfid=fopen("output.dbg", "w");
            for (n=0; n<nt; n++){
               fprintf(debugfid, "%.3e %.3e %.3e .3e\n", n*dt2, accx[n], accy[n], accz[n]);
            }
            fclose(debugfid);
         }
         #endif
         for (f=0; f<nfreq, comp_sa > 0; f++) {
             GM[f][l]=gmrotD50(accx, accy, nt, dt2, xid, freq[f], nang);
             #ifdef DEBUG
             if (l==16) fprintf(stderr, ">>> gmrotD50 @ f=%5.2f Hz: %e\n", freq[f], GM[f][l]);
             #endif
         }

         /* filter for PGVs and duration, after GMROT was computed */
         if (dofilt == 1) {
            xapiir_(x[l], &nt, aproto, &trbndw, &a, &iord, ftype, &lp, &hp, &dt2, &npas);
            xapiir_(y[l], &nt, aproto, &trbndw, &a, &iord, ftype, &lp, &hp, &dt2, &npas);
            xapiir_(z[l], &nt, aproto, &trbndw, &a, &iord, ftype, &lp, &hp, &dt2, &npas);
         }

         comp_cum(x[l], y[l], z[l], nt, dt2, low, high, DUR+l, ENER+l);

         for (n=0; n<nt; n++){
             vh=sqrtf(powf(x[l][n], 2.) + powf(y[l][n], 2.));
             if (PH[l] < vh) PH[l] = vh;
             vz=fabsf(z[l][n]);
             if (PZ[l] < vz) PZ[l] = vz;
             vel = sqrtf(powf(x[l][n], 2.) + powf(y[l][n], 2.) + powf(z[l][n], 2.));
             if (PGV[l] < vel) PGV[l] = vel;
             acc = sqrtf(powf(accx[n], 2.) + powf(accy[n], 2.) + powf(accz[n], 2.));
             if (PGA[l] < acc) PGA[l] = acc;
             AI[l] += powf(acc, 2.);
             // DX[l]+=x[l][n] * dt2;
             // DY[l]+=y[l][n] * dt2;
             // DZ[l]+=z[l][n] * dt2;
         }
         AI[l] = AI[l] * PI * dt2 / 2. / G / 3.;

      }

      if (rank==0) fprintf(stdout, "writing data to disk ...\n");
      MPI_Barrier(MPI_COMM_WORLD);
      off = (MPI_Offset) s0 * sizeof(float);
      MPI_File_write_at(hfid, off, PH, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at(zfid, off, PZ, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at(pgafid, off, PGA, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at(pgvfid, off, PGV, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at(aifid, off, AI, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at(durfid, off, DUR, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at(enerfid, off, ENER, csize, MPI_FLOAT, MPI_STATUS_IGNORE);

      for (f=0; f<nfreq, comp_sa > 0; f++)
         MPI_File_write_at(gmfid[f], off, GM[f], csize, MPI_FLOAT, MPI_STATUS_IGNORE);
    
   }
   free(accx);
   free(accy);
   free(accz);
   free(PH);
   free(PZ);
   free(PGA);
   free(PGV);
   free(AI);
   free(DUR);
   free(ENER);

   free(xfile);
   free(yfile);
   free(zfile);

   MPI_File_close(&hfid); 
   MPI_File_close(&zfid); 
   MPI_File_close(&pgafid); 
   MPI_File_close(&pgvfid); 
   MPI_File_close(&aifid); 
   MPI_File_close(&durfid); 
   MPI_File_close(&enerfid); 


   for (f=0; f<nfreq, comp_sa > 0; f++) {
     free(GM[f]);
     MPI_File_close(gmfid+f); 
   }
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(0);
}

