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

int index_above(float *in, int nt, float thres)
{
   int l = 0, r = nt - 1;
   while (l < r)
   {
      int mid = l + (r - l) / 2;
      if (in[mid] < thres)
      {
         l = mid + 1;
      }
      else
      {
         r = mid;
      }
   }
   return l;
}

void ds_5_95(float *accX, float *accY, int nt, float dt, float *durX, float *durY, float *dur2c){
   float *cumengX, *cumengY, *cumeng2c;
   int l, p, n0[3], n1[3];
   float dur[3];
   float e5[3], e95[3];

   cumengX = (float*) calloc(nt, sizeof(float));
   cumengY = (float*) calloc(nt, sizeof(float));
   cumeng2c = (float*) calloc(nt, sizeof(float));

   cumengX[0] = pow(accX[0], 2.);
   cumengY[0] = pow(accY[0], 2.);
   cumeng2c[0] = pow(accX[0], 2.) + pow(accY[0], 2.);
   for (l=1; l<nt; l++){
      cumengX[l] = cumengX[l-1] + pow(accX[l], 2.);
      cumengY[l] = cumengY[l-1] + pow(accY[l], 2.);
      cumeng2c[l] = cumeng2c[l-1] + pow(accX[l], 2.) + pow(accY[l], 2.);
   }

   e5[0] = cumengX[nt-1] * 0.05;
   e5[1] = cumengY[nt-1] * 0.05;
   e5[2] = cumeng2c[nt-1] * 0.05;

   e95[0]=cumengX[nt-1] * 0.95;
   e95[1]=cumengY[nt-1] * 0.95;
   e95[2]=cumeng2c[nt-1] * 0.95;

   n0[0] = index_above(cumengX, nt, e5[0]);
   n0[1] = index_above(cumengY, nt, e5[1]);
   n0[2] = index_above(cumeng2c, nt, e5[2]);

   n1[0] = index_above(cumengX, nt, e95[0]);
   n1[1] = index_above(cumengY, nt, e95[1]);
   n1[2] = index_above(cumeng2c, nt, e95[2]);

   free(cumengX);
   free(cumengY);
   free(cumeng2c);

   dur[0]=(n1[0]-n0[0]) * dt;
   dur[1]=(n1[1]-n0[1]) * dt;
   dur[2]=(n1[2]-n0[2]) * dt;

   (*durX) = dur[0];
   (*durY) = dur[1];
   (*dur2c) = dur[2];

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
   float *bufx, *bufy, *bufz;
   float *PH, *PZ, vh, vz;
   float *DX, *DY, *DZ;
   MPI_File hfid, zfid;
   MPI_File dxfid, dyfid, dzfid;
   int nchunks = 4, csize; //nchunks=15
   int rank, nprocs;
   int k, l, n;
   MPI_Offset off;
   long int s0;
   int xi, yi, zi;
   float dt2;

   /*parameters for xapiir*/
   int dofilt=0;
   int iord=3, npas=1;
   float trbndw=0., a=0.; /* chebyshev parameters */
   char *aproto="BU";
   char ftype[3];
   float hp=0., lp=1.00; 

   /* parameters for getopt */
   int c;
   extern char *optarg;

   char phfile[200], pzfile[200], dxfile[200], dyfile[200], dzfile[200], durfile[200], gmfile[200];
   int err;

   float freq[]={0.1, 0.2, 0.25, 0.333333, 0.5, 0.6666666, 0.75, 1.0, 1.333333, 1.5};
   int nfreq=10;

   float **GM, *DUR2C;  //gmrotD50 and duration
   int f;
   float xid = 0.05; // damping constant for gmrotD50
   int nang=91;  // number of angles for gmrotD50
   float *accx, *accy;

   float dummy1, dummy2; //for duration in X and Y direction (not saved)
   MPI_File durfid, *gmfid;
 
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  while ((c = getopt (argc, argv, "fl:h:")) != -1) {
     switch(c) {
        case 'f':
           dofilt = 1;
           break;
        case 'l':
           sscanf(optarg, "%f", &lp);
           break;
        case 'h':
           sscanf(optarg, "%f", &hp);
           break;
        case '?':
           fprintf(stdout, "usage: %s [-f] [-l lp] [-h hp]\n", argv[0]);
           exit(0);
        default:
           abort();
     }
  }
  if (hp > 0.) strncpy(ftype, "BP\0", 3);
  else strncpy(ftype, "LP\0", 3);
  if (rank==0) fprintf(stdout, "dofilt=%d, lp=%f, hp=%f, ftype=%s\n", 
          dofilt, lp, hp, ftype);

   /* determine dimensions from these parameters */
   nx = (int)floorf((par.nedx - par.nbgx) / par.nskpx + 1.);
   ny = (int)floorf((par.nedy - par.nbgy) / par.nskpy + 1.);
   nz = (int)floorf((par.nedz - par.nbgz) / par.nskpz + 1.);
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
   DX=(float*) calloc(csize, sizeof(float));
   DY=(float*) calloc(csize, sizeof(float));
   DZ=(float*) calloc(csize, sizeof(float));

   GM=(float**) calloc(nfreq, sizeof(float*)); 
   for (f=0; f < nfreq; f++) GM[f] = (float*) calloc(csize, sizeof(float));

   DUR2C=(float*) calloc(csize, sizeof(float));

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
      sprintf(dxfile, "displacement_X.bin");
      sprintf(dyfile, "displacement_Y.bin");
      sprintf(dzfile, "displacement_Z.bin");
      sprintf(durfile, "duration_2c.bin");
   }
   else {
      sprintf(phfile, "peak_velocity_H_%04.1fHz.bin", lp);
      sprintf(pzfile, "peak_velocity_Z_%04.1fHz.bin", lp);
      sprintf(dxfile, "displacement_X_%04.1fHz.bin", lp);
      sprintf(dyfile, "displacement_Y_%04.1fHz.bin", lp);
      sprintf(dzfile, "displacement_Z_%04.1fHz.bin", lp);
      sprintf(durfile, "duration_2c_%04.1fHz.bin", lp);
   }
   MPI_Barrier(MPI_COMM_WORLD);

   if (rank==0) fprintf(stdout, "%d: opening output files ...\n", rank);

   err=MPI_File_open(MPI_COMM_SELF, phfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &hfid);
   error_check(err, "MPI_File_open phfile");

   err=MPI_File_open(MPI_COMM_SELF, pzfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &zfid);
   error_check(err, "MPI_File_open pzfile");

   err=MPI_File_open(MPI_COMM_SELF, dxfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &dxfid);
   error_check(err, "MPI_File_open dxfile");

   err=MPI_File_open(MPI_COMM_SELF, dyfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &dyfid);
   error_check(err, "MPI_File_open dyfile");

   err=MPI_File_open(MPI_COMM_SELF, dzfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &dzfid);
   error_check(err, "MPI_File_open dzfile");

   err=MPI_File_open(MPI_COMM_SELF, durfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &durfid);
   error_check(err, "MPI_File_open durfile");

   gmfid=(MPI_File*) calloc(nfreq, sizeof(MPI_File));
   for (f=0; f<nfreq; f++){
      sprintf(gmfile, "gmrotD50_%05.2fHz.bin", freq[f]);
      err=MPI_File_open(MPI_COMM_SELF, gmfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
	  MPI_INFO_NULL, gmfid+f);
      error_check(err, "MPI_File_open gmfile");
   }

   MPI_Barrier(MPI_COMM_WORLD);

   accx=(float*) calloc(nt, sizeof(float));
   accy=(float*) calloc(nt, sizeof(float));

   for (k=0; k<nchunks; k++){
      if (rank==0) fprintf(stdout, "processing part %02d of %02d\n", k+1, nchunks);
      s0=(long int) rank* (long int) nchunks* (long int) csize + (long int) k*csize;
      zi=s0 / (nx*ny);
      yi=(s0 % (nx*ny)) / nx;
      xi=s0 % nx;
      fprintf(stdout, "(%d): s0=%d, xi=%d, yi=%d, zi=%d\n", rank, s0, xi, yi, zi);
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
         PH[l] = PZ[l] = DX[l] = DY[l] = DZ[l] = DUR2C[l]= 0.f;

         vel2acc(x[l], nt, dt2, accx);
         vel2acc(y[l], nt, dt2, accy);

         #ifdef DEBUG
         if (l==16){
            FILE *debugfid;
            debugfid=fopen("output.dbg", "w");
            for (n=0; n<nt; n++){
               fprintf(debugfid, "%e %e %e\n", n*dt2, accx[n], accy[n]);
            }
            fclose(debugfid);
         }
         #endif
         for (f=0; f<nfreq; f++) {
             GM[f][l]=gmrotD50(accx, accy, nt, dt2, xid, freq[f], nang);
             #ifdef DEBUG
             if (l==16) fprintf(stderr, ">>> gmrotD50 @ f=%5.2f Hz: %e\n", freq[f], GM[f][l]);
             #endif
         }

         /* filter for PGVs and duration, after GMROT was computed */
         if (dofilt == 1) {
            xapiir_(x[l], &nt, aproto, &trbndw, &a, &iord, ftype, &hp, &lp, &dt2, &npas);
            xapiir_(y[l], &nt, aproto, &trbndw, &a, &iord, ftype, &hp, &lp, &dt2, &npas);
            xapiir_(z[l], &nt, aproto, &trbndw, &a, &iord, ftype, &hp, &lp, &dt2, &npas);
         }

         ds_5_95(x[l], y[l], nt, dt2, &dummy1, &dummy2, DUR2C+l);

         for (n=0; n<nt; n++){
             vh=sqrtf(powf(x[l][n], 2.) + powf(y[l][n], 2.));
             if (PH[l] < vh) PH[l] = vh;
             vz=fabsf(z[l][n]);
             if (PZ[l] < vz) PZ[l] = vz;
             DX[l]+=x[l][n]*par.dt*par.ntiskp;
             DY[l]+=y[l][n]*par.dt*par.ntiskp;
             DZ[l]+=z[l][n]*par.dt*par.ntiskp;
         }
      }

      if (rank==0) fprintf(stdout, "writing data to disk ...\n");
      MPI_Barrier(MPI_COMM_WORLD);
      off = (MPI_Offset) s0 * sizeof(float);
      MPI_File_write_at(hfid, off, PH, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at(zfid, off, PZ, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at(dxfid, off, DX, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at(dyfid, off, DY, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at(dzfid, off, DZ, csize, MPI_FLOAT, MPI_STATUS_IGNORE);

      MPI_File_write_at(durfid, off, DUR2C, csize, MPI_FLOAT, MPI_STATUS_IGNORE);

      for (f=0; f<nfreq; f++)
         MPI_File_write_at(gmfid[f], off, GM[f], csize, MPI_FLOAT, MPI_STATUS_IGNORE);
    
   }
   MPI_File_close(&hfid); 
   MPI_File_close(&zfid); 
   MPI_File_close(&dxfid); 
   MPI_File_close(&dyfid); 
   MPI_File_close(&dzfid); 

   MPI_File_close(&durfid); 

   for (f=0; f<nfreq; f++)
      MPI_File_close(gmfid+f); 

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(0);
}

