/* Computes peak ground velocity from AWP-ODC output.
   Improved version using MPI-IO.
   Daniel Roten, March 2015 <droten@sdsc.edu> */

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

#define PI 3.1415927
#define G 9.8

int main (int argc, char *argv[] ){
   FD3D_param par;
   read_settings(&par);
   char *xfile, *yfile, *zfile;
   int nx=1, ny=1, nt=1, nz=1;
   float **x, **y, **z;
   float *PH, *PZ, *PGV, valh, valz, val;
   float vx, vy, vz;
   float tmpvx, tmpvy, tmpvz;
   float *bufx, *bufy, *bufz;
   float *DUR, durx, dury, durz;
   float *AI, *PGA, acc;
   MPI_File hfid, zfid, pgvfid, pgafid;
   MPI_File aifid, durfid;
   int nchunks = 8, csize;
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
   char *ftype="BP";
   float lp=0.15, hp=5.00; 

   /* parameters for getopt */
   int c;
   extern char *optarg;

   char phfile[200], pzfile[200], pgvfile[200], pgafile[200];
   char aifile[200], durfile[200];
   int err;
 
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  while ((c = getopt (argc, argv, "fl:h:")) != -1) {
     switch(c) {
        case 'f':
           dofilt = 1;
           break;
        case 'h':
           sscanf(optarg, "%f", &hp);
           break;
        case 'l':
           sscanf(optarg, "%f", &lp);
           break;
        case '?':
           fprintf(stdout, "usage: %s [-f] [-l lp] [-h hp]\n", argv[0]);
           exit(0);
        default:
           abort();
     }
  }
  if (rank==0) fprintf(stdout, "dofilt=%d, lp = %f, hp=%f\n", dofilt, lp, hp);

   /* determine dimensions from these parameters */
   nx=(int) floorf( (par.nedx-par.nbgx)/par.nskpx + 1.); 
   ny=(int) floorf( (par.nedy-par.nbgy)/par.nskpy + 1.); 
   nz=(int) floorf( (par.nedz-par.nbgz)/par.nskpz + 1.); 
   nt=(int) floorf( (par.tmax/par.dt/par.ntiskp));

   
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
      fprintf(stdout, "number of points not divisible by number of cpus * buffer size\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
   }
 
   csize = nx*ny*nz / nprocs / nchunks;
   fprintf(stdout, "%d: csize=%ld\n", rank, csize);

   /* arrays for peak velocities, displacements */
   PH = (float*) calloc(csize, sizeof(float));
   PZ = (float*) calloc(csize, sizeof(float));
   PGV = (float*) calloc(csize, sizeof(float));
   PGA = (float*) calloc(csize, sizeof(float));
   AI=(float*) calloc(csize, sizeof(float));
   DUR=(float*) calloc(csize, sizeof(float));
   fprintf(stdout, "%d: Allocating buffers...\n", rank);
   /* buffers for unsorted velocities */
   bufx = (float*) calloc(csize*nt, sizeof(float));
   bufy = (float*) calloc(csize*nt, sizeof(float));
   bufz = (float*) calloc(csize*nt, sizeof(float));

   /* 2D arrays for velocities */
   x = (float**) calloc (csize, sizeof(float*));
   y = (float**) calloc (csize, sizeof(float*));
   z = (float**) calloc (csize, sizeof(float*));

   fprintf(stdout, "%d: Allocating 2D arrays ...\n", rank);
   for (l=0; l<csize; l++){
      x[l] = (float*) calloc(nt, sizeof(float));
      y[l] = (float*) calloc(nt, sizeof(float));
      z[l] = (float*) calloc(nt, sizeof(float));
   }
   if (z[csize-1]==NULL) perror("error using calloc()");

   fprintf(stdout, "%d: defining output files ...\n", rank);
   if (dofilt == 0){
      sprintf(phfile, "peak_velocity_H.bin");
      sprintf(pzfile, "peak_velocity_Z.bin");
      sprintf(pgvfile, "pgv.bin");
      sprintf(pgafile, "pga.bin");
      sprintf(aifile, "arias.bin");
      sprintf(durfile, "dur.bin");
   }
   else {
      sprintf(phfile, "peak_velocity_H_%05.2fHz.bin", hp);
      sprintf(pzfile, "peak_velocity_Z_%05.2fHz.bin", hp);
      sprintf(pgvfile, "pgv_%05.2fHz.bin", hp);
      sprintf(pgafile, "pga_%05.2fHz.bin", hp);
      sprintf(aifile, "arias_%05.2fHz.bin", hp);
      sprintf(durfile, "dur_%05.2fHz.bin", hp);
   }
   MPI_Barrier(MPI_COMM_WORLD);

   fprintf(stdout, "%d: opening output files ...\n", rank);

   err=MPI_File_open(MPI_COMM_WORLD, phfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &hfid);
   error_check(err, "MPI_File_open phfile");

   err=MPI_File_open(MPI_COMM_WORLD, pzfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &zfid);
   error_check(err, "MPI_File_open pzfile");

   err=MPI_File_open(MPI_COMM_WORLD, pgvfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &pgvfid);
   error_check(err, "MPI_File_open pgvfile");

   err=MPI_File_open(MPI_COMM_WORLD, pgafile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &pgafid);
   error_check(err, "MPI_File_open pgafile");

   err=MPI_File_open(MPI_COMM_WORLD, aifile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &aifid);
   error_check(err, "MPI_File_open dxfile");

   err=MPI_File_open(MPI_COMM_WORLD, durfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &durfid);
   error_check(err, "MPI_File_open dyfile");


   MPI_Barrier(MPI_COMM_WORLD);

   for (k=0; k<nchunks; k++){
      if (rank==0) fprintf(stdout, "processing part %02d of %02d\n", k+1, nchunks);
      s0=(long int) rank* (long int) nchunks* (long int) csize + (long int) k*csize;
      zi=s0 / (nx*ny);
      yi=(s0 % (nx*ny)) / nx;
      xi=s0 % nx;
      //fprintf(stdout, "(%d): s0=%d, xi=%d, yi=%d, zi=%d\n", rank, s0, xi, yi, zi);
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
         if (dofilt == 1) {
            dt2=par.dt*par.ntiskp;
            xapiir_(x[l], &nt, aproto, &trbndw, &a, &iord, ftype, &lp, &hp, &dt2, &npas);
            xapiir_(y[l], &nt, aproto, &trbndw, &a, &iord, ftype, &lp, &hp, &dt2, &npas);
            xapiir_(z[l], &nt, aproto, &trbndw, &a, &iord, ftype, &lp, &hp, &dt2, &npas);
         }

         PH[l] = PZ[l] = PGV[l] = PGA[l] = AI[l] = 0;
         vx = vy = vz = 0;
         tmpvx = tmpvy = tmpvz = 0;
         durx = dury = durz = 0;
         x[l][0] = y[l][0] = z[l][0] = 0;
         for (n=1; n<nt; n++){
             valh = sqrtf(powf(x[l][n], 2.) + powf(y[l][n], 2.));
             if (PH[l] < valh) PH[l] = valh;
             valz = fabsf(z[l][n]);
             if (PZ[l] < valz) PZ[l] = valz;
             val = sqrtf(powf(x[l][n], 2.) + powf(y[l][n], 2.) + powf(z[l][n], 2.));
             if (PGV[l] < val) PGV[l] = val;
             acc = sqrtf(powf((x[l][n] - x[l][n - 1]) / (par.dt*par.ntiskp), 2.)
                       + powf((y[l][n] - y[l][n - 1]) / (par.dt*par.ntiskp), 2.)
                       + powf((z[l][n] - z[l][n - 1]) / (par.dt*par.ntiskp), 2.));
             if (PGA[l] < acc) PGA[l] = acc;
             AI[l] += powf(acc, 2.);
             vx += powf(x[l][n], 2.);
             vy += powf(y[l][n], 2.);
             vz += powf(z[l][n], 2.);
             if (rank==0 && n % 100 == 0 && l % 100 == 0) {
               fprintf(stdout, "Processing %d / %d, l: %d / %d\n", n, nt, l, csize);
             }
         }
         // Cumulative energy, define as vel**2 from 5% to 75%
         if (vx > 0) {
             for (n=1; n<nt; n++){
                 tmpvx += powf(x[l][n], 2.);
                 if (vx * 0.75 > tmpvx && vx * 0.05 <= tmpvx) {
                     durx += 1;
                 }
                 tmpvy += powf(y[l][n], 2.);
                 if (vy * 0.75 > tmpvy && vy * 0.05 <= tmpvy) {
                     dury += 1;
                 }
                 tmpvz += powf(z[l][n], 2.);
                 if (vz * 0.75 > tmpvz && vz * 0.05 <= tmpvz) {
                     durz += 1;
                 }
             }
         }
         AI[l] = AI[l] * PI * par.dt * par.ntiskp / 2. / G / 3.;
         DUR[l] = (durx + dury + durz) * par.dt * par.ntiskp / 3.;

      }

      if (rank==0) fprintf(stdout, "writing data to disk ...\n");
      MPI_Barrier(MPI_COMM_WORLD);
      off = (MPI_Offset) s0 * sizeof(float);
      MPI_File_write_at_all(hfid, off, PH, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at_all(zfid, off, PZ, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at_all(pgvfid, off, PGV, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at_all(pgafid, off, PGA, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at_all(aifid, off, AI, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at_all(durfid, off, DUR, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
   }
   MPI_File_close(&hfid); 
   MPI_File_close(&zfid); 
   MPI_File_close(&pgvfid); 
   MPI_File_close(&pgafid); 
   MPI_File_close(&aifid); 
   MPI_File_close(&durfid); 

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(0);
}

