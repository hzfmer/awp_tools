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
   int nchunks = 3, csize; //nchunks=15
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

   char phfile[200], pzfile[200], dxfile[200], dyfile[200], dzfile[200];
   int err;
 
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
   nx=(int) floorf( (par.nedx-par.nbgx)/par.nskpx + 1);
   ny=(int) floorf( (par.nedy-par.nbgy)/par.nskpy + 1);
   nz=(int) floorf( (par.nedz-par.nbgz)/par.nskpz + 1);
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
   }
   else {
      sprintf(phfile, "peak_velocity_H_%04.1fHz.bin", lp);
      sprintf(pzfile, "peak_velocity_Z_%04.1fHz.bin", lp);
      sprintf(dxfile, "displacement_X_%04.1fHz.bin", lp);
      sprintf(dyfile, "displacement_Y_%04.1fHz.bin", lp);
      sprintf(dzfile, "displacement_Z_%04.1fHz.bin", lp);
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

   MPI_Barrier(MPI_COMM_WORLD);

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
         if (dofilt == 1) {
            dt2=par.dt*par.ntiskp;
            xapiir_(x[l], &nt, aproto, &trbndw, &a, &iord, ftype, &hp, &lp, &dt2, &npas);
            xapiir_(y[l], &nt, aproto, &trbndw, &a, &iord, ftype, &hp, &lp, &dt2, &npas);
            xapiir_(z[l], &nt, aproto, &trbndw, &a, &iord, ftype, &hp, &lp, &dt2, &npas);
         }

         PH[l] = PZ[l] = DX[l] = DY[l] = DZ[l] = 0;
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
   }
   MPI_File_close(&hfid); 
   MPI_File_close(&zfid); 
   MPI_File_close(&dxfid); 
   MPI_File_close(&dyfid); 
   MPI_File_close(&dzfid); 

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(0);
}

