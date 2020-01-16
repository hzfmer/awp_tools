/* Computes peak ground velocity from AWP-ODC output.
   Improved version using MPI-IO.
   Daniel Roten, March 2015 <droten@sdsc.edu> */

#include "awp_data.h"
#include "xapiir.h"
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdlib.h>

int main (int argc, char *argv[] ){
   char *xfile, *yfile, *zfile;
   //int nx=1180, ny=1, nt=800, nz=200;
   int nx, ny=1, nt=1800, nz;
   int x0, x1, z0, z1, y0;
   FILE *fid;
   float **x, **z;
   float *bufx, *bufz;
   float *DX, *DZ;
   MPI_File dxfid, dzfid;
   int nchunks = 10, csize;
   int rank, nprocs;
   int k, l, n;
   MPI_Offset off;
   long int s0;
   int xi, yi, zi;
   float dt2=0.05;
   int ntiskp=10;
   int itype=0, writestep=100;

   /*parameters for xapiir*/
   int dofilt=0;
   int iord=3, npas=1;
   float trbndw=0., a=0.; /* chebyshev parameters */
   char *aproto="BU";
   char *ftype="LP";
   float hp=0., lp=1.00; 

   /* parameters for getopt */
   int c;
   extern char *optarg;

   char dxfile[200], dzfile[200];
 
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  while ((c = getopt (argc, argv, "fl:")) != -1) {
     switch(c) {
        case 'f':
           dofilt = 1;
           break;
        case 'l':
           sscanf(optarg, "%f", &lp);
           break;
        case '?':
           fprintf(stdout, "usage: %s [-f] [-l lp]\n", argv[0]);
           exit(0);
        default:
           abort();
     }
  }
  if (rank==0) fprintf(stdout, "dofilt=%d, lp=%f\n", dofilt, lp);

  fid=fopen("fault_output.dat", "r");
  if (fid==NULL) {
     perror("could not open fault_output.dat");
     exit(-1);
  }

  fscanf(fid, "%d %d\n", &x0, &x1);
  fscanf(fid, "%d\n", &y0);
  fscanf(fid, "%d %d\n", &z0, &z1);

  nx=x1-x0+1;
  nz=z1-z0+1;

   /* output parameters for debug */
   fprintf(stdout, "  >> nx=%d, ny=%d, nz=%d, nt=%d," 
                   " ivelocity=%d, ntiskp=%d, writestep=%d\n",
		    nx, ny, nz, nt, itype, ntiskp, writestep);

   xfile=(char*) calloc(100, sizeof(char) );
   zfile=(char*) calloc(100, sizeof(char) );
   /* this removes the remaining ' at the end of the string, if any */
   sprintf(xfile, "output_dyn/ratu");
   sprintf(zfile, "output_dyn/ratw");

   /* number of points to process per read operation */
   if (((nx*ny*nz) % (nprocs * nchunks)) != 0) {
      fprintf(stdout, "number of points not divisible by number of cpus * buffer size\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
   }
 
   csize = nx*ny*nz / nprocs / nchunks;

   /* arrays for peak velocities, displacements */
   DX=(float*) calloc(csize, sizeof(float));
   DZ=(float*) calloc(csize, sizeof(float));

   /* buffers for unsorted velocities */
   bufx = (float*) calloc(csize*nt, sizeof(float));
   bufz = (float*) calloc(csize*nt, sizeof(float));

   x = (float**) calloc(csize, sizeof(float*));
   z = (float**) calloc(csize, sizeof(float*));
   for (l=0; l<csize; l++){
      x[l] = (float*) calloc(nt, sizeof(float));
      z[l] = (float*) calloc(nt, sizeof(float));
   }
   if (z[csize-1]==NULL) perror("error using calloc()");

   if (dofilt == 0){
      sprintf(dxfile, "output_dyn/slipu.bin");
      sprintf(dzfile, "output_dyn/slipw.bin");
   }
   else {
      sprintf(dxfile, "output_dyn/slipu_%04.1fHz.bin", lp);
      sprintf(dzfile, "output_dyn/slipw_%04.1fHz.bin", lp);
   }

   MPI_File_open(MPI_COMM_WORLD, dxfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &dxfid);
   MPI_File_open(MPI_COMM_WORLD, dzfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &dzfid);

   for (k=0; k<nchunks; k++){
      if (rank==0) fprintf(stdout, "processing part %02d of %02d\n", k+1, nchunks);
      s0=(long int) rank* (long int) nchunks* (long int) csize + (long int) k*csize;
      zi=s0 / (nx*ny);
      yi=(s0 % (nx*ny)) / nx;
      xi=s0 % nx;
      //fprintf(stdout, "(%d): s0=%d, xi=%d, yi=%d, zi=%d\n", rank, s0, xi, yi, zi);
      if (itype == 0) {
         read_awp_timeseries(xfile, nx, ny, nz, xi, yi, zi, nt, 
             csize, bufx);
         read_awp_timeseries(zfile, nx, ny, nz, xi, yi, zi, nt, 
             csize, bufz);
      }
      else {
         read_awp_timeseries_multi(xfile, nx, ny, nz, xi, yi, zi, nt, 
            writestep, ntiskp, csize, bufx);
         read_awp_timeseries_multi(zfile, nx, ny, nz, xi, yi, zi, nt, 
            writestep, ntiskp, csize, bufz);
      }

      MPI_Barrier(MPI_COMM_WORLD);

      for (n=0; n<nt; n++){
         for (l=0; l<csize; l++){
            x[l][n]=bufx[n*csize+l];
            z[l][n]=bufz[n*csize+l];
         }
         //if (k==0) fprintf(stdout, "%e\n", x[3851][n]);
         //if (k==1) fprintf(stdout, "%e\n", x[3871][n]);
         //if (k==87) fprintf(stdout, "%e\n", x[87][n]);
      }

      for (l=0; l<csize; l++){
         if (dofilt == 1) {
            xapiir_(x[l], &nt, aproto, &trbndw, &a, &iord, ftype, &hp, &lp, &dt2, &npas);
            xapiir_(z[l], &nt, aproto, &trbndw, &a, &iord, ftype, &hp, &lp, &dt2, &npas);
         }

         DX[l] = DZ[l] = 0;
         for (n=0; n<nt; n++){
             DX[l]+=x[l][n]*dt2;
             DZ[l]+=z[l][n]*dt2;
         }
      }

      MPI_Barrier(MPI_COMM_WORLD);
      off = (MPI_Offset) s0 * sizeof(float);
      MPI_File_write_at_all(dxfid, off, DX, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at_all(dzfid, off, DZ, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
   }
   MPI_File_close(&dxfid); 
   MPI_File_close(&dzfid); 

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(0);
}

