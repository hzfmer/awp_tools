#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <mpi.h>

#include "awp_data.h"
#include "xapiir.h"

int main (int argc, char *argv[]) {
   int nx, ny=1, nz, nt;
   int x0, x1, z0, z1, y0;
   FILE *ffid;
   int rank, nprocs;
   int nstat;
   int nchunks=5, csize;
   char *xfile="output_dyn/ratu", *zfile="output_dyn/ratw";
   long fsize;
   int itype=0; 
   int k, l, n;
   float *PH, *PZ, *PX, vh, vx, vz;
   float *bufx, *bufz, **x, **z;
   char pxfile[100], pzfile[100], pfile[100];
   long int s0;
   int xi, yi, zi;
   MPI_Offset off;
   MPI_File pxfid, pzfid, pfid;

   /*hardcoded parameters */
   float dt2=0.; 
   int writestep = 100, ntiskp=10;

   /*parameters for xapiir*/
   int dofilt=0;
   int iord=3, npas=1;
   float trbndw=0., a=0.; /* chebyshev parameters */
   char *aproto="BU";
   char *ftype="LP";
   float hp=0., lp=1.;

   /* parameters for getopt */
   int c;
   extern char *optarg;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   while ((c = getopt (argc, argv, "fl:d:")) != -1) {
     switch(c) {
        case 'f':
           dofilt = 1;
           break;
        case 'l':
           sscanf(optarg, "%f", &lp);
           break;
        case 'd':
           sscanf(optarg, "%f", &dt2);
           break;
        case '?':
           fprintf(stdout, "usage: %s [-f] [-l lp] [-d dt]\n", argv[0]);
           exit(0);
        default:
           abort();
     }
  }
  if ((dt2 == 0.) && (dofilt == 1)){
     if (rank==0) fprintf(stderr, "Error: need to define dt [-d] to apply filter\n");
     MPI_Finalize();
     exit(0); 
  }
  if (rank==0) fprintf(stdout, "dofilt=%d, lp=%f\n", dofilt, lp);


   if (rank==0) {
      ffid=fopen("fault_output.dat", "r");
      fscanf(ffid, "%d %d\n", &x0, &x1);
      fscanf(ffid, "%d\n", &y0);
      fscanf(ffid, "%d %d\n", &z0, &z1);
      nx=x1-x0+1;
      nz=z1-z0+1;
      fclose(ffid);

      ffid=fopen(xfile, "r");
      fseek(ffid, 0L, SEEK_END);
      fsize=ftell(ffid); 
      fclose(ffid);
      nt=fsize / sizeof(float) / nx / nz;
   } 
   MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&nz, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&nt, 1, MPI_INT, 0, MPI_COMM_WORLD);

   nstat=nx*nz;

   csize = nx*ny*nz / nprocs / nchunks;
   if (rank==0) fprintf(stdout, "nx=%d, nz=%d, nt=%d, csize=%d\n", nx, nz, nt, csize);

   /* arrays for peak slip rates */
   PH=(float*) calloc(csize, sizeof(float));
   PZ=(float*) calloc(csize, sizeof(float));
   PX=(float*) calloc(csize, sizeof(float));

   /* 2D arrays for slip rates */
   x = (float**) calloc (csize, sizeof(float*));
   z = (float**) calloc (csize, sizeof(float*));

   for (l=0; l<csize; l++){
      x[l] = (float*) calloc(nt, sizeof(float));
      z[l] = (float*) calloc(nt, sizeof(float));
   }

   bufx=(float*) calloc(csize*nt, sizeof(float));
   bufz=(float*) calloc(csize*nt, sizeof(float));

   if (dofilt == 0){
      sprintf(pxfile, "output_dyn/maxrate_x.bin");
      sprintf(pzfile, "output_dyn/maxrate_z.bin");
      sprintf(pfile, "output_dyn/maxrate.bin");
   }
   else {
      sprintf(pxfile, "output_dyn/maxrate_x_%04.1fHz.bin", lp);
      sprintf(pzfile, "output_dyn/maxrate_z_%04.1fHz.bin", lp);
      sprintf(pfile, "output_dyn/maxrate_%04.1fHz.bin", lp);
   }

   MPI_File_open(MPI_COMM_WORLD, pxfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &pxfid);
   MPI_File_open(MPI_COMM_WORLD, pzfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &pzfid);
   MPI_File_open(MPI_COMM_WORLD, pfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &pfid);

   for (k=0; k<nchunks; k++){
      fprintf(stdout, "k=%d\n", k);
      s0=(long int) rank* (long int) nchunks* (long int) csize + (long int) k*csize;
      zi=s0 / nx;
      yi=0;
      xi=s0 % nx;

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
      }

      for (l=0; l<csize; l++){
         if (dofilt == 1) {
            xapiir_(x[l], &nt, aproto, &trbndw, &a, &iord, ftype, &hp, &lp, &dt2, &npas);
            xapiir_(z[l], &nt, aproto, &trbndw, &a, &iord, ftype, &hp, &lp, &dt2, &npas);
         }

         PH[l] = PZ[l] = PX[l] = 0;
         for (n=0; n<nt; n++){
             vh=sqrtf(powf(x[l][n], 2.) + powf(z[l][n], 2.));
             if (PH[l] < vh) PH[l] = vh;
             vx=fabsf(x[l][n]);
             if (PX[l] < vx) PX[l] = vx;
             vz=fabsf(z[l][n]);
             if (PZ[l] < vz) PZ[l] = vz;
         }
      }

      MPI_Barrier(MPI_COMM_WORLD);
      off = (MPI_Offset) s0 * sizeof(float);
      MPI_File_write_at_all(pfid, off, PH, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at_all(pxfid, off, PX, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      MPI_File_write_at_all(pzfid, off, PZ, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
   }

   MPI_File_close(&pxfid); 
   MPI_File_close(&pzfid); 
   MPI_File_close(&pfid); 

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   return(0);

}

