/* This code reads the surface velocity produced by the AWM and
   calculates the orientation-independent horizontal response spectra
   gmrotD50 (Boore et all., 2006, BSSA 96(4) 1502-1511) 

   Daniel Roten, Apr 2010 <droten@sciences.sdsu.edu> 
*/

#include "srfdata2.h"
#include "fd3dparam.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include "gmrotD50.h"

int index_above(float *in, int nt, float thres){
   int l;
   int idx;
   for (l=0; l<nt; l++){
      if (in[l] > thres) {
         idx=l;
         break;
      }
   }
   return idx;
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
void bracket_duration(float *accX, float *accY, int nt, float dt, float thres, 
   float *durX, float *durY, float *dur2c){
   float maxacc;
   int l, p, n0[3], n1[3];
   float dur[3];

   for (p=0; p<3; p++) {
     n0[p]=-1;
     n1[p]=-1;
   }

   for (l=0; l<nt; l++){
      if (fabs(accX[l]) > thres){
         n0[0] = l;
         break;
      }
   }

   for (l=0; l<nt; l++){
      if (fabs(accY[l]) > thres){
         n0[1] = l;
         break;
      }
   }

   for (l=0; l<nt; l++){
      maxacc=sqrt(powf(accX[l], 2.) + pow(accY[l], 2.));
      if (maxacc > thres){
         n0[2] = l;
         break;
      }
   }

   for (l=nt-1; l>=0; l--){
      if (fabs(accX[l]) > thres){
         n1[0] = l;
         break;
      }
   }

   for (l=nt-1; l>=0; l--){
      if (fabs(accY[l]) > thres){
         n1[1] = l;
         break;
      }
   }

   for (l=nt-1; l>=0; l--){
      maxacc=sqrt(powf(accX[l], 2.) + pow(accY[l], 2.));
      if (maxacc > thres){
         n1[2] = l;
         break;
      }
   }

   dur[0]=(n1[0]-n0[0]) * dt;
   dur[1]=(n1[1]-n0[1]) * dt;
   dur[2]=(n1[2]-n0[2]) * dt;

   (*durX) = dur[0];
   (*durY) = dur[1];
   (*dur2c) = dur[2];

}

void vel2acc(float ***D, int nrows, int nx, int nt, float dt){
   int k, l, n;
   for (l=0; l<nrows; l++){
      for (k=0; k<nx; k++){
         for (n=0; n<nt-1; n++){
            D[l][k][n] = (D[l][k][n+1] - D[l][k][n]) / dt;
         }
         D[l][k][nt-1] = 0.;
      }
   }
}

int main (int argc, char *argv[] ){
   FD3D_param par;
   read_settings(&par);
   char *xfile, *yfile, *zfile;
   int ix, iy;
   int nx=1, ny=1, nt=1, nz=1;
   FILE *hfid;
   struct srfdata X, Y, Z;
   float ***PH;
   float ***peakh;
   float ***x, ***y;
   int nt2;
   float dt2;
   clock_t c0, c1, c2, c3;
   int nrows=1, k, l;
   float default_freq[]={0.1, 0.2, 0.25, 0.333333, 0.5, 0.6666666, 0.75, 1.0, 1.5};
   float *freq;
   int nfreq=9;
   float xi=0.05;
   char rfileh[200];
   float *tmprsp;
   int nang=91;
   #ifdef I50  
   float *fa, *iresp, minang;
   int nf=100, fidx[5], a;
   #endif

   #ifdef DUR
   char durfilex[100], durfiley[100], durfile2c[100];
   FILE *durfidx, *durfidy, *durfid2c;
   float **durX, **durY, **dur2C;
   float **durx, **dury, **dur2c;
   #endif
   #ifdef DBRACK
   float dthres = 0.05 * 9.81; /*0.05g acc threshold for duration*/
   #endif

   tmprsp=(float*) calloc(1, sizeof(float));

  int rank, ncpus, cpn;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&ncpus);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if (argc != 1){
     nfreq=argc-1;
     freq= (float*) calloc(nfreq, sizeof(float));
     for (l=0; l<nfreq;l++) freq[l] = (float) atof(argv[l+1]);
  }
  else {
     nfreq=9;
     freq= (float*) calloc(nfreq, sizeof(float));
     for (l=0; l<nfreq; l++) freq[l]=default_freq[l];
  }

  if (rank==0) {
    fprintf(stdout, "SAs to compute (Hz): ");
    for (l=0; l<nfreq; l++) fprintf(stdout, " %6.3f", freq[l]);
    fprintf(stdout, "\n");
  }

  /* determine dimensions from these parameters */
  nx=(int) floorf( (par.nedx-par.nbgx)/par.nskpx + 1.); 
  ny=(int) floorf( (par.nedy-par.nbgy)/par.nskpy + 1.); 
  nz=(int) floorf( (par.nedz-par.nbgz)/par.nskpz + 1.); 
  nt=(int) floorf( (par.tmax/par.dt));

   
   /* output parameters for debug */
   if (rank==0){
      fprintf(stdout, "  >> nx=%d, ny=%d, nz=%d, nt=%d, ivelocity=%d, ntiskp=%d, wstep=%d\n",
		nx, ny, nz, nt, par.itype, par.ntiskp, par.writestep);
   }

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

   if (rank==0){
      PH=(float***) calloc(nfreq, sizeof(float**));
      for (l=0; l<nfreq; l++){
         PH[l]=(float**) calloc(ny, sizeof(float*));
         if (PH==NULL) perror("error using calloc()");
         for (iy=0; iy<ny; iy++){
            PH[l][iy]=(float*) calloc(nx, sizeof(float));
            if (PH[l][iy]==NULL) perror("error using calloc()");
         }
      }
   }

   /*3-D array: size (nfreq, nrows, nx) */
   peakh=(float***) calloc(nfreq, sizeof(float**));
   for (l=0; l<nfreq; l++){
      peakh[l]=(float**) calloc(nrows, sizeof(float*));
      for (k=0; k<nrows; k++){
         peakh[l][k]=(float*) calloc(nx, sizeof(float));
      }
   }

   #ifdef DUR
   /*duration: arrays with global values */
   if (rank==0){
      durX=(float**) calloc(ny, sizeof(float*));
      durY=(float**) calloc(ny, sizeof(float*));
      dur2C=(float**) calloc(ny, sizeof(float*));
      if (dur2C==NULL) perror("error using calloc()");
      for (iy=0; iy<ny; iy++){
	 durX[iy]=(float*) calloc(nx, sizeof(float));
	 durY[iy]=(float*) calloc(nx, sizeof(float));
	 dur2C[iy]=(float*) calloc(nx, sizeof(float));
	 if (dur2C[iy]==NULL) perror("error using calloc()");
      }
   }

   /*duration: arrays with local values */
   durx=(float**) calloc(nrows, sizeof(float*));
   dury=(float**) calloc(nrows, sizeof(float*));
   dur2c=(float**) calloc(nrows, sizeof(float*));
   if (dur2c==NULL) perror("error using calloc()");
   for (iy=0; iy<nrows; iy++){
      durx[iy]=(float*) calloc(nx, sizeof(float));
      dury[iy]=(float*) calloc(nx, sizeof(float));
      dur2c[iy]=(float*) calloc(nx, sizeof(float));
      if (dur2c[iy]==NULL) perror("error using calloc()");
   }
   #endif

   X = srfdata_init(nx, ny, nt, par.dt, par.itype, par.writestep, par.ntiskp, xfile);
   Y = srfdata_init(nx, ny, nt, par.dt, par.itype, par.writestep, par.ntiskp, yfile);
   Z = srfdata_init(nx, ny, nt, par.dt, par.itype, par.writestep, par.ntiskp, zfile);

   /* true timestep and nt in skipped data */
   nt2= nt / par.ntiskp;
   dt2= par.dt * par.ntiskp;
   if (rank==0) fprintf(stdout, "dt2=%f\n", dt2);

   /* allocates 3 (nrows, nx, nt) arrays to contain serveral rows 
     of time series */
   x=(float***) calloc(nrows, sizeof(float**));
   y=(float***) calloc(nrows, sizeof(float**));
   for (k=0; k<nrows; k++){
      x[k]=(float**) calloc(nx, sizeof(float*));
      y[k]=(float**) calloc(nx, sizeof(float*));
      for (ix=0; ix<nx; ix++){
         x[k][ix]=(float*) calloc(nt2, sizeof(float));
         y[k][ix]=(float*) calloc(nt2, sizeof(float));
      }
   }

   #ifdef I50
   fa=logspace(0.10, 1.0, nf);
   for (l=0; l<nfreq; l++){
      for (a=1; a < nf; a++){
         if ((fa[a-1] < freq[l] ) && (fa[a] >= freq[l])) 
             fidx[l]=a;
      }
   }
   if (rank==0) {
      for (l=0; l<nfreq; l++) {
         fprintf(stdout, "fidx[%d]=%d, fa[%d]=%e\n", l, fidx[l], fidx[l], fa[fidx[l]]);
      }
   }
   #endif

   MPI_Barrier(MPI_COMM_WORLD);
   for (iy=(rank*nrows); iy < ny; iy+=(ncpus * nrows)){
      /*fprintf(stdout, "\rcomputing shock response (%3d%%)",
          (int) (((float)iy + 1) / (float) ny * 100.));
      fflush(stdout);*/
      //if (rank==0)
         fprintf(stdout, "processing rows %d to %d on cpu %d\n", iy, iy+nrows-1, rank);

      MPI_Barrier(MPI_COMM_WORLD);
      c0=clock();
      srfdata_getrows_wstep(&X, iy, nrows, x);
      srfdata_getrows_wstep(&Y, iy, nrows, y);
      MPI_Barrier(MPI_COMM_WORLD);
      c1=clock();
      if (rank==0){
         fprintf(stdout, "  time for I/O: %f seconds\n", (float) (c1-c0)/CLOCKS_PER_SEC);
      }

      /* convert velocities to acceleration */
      vel2acc(x, nrows, nx, nt2, dt2);
      vel2acc(y, nrows, nx, nt2, dt2);

      for (k=0; k<nrows; k++){
         for (ix=0; ix < nx; ix++){
            if (rank==0)
               if ((ix % 100) == 0){
                  fprintf(stdout, "processing col %d / %d, row %d / %d\n", k, nrows, ix, nx);
                  fflush(stdout);
               }
            #ifdef I50
            iresp=gmrotIpp(x[k][ix], y[k][ix], nt2, dt2, fa, nf, xi, 
                 nang, &minang);
            #endif
            for (l=0; l<nfreq; l++){
               #ifdef I50
               peakh[l][k][ix]=iresp[fidx[l]];
               #else
               peakh[l][k][ix]=gmrotD50(x[k][ix], y[k][ix], nt2, dt2, 
                  xi, freq[l], nang);
               #endif
               #ifdef DEBUG2
               fprintf(stdout, "CPU: %d @ (%d/%d): resp=%e\n", 
                    rank, k, ix, peakh[l][k][ix]);
               #endif
            }
            #ifdef I50
            free(iresp);
            #endif
            #ifdef D595
            ds_5_95(x[k][ix], y[k][ix], nt2, dt2, durx[k]+ix, dury[k]+ix, dur2c[k]+ix);
            #endif
            #ifdef DBRACK
            bracket_duration(x[k][ix], y[k][ix], nt2, dt2, dthres, durx[k]+ix, dury[k]+ix, dur2c[k]+ix);
            #endif
         }
      }
      MPI_Barrier(MPI_COMM_WORLD);
      c2=clock();
      if (rank==0){
         fprintf(stdout, "  time for computing response: %f seconds\n", 
              (float) (c2-c1)/CLOCKS_PER_SEC);
         fflush(stdout);
      }
      
      if (rank==0){
         for (l=0; l<nfreq; l++){
            for (k=0; k<nrows; k++){
               for (ix=0; ix < nx; ix++){
                  PH[l][iy+k][ix]=peakh[l][k][ix];
               }
            }
         }
         #ifdef DUR
	 for (k=0; k<nrows; k++){
	    for (ix=0; ix < nx; ix++){
	       durX[iy+k][ix]=durx[k][ix];
	       durY[iy+k][ix]=dury[k][ix];
	       dur2C[iy+k][ix]=dur2c[k][ix];
	    }
	 }
         #endif
         for (cpn=1; cpn<ncpus; cpn++){
            fprintf(stdout, "receiving SAs from CPU %d...\n", cpn);
            fflush(stdout);
            for (l=0; l<nfreq; l++){
               for (k=0; k<nrows; k++){
                  MPI_Recv(PH[l][iy+cpn*nrows+k], nx, MPI_FLOAT, cpn, 90, MPI_COMM_WORLD, &status);
                  #ifdef DEBUG2
                  fprintf(stdout, "CPU %d: PH[%d][%d][%d]=%e\n", 
                       rank, l, iy+cpn*nrows+k, 500, PH[l][iy+cpn*nrows+k][500]);
                  #endif
               }
            }
            #ifdef DUR
	    for (k=0; k<nrows; k++){
	       MPI_Recv(durX[iy+cpn*nrows+k], nx, MPI_FLOAT, cpn, 91, MPI_COMM_WORLD, &status);
	       MPI_Recv(durY[iy+cpn*nrows+k], nx, MPI_FLOAT, cpn, 92, MPI_COMM_WORLD, &status);
	       MPI_Recv(dur2C[iy+cpn*nrows+k], nx, MPI_FLOAT, cpn, 93, MPI_COMM_WORLD, &status);
	       #ifdef DEBUG2
	       fprintf(stdout, "CPU %d: durX[%d][%d][%d]=%e\n", 
		    rank, l, iy+cpn*nrows+k, 500, durX[iy+cpn*nrows+k][500]);
	       #endif
	    }
            #endif
         }
      }
      else {
          for (l=0; l<nfreq; l++){
             for (k=0; k<nrows; k++){
               MPI_Send(peakh[l][k], nx, MPI_FLOAT, 0, 90, MPI_COMM_WORLD);
             }
          }
          #ifdef DUR
	  for (k=0; k<nrows; k++){
	     MPI_Send(durx[k], nx, MPI_FLOAT, 0, 91, MPI_COMM_WORLD);
	     MPI_Send(dury[k], nx, MPI_FLOAT, 0, 92, MPI_COMM_WORLD);
	     MPI_Send(dur2c[k], nx, MPI_FLOAT, 0, 93, MPI_COMM_WORLD);
	  }
          #endif
      }
      MPI_Barrier(MPI_COMM_WORLD);
      c3=clock();
      if (rank==0){
         fprintf(stdout, "  time for comm: %f seconds\n", (float) (c3-c2)/CLOCKS_PER_SEC);
         fflush(stdout);
      }
   }
   MPI_Barrier(MPI_COMM_WORLD);

   if (rank==0){
      for (l=0; l<nfreq; l++){
      #ifdef I50
      //sprintf(rfileh, "gmrotI50_%05.2fHz.dat", freq[l]);
      sprintf(rfileh, "gmrotI50_%05.2fHz.bin", freq[l]);
      #else
      //sprintf(rfileh, "gmrotD50_%05.2fHz.dat", freq[l]);
      sprintf(rfileh, "gmrotD50_%05.2fHz.bin", freq[l]);
      #endif
      hfid=fopen(rfileh, "w");
      /*for (ix=0; ix < nx; ix++){
         for (iy=0; iy < ny; iy++){
           fprintf(hfid, "%e ", PH[l][iy][ix]);
         }
         fprintf(hfid, "\n");
      }*/
      for (iy=0; iy < ny; iy++){
         fwrite(PH[l][iy], nx, sizeof(float), hfid);
      }
      fclose(hfid);
      }

      #ifdef D595
      sprintf(durfilex, "duration_X.bin");
      sprintf(durfiley, "duration_Y.bin");
      sprintf(durfile2c, "duration_2c.bin");
      #endif
      #ifdef DBRACK
      sprintf(durfilex, "duration_bracket_X.bin");
      sprintf(durfiley, "duration_bracket_Y.bin");
      sprintf(durfile2c, "duration_bracket_2c.bin");
      #endif

      #ifdef DUR
      durfidx=fopen(durfilex, "w");
      durfidy=fopen(durfiley, "w");
      durfid2c=fopen(durfile2c, "w");

      for (iy=0; iy < ny; iy++){
         fwrite(durX[iy], nx, sizeof(float), durfidx);
         fwrite(durY[iy], nx, sizeof(float), durfidy);
         fwrite(dur2C[iy], nx, sizeof(float), durfid2c);
      }

      fclose(durfidx);
      fclose(durfidy);
      fclose(durfid2c);
      #endif
   }

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return(0);
}
