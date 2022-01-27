/* Computes peak ground velocity from AWP-ODC output.
   Improved version using MPI-IO.
   Daniel Roten, March 2015 <droten@sdsc.edu> 
   Zhifeng Hu, Feb 2021 <zhifeng94.hu@gmail.com> */

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
#define MAX(x, y) (((x)>(y))?(x):(y))
const int MAX_NUM_FILTER = 10;

// intercept MPI errors. Same usage as for CUCHK
#ifndef NDEBUG
#define MPICHK(err)                                                         \
    {                                                                       \
        if (err != MPI_SUCCESS)                                             \
        {                                                                   \
            char error_string[2048];                                        \
            int length_of_error_string;                                     \
            MPI_Error_string((err), error_string, &length_of_error_string); \
            fprintf(stderr, "MPI error: %s:%i %s(): %s\n",                  \
                    __FILE__, __LINE__, __func__, error_string);            \
            MPI_Abort(MPI_COMM_WORLD, err);                                 \
            fflush(stderr);                                                 \
            exit(EXIT_FAILURE);                                             \
        }                                                                   \
    }
#else
#define MPICHK(err) \
    {               \
    }
#endif

int index_above(float *in, int nt, float thres)
{
    int i = 0;
    for (i = 0; i < nt; ++i) {
        if (in[i] > thres) break;
    }
    return i;
}

int index_above2(float *in, int nt, float thres)
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

void comp_cum(float *accX, float *accY, float *accZ, int nt, float dt, float start, float end, float *dur, float *ener)
{
   float *cumeng;
   float start_bound, end_bound;
   int l;

   cumeng = (float *)calloc(nt, sizeof(float));
   cumeng[0] = pow(accX[0], 2.) + pow(accY[0], 2.) + pow(accZ[0], 2.);

   for (l = 1; l < nt; l++)
   {
      cumeng[l] = cumeng[l - 1] + pow(accX[l], 2.) + pow(accY[l], 2.) + pow(accZ[l], 2.);
   }
   start_bound = cumeng[nt - 1] * start;
   end_bound = cumeng[nt - 1] * end;

   int start_idx = index_above(cumeng, nt, start_bound);
   int end_idx = index_above(cumeng, nt, end_bound);

   *dur = (end_idx - start_idx) * dt;
   *ener = cumeng[nt - 1] / 3 * dt;
   free(cumeng);
}

void vel2acc(float *vel, int nt, float dt, float * acc)
{
   int n;

   for (n = 0; n < (nt - 1); n++)
   {
      acc[n] = (vel[n + 1] - vel[n]) / dt;
   }
   acc[nt - 1] = 0.;
}

int parsemultiple(char *optarg, float * val) {
    int k = 0;
    char *token;
    const char s[3] = ", ";
    token = strtok(optarg, s);
    while (token != NULL && k < MAX_NUM_FILTER) {
        val[k] = atof(token);
        token = strtok(NULL, s);
        k++;
    }
    return k;
}


int main(int argc, char *argv[])
{
   FD3D_param par;
   read_settings(&par, "IN3D.out");
   char *xfile, *yfile, *zfile;
   int nx = 1, ny = 1, nt = 1, nz = 1;
   float **x, **y, **z;
   float *velx, *vely, *velz;
   float *accx, *accy, *accz, acc;
   float *bufx, *bufy, *bufz;
   float vh, vz, vel;
   MPI_File *hfid, *zfid, *pgvfid, *pgafid;
   MPI_File *aifid, *durfid, *enerfid;
   MPI_Status filestatus;
   int nchunks = 18, csize; //nchunks=15
   int rank, nprocs;
   int i, j, k, l, n, t;
   MPI_Offset off;
   long int s0;
   int xi, yi, zi;
   float dt2;

   /*parameters for xapiir*/
   int dofilt = 0;
   int comp_sa = 0;
   int iord = 3, npas = 1;
   float trbndw = 0., a = 0.; /* chebyshev parameters */
   char *aproto = "BU";
   char ftype[MAX_NUM_FILTER][3];
   int nfilter = 1;;
   float lp = 0.15, hp = 0.0;
   float lps[MAX_NUM_FILTER], hps[MAX_NUM_FILTER];
   for (i = 0; i < MAX_NUM_FILTER; ++i) {
       lps[i] = lp;
       hps[i] = hp;
   }
   float start = 0.05, end = 0.75;

   /* parameters for getopt */
   int c;
   extern char *optarg;

   char phfile[200], pzfile[200], pgvfile[200], pgafile[200], aifile[200], durfile[200], enerfile[200], gmfile[200];
   int err;

   float freq[] = {0.1, 0.333333, 0.5, 0.6666666, 1.0, 1.5, 2.0, 3.0, 3.5};
   int nfreq = (int)sizeof(freq) / sizeof(freq[0]);

   int f;
   float xid = 0.05; // damping constant for gmrotD50
   int nang = 91;    // number of angles for gmrotD50

   MPI_File *gmfid;

   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   while ((c = getopt(argc, argv, "fgl:h:n:s:e:")) != -1)
   {
      switch (c)
      {
      case 'f':
         dofilt = 1;
         break;
      case 'g':
         comp_sa = 1;
         break;
      case 'l':
         parsemultiple(optarg, lps);
         break;
      case 'h':
         nfilter = parsemultiple(optarg, hps);
         break;
      case 'n':
         sscanf(optarg, "%d", &nchunks);
         break;
      case 's':
         sscanf(optarg, "%f", &start);
         break;
      case 'e':
         sscanf(optarg, "%f", &end);
         break;
      case '?':
         fprintf(stdout, "usage: %s [-f] [-l lp0, lp1, lp2...] [-h hp0, hp1, hp2...] [-n nchunks] [-s start] [-e end]\n", argv[0]);
         exit(0);
      default:
         abort();
      }
   }
    for (i = 0; i < nfilter; ++i) {
        if (hps[i] > 0.)
            strncpy(ftype[i], "BP\0", 3);
        else
            strncpy(ftype[i], "LP\0", 3);
    }
   /* determine dimensions from these parameters */
   nx = (int)floorf((par.nedx - par.nbgx) / par.nskpx + 1.);
   ny = (int)floorf((par.nedy - par.nbgy) / par.nskpy + 1.);
   nz = (int)floorf((par.nedz - par.nbgz) / par.nskpz + 1.);
   nt = (int)floorf((par.tmax / par.dt / par.ntiskp));
   dt2 = par.dt * par.ntiskp;

   /* output parameters for debug */
   if (rank == 0)
   {
    fprintf(stdout, "  >> nx=%d, ny=%d, nz=%d, nt=%d,"
                      " ivelocity=%d, ntiskp=%d, writestep=%d\n",
              nx, ny, nz, nt, par.itype, par.ntiskp, par.writestep);
    fprintf(stdout, "nchunks=%d, start=%f, end=%f, dofilt=%d\n", nchunks, start, end, dofilt);
    for (i = 0; i < nfilter && dofilt; ++i) {
        fprintf(stdout, "Filter %d (%s): low=%f, high=%f\n", i, ftype[i], lps[i], hps[i]);
    }
    fflush(stdout);
   }

   xfile = (char *)calloc(100, sizeof(char));
   yfile = (char *)calloc(100, sizeof(char));
   zfile = (char *)calloc(100, sizeof(char));
   sscanf(par.sxrgo, "%s", xfile);
   sscanf(par.syrgo, "%s", yfile);
   sscanf(par.szrgo, "%s", zfile);
   /* this removes the remaining ' at the end of the string, if any */
   xfile = strsep(&xfile, "\'");
   yfile = strsep(&yfile, "\'");
   zfile = strsep(&zfile, "\'");

   /* number of points to process per read operation */
   if (((nx * ny * nz) % (nprocs * nchunks)) != 0)
   {
      if (rank == 0)
         fprintf(stdout, "number of points not divisible by number of cpus * buffer size\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
   }

   csize = nx * ny * nz / nprocs / nchunks;
   if (rank == 0) {
      fprintf(stdout, "%d: csize=%ld\n", rank, csize);
      fflush(stdout);
   }

   /* arrays for peak velocities, displacements */
   float (*PH)[csize] = calloc(nfilter, sizeof(*PH)); // [nfilter, csize]
   float (*PZ)[csize] = calloc(nfilter, sizeof(*PZ));
   float (*PGV)[csize] = calloc(nfilter, sizeof(*PGV));
   float (*PGA)[csize] = calloc(nfilter, sizeof(*PGA));
   float (*AI)[csize] = calloc(nfilter, sizeof(*AI));
   float (*ENER)[csize] = calloc(nfilter, sizeof(*ENER));
   float (*DUR)[csize] = calloc(nfilter, sizeof(*DUR));
   float (*GM)[csize] = calloc(1, sizeof(float[nfreq][csize]));

   /* 2D arrays for velocities */
   x = (float **)calloc(csize, sizeof(float *));
   y = (float **)calloc(csize, sizeof(float *));
   z = (float **)calloc(csize, sizeof(float *));

   if (rank == 0)
      fprintf(stdout, "%d: Allocating 2D arrays ...\n", rank);
   for (l = 0; l < csize; l++)
   {
      x[l] = (float *)calloc(nt, sizeof(float));
      y[l] = (float *)calloc(nt, sizeof(float));
      z[l] = (float *)calloc(nt, sizeof(float));
   }
   if (z[csize - 1] == NULL)
      perror("error using calloc()");
   /* buffers for unsorted velocities */
   bufx = (float *)calloc(csize * nt, sizeof(float));
   bufy = (float *)calloc(csize * nt, sizeof(float));
   bufz = (float *)calloc(csize * nt, sizeof(float));
   accx = (float *)calloc(nt, sizeof(float));
   accy = (float *)calloc(nt, sizeof(float));
   accz = (float *)calloc(nt, sizeof(float));
   velx = (float *)calloc(nt, sizeof(float));
   vely = (float *)calloc(nt, sizeof(float));
   velz = (float *)calloc(nt, sizeof(float));

   if (rank == 0)
      fprintf(stdout, "%d: defining output files ...\n", rank);
    hfid = (MPI_File *)calloc(nfilter, sizeof(MPI_File));
    zfid = (MPI_File *)calloc(nfilter, sizeof(MPI_File));
    pgvfid = (MPI_File *)calloc(nfilter, sizeof(MPI_File));
    pgafid = (MPI_File *)calloc(nfilter, sizeof(MPI_File));
    aifid = (MPI_File *)calloc(nfilter, sizeof(MPI_File));
    durfid = (MPI_File *)calloc(nfilter, sizeof(MPI_File));
    enerfid = (MPI_File *)calloc(nfilter, sizeof(MPI_File));
    for (i = 0; i < nfilter; ++i) {
        if (dofilt == 0) {
           sprintf(phfile, "pgv_H.bin");
           sprintf(pzfile, "pgv_Z.bin");
           sprintf(pgvfile, "pgv.bin");
           sprintf(pgafile, "pga.bin");
           sprintf(aifile, "arias.bin");
           sprintf(enerfile, "ener.bin");
           sprintf(durfile, "dur.bin");
        }
        else
        {
            sprintf(phfile, "pgv_H_%05.2f_%05.2fHz.bin", lps[i], hps[i]);
            sprintf(pzfile, "pgv_Z_%05.2f_%05.2fHz.bin", lps[i], hps[i]);
            sprintf(pgvfile, "pgv_%05.2f_%05.2fHz.bin", lps[i], hps[i]);
            sprintf(pgafile, "pga_%05.2f_%05.2fHz.bin", lps[i], hps[i]);
            sprintf(aifile, "arias_%05.2f_%05.2fHz.bin", lps[i], hps[i]);
            sprintf(durfile, "dur_%05.2f_%05.2fHz.bin", lps[i], hps[i]);
            sprintf(enerfile, "ener_%05.2f_%05.2fHz.bin", lps[i], hps[i]);
        }
        MPICHK(MPI_File_open(MPI_COMM_SELF, phfile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                             MPI_INFO_NULL, hfid + i));

        MPICHK(MPI_File_open(MPI_COMM_SELF, pzfile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                             MPI_INFO_NULL, zfid + i));

        MPICHK(MPI_File_open(MPI_COMM_SELF, pgvfile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                             MPI_INFO_NULL, pgvfid + i));

        MPICHK(MPI_File_open(MPI_COMM_SELF, pgafile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                             MPI_INFO_NULL, pgafid + i));

        MPICHK(MPI_File_open(MPI_COMM_SELF, aifile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                             MPI_INFO_NULL, aifid + i));

        MPICHK(MPI_File_open(MPI_COMM_SELF, durfile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                             MPI_INFO_NULL, durfid + i));

        MPICHK(MPI_File_open(MPI_COMM_SELF, enerfile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                             MPI_INFO_NULL, enerfid + i));
    }

   gmfid = (MPI_File *)calloc(nfreq, sizeof(MPI_File));

   for (f = 0; f<nfreq && comp_sa> 0; f++)
   {
      sprintf(gmfile, "gmrotD50_%05.2fHz.bin", freq[f]);
      MPICHK(MPI_File_open(MPI_COMM_SELF, gmfile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                    MPI_INFO_NULL, gmfid + f));
   }

   MPI_Barrier(MPI_COMM_WORLD);


   for (k = 0; k < nchunks; k++)
   {
      if (rank == 0)
         fprintf(stdout, "processing part %02d of %02d\n", k + 1, nchunks);
      fflush(stdout);
      s0 = (long int)rank * (long int)nchunks * (long int)csize + (long int)k * csize;
      zi = s0 / (nx * ny);
      yi = (s0 % (nx * ny)) / nx;
      xi = s0 % nx;
      if (k == 0 && rank % 10 == 0)
      {
         fprintf(stdout, "(%d): s0=%d, xi=%d, yi=%d, zi=%d\n", rank, s0, xi, yi, zi);
         fflush(stdout);
      }
      if (par.itype == 0)
      {
         read_awp_timeseries(xfile, nx, ny, nz, xi, yi, zi, nt,
                             csize, bufx);
         read_awp_timeseries(yfile, nx, ny, nz, xi, yi, zi, nt,
                             csize, bufy);
         read_awp_timeseries(zfile, nx, ny, nz, xi, yi, zi, nt,
                             csize, bufz);
         if (rank == 0)
            fprintf(stdout, "done reading.\n");
      }
      else
      {
         read_awp_timeseries_multi(xfile, nx, ny, nz, xi, yi, zi, nt,
                                   par.writestep, par.ntiskp, csize, bufx);
         read_awp_timeseries_multi(yfile, nx, ny, nz, xi, yi, zi, nt,
                                   par.writestep, par.ntiskp, csize, bufy);
         read_awp_timeseries_multi(zfile, nx, ny, nz, xi, yi, zi, nt,
                                   par.writestep, par.ntiskp, csize, bufz);
      }

      MPI_Barrier(MPI_COMM_WORLD);

      if (rank == 0)
      {
         fprintf(stdout, "Done. Reshaping data ...\n");
         fflush(stdout);
      }

      for (t = 0; t < nt; t++)
      {
         for (l = 0; l < csize; l++)
         {
            x[l][t] = bufx[t * csize + l];
            y[l][t] = bufy[t * csize + l];
            z[l][t] = bufz[t * csize + l];
         }
      }

      if (rank == 0)
         fprintf(stdout, "filtering time series ...\n");

      for (l = 0; l < csize; l++)
      {

         vel2acc(x[l], nt, dt2, accx);
         vel2acc(y[l], nt, dt2, accy);
         vel2acc(z[l], nt, dt2, accz);

#ifdef DEBUG
         if (l == 16)
         {
            FILE *debugfid;
            debugfid = fopen("output.dbg", "w");
            for (t = 0; t < nt; t++)
            {
               fprintf(debugfid, "%.3e %.3e %.3e .3e\n", t * dt2, accx[n], accy[n], accz[n]);
            }
            fclose(debugfid);
         }
#endif
         for (f = 0; f<nfreq && comp_sa> 0; f++)
         {
            GM[f][l] = gmrotD50(accx, accy, nt, dt2, xid, freq[f], nang);
#ifdef DEBUG
            if (l == 16)
               fprintf(stderr, ">>> gmrotD50 @ f=%5.2f Hz: %e\n", freq[f], GM[f][l]);
#endif
         }

         /* filter for PGVs and duration, after GMROT was computed */
         for (i = 0; i < nfilter; ++i) {
             PH[i][l] = PZ[i][l] = PGV[i][l] = PGA[i][l] = AI[i][l] = DUR[i][l] = ENER[i][l] = 0.f;
             vel2acc(x[l], nt, dt2, accx);
             vel2acc(y[l], nt, dt2, accy);
             vel2acc(z[l], nt, dt2, accz);
             memcpy(velx, x[l], nt * sizeof(x[l][0]));
             memcpy(vely, y[l], nt * sizeof(y[l][0]));
             memcpy(velz, z[l], nt * sizeof(z[l][0]));
            if (dofilt == 1)
            {
                xapiir_(velx, &nt, aproto, &trbndw, &a, &iord, ftype[i], &lps[i], &hps[i], &dt2, &npas);
                xapiir_(vely, &nt, aproto, &trbndw, &a, &iord, ftype[i], &lps[i], &hps[i], &dt2, &npas);
                xapiir_(velz, &nt, aproto, &trbndw, &a, &iord, ftype[i], &lps[i], &hps[i], &dt2, &npas);
                xapiir_(accx, &nt, aproto, &trbndw, &a, &iord, ftype[i], &lps[i], &hps[i], &dt2, &npas);
                xapiir_(accy, &nt, aproto, &trbndw, &a, &iord, ftype[i], &lps[i], &hps[i], &dt2, &npas);
                xapiir_(accz, &nt, aproto, &trbndw, &a, &iord, ftype[i], &lps[i], &hps[i], &dt2, &npas);
            }

            comp_cum(velx, vely, velz, nt, dt2, start, end, DUR[i] + l, ENER[i] + l);

            for (t = 0; t < nt; t++)
            {
                vh = sqrtf(powf(velx[t], 2.) + powf(vely[t], 2.));
                PH[i][l] = MAX(PH[i][l], vh);
                vz = fabsf(z[l][t]);
                PZ[i][l] = MAX(PZ[i][l], vz);
                vel = sqrtf(powf(velx[t], 2.) + powf(vely[t], 2.) + powf(velz[t], 2.));
                PGV[i][l] = MAX(PGV[i][l], vel);
                acc = sqrtf(powf(accx[t], 2.) + powf(accy[t], 2.) + powf(accz[t], 2.));
                PGA[i][l] = MAX(PGA[i][l], acc);
                AI[i][l] += powf(acc, 2.);
            }
            AI[i][l] = AI[i][l] * PI * dt2 / 2. / G / 3.;
        }
      }

    off = (MPI_Offset)s0 * sizeof(float);
    if (rank == 0) {
        fprintf(stdout, "writing data to disk ...\n");
        fflush(stdout);
    }
    for (i = 0; i < nfilter; ++i) {
        MPI_File_write_at(hfid[i], off, PH[i], csize, MPI_FLOAT, MPI_STATUS_IGNORE);
        MPI_File_write_at(zfid[i], off, PZ[i], csize, MPI_FLOAT, MPI_STATUS_IGNORE);
        MPI_File_write_at(pgafid[i], off, PGA[i], csize, MPI_FLOAT, MPI_STATUS_IGNORE);
        MPI_File_write_at(pgvfid[i], off, PGV[i], csize, MPI_FLOAT, MPI_STATUS_IGNORE);
        MPI_File_write_at(aifid[i], off, AI[i], csize, MPI_FLOAT, MPI_STATUS_IGNORE);
        MPI_File_write_at(durfid[i], off, DUR[i], csize, MPI_FLOAT, MPI_STATUS_IGNORE);
        MPI_File_write_at(enerfid[i], off, ENER[i], csize, MPI_FLOAT, MPI_STATUS_IGNORE);
    }
    for (f = 0; f < nfreq && comp_sa > 0; f++)
    {
        MPI_File_write_at(gmfid[f], off, GM[f], csize, MPI_FLOAT, MPI_STATUS_IGNORE);
    }


    MPI_Barrier(MPI_COMM_WORLD);
   }
    for (i = 0; i < csize; ++i)
    {
       free(x[i]); free(y[i]); free(z[i]);
    }
    free(x); free(y); free(z);
   free(velx); free(vely); free(velz);
   free(accx); free(accy); free(accz);
   free(bufx); free(bufy); free(bufz);
   free(PH);
   free(PZ);
   free(PGA);
   free(PGV);
   free(AI);
   free(DUR);
   free(ENER);
   free(GM);

   free(xfile); free(yfile); free(zfile);
    for (i = 0; i < nfilter; ++i) {
        MPI_File_close(hfid + i);
        MPI_File_close(zfid + i);
        MPI_File_close(pgafid + i);
        MPI_File_close(pgvfid + i);
        MPI_File_close(aifid + i);
        MPI_File_close(durfid + i);
        MPI_File_close(enerfid + i);
    }
   for (f = 0; f<nfreq && comp_sa> 0; f++)
   {
      MPI_File_close(gmfid + f);
   }
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(0);
}
