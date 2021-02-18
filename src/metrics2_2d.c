/* Computes peak ground velocity from AWP-ODC output.
   Improved version using MPI-IO.
   Daniel Roten, March 2015 <droten@sdsc.edu> */

/* This version also computes gmrotD50 and duration, in addition to 
 * PGVs. */

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

// intercept MPI errors. Same usage as for CUCHK
#ifndef NDEBUG
#define MPICHK(err) {                                                         \
    if (err != MPI_SUCCESS) {                                                    \
        char error_string[2048];                                                     \
        int length_of_error_string;                                                  \
        MPI_Error_string((err), error_string, &length_of_error_string);              \
        fprintf(stderr, "MPI error: %s:%i %s(): %s\n",                               \
                    __FILE__, __LINE__, __func__, error_string);                         \
        MPI_Abort(MPI_COMM_WORLD, err);                                              \
        fflush(stderr);                                                              \
        exit(EXIT_FAILURE);                                                          \
    }                                                                             \
}
#else
#define MPICHK(err) {}
#endif

struct Mpi
{
    // MPI vars
    int rank;
    int size;
    int err;
    int dim[2];
    int period[2];
    int reorder;
    int nxt, nyt, ntt;
    int coord[2];
    MPI_Comm MCW, MC1;
};

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

void comp_cum(float *accX, float *accY, float *accZ, int nt, float dt, float low, float high, float *dur, float *ener)
{
    float *cumeng;
    float low_bound, high_bound;
    int l;

    cumeng = (float *)calloc(nt, sizeof(float));
    cumeng[0] = pow(accX[0], 2.) + pow(accY[0], 2.) + pow(accZ[0], 2.);

    for (l = 1; l < nt; l++)
    {
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

void vel2acc(float *vel, int nt, float dt, float *acc)
{
    int n;

    for (n = 0; n < (nt - 1); n++)
    {
        acc[n] = (vel[n + 1] - vel[n]) / dt;
    }
    acc[nt - 1] = 0.;
}

void read_timeseries(char *fbase, int step, struct Mpi *m, MPI_Datatype dtype, float* output, int buffer_size) { 
    char fname[200];
    sprintf(fname, "%s%07d", fbase, step);
    MPI_File fp;
    MPI_Status   filestatus;
    char mpiErrStr[100];
    int mpiErrStrLen;
    MPICHK(MPI_File_open(m->MCW, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp));
    MPICHK(MPI_File_set_view(fp, 0, MPI_FLOAT, dtype, "native", MPI_INFO_NULL));
    fprintf(stdout, "%s,rank=%d, buffer_size=%d \n", fname, m->rank, buffer_size);
    fflush(stdout);
    MPICHK(MPI_File_read_all(fp, output, buffer_size, MPI_FLOAT, &filestatus));
}

void mpi_init(struct Mpi *m, int argc, char **argv);
void mpi_cart(struct Mpi *m, const int *size, const int *part);
MPI_Datatype data_type(const struct Mpi *mpi, int nx, int ny, int nt);


int main(int argc, char *argv[])
{
    struct Mpi m;
    mpi_init(&m, argc, argv);
    FD3D_param par;
    read_settings(&par, "IN3D.out");
    char *xfile, *yfile, *zfile;
    int nx = 1, ny = 1, nz = 1, nt = 1;
    float **x, **y, **z;
    float *accx, *accy, *accz, acc;
    float *bufx, *bufy, *bufz;
    float *PH, *PZ, *PGV, vh, vz, vel;
    float *AI, *PGA, *ENER, *DUR;
    MPI_File hfid, zfid, pgvfid, pgafid;
    MPI_File aifid, durfid, enerfid, *gmfid;
    MPI_Status   filestatus;
    int csize;
    int px = 1, py = 1, pz = 1;
    int rank, nprocs;
    int i, j, k, l, n;
    MPI_Offset off;
    long int offset;
    long int s0;
    int xi, yi, zi;
    float dt2;

    /*parameters for xapiir*/
    int dofilt = 0;
    int comp_sa = 0;
    int iord = 3, npas = 1;
    float trbndw = 0., a = 0.; /* chebyshev parameters */
    char *aproto = "BU";
    char ftype[3];
    float lp = 0.15, hp = 5.00;
    float low = 0.05, high = 0.75;

    /* parameters for getopt */
    int c;
    extern char *optarg;

    char phfile[200], pzfile[200], pgvfile[200], pgafile[200], aifile[200], durfile[200], enerfile[200], gmfile[200];
    int err;

    float freq[] = {0.1, 0.333333, 0.5, 0.6666666, 1.0, 1.5, 2.0, 3.0, 3.5};
    int nfreq = (int)sizeof(freq) / sizeof(freq[0]);

    float **GM; //gmrotD50 and duration
    int f;
    float xid = 0.05; // damping constant for gmrotD50
    int nang = 91;    // number of angles for gmrotD50


    //MPI_Init(&argc, &argv);
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    while ((c = getopt(argc, argv, "fgl:h:s:e:x:y:z:")) != -1)
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
                sscanf(optarg, "%f", &lp);
                break;
            case 'h':
                sscanf(optarg, "%f", &hp);
                break;
            case 's':
                sscanf(optarg, "%f", &low);
                break;
            case 'e':
                sscanf(optarg, "%f", &high);
                break;
            case 'x':
                sscanf(optarg, "%d", &px);
                break;
            case 'y':
                sscanf(optarg, "%d", &py);
                break;
            case 'z':
                sscanf(optarg, "%d", &pz);
                break;
            case '?':
                fprintf(stdout, "usage: %s [-f] [-l lp] [-h hp] [-s low] [-e high] [-x px] [-y py] [-z pz]\n", argv[0]);
                exit(0);
            default:
                abort();
        }
    }

    if (hp > 0.)
      strncpy(ftype, "BP\0", 3);
    else
      strncpy(ftype, "LP\0", 3);
    /* determine dimensions from these parameters */
    nx = (int)floorf((par.nedx - par.nbgx) / par.nskpx + 1.);
    ny = (int)floorf((par.nedy - par.nbgy) / par.nskpy + 1.);
    nz = (int)floorf((par.nedz - par.nbgz) / par.nskpz + 1.);
    nt = (int)floorf((par.tmax / par.dt / par.ntiskp));
    dt2 = par.dt * par.ntiskp;

    int size[3] = {nx, ny, par.itype == 0 ? nt : par.writestep};
    int part[2] = {px, py};
    mpi_cart(&m, size, part);

    /* output parameters for debug */
    if (m.rank == 0)
    {
        fprintf(stdout, "  >> nx=%d, ny=%d, nz=%d, nt=%d,"
                    " ivelocity=%d, ntiskp=%d, writestep=%d\n",
                    nx, ny, nz, nt, par.itype, par.ntiskp, par.writestep);
        fprintf(stdout, "MPI Comm m: nxt=%d, nyt=%d, ntt=%d, dim=(%d, %d), period=(%d, %d), reorder=%d\n", m.nxt, m.nyt, m.ntt, m.dim[0], m.dim[1], m.period[0], m.period[1], m.reorder); 
        fflush(stdout);
        /* number of points to process per read operation */
        fprintf(stdout, "dofilt=%d, lp=%f, hp=%f, filter type=%s, low=%f, high=%f, px=%d, py=%d, pz=%d\n", dofilt, lp, hp, ftype, low, high, px, py, pz);
        if ((nx / px != m.nxt) || (ny / py != m.nyt)) {
            MPI_Finalize();
            exit(1);
        }
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


    csize = m.nxt * m.nyt;
    if (m.rank == 0) {
        fprintf(stdout, "rank=%d: csize=%ld\n", m.rank, csize);
        fflush(stdout);
    }

    /* arrays for peak velocities, displacements */
    PH = (float *)calloc(csize, sizeof(float));
    PZ = (float *)calloc(csize, sizeof(float));
    PGV = (float *)calloc(csize, sizeof(float));
    PGA = (float *)calloc(csize, sizeof(float));
    AI = (float *)calloc(csize, sizeof(float));
    ENER = (float *)calloc(csize, sizeof(float));
    DUR = (float *)calloc(csize, sizeof(float));

    GM = (float **)calloc(nfreq, sizeof(float *));
    for (f = 0; f<nfreq && comp_sa> 0; f++)
      GM[f] = (float *)calloc(csize, sizeof(float));

    if (m.rank == 0) {
        fprintf(stdout, "%d: Allocating buffers...\n", m.rank);
        fflush(stdout);
    }

    /* 2D arrays for velocities */
    x = (float **)calloc(csize, sizeof(float *));
    y = (float **)calloc(csize, sizeof(float *));
    z = (float **)calloc(csize, sizeof(float *));

    if (m.rank == 0) {
      fprintf(stdout, "%d: Allocating 2D arrays ...\n", m.rank);
      fflush(stdout);
    }
    for (l = 0; l < csize; l++)
    {
        x[l] = (float *)calloc(nt, sizeof(float));
        y[l] = (float *)calloc(nt, sizeof(float));
        z[l] = (float *)calloc(nt, sizeof(float));
    }
    if (z[csize - 1] == NULL)
      perror("error using calloc()");

    if (dofilt == 0)
    {
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
        sprintf(phfile, "pgv_H_%05.2f_%05.2fHz.bin", lp, hp);
        sprintf(pzfile, "pgv_Z_%05.2f_%05.2fHz.bin", lp, hp);
        sprintf(pgvfile, "pgv_%05.2f_%05.2fHz.bin", lp, hp);
        sprintf(pgafile, "pga_%05.2f_%05.2fHz.bin", lp, hp);
        sprintf(aifile, "arias_%05.2f_%05.2fHz.bin", lp, hp);
        sprintf(durfile, "dur_%05.2f_%05.2fHz.bin", lp, hp);
        sprintf(enerfile, "ener_%05.2f_%05.2fHz.bin", lp, hp);
    }
    MPI_Barrier(m.MCW);

    if (m.rank == 0) {
      fprintf(stdout, "%d: opening output files ...\n", m.rank);
      fflush(stdout);
    }

    MPI_Datatype write_dtype = data_type(&m, nx, ny, 1);
    MPI_Datatype read_dtype = data_type(&m, nx, ny, size[2]);
    MPICHK(MPI_File_open(MPI_COMM_SELF, phfile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                    MPI_INFO_NULL, &hfid));
    MPICHK(MPI_File_set_view(hfid, 0, MPI_FLOAT, write_dtype, "native",
                    MPI_INFO_NULL));

    MPICHK(MPI_File_open(MPI_COMM_SELF, pzfile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                    MPI_INFO_NULL, &zfid));
    MPICHK(MPI_File_set_view(zfid, 0, MPI_FLOAT, write_dtype, "native",
                    MPI_INFO_NULL));

    MPICHK(MPI_File_open(MPI_COMM_SELF, pgvfile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                    MPI_INFO_NULL, &pgvfid));
    MPICHK(MPI_File_set_view(pgvfid, 0, MPI_FLOAT, write_dtype, "native",
                    MPI_INFO_NULL));

    MPICHK(MPI_File_open(MPI_COMM_SELF, pgafile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                    MPI_INFO_NULL, &pgafid));
    MPICHK(MPI_File_set_view(pgafid, 0, MPI_FLOAT, write_dtype, "native",
                    MPI_INFO_NULL));

    MPICHK(MPI_File_open(MPI_COMM_SELF, aifile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                    MPI_INFO_NULL, &aifid));
    MPICHK(MPI_File_set_view(aifid, 0, MPI_FLOAT, write_dtype, "native",
                    MPI_INFO_NULL));

    MPICHK(MPI_File_open(MPI_COMM_SELF, durfile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                    MPI_INFO_NULL, &durfid));
    MPICHK(MPI_File_set_view(durfid, 0, MPI_FLOAT, write_dtype, "native",
                    MPI_INFO_NULL));

    MPICHK(MPI_File_open(MPI_COMM_SELF, enerfile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                    MPI_INFO_NULL, &enerfid));
    MPICHK(MPI_File_set_view(enerfid, 0, MPI_FLOAT, write_dtype, "native",
                    MPI_INFO_NULL));

    gmfid = (MPI_File *)calloc(nfreq, sizeof(MPI_File));

    for (f = 0; f<nfreq && comp_sa> 0; f++)
    {
        sprintf(gmfile, "gmrotD50_%05.2fHz.bin", freq[f]);
        MPICHK(MPI_File_open(MPI_COMM_SELF, gmfile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                        MPI_INFO_NULL, gmfid + f));
        MPICHK(MPI_File_set_view(*(gmfid + f), 0, MPI_FLOAT, write_dtype, "native",
                        MPI_INFO_NULL));
    }
    MPI_Barrier(MPI_COMM_WORLD);

    accx = (float *)calloc(nt, sizeof(float));
    accy = (float *)calloc(nt, sizeof(float));
    accz = (float *)calloc(nt, sizeof(float));

    int step, nfiles;
    if (par.itype == 0) {
        step = nt;
        nfiles = 1;
    } else {
        step = par.writestep * par.ntiskp;
        nfiles = nt / par.writestep;
    }
    /* buffers for unsorted velocities */
    int buffer_size = csize * par.writestep;
    bufx = (float *)calloc(buffer_size, sizeof(float));
    bufy = (float *)calloc(buffer_size, sizeof(float));
    bufz = (float *)calloc(buffer_size, sizeof(float));
    for (k = 0; k < nfiles; ++k)
    {
        if (m.rank == 0) {
            fprintf(stdout, "processing part %02d of %02d\n", k + 1, nfiles);
            read_timeseries(xfile, (k + 1) * step, &m, read_dtype, bufx, buffer_size); 
            read_timeseries(yfile, (k + 1) * step, &m, read_dtype, bufy, buffer_size); 
            read_timeseries(zfile, (k + 1) * step, &m, read_dtype, bufz, buffer_size); 
            fflush(stdout);
        }
        // buffer[writestep, nzt, nyt, nxt]
        // x[nz, ny, nx, nt]
        for (i = 0; i < par.writestep; ++i) {
            for (j = 0; j < csize; ++j) {
                x[j][i + k * par.writestep] = bufx[i * csize + j];
                y[j][i + k * par.writestep] = bufy[i * csize + j];
                z[j][i + k * par.writestep] = bufz[i * csize + j];
            }
        }
        MPI_Barrier(m.MCW);
    }
    if (m.rank == 0) {
        fprintf(stdout, "Done reading time series\n");
        fflush(stdout);
    }

    for (l = 0; l < csize; l++)
    {
        if (m.rank == 0 && l % 2000 == 0) {
            fprintf(stdout, "Processing time series at idx =%d\n", l);
            fflush(stdout);
        }
        PH[l] = PZ[l] = PGV[l] = PGA[l] = AI[l] = DUR[l] = ENER[l] = 0.f;

        vel2acc(x[l], nt, dt2, accx);
        vel2acc(y[l], nt, dt2, accy);
        vel2acc(z[l], nt, dt2, accz);

#ifdef DEBUG
        if (l == 16)
        {
            FILE *debugfid;
            debugfid = fopen("output.dbg", "w");
            for (n = 0; n < nt; n++)
            {
                fprintf(debugfid, "%.3e %.3e %.3e .3e\n", n * dt2, accx[n], accy[n], accz[n]);
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
        if (dofilt == 1)
        {
            if (m.rank == 0) {
                fprintf(stdout, "Filtering time series\n");
                fflush(stdout);
            }
            xapiir_(x[l], &nt, aproto, &trbndw, &a, &iord, ftype, &lp, &hp, &dt2, &npas);
            xapiir_(y[l], &nt, aproto, &trbndw, &a, &iord, ftype, &lp, &hp, &dt2, &npas);
            xapiir_(z[l], &nt, aproto, &trbndw, &a, &iord, ftype, &lp, &hp, &dt2, &npas);
            xapiir_(accx, &nt, aproto, &trbndw, &a, &iord, ftype, &lp, &hp, &dt2, &npas);
            xapiir_(accy, &nt, aproto, &trbndw, &a, &iord, ftype, &lp, &hp, &dt2, &npas);
            xapiir_(accz, &nt, aproto, &trbndw, &a, &iord, ftype, &lp, &hp, &dt2, &npas);
        }
        if (m.rank == 0) {
            fprintf(stdout, "Done Filtering time series\n");
            fflush(stdout);
        }

        comp_cum(x[l], y[l], z[l], nt, dt2, low, high, DUR + l, ENER + l);

        for (n = 0; n < nt; n++)
        {
            vh = sqrtf(powf(x[l][n], 2.) + powf(y[l][n], 2.));
            if (PH[l] < vh)
              PH[l] = vh;
            vz = fabsf(z[l][n]);
            if (PZ[l] < vz)
              PZ[l] = vz;
            vel = sqrtf(powf(x[l][n], 2.) + powf(y[l][n], 2.) + powf(z[l][n], 2.));
            if (PGV[l] < vel)
              PGV[l] = vel;
            acc = sqrtf(powf(accx[n], 2.) + powf(accy[n], 2.) + powf(accz[n], 2.));
            if (PGA[l] < acc)
              PGA[l] = acc;
            AI[l] += powf(acc, 2.);
        }
        AI[l] = AI[l] * PI * dt2 / 2. / G / 3.;
    }

    if (m.rank == 0)
      fprintf(stdout, "writing data to disk ...\n");
    fflush(stdout);
    MPI_Barrier(m.MCW);
    MPICHK(MPI_File_write_all(hfid, PH, csize, MPI_FLOAT, &filestatus));
    MPICHK(MPI_File_write_all(zfid, PZ, csize, MPI_FLOAT, &filestatus));
    MPICHK(MPI_File_write_all(pgvfid, PGV, csize, MPI_FLOAT, &filestatus));
    MPICHK(MPI_File_write_all(pgafid, PGA, csize, MPI_FLOAT, &filestatus));
    MPICHK(MPI_File_write_all(aifid, AI, csize, MPI_FLOAT, &filestatus));
    MPICHK(MPI_File_write_all(durfid, DUR, csize, MPI_FLOAT, &filestatus));
    MPICHK(MPI_File_write_all(enerfid, ENER, csize, MPI_FLOAT, &filestatus));

    for (f = 0; f<nfreq && comp_sa> 0; f++)
      MPICHK(MPI_File_write_all(gmfid[f], GM[f], csize, MPI_FLOAT, &filestatus));
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

    for (f = 0; f<nfreq && comp_sa> 0; f++)
    {
        free(GM[f]);
        MPI_File_close(gmfid + f);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(0);
}

void mpi_init(struct Mpi *m, int argc, char **argv)
{

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&m->rank);
    MPI_Comm_size(MPI_COMM_WORLD,&m->size);
    MPI_Comm_dup(MPI_COMM_WORLD, &m->MCW);
    m->coord[0] = 0;
    m->coord[1] = 0;

}

void mpi_cart(struct Mpi *m, const int *size, const int *part)
{
    m->nxt       = size[0] / part[0];
    m->nyt       = size[1] / part[1];
    m->ntt       = size[2];
    m->dim[0]    = part[0];
    m->dim[1]    = part[1];
    m->period[0] = 0;
    m->period[1] = 0;
    m->reorder   = 0;
    MPICHK(MPI_Cart_create(m->MCW, 2, m->dim, m->period, m->reorder,
                    &m->MC1));
    MPICHK(MPI_Cart_coords(m->MC1, m->rank, 2, m->coord));
}

MPI_Datatype data_type(const struct Mpi *mpi, int nx, int ny, int nt)
{
    int old[3], new[3], offset[3];
    old[0] = nt;
    old[1] = ny;
    old[2] = nx;
    new[0] = nt;
    new[1] = mpi->nyt;
    new[2] = mpi->nxt;
    offset[0] = 0;
    offset[1] = mpi->nyt * mpi->coord[1];
    offset[2] = mpi->nxt * mpi->coord[0];

    MPI_Datatype readtype;
    MPICHK(MPI_Type_create_subarray(3, old, new, offset, MPI_ORDER_C,
                    MPI_FLOAT, &readtype));
    MPICHK(MPI_Type_commit(&readtype));
    return readtype;
}
