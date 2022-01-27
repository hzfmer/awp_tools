/* Computes peak ground velocity from AWP-ODC output.
   Improved version using MPI-IO.
   Daniel Roten, March 2015 <droten@sdsc.edu> 
   Zhifeng Hu, Feb 2021 <zhifeng94.hu@gmail.com> */

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

struct Mpi
{
    // MPI vars
    int rank;
    int size;
    int err;
    int dim[3];
    int period[3];
    int reorder;
    int nxt, nyt, nzt, ntt;
    int coord[3];
    MPI_Comm MCW, MC1;
};

int index_above(float *in, int nt, float thres)
{
    int i = 0;
    for (i = 0; i < nt; ++i) {
        if (in[i] > thres) break;
    }
    return i;
}

int index_above1(float *in, int nt, float thres)
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

void comp_cum(float *accX, float *accY, float *accZ, int nt, float dt, float start, float end, float *dur, float *ener, int p)
{
    float *cumeng;
    float start_bound, end_bound;
    int l;

    cumeng = (float *)calloc(nt, sizeof(float));
    cumeng[0] = powf(accX[0], 2.) + powf(accY[0], 2.) + powf(accZ[0], 2.);

    for (l = 1; l < nt; l++)
    {
//        if (p)
//            fprintf(stderr, "cumeng: %f %f %f %f ", cumeng[l], accX[l], accY[l], accZ[l]);
        cumeng[l] = cumeng[l - 1] + powf(accX[l], 2.) + powf(accY[l], 2.) + powf(accZ[l], 2.);
//        if (p) fprintf(stderr, "%f\n", cumeng[l]);
    }
    start_bound = cumeng[nt - 1] * start;
    end_bound = cumeng[nt - 1] * end;

    int start_idx = index_above(cumeng, nt, start_bound);
    int end_idx = index_above(cumeng, nt, end_bound);
    *dur = (end_idx - start_idx) * dt;
    *ener = cumeng[nt - 1] / 3 * dt;
//    if (p)
//    fprintf(stderr, "res: %d %d %f %f %f %f\n", start_idx, end_idx, start_bound, end_bound, *dur, *ener);
    free(cumeng);
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


void read_timeseries(char *fbase, int step, struct Mpi *m, MPI_Datatype dtype, float *output, int buffer_size)
{
    char fname[200];
    sprintf(fname, "%s%07d", fbase, step);
    MPI_File fp;
    MPI_Status filestatus;
    MPICHK(MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp));
    MPICHK(MPI_File_set_view(fp, 0, MPI_FLOAT, dtype, "native", MPI_INFO_NULL));
    MPICHK(MPI_File_read_all(fp, output, buffer_size, MPI_FLOAT, &filestatus));
}

void mpi_init(struct Mpi *m, int argc, char **argv);
void mpi_cart(struct Mpi *m, const int *size, const int *part, int ysplit);
MPI_Datatype data_type(const struct Mpi *mpi, int nx, int ny, int nz, int nt, int yoff);

int main(int argc, char *argv[])
{
    struct Mpi m;
    mpi_init(&m, argc, argv);
    FD3D_param par;
    read_settings(&par, "IN3D.out");
    char *xfile, *yfile, *zfile;
    int nx = 1, ny = 1, nz = 1, nt = 1;
    float **x, **y, **z;
    float *velx, *vely, *velz;
    float *accx, *accy, *accz, acc;
    float *bufx, *bufy, *bufz;
    float vh, vz, vel;
    MPI_File *hfid, *zfid, *pgvfid, *pgafid;
    MPI_File *aifid, *durfid, *enerfid, *gmfid;
    MPI_Status filestatus;
    int csize, buffer_size;
    int px = 1, py = 1, pz = 1, ysplit = 1;
    int rank, nprocs;
    int i, j, k, l, t, n;
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
    char ftype[MAX_NUM_FILTER][3];
    int nfilter = 1;
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

    int err;

    float freq[] = {0.1, 0.333333, 0.5, 0.6666666, 1.0, 1.5, 2.0, 3.0, 3.5};
    int nfreq = (int)sizeof(freq) / sizeof(freq[0]);

    int f;
    float xid = 0.05; // damping constant for gmrotD50
    int nang = 91;    // number of angles for gmrotD50

    while ((c = getopt(argc, argv, "fgl:h:s:e:p:x:y:z:")) != -1)
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
        case 's':
            sscanf(optarg, "%f", &start);
            break;
        case 'e':
            sscanf(optarg, "%f", &end);
            break;
        case 'x':
            sscanf(optarg, "%d", &px);
            break;
        case 'y':
            sscanf(optarg, "%d", &py);
            break;
        case 'p':
            sscanf(optarg, "%d", &ysplit);
            break;
        case 'z':
            sscanf(optarg, "%d", &pz);
            break;
        case '?':
            fprintf(stdout, "usage: %s [-f] [-l lp0, lp1, lp2...] [-h hp0, hp1, hp2...]"
                            "[-s start] [-e end] [-x px] [-y py] [-z pz] [-p ysplit]\n",
                    argv[0]);
            exit(0);
        default:
            abort();
        }
    }

    char phfile[200], pzfile[200], pgvfile[200], pgafile[200], aifile[200], durfile[200], enerfile[200], gmfile[200];
    if (dofilt == 0) nfilter = 1;
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

    int size[4] = {nx, ny, nz, par.itype == 0 ? nt : par.writestep};
    int part[3] = {px, py, pz};
    mpi_cart(&m, size, part, ysplit);

    /* output parameters for debug */
    if (m.rank == 0)
    {
        fprintf(stdout, "  >> nx=%d, ny=%d, nz=%d, nt=%d,"
                        " ivelocity=%d, ntiskp=%d, writestep=%d\n",
                nx, ny, nz, nt, par.itype, par.ntiskp, par.writestep);
        fprintf(stdout, "MPI Comm m: nxt=%d, nyt=%d, nzt=%d, ntt=%d, "
                        "dim=(%d, %d, %d), period=(%d, %d, %d), reorder=%d\n",
                m.nxt, m.nyt, m.nzt, m.ntt, m.dim[0], m.dim[1], m.dim[2],
                m.period[0], m.period[1], m.period[2], m.reorder);
        /* number of points to process per read operation */
        fprintf(stdout, "start=%f, end=%f, px=%d, py=%d, pz=%d, ysplit=%d, dofilt=%d\n",
                start, end, px, py, pz, ysplit, dofilt);
        for (i = 0; i < nfilter && dofilt; ++i) {
            fprintf(stdout, "Filter %d (%s): low=%f, high=%f\n", i, ftype[i], lps[i], hps[i]);
        }
        fflush(stdout);
    }
    if ((nx / px != m.nxt) || (ny / py /ysplit != m.nyt) || (nz / pz != m.nzt))
    {
        fprintf(stdout, "Mismatch in computation nodes allocation!\n");
        fflush(stdout);
        MPI_Finalize();
        abort();
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

    csize = m.nxt * m.nyt * m.nzt;
    buffer_size = csize * par.writestep;
    if (m.rank == 0)
    {
        fprintf(stdout, "rank=%d: csize=%d, buffer_size=%d\n", m.rank, csize, buffer_size);
        fflush(stdout);
    }


    int step, nfiles;
    if (par.itype == 0)
    {
        step = nt;
        nfiles = 1;
    }
    else
    {
        step = par.writestep * par.ntiskp;
        nfiles = nt / par.writestep;
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

    for (l = 0; l < csize; l++)
    {
        x[l] = (float *)calloc(nt, sizeof(float));
        y[l] = (float *)calloc(nt, sizeof(float));
        z[l] = (float *)calloc(nt, sizeof(float));
    }
    if (z[csize - 1] == NULL)
        perror("error using calloc()");

    /* buffers for unsorted velocities */
    bufx = (float *)calloc(buffer_size, sizeof(float));
    bufy = (float *)calloc(buffer_size, sizeof(float));
    bufz = (float *)calloc(buffer_size, sizeof(float));
    accx = (float *)calloc(nt, sizeof(float));
    accy = (float *)calloc(nt, sizeof(float));
    accz = (float *)calloc(nt, sizeof(float));
    velx = (float *)calloc(nt, sizeof(float));
    vely = (float *)calloc(nt, sizeof(float));
    velz = (float *)calloc(nt, sizeof(float));

    gmfid = (MPI_File *)calloc(nfreq, sizeof(MPI_File));
    hfid = (MPI_File*)calloc(nfilter, sizeof(MPI_File));
    zfid = (MPI_File*)calloc(nfilter, sizeof(MPI_File));
    pgvfid = (MPI_File*)calloc(nfilter, sizeof(MPI_File));
    pgafid = (MPI_File*)calloc(nfilter, sizeof(MPI_File));
    aifid = (MPI_File*)calloc(nfilter, sizeof(MPI_File));
    enerfid = (MPI_File*)calloc(nfilter, sizeof(MPI_File));
    durfid = (MPI_File*)calloc(nfilter, sizeof(MPI_File));

    for (i = 0; i < nfilter; ++i) {
        if (dofilt == 0) {
           sprintf(phfile, "pgv_H.bin");
           sprintf(pzfile, "pgv_Z.bin");
           sprintf(pgvfile, "pgv.bin");
           sprintf(pgafile, "pga.bin");
           sprintf(aifile, "ai.bin");
           sprintf(enerfile, "ener.bin");
           sprintf(durfile, "dur.bin");
        }
        else
        {
            sprintf(phfile, "pgv_H_%05.2f_%05.2fHz.bin", lps[i], hps[i]);
            sprintf(pzfile, "pgv_Z_%05.2f_%05.2fHz.bin", lps[i], hps[i]);
            sprintf(pgvfile, "pgv_%05.2f_%05.2fHz.bin", lps[i], hps[i]);
            sprintf(pgafile, "pga_%05.2f_%05.2fHz.bin", lps[i], hps[i]);
            sprintf(aifile, "ai_%05.2f_%05.2fHz.bin", lps[i], hps[i]);
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
    for (f = 0; f < nfreq && comp_sa > 0; f++)
    {
        sprintf(gmfile, "gmrotD50_%05.2fHz.bin", freq[f]);
        MPICHK(MPI_File_open(MPI_COMM_SELF, gmfile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                        MPI_INFO_NULL, gmfid + f));
    }

    for (n = 1; n <= ysplit; ++n)
    {
        MPI_Datatype read_dtype = data_type(&m, nx, ny, nz, nt, (n - 1) * ny / ysplit);
        MPI_Datatype write_dtype = data_type(&m, nx, ny, nz, 1, (n - 1) * ny / ysplit);
        for (i = 0; i < nfilter; ++i) {
            MPICHK(MPI_File_set_view(hfid[i], 0, MPI_FLOAT, write_dtype, "native",
                                    MPI_INFO_NULL));

            MPICHK(MPI_File_set_view(zfid[i], 0, MPI_FLOAT, write_dtype, "native",
                                    MPI_INFO_NULL));

            MPICHK(MPI_File_set_view(pgvfid[i], 0, MPI_FLOAT, write_dtype, "native",
                                    MPI_INFO_NULL));

            MPICHK(MPI_File_set_view(pgafid[i], 0, MPI_FLOAT, write_dtype, "native",
                                    MPI_INFO_NULL));

            MPICHK(MPI_File_set_view(aifid[i], 0, MPI_FLOAT, write_dtype, "native",
                                    MPI_INFO_NULL));

            MPICHK(MPI_File_set_view(durfid[i], 0, MPI_FLOAT, write_dtype, "native",
                                    MPI_INFO_NULL));

            MPICHK(MPI_File_set_view(enerfid[i], 0, MPI_FLOAT, write_dtype, "native",
                                    MPI_INFO_NULL));
        }
    
        for (f = 0; f < nfreq && comp_sa > 0; f++)
        {
            MPICHK(MPI_File_set_view(*(gmfid + f), 0, MPI_FLOAT, write_dtype, "native",
                                     MPI_INFO_NULL));
        }

        // Read and processing each ysplit chunk
        for (k = 0; k < nfiles; ++k)
        {
            if (m.rank == 0)
            {
                fprintf(stdout, "processing part %d of %d, chunk %d of %d\n", k + 1, nfiles, n, ysplit);
                fflush(stdout);
            }
            read_timeseries(xfile, (k + 1) * step, &m, read_dtype, bufx, buffer_size);
            if (m.rank == 0)
            {
                fprintf(stdout, "processing component x\n");
                fflush(stdout);
            }
            read_timeseries(yfile, (k + 1) * step, &m, read_dtype, bufy, buffer_size);
            if (m.rank == 0)
            {
                fprintf(stdout, "processing component y\n");
                fflush(stdout);
            }
            read_timeseries(zfile, (k + 1) * step, &m, read_dtype, bufz, buffer_size);
            if (m.rank == 0)
            {
                fprintf(stdout, "processing component z\n");
                fflush(stdout);
            }

            // buffer[writestep, nzt, nyt, nxt]
            // x[nzt, nyt, nxt, nt]
            for (i = 0; i < par.writestep; ++i)
            {
                for (j = 0; j < csize; ++j)
                {
                    x[j][i + k * par.writestep] = bufx[i * csize + j];
                    y[j][i + k * par.writestep] = bufy[i * csize + j];
                    z[j][i + k * par.writestep] = bufz[i * csize + j];
                }
            }
        }
        
        for (l = 0; l < csize; l++)
        {
            vel2acc(x[l], nt, dt2, accx);
            vel2acc(y[l], nt, dt2, accy);
            vel2acc(z[l], nt, dt2, accz);

            for (f = 0; f < nfreq && comp_sa > 0; f++)
            {
                GM[f][l] = gmrotD50(accx, accy, nt, dt2, xid, freq[f], nang);
#ifdef DEBUG
                if (l == 16)
                    fprintf(stderr, ">>> gmrotD50 @ f=%5.2f Hz: %e\n", freq[f], GM[f][l]);
#endif
            }

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

//                if (m.rank == 0 && l == 1000 && i == 0 ) {
//                    fprintf(stderr, ">>> Before cumsum: vel @ f=%5.2f-%5.2f Hz: nt=%d\n", lps[i], hps[i], nt);
//                    for (j = 0; j < nt; ++j) {
//                        fprintf(stderr, "%f %f %f\n", velx[j], vely[j], velz[j]);
//                    }
//                    fprintf(stderr, "dur=%f,%zu, %zu, ener=%f,%zu, %zu\n", DUR[i][l], &DUR[i][l], DUR[i] + l, ENER[i][l], &ENER[i][l], ENER[i] + l);
//                    fflush(stderr);
//                }
                comp_cum(velx, vely, velz, nt, dt2, start, end, DUR[i] + l, ENER[i] + l, m.rank == 0 && l == 1000 && i == 0);
//                if (m.rank == 0 && l == 1000 && i == 0 ) {
//                    fprintf(stderr, ">>> vel @ f=%5.2f-%5.2f Hz: nt=%d\n", lps[i], hps[i], nt);
//                    for (j = 0; j < nt; ++j) {
//                        fprintf(stderr, "%f %f %f\n", velx[j], vely[j], velz[j]);
//                    }
//                    fprintf(stderr, "dur=%f,%zu, %zu, ener=%f,%zu, %zu\n", DUR[i][l], &DUR[i][l], DUR[i] + l, ENER[i][l], &ENER[i][l], ENER[i] + l);
//                    fflush(stderr);
//                }

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

        if (m.rank == 0) {
            fprintf(stdout, "writing data to disk ...\n");
            fflush(stdout);
        }
        for (i = 0; i < nfilter; ++i) {
            MPICHK(MPI_File_write_all(hfid[i], PH[i], csize, MPI_FLOAT, &filestatus));
            MPICHK(MPI_File_write_all(zfid[i], PZ[i], csize, MPI_FLOAT, &filestatus));
            MPICHK(MPI_File_write_all(pgvfid[i], PGV[i], csize, MPI_FLOAT, &filestatus));
            MPICHK(MPI_File_write_all(pgafid[i], PGA[i], csize, MPI_FLOAT, &filestatus));
            MPICHK(MPI_File_write_all(aifid[i], AI[i], csize, MPI_FLOAT, &filestatus));
            MPICHK(MPI_File_write_all(durfid[i], DUR[i], csize, MPI_FLOAT, &filestatus));
            MPICHK(MPI_File_write_all(enerfid[i], ENER[i], csize, MPI_FLOAT, &filestatus));
        }

        for (f = 0; f < nfreq && comp_sa > 0; f++)
        {
            MPICHK(MPI_File_write_all(gmfid[f], GM[f], csize, MPI_FLOAT, &filestatus));
        }

        MPI_Barrier(m.MCW);
    }
    for (i = 0; i < csize; ++i)
    {
        free(x[i]); free(y[i]); free(z[i]);
    }
    free(x); free(y); free(z);
    free(accx); free(accy); free(accz);
    free(bufx); free(bufy); free(bufz);

    for (i = 0; i < nfilter; ++i) {
        MPI_File_close(hfid + i);
        MPI_File_close(zfid + i);
        MPI_File_close(pgafid + i);
        MPI_File_close(pgvfid + i);
        MPI_File_close(aifid + i);
        MPI_File_close(durfid + i);
        MPI_File_close(enerfid + i);
    }
    for (f = 0; f < nfreq && comp_sa > 0; f++)
    {
        MPI_File_close(gmfid + f);
    }
    free(xfile); free(yfile); free(zfile);
    free(GM);
    free(PH);
    free(PZ);
    free(PGA);
    free(PGV);
    free(AI);
    free(DUR);
    free(ENER);

    MPI_Barrier(m.MCW);
    MPI_Finalize();
    exit(0);
}

void mpi_init(struct Mpi *m, int argc, char **argv)
{

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &m->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &m->size);
    MPI_Comm_dup(MPI_COMM_WORLD, &m->MCW);
    m->coord[0] = 0;
    m->coord[1] = 0;
    m->coord[2] = 0;
}

void mpi_cart(struct Mpi *m, const int *size, const int *part, int ysplit)
{
    m->nxt = size[0] / part[0];
    m->nyt = size[1] / part[1] / ysplit;
    m->nzt = size[2] / part[2];
    m->ntt = size[3];
    m->dim[0] = part[0];
    m->dim[1] = part[1];
    m->dim[2] = part[2];
    m->period[0] = 0;
    m->period[1] = 0;
    m->period[2] = 0;
    m->reorder = 0;
    MPICHK(MPI_Cart_create(m->MCW, 3, m->dim, m->period, m->reorder,
                           &m->MC1));
    MPICHK(MPI_Cart_coords(m->MC1, m->rank, 3, m->coord));
}

MPI_Datatype data_type(const struct Mpi *mpi, int nx, int ny, int nz, int nt, int yoff)
{
    int old[4] = {nt, nz, ny, nx};
    int new[4] = {nt, mpi->nzt, mpi->nyt, mpi->nxt};
    int offset[4];
    offset[0] = 0;
    offset[1] = new[1] * mpi->coord[2];
    offset[2] = new[2] * mpi->coord[1] + yoff;
    offset[3] = new[3] * mpi->coord[0];

    MPI_Datatype dtype;
    MPICHK(MPI_Type_create_subarray(4, old, new, offset, MPI_ORDER_C,
                                    MPI_FLOAT, &dtype));
    MPICHK(MPI_Type_commit(&dtype));
    return dtype;
}
