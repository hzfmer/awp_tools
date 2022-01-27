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
    MPI_Status   filestatus;

    int px, py;
    char fname[1024];
    int c;
    extern char *optarg;
    while ((c = getopt(argc, argv, "s:x:y:")) != -1)
    {
        switch (c)
        {
            case 's':
                //sscanf(optarg, "%s", fname);
                strcpy(fname, optarg);
                break;
            case 'x':
                sscanf(optarg, "%d", &px);
                break;
            case 'y':
                sscanf(optarg, "%d", &py);
                break;
            case '?':
                fprintf(stdout, "usage: %s [-f] [-l lp] [-h hp] [-s low] [-e high] [-x px] [-y py] [-z pz]\n", argv[0]);
                exit(0);
            default:
                abort();
        }
    }

    /* determine dimensions from these parameters */
    int nx=9504, ny=7020, nt=800;
    int i, nio=2;
    for (i = 0; i < nio; ++i) {
        int size[3] = {nx, ny / nio, nt};
        int part[2] = {px, py};
        mpi_cart(&m, size, part);
        MPI Offset off = 


        /* this removes the remaining ' at the end of the string, if any */

        MPI_Datatype read_dtype = data_type(&m, nx, ny, size[2]);

        int csize = m.nxt * m.nyt;
        int buffer_size = csize * nt;
        float *bufx; 
        bufx = (float *)calloc(buffer_size, sizeof(float));
        if (m.rank == 0) {
            fprintf(stdout, "rank=%ld: fname=%s, csize=%d, buffer_size=%ld\n", m.rank, fname, csize, buffer_size);
            fflush(stdout);
        }

        MPI_File fp;
        MPICHK(MPI_File_open(m.MCW, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp));
        MPICHK(MPI_File_set_view(fp, 0, MPI_FLOAT, read_dtype, "native", MPI_INFO_NULL));
        MPICHK(MPI_File_read_all(fp, bufx, buffer_size, MPI_FLOAT, &filestatus));

        if (m.rank == 0) {
            fprintf(stdout, "Done reading\n");
            fflush(stdout);
        }

        MPI_File fh;
        MPICHK(MPI_File_open(m.MCW, "test.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh));
        MPICHK(MPI_File_set_view(fh, 0, MPI_FLOAT, read_dtype, "native", MPI_INFO_NULL));
        MPICHK(MPI_File_write_all(fh, bufx, buffer_size, MPI_FLOAT, &filestatus));
        MPI_Barrier(m.MCW);
    }
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
