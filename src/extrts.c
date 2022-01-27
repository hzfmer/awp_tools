#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>
#include "awp_data.h"
#include "fd3dparam.h"

void writeasc(char *fname, int nt, float *buf)
{
    FILE *fid;
    int k;

    fid = fopen(fname, "w");
    for (k = 0; k < nt; k++)
      fprintf(fid, "%e\n", buf[k]);
    fclose(fid);
}

int main(int argc, char *argv[])
{
#ifdef SAF50
    int nx = 6000, ny = 2744, nz = 1;
    int nt = 54000, wstep = 100, ntiskp = 20;
    int multiplexed = 1;
#endif
#ifdef SAF50SM
    int nx = 5400, ny = 400, nz = 1;
    int nt = 2400, wstep = 300, ntiskp = 20;
    int multiplexed = 0;
#endif
#ifdef SAF100
    int nx = 3000, ny = 1372, nz = 1;
    int nt = 24000, wstep = 50, ntiskp = 10;
    int multiplexed = 0;
#endif
#ifdef USE_IN3D
    int nx, ny, nz, nt, wstep, ntiskp, multiplexed;
    FD3D_param par;
    read_settings(&par, "IN3D.out");
    nx = (int)floorf((par.nedx - par.nbgx) / par.nskpx + 1);
    ny = (int)floorf((par.nedy - par.nbgy) / par.nskpy + 1);
    nz = (int)floorf((par.nedz - par.nbgz) / par.nskpz + 1);
    nt = (int)floorf((par.tmax / par.dt));
    wstep = par.writestep;
    ntiskp = par.ntiskp;
    multiplexed = par.itype;
#endif

    /*char xfile[]="output_sfc/SX";
      char yfile[]="output_sfc/SY";
      char zfile[]="output_sfc/SZ";*/
    int rank, nprocs;
    FILE *fid;
    int nstat, *xi, *yi, *zi, *xt, *yt, *zt;
    int nt2, k;
    float *bufx, *bufy, *bufz;
    char fnamex[200], fnamey[200], fnamez[200];

    char *xfile, *yfile, *zfile;
    xfile = (char *)calloc(100, sizeof(char));
    yfile = (char *)calloc(100, sizeof(char));
    zfile = (char *)calloc(100, sizeof(char));
    sscanf(par.sxrgo, "%s", xfile);
    sscanf(par.syrgo, "%s", yfile);
    sscanf(par.szrgo, "%s", zfile);
    // this removes the remaining ' at the end of the stri
    xfile = strsep(&xfile, "\'");
    yfile = strsep(&yfile, "\'");
    zfile = strsep(&zfile, "\'");

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (rank == 0)
    {
        fid = fopen("stat.txt", "r");
        fscanf(fid, "%d\n", &nstat);
        fprintf(stdout, "  >> nx=%d, ny=%d, nz=%d, nt=%d,"
                    " ivelocity=%d, ntiskp=%d, writestep=%d\n",
                    nx, ny, nz, nt, par.itype, par.ntiskp, par.writestep);
    }
    MPI_Bcast(&nstat, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if ((nstat % nprocs) != 0)
    {
        if (rank == 0)
          fprintf(stderr, "Number of stations %d not divisible by number of CPUs %d\n", nstat, nprocs);
        MPI_Finalize();
        exit(0);
    }
    xi = (int *)calloc(nstat, sizeof(int));
    yi = (int *)calloc(nstat, sizeof(int));
    zi = (int *)calloc(nstat, sizeof(int));
    xt = (int *)calloc(nstat, sizeof(int));
    yt = (int *)calloc(nstat, sizeof(int));
    zt = (int *)calloc(nstat, sizeof(int));

    if (rank == 0)
    {
        for (k = 0; k < nstat; k++)
        {
            fscanf(fid, "%d %d %d\n", xi + k, yi + k, zi + k);
            xt[k] = (xi[k] - par.nbgx) / par.nskpx;
            yt[k] = (yi[k] - par.nbgy) / par.nskpy;
            zt[k] = (zi[k] - par.nbgz) / par.nskpz;
        }
        fflush(stdout);
        fclose(fid);
    }

    MPI_Bcast(xi, nstat, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(yi, nstat, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(zi, nstat, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(xt, nstat, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(yt, nstat, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(zt, nstat, MPI_INT, 0, MPI_COMM_WORLD);

    nt2 = nt / ntiskp;

    bufx = (float *)calloc(nt2, sizeof(float));
    bufy = (float *)calloc(nt2, sizeof(float));
    bufz = (float *)calloc(nt2, sizeof(float));

    for (k = rank; k < nstat; k += nprocs)
    {
        printf("Processing site %d / %d: (%d, %d, %d)\n", k, nstat, xt[k], yt[k], zt[k]);
        if (multiplexed)
        {
            read_awp_timeseries_multi(xfile, nx, ny, nz, xt[k], yt[k], zt[k], nt2,
                        wstep, ntiskp, 1, bufx);
            read_awp_timeseries_multi(yfile, nx, ny, nz, xt[k], yt[k], zt[k], nt2,
                        wstep, ntiskp, 1, bufy);
            read_awp_timeseries_multi(zfile, nx, ny, nz, xt[k], yt[k], zt[k], nt2,
                        wstep, ntiskp, 1, bufz);
        }
        else
        {
            read_awp_timeseries(xfile, nx, ny, nz, xt[k], yt[k], zt[k], nt2, 1, bufx);
            read_awp_timeseries(yfile, nx, ny, nz, xt[k], yt[k], zt[k], nt2, 1, bufy);
            read_awp_timeseries(zfile, nx, ny, nz, xt[k], yt[k], zt[k], nt2, 1, bufz);
        }

        sprintf(fnamex, "%s%04d_%04d_%04d.dat", xfile, xi[k], yi[k], zi[k]);
        sprintf(fnamey, "%s%04d_%04d_%04d.dat", yfile, xi[k], yi[k], zi[k]);
        sprintf(fnamez, "%s%04d_%04d_%04d.dat", zfile, xi[k], yi[k], zi[k]);

        writeasc(fnamex, nt2, bufx);
        writeasc(fnamey, nt2, bufy);
        writeasc(fnamez, nt2, bufz);
        fprintf(stdout, "Done extracting site %d / %d\n", k, nstat);
        fflush(stdout);
    }
    free(xi);
    free(yi);
    free(zi);
    free(xt);
    free(yt);
    free(zt);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return (0);
}
