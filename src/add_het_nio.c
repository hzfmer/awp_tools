#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <getopt.h>
#include <string.h>
#include <limits.h>

#define DEBUG 1
#ifdef DEBUG
#define MPICHK(err) {                                                              \
  if (err != MPI_SUCCESS) {                                                        \
    char error_string[2048];                                                       \
    int length_of_error_string;                                                    \
    MPI_Error_string((err), error_string, &length_of_error_string);                \
    fprintf(stderr, "MPI error: %s:%i %s(): %s\n",                                 \
          __FILE__, __LINE__, __func__, error_string);                             \
    MPI_Abort(MPI_COMM_WORLD, err);                                                \
    fflush(stderr);                                                                \
    exit(-1);                                                                      \
    }                                                                              \
}                                                                                  
#else
#define MPICHK(err) {}
#endif

#define COUNT 0
#define ABS(x, y) ((x) > (y)) ? (x) - (y) : (y) - (x)


struct Mpi
{
    int rank;
    int size;
    int err;
    int nxt, nyt, nzt;
    int dim[2], period[2], coord[2];
    int reorder;
    MPI_Comm MCW, MC1;
};

void mpi_init(struct Mpi * m, int argc, char ** argv);
void mpi_cart(struct Mpi *m, const int *size, const int *part);
MPI_Datatype data_type(const struct Mpi *m, int num_z, int nvar);


int main(int argc, char *argv []) {
    MPI_Status status;

    float std = 0.1;
    int nx = 19440, ny = 14904, nz = 1250;
    int px = 1, py = 1;
    int nz_ssh = nz, nz_tap = nz * 3 / 4;
    int nvar = 3;
    int nio = nz;
    float vpmax = 9500., vsmax=6000.;  // Note vpmax should meet CFL criterion, vpmax * dt / dx < 0.5
    char fh_ssh[1024], fh_hom[1024], fh_het[1024];


    int c = 0;
    static struct option long_options[] = {
        {"std",     required_argument,   0,  'a'},
        {"nx",      required_argument,   0,  'x'},
        {"ny",      required_argument,   0,  'y'},
        {"nz",      required_argument,   0,  'z'},
        {"nz_ssh",  required_argument,   0,  'b'},
        {"nz_tap",  optional_argument,   0,  'c'},
        {"vpmax",   required_argument,   0,  'p'},
        {"vsmax",   required_argument,   0,  's'},
        {"fh_ssh",  required_argument,   0,  'd'},
        {"fh_hom",  required_argument,   0,  'e'},
        {"fh_het",  required_argument,   0,  'f'},
        {"px",      required_argument,   0,  'g'},
        {"py",      required_argument,   0,  'i'},
        {"nio",     required_argument,   0,  'j'},
        {NULL, 0, NULL, 0}
    };
    while (1) {
        c = getopt_long_only(argc, argv, "a:x:y:z:b:c::p:s:d:e:f:g:i:",
                long_options, NULL);
        if (c == -1) break;
        switch (c) {
            case 'a':
                std = atof(optarg); break;
            case 'x':
                nx = atoi(optarg); break;
            case 'y':
                ny = atoi(optarg); break;
            case 'z':
                nz = atoi(optarg); break;
            case 'b':
                nz_ssh = atoi(optarg); break;
            case 'c':
                if (optarg == NULL && argv[optind] != NULL
                        && argv[optind][0] != '-') {
                    nz_tap = atoi(argv[optind]);
                    optind++;
                } else if (optarg != NULL)  {
                    nz_tap = atoi(optarg);
                } else {
                    nz_tap = nz_ssh + 1;
                }
                break;
            case 'p':
                vpmax = atof(optarg); break;
            case 's':
                vsmax = atof(optarg); break;
            case 'd':
                strcpy(fh_ssh, optarg); break;
            case 'e':
                strcpy(fh_hom, optarg); break;
            case 'f':
                strcpy(fh_het, optarg); break;
            case 'g':
                px = atoi(optarg); break;
            case 'i':
                py = atoi(optarg); break;
            case 'j':
                nio = atoi(optarg); break;
            default:
                printf("Usage %s\nOptions:\n--std\t--nx\t--ny\t--nz\t--nz_ssh\t--nz_tap\t--vpmax\t--vsmax--nio\t--fh_ssh\t--fh_hom\t--fh_het\n", argv[0]);
                printf("nz_tap controls how deep (# of layers) to start tapering, default is no tapering.\n");
                abort();
        }
    }

    if (optind < argc) {
        printf("Non-option ARGV-elements: ");
        while (optind < argc) 
          printf("%s\t", argv[optind++]);
        printf("\n");
    }

    
    struct Mpi m;
    mpi_init(&m, argc, argv);

    int size[3] = {nx, ny, nz};
    int part[2] = {px, py};
    int err = 0;
    char mpiErrStr[1024];
    int mpiErrStrLen;
    mpi_cart(&m, size, part);
        
    if (m.rank == 0) {
        fprintf(stdout, "nprocs = %d\n", m.size);
        fprintf(stdout, "nx=%d, ny=%d, nz=%d, std=%f, nz_ssh=%d, nz_tap=%d, vpmax=%f, vsmax=%f, px=%d, py=%d, nio=%d\n", nx, ny, nz, std, nz_ssh, nz_tap, vpmax, vsmax, px, py, nio);
        fprintf(stdout, "fh_ssh=%s\nfh_hom=%s\nfh_het=%s\n", fh_ssh, fh_hom, fh_het);
        fflush(stdout);
        if (m.nxt * px != nx || m.nyt * py != ny) {
            fprintf(stderr, "Number of grid points nx/ny should be divisible by px/py.\n");
            MPI_Finalize();
            return -1;
        }
        if ((nx / px) * (ny / py)  >= INT_MAX / nvar / nz) {
            fprintf(stdout, "(nx / px) * (ny / py) * nvar * nz = %d\n", (nx / px) * (ny / py) * nvar * nz);
            fprintf(stderr, "Number of grids in each processor is larger than INT_MAX, may overflow!\n");
            MPI_Finalize();
            return -1;
        }
    }

    // Determine the residual of z-layers
    int resi_hom = nz % nio;  // residual of reading hom file
    int kz_hom = nz / nio;  // 
    if (resi_hom > 0) kz_hom++;

    int resi_ssh = nz_ssh % nio;
    int kz_ssh = nz_ssh / nio;
    if (resi_ssh > 0) kz_ssh++;

    int readz_hom, readz_ssh;

    MPI_File fm, fs, ft;  // fm: fh_hom; fs: fh_ssh; ft: fh_het
    MPI_Status filestatus;
    MPICHK(MPI_File_open(m.MCW, fh_ssh, MPI_MODE_RDONLY, MPI_INFO_NULL, &fs));
    MPICHK(MPI_File_open(m.MCW, fh_hom, MPI_MODE_RDONLY, MPI_INFO_NULL, &fm));
    MPICHK(MPI_File_open(m.MCW, fh_het, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                        MPI_INFO_NULL, &ft));
    // Interestingly, the size_z of hom_type in MPI_File_set_view doesn't really matter, that is
    // if the dtype::nxt, dtype::nyt are fixed, then we can read/write more z levels than given in dtype::nzt
    MPI_Datatype ssh_type = data_type(&m, nio, 1);
    MPI_Datatype hom_type = data_type(&m, nio, nvar);  // shared as het_type
    // ssh file
    MPICHK(MPI_File_set_view(fs, 0, MPI_FLOAT, ssh_type, "native", MPI_INFO_NULL));
    // homogeneous file
    MPICHK(MPI_File_set_view(fm, 0, MPI_FLOAT, hom_type, "native", MPI_INFO_NULL));
    // heterogenous file
    MPICHK(MPI_File_set_view(ft, 0, MPI_FLOAT, hom_type, "native", MPI_INFO_NULL));
    
    for (c = 0; c < kz_hom; c++) {
        
        readz_hom = (c < kz_hom - 1) ? nio : nz - c * nio;
        readz_ssh = (c < kz_ssh - 1) ? nio : (c == kz_ssh - 1) ? nz_ssh - c * nio : 0;

        int bufsize_ssh = m.nxt * m.nyt * readz_ssh;
        int bufsize_hom = m.nxt * m.nyt * nvar * readz_hom;
        float *buf_ssh = (float *) calloc(bufsize_ssh, sizeof(float));
        float *buf_hom = (float *) calloc(bufsize_hom, sizeof(float));
        
        if (m.rank==0) printf("c=%d/%d, kz_hom=%d, kz_ssh=%d, sizeof(buf_ssh)=%ld, bufsize_ssh=%d, bufsize_hom=%d, readz_ssh=%d, readz_hom=%d\n", c, kz_hom, kz_hom, kz_ssh, sizeof(buf_ssh), bufsize_ssh, bufsize_hom, readz_ssh, readz_hom);
        fflush(stdout);

        // Read SSH file
        if (readz_ssh > 0) {
            err = MPI_File_read_all(fs, buf_ssh, bufsize_ssh, MPI_FLOAT, &filestatus);
            if (err != MPI_SUCCESS) {
                MPICHK(MPI_Error_string(err, mpiErrStr, &mpiErrStrLen));
                printf("Rank=%d, ERROR=%d! MPI-IO reading SSH file, %s\n", m.rank, err, mpiErrStr);
            }
        }

        if (m.rank == 0) printf("c=%d/%d, Done reading homogeneous files.\n", c, kz_ssh);
        fflush(stdout);
        // Read homogeneous file
        err = MPI_File_read_all(fm, buf_hom, bufsize_hom, MPI_FLOAT, &filestatus);
        if (err != MPI_SUCCESS) {
            MPICHK(MPI_Error_string(err, mpiErrStr, &mpiErrStrLen));
            printf("Rank=%d, ERROR=%d! MPI-IO reading homogeneous file, %s\n", m.rank, err, mpiErrStr);
        }
        if (m.rank == 0) printf("c=%d/%d, Done reading input files.\n", c, kz_hom);
        fflush(stdout);
        
        float tapval = 1.0; //taper ratio
        int ind = 0;  // Be careful if index can be too large
        int i, j, k, z;
        float tmp_vp, tmp_vs, ssh;
        for (k = 0; k < readz_hom; k++) {
            if (m.rank == 0 && k == 0) {
                printf("Start\n");
                fflush(stdout);
            }
            z = c * nio + k;
            if (z >= nz_ssh) break;
            tapval = 1.0 - (z >= nz_tap ? (float)(z + 1 - nz_tap) / (nz_ssh - nz_tap) : 0);
            #if COUNT
            int count = 0;
            float vvp_max = 0., ssh_max = 0.;
            #endif
            for (j = 0; j < m.nyt; ++j) {
                for (i = 0; i < m.nxt; ++i) { 
                    ind = k * m.nxt * m.nyt + j * m.nxt + i;
                    ssh = std * buf_ssh[ind] * tapval;
                    if (buf_hom[nvar * ind + 1] <= 0.0001) continue; // No SSH for water
                    // Vp
                    tmp_vp = buf_hom[nvar * ind] / (1 - ssh); // David Gill use (1 + ssh), should be equivalent
                    buf_hom[nvar * ind] = (tmp_vp > vpmax) ? vpmax : tmp_vp;
                    // Vs
                    tmp_vs = buf_hom[nvar * ind + 1] / (1 - ssh);
                    buf_hom[nvar * ind + 1] = (tmp_vs > vsmax) ? vsmax : tmp_vs;
                    // Rho
                    buf_hom[nvar * ind + 2] = buf_hom[nvar * ind + 2] * (1 + ssh);

                    #if COUNT
                    // how many points having vp exceeding the given vpmax
                    if (tmp_vp > vpmax) {
                        count += 1;
                    }
                    vvp_max = vvp_max > tmp_vp ? vvp_max : tmp_vp;
                    ssh_max = ssh_max > buf_ssh[ind] ? ssh_max : buf_ssh[ind];
                    #endif
                }
            }
            #if COUNT
            if (count > 1000) {
                fprintf(stdout, "Rank=%d, Processing layer %d / %d, tapval=%f.\n", m.rank, k, nz_ssh, tapval);
                fprintf(stdout, "count = %d, vvp_max=%f, vpmax=%f, ssh_max=%f\n", count, vvp_max, vpmax, ssh_max);
                fflush(stdout);
            }
            #endif
        }
        fprintf(stdout, "Rank=%d, c=%d/%d, writing to the heterogenous file.\n", m.rank, c, kz_hom);
        fflush(stdout);
        MPICHK(MPI_File_write_all(ft, buf_hom, bufsize_hom, MPI_FLOAT, &filestatus));
        fprintf(stdout, "Rank=%d, written done to the heterogenous file.\n", m.rank);
        fflush(stdout);

        free(buf_ssh);
        free(buf_hom);
        MPI_Barrier(m.MCW);
    }
    if (m.rank == 0) fprintf(stdout, "Done\n");
    MPICHK(MPI_File_close(&fs));
    MPICHK(MPI_File_close(&fm));
    MPICHK(MPI_File_close(&ft));
    MPI_Barrier(m.MCW);
    MPI_Finalize();

    return err;
}

void mpi_init(struct Mpi * m, int argc, char ** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &m->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &m->size);
    MPI_Comm_dup(MPI_COMM_WORLD, &m->MCW);
}


void mpi_cart(struct Mpi *m, const int * size, const int * part)
{
    m->nxt = size[0] / part[0];
    m->nyt = size[1] / part[1];
    m->nzt = size[2];
    m->dim[0] = part[0];
    m->dim[1] = part[1];
    m->period[0] = 0;
    m->period[1] = 0;
    m->reorder = 0;
    MPICHK(MPI_Cart_create(m->MCW, 2, m->dim, m->period, m->reorder, 
                    &m->MC1));
    MPICHK(MPI_Cart_coords(m->MC1, m->rank, 2, m->coord));
}


MPI_Datatype data_type(const struct Mpi *m, int num_z, int nvar)
{
    int old[3] = {num_z, m->nyt * m->dim[1], m->nxt * m->dim[0] * nvar};
    int new[3] = {num_z, m->nyt, m->nxt * nvar};
    int offset[3] = {0, m->nyt * m->coord[1], m->nxt * m->coord[0] * nvar};
    MPI_Datatype dtype;
    MPICHK(MPI_Type_create_subarray(3, old, new, offset, MPI_ORDER_C,
                    MPI_FLOAT, &dtype));
    MPICHK(MPI_Type_commit(&dtype));
    return dtype;
}
