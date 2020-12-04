#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

void error_check(int ierr, char *message){
   char errmsg[500];
   int errlen;
   if (ierr != MPI_SUCCESS) {
      fprintf(stderr, "%d: Error in %s\n", ierr, message);
      MPI_Error_string(ierr, errmsg, &errlen);
      fprintf(stderr, errmsg);
   }
}


/* reads several time series from contiguous sites (non-multiplexed) */
void read_awp_timeseries(char *fname, int nx, int ny, int nz, int xi, int yi, int zi, 
    int nt, int nsites, float *buf){
    MPI_File fh;
    MPI_Datatype filetype;
    MPI_Offset disp;
    int ierr;
    MPI_Info info;
    //MPI_Aint extent1, extent2;

    ierr=MPI_Info_create(&info);
    error_check(ierr, "MPI_Info_create()");
    /* Enable collective buffering optimization */
    ierr=MPI_Info_set(info, "romio_cb_write", "enable");
    error_check(ierr, "MPI_Info_set()");

    #ifdef NONCOLLECTIVE
    ierr=MPI_File_open(MPI_COMM_SELF, fname, MPI_MODE_RDONLY, info, &fh);
    #else
    ierr=MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, info, &fh);
    #endif
    error_check(ierr, "MPI_File_open()");

    ierr=MPI_Type_vector(nt, nsites, nx*ny*nz, MPI_FLOAT, &filetype);
    error_check(ierr, "MPI_Type_vector()");
    ierr=MPI_Type_commit(&filetype);
    error_check(ierr, "MPI_Type_commit()");

    /*MPI_File_get_type_extent(fh, MPI_FLOAT, &extent1);   
    MPI_File_get_type_extent(fh, filetype, &extent2);   
    fprintf(stdout, "extent of MPI_FLOAT = %ld, extent of filetype = %ld\n", extent1, extent2);*/

    disp = (MPI_Offset) sizeof(float) * ((MPI_Offset) nx*ny*zi + (MPI_Offset) nx*yi + (MPI_Offset) xi);
    if (disp==0) fprintf(stdout, "[Single]reading %s ...\n", fname);
    ierr=MPI_File_set_view(fh, disp, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    //fprintf(stderr, "disp=%ld, %s\n", disp, fname);
    error_check(ierr, "MPI_File_set_view()");

    #ifdef NONCOLLECTIVE
    ierr=MPI_File_read(fh, buf, nt*nsites, MPI_FLOAT, MPI_STATUS_IGNORE);
    #else
    ierr=MPI_File_read_all(fh, buf, nt*nsites, MPI_FLOAT, MPI_STATUS_IGNORE);
    #endif
    error_check(ierr, "MPI_File_read_all()");
    MPI_File_close(&fh);

    ierr=MPI_Type_free(&filetype);
    error_check(ierr, "MPI_Type_free()");

    ierr=MPI_Info_free(&info);
    error_check(ierr, "MPI_Info_free()");
}

/* reads several time series from contiguous sites (multiplexed) */
void read_awp_timeseries_multi(char *fbase, int nx, int ny, int nz, int xi, int yi, int zi, 
    int nt, int wstep, int nskip, int nsites, float *buf){
   int nfiles, k;
   MPI_Offset bpos;
   char fname[200];

   nfiles = nt / wstep;
   for (k=0; k<nfiles; k++){
      sprintf(fname, "%s%07d", fbase, (k+1)*wstep*nskip);
      fprintf(stdout, "[multiplexed] reading %d / %d: %s ...\n", k + 1, nfiles, fname);
      fflush(stdout);
      bpos=k*wstep*nsites;
      read_awp_timeseries(fname, nx, ny, nz, xi, yi, zi, wstep, nsites, buf+bpos);
   }
};

/*int main (int argc, char **argv){
    int rank, nprocs;
    float *buf, **x;
    int k, l;
    int nsites=701724, nt=8700;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    buf=(float*) calloc(nt*nsites, sizeof(float));

    read_awp_timeseries("output_sfc_demulti/SX", 3270, 5365, 1, 0, 0, 0, nt, nsites, buf);

    x=(float**) calloc(nsites, sizeof(float*));
    for (l=0; l<nsites; l++) x[l] = (float*) calloc(nt, sizeof(float));

    for (k=0; k<nt; k++){
       for (l=0; l<nsites; l++){
          x[l][k]=buf[k*nsites+l];
       }
    }

    free(buf);
    for (k=0; k<nt; k++) fprintf(stdout, "%e %e\n", x[0][k], x[19][k]);
    return(0);

}*/
