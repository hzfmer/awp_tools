#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

MPI_Aint compute_odc_offset(int nx, int ny, int nz, int nvar, 
    int xi, int yi, int zi){
   MPI_Aint offset;
   offset = (MPI_Aint) ( (MPI_Aint) zi*nx*ny + yi*nx + xi) * nvar * sizeof(float);
   if (offset < 0) fprintf(stderr, "Error: offset is %ld.\n", offset);
   return offset;
}

void errhandle(int ierr, char *where){
   int errlen;
   char *errstr;

   if (ierr != 0) {
      fprintf(stderr, "error in %s\n", where);
      errstr=calloc(500, sizeof(char));
      MPI_Error_string(ierr, errstr, &errlen);
      fprintf(stderr, "%s", errstr);
      MPI_Finalize();
   }
}

void downsample_odc_mod(int rank, int nprocs, int nvar, int nx, int ny, int nz,
   int step, char *infile, char *outfile){

   int i, j, k, np, p;
   MPI_File infid, outfid;
   int *blocklen;
   MPI_Aint *map;
   MPI_Datatype filetype;
   MPI_Offset disp=0, disp2;
   float *buf;
   int ierr;
   int nx2, ny2, nz2;

   nx2=nx / 2;
   ny2=ny / 2;
   nz2=nz / 2;

   np=nx2*ny2;

   if ((nz % nprocs) != 0) {
        fprintf(stderr, "NZ of output is not divisible by NCPUS. Quitting\n");
        MPI_Finalize();
   }
   
   blocklen=(int*) calloc(np, sizeof(int));
   map=(MPI_Aint*) calloc (np, sizeof(MPI_Aint));
   buf=(float*) calloc(np*nvar, sizeof(float));

   MPI_File_open(MPI_COMM_WORLD, infile, MPI_MODE_RDONLY, 
       MPI_INFO_NULL, &infid);
   MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &outfid);

   for (k=0+rank*step; k<nz; k+=nprocs*step){
      if (rank==0) fprintf(stdout, "reading slice %d\n", k);
      p=0;
      for (j=0; j<ny; j+=step){
         for (i=0; i<nx; i+=step){
             blocklen[p]=nvar;
             map[p]=compute_odc_offset(nx, ny, nz, nvar, i, j, k);
             p++;
         }
      }
      //this won't work for large files
      //ierr=MPI_Type_create_indexed_block(np, nvar, map, MPI_FLOAT, &filetype);
      ierr=MPI_Type_create_hindexed(np, blocklen, map, MPI_FLOAT, &filetype);
      errhandle(ierr, "create_indexed_block");
      ierr=MPI_Type_commit(&filetype);
      errhandle(ierr, "commit");
      ierr=MPI_File_set_view(infid, disp, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
      errhandle(ierr, "set_view");
      ierr=MPI_File_read_all(infid, buf, np*nvar, MPI_FLOAT, MPI_STATUS_IGNORE);
      errhandle(ierr, "read_all");

      p=0;
      for (j=0; j<ny2; j++){
         for (i=0; i<nx2; i++){
             int xi = (int) buf[p*8] - 1;
             if ( (xi % step) != 0) {
               fprintf(stderr, "index %d not divisible by step (%d,%d).  stopping\n", xi, j, i);
               MPI_Finalize();
             }
             buf[p*8] = (buf[p*8] - 1) / 2 + 1;
             buf[p*8+1] = (buf[p*8+1] - 1) / 2 + 1;
             buf[p*8+2] = (buf[p*8+2] - 1) / 2 + 1;
             p++;
         }
      }
      disp2=(MPI_Aint) nx2* (MPI_Aint) ny2 * k/2 *sizeof(float) * (MPI_Aint) nvar;
      MPI_File_write_at_all(outfid, disp2, buf, np*nvar, MPI_FLOAT, MPI_STATUS_IGNORE);
   }

   MPI_File_close(&infid);
   MPI_File_close(&outfid);
   MPI_Barrier(MPI_COMM_WORLD);
   }

int main (int argc, char *argv[] ) {
   int rank, nprocs;
   #ifdef SAFDYN50SM
   int nx=10800, ny=800, nz=800;
   int nvar=8;
   int skip=2;
   char *infile="SAF_dyn_25m_small_het+idx";
   char *outfile="SAF_dyn_50m_small_het+idx";
   #endif
   #ifdef SAFDYN50
   int nx=12000, ny=5488, nz=2048;
   int nvar=8;
   int skip=2;
   char *infile="SAF_dyn_25m_full_het+idx";
   char *outfile="SAF_dyn_50m_full_het+idx";
   #endif
   #ifdef TEST
   int nx=800, ny=600, nz=400;
   int nvar=8;
   int skip=2;
   char *infile="/Users/daniel/Work.4//tpv31-32/tpv31/50m/velmod/tpv31_50m.odc";
   char *outfile="test.odc";
   #endif

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   downsample_odc_mod(rank, nprocs, nvar, nx, ny, nz, skip, infile, outfile);

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   return(0);

}
