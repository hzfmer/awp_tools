#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

MPI_Aint compute_odc_offset(int nx, int ny, int nvar, int xi, int yi, int zi){
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

void downsample_odc_mod(int rank, int nprocs, int nvar, int nx, int ny, int nz, int k0,
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

   nx2=nx / step;
   ny2=ny / step;
   nz2=nz / step;
   
   np=nx2*ny2;
   fprintf(stdout, "nx2, ny2, nz2, np, is %d, %d, %d, %d\n", nx2, ny2, nz2, np);
   if ((nz2 % nprocs) != 0) {
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

   for (k=k0+rank*step; k<k0+nz; k+=nprocs*step){
      if (rank==0) fprintf(stdout, "reading slice %d\n", k);
      p=0;
      for (j=1; j<ny+1; j+=step){
         for (i=0; i<nx; i+=step){
             blocklen[p]=nvar;
             map[p]=compute_odc_offset(nx, ny, nvar, i, j, k);
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
      ierr=MPI_Type_free(&filetype);
      errhandle(ierr, "Type_free");

      if (nvar == 8){
	  p=0;
	  for (j=0; j<ny2; j++){
	     for (i=0; i<nx2; i++){
		 int xi = (int) buf[p*8] - 1;
		 if ( (xi % step) != 0) {
		   fprintf(stderr, "index %d not divisible by step (%d,%d).  stopping\n", xi, j, i);
		   MPI_Finalize();
		 }
		 buf[p*8] = (buf[p*8] - 1) / step + 1;
		 buf[p*8+1] = (buf[p*8+1] - 1) / step + 1;
		 buf[p*8+2] = (buf[p*8+2] - 1) / step + 1;
		 p++;
	     }
	  }
      }
      disp2=(MPI_Aint) nx2* (MPI_Aint) ny2 * (k-k0)/step *sizeof(float) * (MPI_Aint) nvar;
      MPI_File_write_at_all(outfid, disp2, buf, np*nvar, MPI_FLOAT, MPI_STATUS_IGNORE);
   }

   MPI_File_close(&infid);
   MPI_File_close(&outfid);
   MPI_Barrier(MPI_COMM_WORLD);
   }

int main (int argc, char *argv[] ) {
   int rank, nprocs; 
   //int nx=9000, ny=6750, nz=1296, k0=252;
   int nx=2000, ny=960, nz=180, k0=60;
   int nvar=5;
   int skip=3;
   char *infile="mesh_TKCH05_het_081318_h005l200.bin";
   char *outfile="mesh_TKCH05_het_081318_h005l200_p1.bin";

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   downsample_odc_mod(rank, nprocs, nvar, nx, ny, nz, k0, skip, infile, outfile);

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   return(0);

}

