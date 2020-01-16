#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <errno.h>

float merge_multiplexed_data(char *infbase, char *outfile, int rank, int step, long int count, int nwrite){
   MPI_Offset disp;
   char infilename[200];
   int fnum;
   float *buf;
   long int count2;
   int p;

   FILE *fid;
   long int nread;

   MPI_File fh;
   MPI_Offset offset, offset2;
   int err;
   char errmsg[200];
   int errlen;

   if ((count % nwrite) != 0) {
      fprintf(stdout, "total buffer size %ld must be divisible by number of writes %d\n", count, nwrite);
      return(-1);
   }

   count2=count / nwrite;

   fnum=(rank+1)*step;
   sprintf(infilename, "%s%07d", infbase, fnum);
   if (rank==0) fprintf(stdout, "Reading input file %s on rank 0\n", infilename);
   
   fid=fopen(infilename, "r");
   if (fid==NULL) {
      sprintf(errmsg, "could not open %s", infilename);
      perror(errmsg);
      return(-1);
   }

   buf=(float*) calloc(count, sizeof(float));
   if (buf==NULL){
      sprintf(errmsg, "could not allocate %d bytes of memory", count*sizeof(float));
      perror(errmsg);
      return(-1);
   } 

   nread=fread(buf, sizeof(float),count, fid);
   if (nread != count){
      fprintf(stderr, "nread=%d\n", nread);
      perror("could not read input data");
      return(-1);
   }

   MPI_Barrier(MPI_COMM_WORLD);
   if (rank==0) fprintf(stdout, "Done reading input files\nWriting data to single output file %s.\n", outfile);

   err=MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_CREATE | MPI_MODE_WRONLY,
        MPI_INFO_NULL, &fh);

   if (err != MPI_SUCCESS) {
      fprintf(stderr, "Error in MPI_File_open");
      MPI_Error_string(err, errmsg, &errlen);
      fprintf(stderr, errmsg);
      return(-1);
   }

   for (p=0; p<nwrite; p++){
      offset2 = (MPI_Offset) p * (MPI_Offset) count2 * (MPI_Offset) sizeof(float);
      offset = (MPI_Offset) count * (MPI_Offset) rank * (MPI_Offset) sizeof(float) + offset2;
      err=MPI_File_write_at_all(fh, offset, buf+p*count2, count2, MPI_FLOAT, MPI_STATUS_IGNORE);
      //fprintf(stderr, "count2=%ld, offset=%ld\n", count2, offset);
      if (rank==0) fprintf(stdout, "writing chunk %d out of %d\n", p+1, nwrite);
      if (err != MPI_SUCCESS) {
	 fprintf(stderr, "count2=%ld, offset=%ld\n", count2, offset);
	 fprintf(stderr, "Error in MPI_File_write_at_all\n");
	 MPI_Error_string(err, errmsg, &errlen);
	 fprintf(stderr, errmsg);
	 return(-1);
      }
   }

   MPI_Barrier(MPI_COMM_WORLD);
   err=MPI_File_close(&fh);
   if (err != MPI_SUCCESS) {
      fprintf(stderr, "Error in MPI_File_close");
      MPI_Error_string(err, errmsg, &errlen);
      fprintf(stderr, errmsg);
      return(-1);
   }
   if (rank==0) fprintf(stdout, "MPI_File_write_at_all() completed, file closed.\n");

   free(buf);

   return(0);

}
      

int main (int argc, char **argv){
   int rank, ncpus; 
   int nx=6540, ny=10728, wstep=100, ntiskp=10, step;
   //int nx=480, ny=1120, wstep=100, ntiskp=1, step;
   long int count;
   char infbase[100], outfname[100];
   int nwrite=20;

   count=(long int) nx*(long int) ny*(long int) wstep;
   step=ntiskp*wstep;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   merge_multiplexed_data("output_sfc/SX_0_", "output_sfc_demulti/SX", rank, step, count, nwrite);
   merge_multiplexed_data("output_sfc/SY_0_", "output_sfc_demulti/SY", rank, step, count, nwrite);
   merge_multiplexed_data("output_sfc/SZ_0_", "output_sfc_demulti/SZ", rank, step, count, nwrite);

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   return(0);

}
