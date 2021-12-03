#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

void program_abort(char *exec_name, char *message) {
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  if (my_rank == 0) {
    if (message) {
      fprintf(stderr,"%s",message);
    }
    if (exec_name) {
      print_usage(exec_name);
    }
  }
  MPI_Finalize();
  exit(1);
}

// Print the usage information
void print_usage(char *exec_name) {
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

  if (my_rank == 0) {
    fprintf(stderr,"Usage: smpirun --cfg=smpi/bcast:mpich -np <num processes>\n");
    fprintf(stderr,"              -platform <XML platform file> -hostfile <host file>\n");
    fprintf(stderr,"              %s <image path> <angle> [-d <distance>] [-i implementation] \n",exec_name);
    fprintf(stderr,"MPIRUN arguments:\n");
    fprintf(stderr,"\t<num processes>: number of MPI processes\n");
    fprintf(stderr,"\t<XML platform file>: a Simgrid platform description file\n");
    fprintf(stderr,"\t<host file>: MPI host file with host names from the platform file\n");
    fprintf(stderr,"PROGRAM arguments:\n");
    fprintf(stderr,"\t<image path>: path to image to be processed\n");
    fprintf(stderr,"\t<angle>: haralick angle, can be one of the following (0, 90, 45, 135) \n");
    fprintf(stderr,"\t[-d <distance>]: haralick distance (optional, default to 1)\n");
    fprintf(stderr,"\t[-i <implementation>]: haralick implementation (optional)\n");
    fprintf(stderr,"\n");
  }
  return;
}

void init_arr_zero(int* arr, int size) {
	for(int i = 0; i < size; i++){
		arr[i] = 0;
	}
}

