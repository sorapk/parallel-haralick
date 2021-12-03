#ifndef MISC_H
#define MISC_H

typedef struct {
    int rank;
    int num_procs;
} mpi_data;

typedef struct {

} param;

void program_abort(char *exec_name, char *message);
void print_usage();
void init_arr_zero(int* arr, int size);
#endif
