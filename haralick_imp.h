#ifndef HARALICK_IMP_H
#define HARALICK_IMP_H

#include "misc.h"

/*
    Structs
*/
typedef struct {
    unsigned char* arr;
    int width;
    int height;
} img_data;

typedef struct {
    int* arr;
    int size;       //total array length  = size * size
    int angle;
    int dist;
} gclm_data;


/*
    Helpers
*/
unsigned char get_first_neighbor(unsigned char image[], int x, int y, int angle, int distance, int width, int height);
unsigned char get_second_neighbor(unsigned char image[], int x, int y, int angle, int distance, int width, int height);
int num_gray_levels(unsigned char image[], int size);
void check_gclm(mpi_data mpi, img_data img, gclm_data target_gclm);

/*
    Implementation
*/
void sequential(mpi_data mpi, img_data img, gclm_data gclm);
void sync_vertical_split(mpi_data mpi, img_data img, gclm_data gclm);

#endif
