#include "haralick_imp.h"
#include <mpi.h>
#include <math.h>



#define ANGLE_0   0
#define ANGLE_90  90
#define ANGLE_45  45
#define ANGLE_135 135

/*
 Private Prototypes
*/
void get_neighbors(int x, int y, img_data img, gclm_data gclm, int* neighbor1, int* neighbor2);

/* 
 * Deprecated Prototypes
 */
void sequential(mpi_data mpi, img_data img, gclm_data gclm); //Still used in check_gclm for consistencies
unsigned char get_first_neighbor(unsigned char image[], int x, int y, int angle, int distance, int width, int height);
unsigned char get_second_neighbor(unsigned char image[], int x, int y, int angle, int distance, int width, int height);


/*
  Helpers
*/
void get_neighbors(int x, int y, img_data img, gclm_data gclm, int* neighbor1, int* neighbor2) {
	if (gclm.angle == ANGLE_0) {
		*neighbor1 = (x - gclm.dist) >= 0        ? img.arr[ img.width * y + (x - gclm.dist)] : -1;
		*neighbor2 = (x + gclm.dist) < img.width ? img.arr[ img.width * y + (x + gclm.dist)] : -1;
	} else if (gclm.angle == ANGLE_45) {
		*neighbor1 = (x + gclm.dist) < img.width && (y - gclm.dist) >= 0         ? img.arr[ img.width * (y - gclm.dist) + (x + gclm.dist)] : -1;
		*neighbor2 = (x - gclm.dist) >= 0        && (y + gclm.dist) < img.height ? img.arr[ img.width * (y + gclm.dist) + (x - gclm.dist)] : -1;
	} else if (gclm.angle == ANGLE_90) {
		*neighbor1 = (y - gclm.dist) >= 0         ? img.arr[ img.width * (y - gclm.dist) + x] : -1;
		*neighbor2 = (y + gclm.dist) < img.height ? img.arr[ img.width * (y + gclm.dist) + x] : -1;
	} else if (gclm.angle == ANGLE_135) {
		*neighbor1 = (x - gclm.dist) >= 0        && (y - gclm.dist) >= 0         ? img.arr[ img.width * (y - gclm.dist) + (x - gclm.dist)] : -1;
		*neighbor2 = (x + gclm.dist) < img.width && (y + gclm.dist) < img.height ? img.arr[ img.width * (y + gclm.dist) + (x + gclm.dist)] : -1;
	} else {
		program_abort(NULL, "Invalid Angle!");
	}
}
int num_gray_levels(unsigned char image[], int size){
     int max = 0;
     for (int i = 1; i < size; ++i)
     {
        if (image[i] > max){
			max = image[i];
		}
     }
     return max+1;
}
int check_gclm(mpi_data mpi, img_data img, gclm_data target_gclm) {

	int good = 1;
	int* ref_glcm;
	
	//Perform Reference Sequential Calculation
	#ifdef DEBUG
		printf("checking calculation...\n");
	#endif
	if ((ref_glcm = (int*) malloc(sizeof(int) * target_gclm.size * target_gclm.size)) == NULL) {
		program_abort(NULL, "Out of memory - 3");
	}  
	init_arr_zero(ref_glcm, target_gclm.size * target_gclm.size);
	sequential(mpi, img, (gclm_data){
		.arr = ref_glcm,
		.dist = target_gclm.dist,
		.size = target_gclm.size,
		.angle = target_gclm.angle, 
	});
	// Compare
	for(int j = 0; j < target_gclm.size * target_gclm.size; j++ ){
		if (ref_glcm[j] != target_gclm.arr[j]) {
			#ifdef DEBUG
				printf("ERROR: mismatch found at index = %d\n",j);
			#endif 
			good = 0;
			break;
		}
	}
	
	#ifdef DEBUG
	if (good){
		printf("all good!\n");
	}
	#endif 
	free(ref_glcm);

	return good;
}

/*
 Implementations
*/
void sequential_v2(mpi_data mpi, img_data img, gclm_data gclm) {
	int x, y, neighbor1, neighbor2;
	if (mpi.rank == 0) {
		for (y = 0; y < img.height; y++){
			for (x = 0; x < img.width; x++){
				int my_value = (int) img.arr[y*img.width + x];
				get_neighbors(x, y, img, gclm, &neighbor1, &neighbor2);
				if (neighbor1 >= 0) 
					gclm.arr[my_value*gclm.size + neighbor1]++;
				if (neighbor2 >= 0)
					gclm.arr[my_value*gclm.size + neighbor2]++;
			}
		}
	}
}
void sync_vertical_split_v2(mpi_data mpi, img_data img, gclm_data gclm){
	int x, y, neighbor1, neighbor2;
	int range = img.height / mpi.num_procs;
	int  *partial_GLCM;
	MPI_Status status;

	// Init Partial GCLM
	if ((partial_GLCM = (int*) malloc(sizeof(int) * gclm.size * gclm.size)) == NULL) {
    	program_abort(NULL, "Out of memory! - 2");
  	}
	init_arr_zero(partial_GLCM, gclm.size * gclm.size);
	
	// Rank 0
	if (mpi.rank == 0) {
		// Compute Partial GCLM
		for (y = 0; y < fmin(range, img.height); y++){
			for (x = 0; x < img.width; x++){
				int my_value = img.arr[y*img.width + x];
				get_neighbors(x, y, img, gclm, &neighbor1, &neighbor2);
				if (neighbor1 >= 0) 
					gclm.arr[my_value*gclm.size + neighbor1]++;
				if (neighbor2 >= 0)
					gclm.arr[my_value*gclm.size + neighbor2]++;
			}
		}
		// Aggregate Partial GCLM into GCLM
		for(int id = 1; id < mpi.num_procs; id++) {
			MPI_Recv(partial_GLCM, gclm.size*gclm.size, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			for(int j = 0; j < gclm.size * gclm.size; j++ ){
					gclm.arr[j] = gclm.arr[j] + partial_GLCM[j];
			}
		}
	} 
	// Other Ranks
	else {
		// Compute Partial GCLM
		for (y = fmax(0, mpi.rank * range); y < fmin((mpi.rank + 1) * range, img.height); y++){
			for (x = 0; x < img.width; x++){
				int my_value = img.arr[y*img.width + x];
				get_neighbors(x, y, img, gclm, &neighbor1, &neighbor2);
				if (neighbor1 >= 0) 
					partial_GLCM[my_value*gclm.size + neighbor1]++;
				if (neighbor2 >= 0)
					partial_GLCM[my_value*gclm.size + neighbor2]++;
			}
		}
		// Send Partial GCLM
		MPI_Send(partial_GLCM, gclm.size*gclm.size, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
}
void sync_tiling_split(mpi_data mpi, img_data img, gclm_data gclm) {
	
	int N = img.width;
	int root_p = (int) sqrt((double) mpi.num_procs);
	int rank_row = mpi.rank / root_p;
	int rank_col = mpi.rank % root_p;
	int block = img.width / root_p;
	int block_size = block * block;

	int i,j;
	int neighbor1, neighbor2;
	int  *partial_GLCM;
	MPI_Status status;

	// Init Partial GCLM
	if ((partial_GLCM = (int*) malloc(sizeof(int) * gclm.size * gclm.size)) == NULL) {
    	program_abort(NULL, "Out of memory! - 2");
  	}
	init_arr_zero(partial_GLCM, gclm.size * gclm.size);
	
	// Rank 0
	if (mpi.rank == 0) {
		for (i = 0; i < block; i++) {
			for (j = 0; j < block; j++) {
				int y = rank_row * block + i;
				int x = rank_col * block + j;
				int my_value = img.arr[y*img.width + x];
				get_neighbors(x, y, img, gclm, &neighbor1, &neighbor2);
				if (neighbor1 >= 0) 
					gclm.arr[my_value*gclm.size + neighbor1]++;
				if (neighbor2 >= 0)
					gclm.arr[my_value*gclm.size + neighbor2]++;
			}
		}
		// Aggregate Partial GCLM into GCLM
		for(int id = 1; id < mpi.num_procs; id++) {
			MPI_Recv(partial_GLCM, gclm.size*gclm.size, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			for(int j = 0; j < gclm.size * gclm.size; j++ ){
				gclm.arr[j] = gclm.arr[j] + partial_GLCM[j];
			}
		}
	} 
	// Other Ranks
	else {
		// Compute Partial GCLM
		for (i = 0; i < block; i++) {
			for (j = 0; j < block; j++) {
				int y = rank_row * block + i;
				int x = rank_col * block + j;
				int my_value = img.arr[y*img.width + x];
				get_neighbors(x, y, img, gclm, &neighbor1, &neighbor2);
				if (neighbor1 >= 0) 
					partial_GLCM[my_value*gclm.size + neighbor1]++;
				if (neighbor2 >= 0)
					partial_GLCM[my_value*gclm.size + neighbor2]++;
			}
		}
		// Send Partial GCLM
		MPI_Send(partial_GLCM, gclm.size*gclm.size, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
}
/*
	Deprecated Implementations
*/
unsigned char get_first_neighbor(unsigned char image[], int x, int y, int angle, int distance, int width, int height) {
	// Get the gray value of the left neighbor for 0 angle 
	if (angle == 0) {
		//printf("myself: %d, %d \n ", x, y);
		//printf("neighbor: %d, %d \n ", x - 1, y);
		if (x == 0) {
			return -1;
		}
		else {
			//printf("value: %u", image[(x - distance) + width*y]);
			return image[(x - distance) + width*y];
		}
	}
	// Get the gray value of the top neighbor for 90 angle 
	else if (angle == 90) {
		if (y == 0) {
			return -1;
		}
		else {
			return image[x + width*(y - distance)];
		}
	}
	// Get the gray value of the top/right neighbor for 45 angle
	else if (angle == 45) {
		if (y == 0 ||x == width - 1) {
			return -1;
		}
		else {
			return image[(x + distance) + width*(y - distance)];
		}
	}
	// Get the gray value of the top/left neighbor for 135 angle 
	else {
		if (y == 0 || x == 0) {
			return -1;
		}
		else {
			return image[(x - distance) + width*(y - distance)];
		}
	}
}
unsigned char get_second_neighbor(unsigned char image[], int x, int y, int angle, int distance, int width, int height) {
	// Get the gray value of the right neighbor for 0 angle 
	if (angle == 0) {
		if (x == width - 1) {
			return -1;
		}
		else {
			return image[(x + distance) + width*y];
		}
	}
	// Get the gray value of the bottom neighbor for 90 angle 
	else if (angle == 90) {
		if (y == height - 1) {
			return -1;
		}
		else {
			return image[x + width*(y + distance)];
		}
	}
	// Get the gray value of the bottom/left neighbor for 45 angle
	else if (angle == 45) {
		if ( x == 0 || y == height - 1) {
			return -1;
		}
		else {
			return image[(x - distance) + width*(y + distance)];
		}
	}
	// Get the gray value of the bottom/right neighbor for 135 angle 
	else {
		if (y == height - 1 || x == width - 1) {
			return -1;
		}
		else {
			return image[(x + distance) + width*(y + distance)];
		}
	}
}
void sequential(mpi_data mpi, img_data img, gclm_data gclm) {
	int x, y;
	if (mpi.rank == 0) {
		for (y = 0; y < img.height; y++){
			for (x = 0; x < img.width; x++){
				int my_value = (int) img.arr[y*img.width + x];
				int first_neighbor_value =  (int) get_first_neighbor(img.arr, x, y, gclm.angle, gclm.dist, img.width, img.height);
				int second_neighbor_value = (int) get_second_neighbor(img.arr, x, y, gclm.angle, gclm.dist, img.width, img.height);
				if ((gclm.angle == 0 && x > 0) || (gclm.angle == 90 && y > 0) || (gclm.angle == 45 && (y > 0 && x < img.width - 1)) || (gclm.angle == 135 && (y > 0 && x > 0))) {
					gclm.arr[my_value*gclm.size + first_neighbor_value]++;
				}
				if ((gclm.angle == 0 && x < img.width - gclm.dist) || (gclm.angle == 90 && y < img.height - gclm.dist) || (gclm.angle == 45 && (y < img.height - gclm.dist && x > 0)) || (gclm.angle == 135 && (y < img.height - gclm.dist && x < img.width - gclm.dist))) {
					gclm.arr[my_value*gclm.size + second_neighbor_value]++;
				}
			}
		}
	}
}
void sync_vertical_split(mpi_data mpi, img_data img, gclm_data gclm){
	int y, x;
	int range = img.height / mpi.num_procs;
	int  *partial_GLCM;
	MPI_Status status;

	// Init Partial GCLM
	if ((partial_GLCM = (int*) malloc(sizeof(int) * gclm.size * gclm.size)) == NULL) {
    	program_abort(NULL, "Out of memory! - 2");
  	}
	init_arr_zero(partial_GLCM, gclm.size * gclm.size);
	
	// Rank 0
	if (mpi.rank == 0) {
		// Compute Partial GCLM
		for (y = 0; y < fmin(range, img.height); y++){
			for (x = 0; x < img.width; x++){
				int my_value = img.arr[y*img.width + x];
				int first_neighbor_value = (int) get_first_neighbor(img.arr, x, y, gclm.angle, 1, img.width, img.height);
				int second_neighbor_value = (int) get_second_neighbor(img.arr, x, y, gclm.angle, 1, img.width, img.height);
				if ((gclm.angle == 0 && x > 0) || (gclm.angle == 90 && y > 0) || (gclm.angle == 45 && (y > 0 && x < img.width - 1)) || (gclm.angle == 135 && (y > 0 && x > 0))) {
					gclm.arr[my_value*gclm.size + (int) first_neighbor_value]++;
				}
				if ((gclm.angle == 0 && x < img.width - gclm.dist) || (gclm.angle == 90 && y < img.height - gclm.dist) || (gclm.angle == 45 && (y < img.height - gclm.dist && x > 0)) || (gclm.angle == 135 && (y < img.height - gclm.dist && x < img.width - gclm.dist))) {
					gclm.arr[((int) my_value)*gclm.size + (int) second_neighbor_value]++;
				}
			}
		}

		// Aggregate Partial GCLM into GCLM
		for(int id = 1; id < mpi.num_procs; id++) {
			MPI_Recv(partial_GLCM, gclm.size*gclm.size, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			for(int j = 0; j < gclm.size * gclm.size; j++ ){
					gclm.arr[j] = gclm.arr[j] + partial_GLCM[j];
			}
		}
	} 
	// Other Ranks
	else {
		// Compute Partial GCLM
		for (y = fmax(0, mpi.rank * range); y < fmin((mpi.rank + 1) * range, img.height); y++){
			for (x = 0; x < img.width; x++){
				int my_value = img.arr[y*img.width + x];
				int first_neighbor_value = (int) get_first_neighbor(img.arr, x, y, gclm.angle, 1, img.width, img.height);
				int second_neighbor_value = (int) get_second_neighbor(img.arr, x, y, gclm.angle, 1, img.width, img.height);
				if ((gclm.angle == 0 && x > 0) || (gclm.angle == 90 && y > 0) || (gclm.angle == 45 && (y > 0 && x < img.width - 1)) || (gclm.angle == 135 && (y > 0 && x > 0))) {
					partial_GLCM[my_value*gclm.size + (int) first_neighbor_value]++;
				}
				if ((gclm.angle == 0 && x < img.width - gclm.dist) || (gclm.angle == 90 && y < img.height - gclm.dist) || (gclm.angle == 45 && (y < img.height - gclm.dist && x > 0)) || (gclm.angle == 135 && (y < img.height - gclm.dist && x < img.width - gclm.dist))) {
					partial_GLCM[((int) my_value)*gclm.size + (int) second_neighbor_value]++;
				}
			}
		}
		// Send Partial GCLM
		MPI_Send(partial_GLCM, gclm.size*gclm.size, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
}
