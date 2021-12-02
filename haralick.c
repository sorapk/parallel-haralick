#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

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

/*
int get_first_neighbor_coords(int x, int y, int angle, int distance) {
	// Get the location of the left neighbor for 0 angle 
	if (angle == 0) {
		return (x - distance) + width*y;
	}
	// Get the location of the top neighbor for 90 angle 
	else if (angle == 90) {
		return x + width*(y - distance);
	}
	// Get the location of the top/right neighbor for 45 angle
	else if (angle == 45) {
		return (x + distance) + width*(y - distance);
	}
	// Get the location of the top/left neighbor for 135 angle 
	else {
		return (x - distance) + width*(y - distance);
	}
}
*/
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
/*
char get_second_neighbor_coords(int x, int y, int angle, int distance) {
	// Get the gray value of the right neighbor for 0 angle 
	if (angle == 0) {
		return (x + distance) + width*y;
	}
	// Get the gray value of the bottom neighbor for 90 angle 
	else if (angle == 90) {
		return x + width*(y + distance);
	}
	// Get the gray value of the bottom/left neighbor for 45 angle
	else if (angle == 45) {
		return (x - distance) + width*(y + distance);
	}
	// Get the gray value of the bottom/right neighbor for 135 angle 
	else {
		return (x + distance) + width*(y + distance);
	}
}
*/
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


/*
 * main() function
 */

int main(int argc, char **argv)
{
	if(argc != 4) {
		printf("Please specify a grayscale image path, an angle in (0, 90, 45, 135), and a number of processes \n");
        exit(1);
	}
	FILE *file;
	int angle = atoi(argv[2]);
	int num_procs = atoi(argv[3]);
	int distance = 1;
	
	// Check that we can open the file path provided
	if (file = fopen(argv[1], "r")) {
        fclose(file);
    }
	else {
		printf("%s is not a valid image path \n", argv[1]);
        exit(1);
	}
	
	// Read in the provided image path as a 1D array 
	int width, height, bpp;
	unsigned char *gray_img;
	int *GLCM, *partial_GLCM;
	
	// Load in grayscale image 
	//gray_img = malloc(sizeof(char) * width * height);
    gray_img = stbi_load(argv[1], &width, &height, &bpp, 1);
	
	// Validate angle argument 
	if (angle != 0 && angle != 90 && angle != 135 && angle != 45) {
		printf("Invalid angle provided, angle must be one of 0, 45, 90, 135 \n");
		exit(1);
	}
	
	if (height % num_procs != 0) {
		printf("Image size is not evenly divisible by the number of processes \n");
		exit(1);
	}
	
	MPI_Init(&argc, &argv);
	
	int rank, range;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Status status;
	
	int GLCM_size = num_gray_levels(gray_img, width*height);
	//int GLCM_size = width
	
	partial_GLCM = (int*) malloc(sizeof(int) * GLCM_size * GLCM_size);
	if (!partial_GLCM) { perror("malloc arr"); exit(EXIT_FAILURE); };
	range = height / num_procs;
	
	//printf("I should compute pixel rows %d to %d, for a total of %d rows \nI need access to rows %d to %d to compute my rows \n\n", rank * range, (rank + 1) * range - 1, range, (int) fmax(0, rank * range - distance), (int) fmin((rank + distance) * range, height - 1));
	
	for(int d = 0; d < GLCM_size * GLCM_size; d++ ){
			partial_GLCM[d] = 0;
	}
	
	int y, x;
	unsigned char first_neighbor_value, second_neighbor_value;
	unsigned char my_value;
	
	  double start_time;
  MPI_Barrier(MPI_COMM_WORLD);
	
	if (rank == 0) {
		start_time = MPI_Wtime();
		GLCM = (int*) malloc(sizeof(int) * GLCM_size * GLCM_size);
		if (!GLCM) { perror("malloc arr"); exit(EXIT_FAILURE); };
		
		for(int j = 0; j < GLCM_size * GLCM_size; j++ ){
			GLCM[j] = 0;
		}
		
		
		for (y = 0; y < fmin(range, height); y++){
			for (x = 0; x < width; x++){
				my_value = gray_img[y*width + x];
				//printf("%u <- %u \n", my_value, gray_img[y*width + x]);
				first_neighbor_value = get_first_neighbor(gray_img, x, y, angle, 1, width, height);
				second_neighbor_value = get_second_neighbor(gray_img, x, y, angle, 1, width, height);
				
				if ((angle == 0 && x > 0) || (angle == 90 && y > 0) || (angle == 45 && (y > 0 && x < width - 1)) || (angle == 135 && (y > 0 && x > 0))) {
					//printf("%u, %d, %u => %d \n", my_value, GLCM_size, first_neighbor_value, ((int) my_value)*GLCM_size + (int) first_neighbor_value);
				//printf("%d \n", my_value*GLCM_size + (int) first_neighbor_value);
					GLCM[((int) my_value)*GLCM_size + (int) first_neighbor_value]++;
				}
				if ((angle == 0 && x < width - distance) || (angle == 90 && y < height - distance) || (angle == 45 && (y < height - distance && x > 0)) || (angle == 135 && (y < height - distance && x < width - distance))) {
					GLCM[((int) my_value)*GLCM_size + (int) second_neighbor_value]++;
				}
			   
			}
		}
		
		for(int id = 1; id < num_procs; id++) {
            MPI_Recv(partial_GLCM, GLCM_size*GLCM_size, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			
        for(int j = 0; j < GLCM_size * GLCM_size; j++ ){
				GLCM[j] = GLCM[j] + partial_GLCM[j];
			}
         }
		 
//	for(int j = 0; j < GLCM_size * GLCM_size; j++ ){
	//	printf("%d ", GLCM[j]);
	//}
	}
	else {
		for (y = fmax(0, rank * range); y < fmin((rank + 1) * range, height); y++){
			for (x = 0; x < width; x++){
				my_value = gray_img[y*width + x];
				first_neighbor_value = get_first_neighbor(gray_img, x, y, angle, 1, width, height);
				second_neighbor_value = get_second_neighbor(gray_img, x, y, angle, 1, width, height);
				
				if ((angle == 0 && x > 0) || (angle == 90 && y > 0) || (angle == 45 && (y > 0 && x < width - distance)) || (angle == 135 && (y > 0 && x > 0))) {
					
					//printf("%u, %d, %u => %d/%d \n", my_value, GLCM_size, first_neighbor_value, ((int) my_value)*GLCM_size + (int) first_neighbor_value, GLCM_size*GLCM_size);
					
					partial_GLCM[((int) my_value)*GLCM_size + (int) first_neighbor_value]++;
					
					//partial_GLCM[0]++;
				}
				if ((angle == 0 && x < width - distance) || (angle == 90 && y < height - distance) || (angle == 45 && (y < height - distance && x > 0)) || (angle == 135 && (y < height - distance && x < width - distance))) {
					partial_GLCM[((int) my_value)*GLCM_size + (int) second_neighbor_value]++;
				}
			   
			}
		}
		
		MPI_Send(partial_GLCM, GLCM_size*GLCM_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if (0 == rank) {
    fprintf(stdout,"image size: %d | number of processes: %d |  time: %.3lf seconds\n",
		    width, 
		    num_procs,
		    MPI_Wtime() - start_time);
  }
	  

  // Clean-up
	free(partial_GLCM);
	//free(GLCM);

    MPI_Finalize();
	
	
}
