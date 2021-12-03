#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "haralick_imp.h"
#include "misc.h"

#define DEFAULT_DISTANCE 1
#define IMPLEMENTATION_SIZE 20
#define NUM_ARG 3

/*
 * main() function
 */
int main(int argc, char **argv)
{	
	int i;
	char implementation[IMPLEMENTATION_SIZE];
	implementation[0] = '\0';

	FILE *file;
	int angle, num_procs;
	int distance = DEFAULT_DISTANCE;	
	char *img_path; 

	/*
	 * Init MPI ENV
	*/
	int rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	
	/*
	 *	Parse input arguments
	 */
	if(argc < NUM_ARG) {
		program_abort(argv[0],"Missing arguments\n");
	}
	angle = atoi(argv[2]);
	img_path = argv[1];
	for (i=1; i < argc; i++) {
		if (!strcmp(argv[i],"-d")) {
			if ((i+1 >= argc) || (sscanf(argv[i+1],"%d",&distance) != 1)) {
				program_abort(argv[0],"Invalid <distance> argument\n");
			}
		} else if (!strcmp(argv[i],"-i")) {
			if ((i+1 >= argc) || (sscanf(argv[i+1],"%s",implementation) != 1)) {
				program_abort(argv[0],"Invalid <implementation> argument\n");
			}
			// Check that the implementation name is valid
			if (strcmp(implementation, "sequential") &&
				strcmp(implementation, "sync_vertical_split")) {
				char message[256];
				sprintf(message, "Unknown  implementation name '%s'\n",implementation);
				program_abort(NULL,message);
			}
		}
	}
	if (rank == 0) {
		printf("=== parsed params === \n");
		printf(" implementation: %s \n", implementation);
		printf(" distance: %d \n", distance);
		printf(" angle: %d \n", angle);
		printf(" image path: %s \n", img_path);
		printf("===================== \n");
	}
	/*
	 *	Check arguments
	 */
	//Check if the file can be open
	if (file = fopen(img_path, "r")) {
        fclose(file);
    } else {
		program_abort(NULL, " :Failed to open the image file.\n");
	}
	// Validate angle argument 
	if (angle != 0 && angle != 90 && angle != 135 && angle != 45) {
		program_abort(NULL, "Invalid angle provided, angle must be one of 0, 45, 90, 135 \n");
	}

	/*
	 *	Load Image
	 */
	int width, height, bpp;
	unsigned char *grey_img;
	int *GLCM;

	// Load in grayscale image 
    grey_img = stbi_load(img_path, &width, &height, &bpp, 1);
	int GLCM_size = num_gray_levels(grey_img, width*height);
	
	if (height % num_procs != 0) {
		program_abort(NULL, "Image size is not evenly divisible by the number of processes \n");
	}
	
	/*
	 *	Process Image
	 */
	//Init Result GLCM
	if (rank == 0) {
		if ((GLCM = (int*) malloc(sizeof(int) * GLCM_size * GLCM_size)) == NULL) {
			program_abort(NULL, "Out of memory! - 1");
		}
		init_arr_zero(GLCM, GLCM_size * GLCM_size);
	}
	
	//Execute
	double start_time;
	start_time = MPI_Wtime();

	if (strcmp(implementation, "sequential") == 0){
		if (rank == 0) printf("running implementation: sequential \n");
		sequential(
			(mpi_data) {.num_procs=num_procs, .rank=rank},
			(img_data) {.arr=grey_img, .height=height, .width=width},
			(gclm_data) {.angle=angle, .arr=GLCM, .dist=distance, .size=GLCM_size}
		);
	} else {
		if (rank == 0) printf("running implementation: sync_vertical_split \n");
		sync_vertical_split(
			(mpi_data) {.num_procs=num_procs, .rank=rank},
			(img_data) {.arr=grey_img, .height=height, .width=width},
			(gclm_data) {.angle=angle, .arr=GLCM, .dist=distance, .size=GLCM_size}
		);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	//Aggregate Time
	if (rank == 0) {
		fprintf(stdout,"image size: %d | number of processes: %d |  time: %.3lf seconds\n",
			width, 
			num_procs,
			MPI_Wtime() - start_time);
	}
	
	//Check Results
	check_gclm(
		(mpi_data){.rank=rank},
		(img_data){.arr=grey_img, .height=height, .width=width},
		(gclm_data){.angle=angle, .arr=GLCM, .dist=distance, .size=GLCM_size}
	);
	
  	// Clean-up
	if (rank == 0)
		free(GLCM);
	stbi_image_free(grey_img); 

    MPI_Finalize();
}
