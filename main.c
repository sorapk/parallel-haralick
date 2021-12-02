#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

#define DEFAULT_DISTANCE 1
#define IMPLEMENTATION_SIZE 20
#define NUM_ARG 3
#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"


#include "haralick_imp.h"

static void program_abort(char *exec_name, char *message);
static void print_usage();

static void program_abort(char *exec_name, char *message) {
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
static void print_usage(char *exec_name) {
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


/*
 * main() function
 */
int main(int argc, char **argv)
{	
	int i;
	char implementation[IMPLEMENTATION_SIZE];

	FILE *file;
	int angle, num_procs;
	int distance = DEFAULT_DISTANCE;	
	char *img_path; 

	/*
	 * Init MPI ENV
	*/
	MPI_Status status;
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
		}else if (!strcmp(argv[i],"-i")) {
			if ((i+1 >= argc) || (sscanf(argv[i+1],"%s",implementation) != 1)) {
				program_abort(argv[0],"Invalid <implementation> argument\n");
			}
		}
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
	unsigned char *gray_img;
	int *GLCM, *partial_GLCM;
	int range;

	// Load in grayscale image 
    gray_img = stbi_load(img_path, &width, &height, &bpp, 1);
	int GLCM_size = num_gray_levels(gray_img, width*height);
	
	if ((partial_GLCM = (int*) malloc(sizeof(int) * GLCM_size * GLCM_size)) == NULL) {
    	program_abort(argv[0], "Out of memory!");
  	}
	range = height / num_procs;
	
	if (height % num_procs != 0) {
		program_abort(NULL, "Image size is not evenly divisible by the number of processes \n");
	}
	
	/*
	 *	Process Image
	 */
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
