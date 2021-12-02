#include "haralick_imp.h"


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
    Implementation
*/

