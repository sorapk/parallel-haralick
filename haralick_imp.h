#ifndef HEADER_FILE
#define HEADER_FILE


/*
    Helpers
*/
unsigned char get_first_neighbor(unsigned char image[], int x, int y, int angle, int distance, int width, int height);
unsigned char get_second_neighbor(unsigned char image[], int x, int y, int angle, int distance, int width, int height);
int num_gray_levels(unsigned char image[], int size);

/*
    Implementation
*/

#endif
