# Parallel Implementation of Gray Level Co-occurrence Matrix Construction using MPI
#### ICS632 Final Project Submission  
Aditi Jaiswal, Arianna Bunnell, and Sorapong Khongnawang
## Docker Container
SMPI Simulations were run via a Docker container. Instructions on how to install and run the SMPI via Docker can be found at https://simgrid.github.io/SMPI_CourseWare/topic_getting_started/installing_simgrid/ 

## Instructions to Run 
The program can be run either on a physical cluster with MPI or using SMPI. To use MPI, replace the SMPI commands below with the MPI equivalent. To run with SMPI, first verify that you are in the Docker container which contains your installation of SMPI. 
1. ``` smpicc main.c haralick_imp.c haralick_imp.h stb_image.h misc.c misc.h -o <name> -lm -Ofast ```
2. smpirun -np <X> -hostfile ./config/hostfile_64.txt -platform ./config/cluster_crossbar_64.xml ./<name> ./data/<image file> <angle> -d <distance> -i <implementation>

```-d``` and ```-i``` are optional parameters. The default value for distance is 1, and any integer up to the size of the image can be provided. The default value for implementation is row-based, but ```sequential``` and ```tiling``` are also valid parameter values. 

## File Structure 
- ```config``` contains files needed to run the code using SMPI 
   - ```calibrate_flops.py``` contains code needed to calibrate the SMPI simulation, taken from https://simgrid.github.io/SMPI_CourseWare/topic_understanding_performance/matrixmultiplication/
   - ```cluster_crossbar_64.txt``` contains code to simulate a 100-processor cluster 
   - ```hostfile_64.txt``` contains code to simulate 100 hosts
- ```data``` contains examples of random png image files which were used to obtain results in the writeup. Images were generated with a simple Python script using Pillow and NumPy
- ```docs```contains relevant documentation for the project 
   - ```reference``` contains papers which served as the inspiration for the project 
   - The base folder contains the LaTeX code to create the final writeup 