#ifndef _TETRAD_INIT_H
#define _TETRAD_INIT_H

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "reader.h"

// #Defines
#define FILE_NAME_SIZE 256
#define CHAR_BUF_SIZE 256
#define ROOT_DIR "."
#define INPUT_DIR "input"
#define OUTPUT_DIR "output"
#define DATA_OUT_DIR "data-part-pair"
#define MAX_THREADS_1D 128
#define MAX_THREADS_DIM 16

#define PERIODIC 0
#define DIRICHLET 1
#define NEUMANN 2

#define ALPHA_MAX 0.74048
#define PI 3.1415926535897932385
#define nDim 3
#define nDim2 nDim*nDim

/**** STRUCTURES ****/
/**** VARIABLES ****/
// Declare global variables
extern double tStart;       // start time
extern double tEnd;         // end time
extern double R0_a;         // tetrad init length
extern double eps_a;        // Elmer Fudge factor
extern int orderL;          // order of lagrange polynomials to use
extern int orderN;          // order of fourier series to use

extern int dev_start;       // cuda device number

extern int nRegular;

extern int tt;

/**** FUNCTIONS ****/
// dom and part push
void cuda_dom_push(void);
void cuda_part_push(void);

// part_pull
void cuda_dev_pull(void);

// allocate device memory
//  - _parts
//  - _dom
//  - _binDom
//  - _neighborList
void cuda_dev_malloc(void);

// bin and find pairs
void cuda_find_pairs(void);

// allocate memory for tetrads and pull from device
void cuda_tetrad_malloc(void);

// Deal with periodicity
void cuda_periodic_flip(void);

// calculate tetrad stats
void cuda_tetrad_stats(void);

// save previous timestep part positions
void cuda_save_parts_prev(void);

// calculate higher order statistics
void cuda_higher_moments(double *_array, int length, double *moments);

// wrap cast to pointer and reduce in one function
double cuda_sum(double *_array, int N);

// free cuda memory
void cuda_dev_free(void);

#endif
