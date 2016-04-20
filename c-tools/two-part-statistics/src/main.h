#ifndef _MAIN_H
#define _MAIN_H

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
#define DATA_OUT_DIR "data-two-part-stats"
#define MAX_THREADS_1D 128
#define MAX_THREADS_DIM 16

#define PERIODIC 0
#define DIRICHLET 1
#define NEUMANN 2

#define ALPHA_MAX 0.74048
#define PI 3.1415926535897932385
#define nDim 3
#define nDim2 nDim*nDim
#define NPOINTS 2

/**** STRUCTURES ****/
/**** VARIABLES ****/
// Declare global variables
extern double tStart;       // start time
extern double tEnd;         // end time
extern double R0_a;         // pair init length
extern double eps_a;        // Elmer Fudge factor
extern int multRuns;        // flag for if multiple runs will be performed

extern int dev_start;       // cuda device number

extern int nPairs;

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
void cuda_dev_malloc(void);

// bin and find pairs
void cuda_find_pairs(void);

// allocate memory for pairs and pull from device
void cuda_pair_malloc(void);

// Deal with periodicity
void cuda_periodic_flip(void);

// calculate pair stats
void cuda_pair_stats(void);

// save previous timestep part positions
void cuda_save_parts_prev(void);

// calculate higher order statistics
void cuda_higher_moments(double *_array, int length, double *moments);

// init cuda threads blocks
void cuda_init_threads_blocks(int *threads, int *blocks, int length);

// wrap cast to pointer and reduce in one function
double cuda_sum(double *_array, int N);

// free cuda memory
void cuda_dev_free(void);

/* MORE EXTERNAL VARIABLES */
// Shape Moments
// m[0] -- mean, m[1] -- sdev
// m[2] -- skew, m[3] -- kurt
extern double m_rSep[4];

// Alignment moments 
extern double m_g1_s1[4];   // Alignment of principal shape axes with initial 
extern double *_g1_s1;   // Alignment of principal shape axes with initial

#endif
