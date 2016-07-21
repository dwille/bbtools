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
#define SIM_ROOT_DIR ".."
#define ROOT_DIR "."
#define INPUT_DIR "input"
#define OUTPUT_DIR "output"
#define DATA_DIR "data-tetrads"
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
extern double EVarCutLow;    // eigval variance cutoff
extern double EVarCutHigh;
extern double shapeCutLow;  // tet shape cutoff
extern double shapeCutHigh;
extern int findTets;        // find tets or no (if it already exists)
extern int multRuns;        // flag for if multiple runs will be performed

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

// bin and find tetrads
// also, find unique and then regular tetrads
void cuda_find_tetrads(void);

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

/* MORE EXTERNAL VARIABLES */
// Shape Moments
// m[0] -- mean, m[1] -- sdev
// m[2] -- skew, m[3] -- kurt
extern double m_RoG[4];
extern double m_EVar[4];
extern double m_Shape[4];
extern double m_I1[4];
extern double m_I2[4];
extern double m_I3[4];
extern double m_S11[4];
extern double m_S22[4];
extern double m_S33[4];

// Alignment moments 
extern double m_g1_s1[4];   // Alignment of principal shape axes with initial 
extern double m_g2_s1[4];   // strain axes
extern double m_g3_s1[4];
extern double m_g1_s2[4];
extern double m_g2_s2[4];
extern double m_g3_s2[4];
extern double m_g1_s3[4];
extern double m_g2_s3[4];
extern double m_g3_s3[4];
extern double m_g1_z[4];    // Alignment of shape, strain, vorticity with 
extern double m_g2_z[4];    // gravity
extern double m_g3_z[4];
extern double m_s1_z[4];
extern double m_s2_z[4];
extern double m_s3_z[4];
extern double m_w_z[4];
extern double m_w_g1[4];    // Alignment of vorticity with initial shape, 
extern double m_w_g2[4];    // strain axes
extern double m_w_g3[4];
extern double m_w_s1[4];
extern double m_w_s2[4];
extern double m_w_s3[4];
extern double m_vortMag[4];

extern double *_g1_s1;   // Alignment of principal shape axes with initial
extern double *_g2_s1;   // strain axes
extern double *_g3_s1;
extern double *_g1_s2;
extern double *_g2_s2;
extern double *_g3_s2;
extern double *_g1_s3;
extern double *_g2_s3;
extern double *_g3_s3;
extern double *_g1_z;    // Alignment of shape, strain, vorticity with gravity
extern double *_g2_z;
extern double *_g3_z;
extern double *_s1_z;
extern double *_s2_z;
extern double *_s3_z;
extern double *_w_z;
extern double *_w_g1;    // Alignment of vorticity with initial shape, 
extern double *_w_g2;    // strain axes
extern double *_w_g3;
extern double *_w_s1;
extern double *_w_s2;
extern double *_w_s3;
extern double *_vortMag;

#endif
