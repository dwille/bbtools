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

// set up input if we need batch job submission
#ifdef BATCH
  extern char *ROOT_DIR;
  extern char *SIM_ROOT_DIR;
#else
  #define ROOT_DIR "."        // analysis
  #define SIM_ROOT_DIR ".."   // simulation
#endif
#define INPUT_DIR "input"
#define OUTPUT_DIR "output"
#define DATA_DIR "data"
#define CONFIG_FILE "flowvel.config"

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

extern int dev_start;       // cuda device number

extern int tt;

/**** FUNCTIONS ****/
// allocate device memory
//  - _dom
//  - _uf, _vf, _wf, _phase
void cuda_dev_malloc(void);

// dom and part push
void cuda_dom_push(void);
void cuda_flow_push(void);

// flow_pull
void cuda_flow_pull(void);

// calcualte phase average velcoity
void cuda_phase_averaged_vel(void);

// free cuda memory
void cuda_dev_free(void);

#endif
