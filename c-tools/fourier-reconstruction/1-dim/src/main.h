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

#include "cgns_reader.h"
#include "reconstruct.h"

// #Defines
#define FILE_NAME_SIZE 256
#define CHAR_BUF_SIZE 256

// take care of batch job submission
#ifdef BATCH
  extern char *ROOT_DIR;
  extern char *SIM_ROOT_DIR;
#else
  #define SIM_ROOT_DIR ".."   // sim
  #define ROOT_DIR "."        // analysis
#endif

#define DATA_DIR "data"
#define INPUT_DIR "input"
#define OUTPUT_DIR "output"
#define CONFIG_FILE "f-rec-1D.config"

#define PERIODIC 0
#define DIRICHLET 1
#define NEUMANN 2

#define ALPHA_MAX 0.74048
#define PI 3.1415926535897932385
#define nDim 3
#define nDim2 nDim*nDim

/**** VARIABLES ****/
// Declare global variables
extern double tStart;       // start time
extern double tEnd;         // end time
extern int order;           // order of fourier expansion
extern int coeffsOut;       // output coeffs or not
extern int npoints;         // number of points to evaluate at
extern int tt;              // time iterator

// constants
extern double pref;         // prefactor constant


#endif
