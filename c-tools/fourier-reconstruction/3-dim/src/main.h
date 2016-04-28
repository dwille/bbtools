#ifndef _MAIN_H
#define _MAIN_H

#include <float.h>
#include <math.h>
#include <complex.h>
#undef I
#define J _Complex_I
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
#define ROOT_DIR "."
#define INPUT_DIR "input"
#define OUTPUT_DIR "output"
#define DATA_DIR "data"
#define DATA_SUBDIR "reconstruct-3D"

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
extern int orderL;          // order of fourier expansion
extern int orderM;
extern int orderN;
extern int tt;              // time iterator

#endif
