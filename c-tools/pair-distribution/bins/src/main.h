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

// Define File Structure
extern char *SIM_ROOT_DIR;      // Simulation root directory
extern char *ANALYSIS_DIR;      // Analysis directory

#define PAIR_DIR "part-pair-distribution"
#define ANALYSIS "bins"
#define DATA_DIR "data"
#define INPUT_DIR "input"
#define OUTPUT_DIR "output"
#define CONFIG_FILE "part-pair-bin.config"

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
extern double R0;           // characteristic length
extern int nPointsR;        // npoints to reconstruct at
extern int nPointsTh;
extern int tt;              // time iterator

extern double r_ab;
extern double th_ab;


#endif
