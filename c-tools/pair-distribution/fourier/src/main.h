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
#define ANALYSIS "fourier"
#define DATA_DIR "data"
#define INPUT_DIR "input"
#define OUTPUT_DIR "output"
#define CONFIG_FILE "part-pair-fourier.config"

#define PERIODIC 0
#define DIRICHLET 1
#define NEUMANN 2

#define PI 3.1415926535897932385

/**** VARIABLES ****/
// Declare global variables
extern double tStart;       // start time
extern double tEnd;         // end time
extern double R0;           // characteristic length
extern int legendreOrder;   // order of legendre reconstruction
extern int fourierOrder;    // order of fourier reconstruction
extern int printAvgFlag;    // 1 prints avg, 0 prints timesteps
extern int printCoeffsFlag;    // 1 prints avg, 0 prints timesteps
extern int nPointsR;        // npoints to reconstruct at
extern int nPointsTh;
extern int tt;              // time iterator

// geometry
extern double r_ab;
extern double mu_ab;
extern double p_ell;


#endif
