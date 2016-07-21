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
// Deal with fftw for Marcc, lucan
#include <fftw3.h>        // if lucan
// #include "fftw3.h"        // if marcc

#include "cgns_reader.h"


// #Defines
#define FILE_NAME_SIZE 256
#define CHAR_BUF_SIZE 256

// take care of batch job submission
extern char *SIM_ROOT_DIR;
extern char *ANALYSIS_DIR;

#define FREC_DIR "fourier-reconstruction"
#define ANALYSIS "1-dim-part"
#define DATA_DIR "data"
#define INPUT_DIR "input"
#define OUTPUT_DIR "output"
#define CONFIG_FILE "f-rec-1D-flow.config"

#define PI 3.1415926535897932385

/**** VARIABLES ****/
// Declare global variables
extern double tStart;       // start time
extern double tEnd;         // end time
extern int order;           // order of fourier expansion
extern int coeffsOut;       // output coeffs or not
extern int npoints;         // number of points to evaluate at
extern int tt;              // time iterator

extern fftw_plan pU;
extern fftw_plan pV;
extern fftw_plan pW;


#endif
