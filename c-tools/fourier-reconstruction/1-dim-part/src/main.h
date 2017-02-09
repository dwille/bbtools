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

#define FREC_DIR "fourier-reconstruction"
#define ANALYSIS "1-dim-part"
#define DATA_DIR "data"
#define INPUT_DIR "input"
#define OUTPUT_DIR "output"
#define CONFIG_FILE "f-rec-1D-part.config"

#define PI 3.1415926535897932385

/**** VARIABLES ****/
// Declare global variables
extern double tStart;       // start time
extern double tEnd;         // end time
extern int order_s;         // starting order of fourier expansion
extern int order_e;         // ending order of fourier expansion
extern int n_order;         // number of orders
extern int coeffsOut;       // output coeffs or not
extern int npoints;         // number of points to evaluate at
extern int tt;              // time iterator

// constants
extern double pref;         // prefactor constant

#endif
