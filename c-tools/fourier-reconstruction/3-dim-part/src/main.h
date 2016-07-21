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

// On lucan: (system lib)
#include <fftw3.h>
// On MARCC: (install from source)
// #include "fftw3.h"

#include "cgns_reader.h"

// #Defines
#define FILE_NAME_SIZE 256
#define CHAR_BUF_SIZE 256

// take care of batch job submission
extern char *SIM_ROOT_DIR;   // Simulation root dir
extern char *ANALYSIS_DIR;   // Analysis dir

#define FREC_DIR "fourier-reconstruction"
#define ANALYSIS "3-dim-part"
#define DATA_DIR "data"
#define INPUT_DIR "input"
#define OUTPUT_DIR "output"
#define CONFIG_FILE "f-rec-3D-part.config"

#define PI 3.1415926535897932385

/**** VARIABLES ****/
// Declare global variables
extern double tStart;       // start time
extern double tEnd;         // end time
extern int orderX;           // order of fourier expansion
extern int orderY;           // order of fourier expansion
extern int orderZ;           // order of fourier expansion
extern int tt;              // time iterator

extern fftw_plan chi2phi_k;
extern fftw_plan phi_k2phi;

#endif
