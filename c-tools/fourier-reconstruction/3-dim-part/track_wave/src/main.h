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

// #Defines
#define FILE_NAME_SIZE 256
#define CHAR_BUF_SIZE 256

// take care of batch job submission
extern char *SIM_ROOT_DIR;   // Simulation root dir
extern char *ANALYSIS_DIR;   // Analysis dir

#define FREC_DIR "fourier-reconstruction"
#define ANALYSIS "3-dim-part"
#define DATA_DIR "track_data"
#define INPUT_DIR "input"
#define OUTPUT_DIR "output"
#define CONFIG_FILE "track_wave.config"

#define PI 3.1415926535897932385

/**** VARIABLES ****/
// Declare global variables
extern double tStart;       // start time
extern double tEnd;         // end time
extern int tt;              // time iterator

// Averaging variables
extern double initial_location;   // location of where to start
extern double averaging_distance;   // TOTAL averaging distance
extern double radius;             // so we don't have to read the damn thing in
extern double wavespeed;          // ditto

#endif
