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

// set up batch job submission
extern char *SIM_ROOT_DIR;
extern char *ANALYSIS_DIR;

#define PHASE_DIR "phase-averages"
#define ANALYSIS "collision"
#define DATA_DIR "data"
#define DATA_FILE "collision_stats"
#define INPUT_DIR "input"
#define OUTPUT_DIR "output"
#define CONFIG_FILE "collision.config"

/**** VARIABLES ****/
// Declare global variables
extern double tStart;       // start time
extern double tEnd;         // end time

#endif
