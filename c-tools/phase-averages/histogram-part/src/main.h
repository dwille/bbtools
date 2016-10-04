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

/* Directory structure
 * SIM_ROOT_DIR/
   * input
   * output
   * analysis/FREC_DIR/ANALYSIS == ANALYSIS_DIR
     * DATA_DIR
        * OUTPUT_FILE
        * INFO_FILE
     * CONFIG_FILE 
*/
extern char *SIM_ROOT_DIR;   // Simulation root dir
#define INPUT_DIR "input"
#define OUTPUT_DIR "output"
#define PHASE_DIR "phase-averages"
#define ANALYSIS "histograms"
extern char *ANALYSIS_DIR;   // Analysis dir
#define DATA_DIR "data"
#define INFO_FILE "info"
#define CONFIG_FILE "part-hist.config"

#define CGNS_FILE "part-" // <*.cgns>

#define PI 3.1415926535897932385

/**** VARIABLES ****/
// Declare global variables
extern double tStart;       // start time
extern double tEnd;         // end time
extern int tt;              // time iterator

/**** FUNCTIONS ****/
void minmax(void);
void bin_init(void);

#endif
