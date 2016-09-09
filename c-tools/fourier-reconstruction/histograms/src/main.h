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
#define FREC_DIR "fourier-reconstruction"
#define ANALYSIS "histograms"
extern char *ANALYSIS_DIR;   // Analysis dir
#define DATA_DIR "data"
#define INFO_FILE "info"
#define CONFIG_FILE "f-rec-3D-hist.config"

#define CGNS_FILE "f-rec-3D-" // <*.cgns>

#define PI 3.1415926535897932385

/**** VARIABLES ****/
// Declare global variables
extern double tStart;       // start time
extern double tEnd;         // end time
extern int tt;              // time iterator

/* Volume Fraction */
extern double min_vf;
extern double mean_vf;
extern double max_vf;
extern double dBin_vf;
extern double binStart_vf;
extern double binEnd_vf;

/* Vertical Velocity */
extern double min_wp;
extern double mean_wp;
extern double max_wp;
extern double dBin_wp;
extern double binStart_wp;
extern double binEnd_wp;

/* Kinetic Energy */
extern double min_ke;
extern double mean_ke;
extern double max_ke;
extern double dBin_ke;
extern double binStart_ke;
extern double binEnd_ke;

/**** FUNCTIONS ****/
void minmax(double *min, double *max, double *array);
void bin_init(double min, double max, double *dBin, double *binStart, 
  double *binEnd);

#endif
