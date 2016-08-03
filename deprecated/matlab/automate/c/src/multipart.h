#ifndef _MULTIPART_H
#define _MULTIPART_H

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "reader.h"
#include "array.h"

// #Defines
#define FILE_NAME_SIZE 256
#define CHAR_BUF_SIZE 256
#define ROOT_DIR "."
#define INPUT_DIR "input"
#define OUTPUT_DIR "output"

// Declare global variables
extern double ts;         // starting time
extern double te;         // ending time
extern double dt;         // timestep
extern int tetrad_flag;   // tetrad flag
extern int triad_flag;    // triad flag
extern int volume_flag;   // volume flag
extern int rog_flag;      // radius of gyration flag
extern int shape_factor_flag;   // shape factor flag
extern int lambda_flag;   // lambda flag
extern int invari_flag;   // invariants flag

extern double *fileTime;  // vector contining times written to file name
extern int nparts;        // number of particles
#endif
