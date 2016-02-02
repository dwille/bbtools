#ifndef _TEMPLATE_H
#define _TEMPLATE_H

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

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
extern double *fileTime;  // vector contining times written to file name
extern int nparts;        // number of particles
#endif
