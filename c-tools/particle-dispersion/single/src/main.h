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

#include "reader.h"

#define FILE_NAME_SIZE 256
#define CHAR_BUF_SIZE 256

#define MAX_THREADS_1D 128
#define MAX_THREADS_DIM 16

/*******************/
/**** VARIABLES ****/
/*******************/
extern double t_start;       // start time
extern double t_end;         // end time
extern double *sim_time;     // value of timesteps
extern int nparts;           // number of particles
extern int tt;               // time series iterator

typedef struct part_struct {
  double r;
  double x;
  double y;
  double z;
  double u;
  double v;
  double w;
  int flip_count_i;
  int flip_count_j;
  int flip_count_k;
} part_struct;

extern part_struct *parts;
extern part_struct *_parts;

typedef struct grid_info {
  int is;
  int js;
  int ks;
  int ie;
  int je;
  int ke;
  int in;
  int jn;
  int kn;
  int s1;
  int s2;
  int s3;
} grid_info;

// dom_struct
typedef struct dom_struct {
  grid_info Gcc;
  double xs;
  double ys;
  double zs;
  double xe;
  double ye; 
  double ze;
  double xl;
  double yl;
  double zl;
  int xn;
  int yn;
  int zn;
  double dx;
  double dy;
  double dz;
} dom_struct;

extern dom_struct dom;
extern dom_struct *_dom;

extern double *r2_total; // mean squared separation time series
extern double *r2_horiz; // mean squared separation time series (horizontal)
extern double *r2_verti; // mean squared separation time series (vertical)

extern double *_x0;      // particle x position at t = t0
extern double *_y0;      // particle y position at t = t0
extern double *_z0;      // particle z position at t = t0
extern double *_x;       // particle x position at t = t
extern double *_y;       // particle y position at t = t
extern double *_z;       // particle z position at t = t
extern double *_r2_total; // mean sqrd sep -- per particle (total)
extern double *_r2_horiz; // mean sqrd sep -- per particle (horiz)
extern double *_r2_verti; // mean sqrd sep -- per particle (verti)


/*******************/
/**** FUNCTIONS ****/
/*******************/
void cuda_set_device(void);
void cuda_dev_malloc(void);
void cuda_dom_push(void);
void cuda_part_push(void);
void cuda_fill_initial_positions(void);
void cuda_periodic_flip(void);
void cuda_find_separations(void);
void cuda_free(void);

#endif // _MAIN_H
