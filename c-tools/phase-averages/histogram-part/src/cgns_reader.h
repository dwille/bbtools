#ifndef _CGNS_READER_H
#define _CGNS_READER_H

#include "main.h"
#include <cgnslib.h>
#include <dirent.h>

/**** STRUCTURES ****/
//part_sttruct
typedef struct part_struct {
  double r;
  double x;
  double y;
  double z;
} part_struct;
//extern part_struct *parts;

// bin_params
typedef struct bin_params {
  double min;
  double max;
  double mean;
  double dBin;
} bin_params;

// bin_struct
typedef struct bin_struct {
  bin_params U;
  bin_params ke_trans;
  bin_params ke_rot;
  bin_params ke;
  bin_params T_z;
  bin_params T_perp;
  bin_params T;
  bin_params ux;
  bin_params vy;
  bin_params wz;
} bin_struct;
extern bin_struct bins;

/**** VARIABLES ****/
// File Variables
extern int nFiles;
extern char **files;
extern double *fileTime;
extern int *fileMap;

// Part stuff
extern int nparts;
extern double mass;
extern double meanR;
extern double meanI;
extern double *up;
extern double *vp;
extern double *wp;
extern double *ox;
extern double *oy;
extern double *oz;
extern double *U;
extern double *ke;
extern double *ke_rot;
extern double *ke_trans;
extern double *T;
extern double *T_perp;
extern double *T_z;

// Solution variables
extern int nBins;
extern int *hist_U;
extern int *hist_ke;
extern int *hist_ke_rot;
extern int *hist_ke_trans;
extern int *hist_T;
extern int *hist_T_perp;
extern int *hist_T_z;
extern int *hist_ux;
extern int *hist_vy;
extern int *hist_wz;

/**** FUNCTIONS ****/
// set up directory structure
void directory_init(int argc, char *argv[]);
 
// read input file
void main_read_input(void);

// read and sort flow files
void init_files(void);

// Merge sort
void merge_sort(double *A, int n, int *A2);
void merge(double *A, int n, int m, int *A2);

// Create directory for output data, init output files
void create_output(void);

// initialize dom_struct binDom
void domain_init(void);

// Read data
void cgns_fill(void);

// show domain
void show_domain(void);

// write reconsructed data
void write_field(void);

// Free parts
void free_vars(void);

#endif
