#ifndef _CGNS_READER_H
#define _CGNS_READER_H

#include "main.h"
#include <cgnslib.h>
#include <dirent.h>

/**** STRUCTURES ****/
// grid_info
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

/**** VARIABLES ****/
// File Variables
extern int nFiles;
extern char **files;
extern double *fileTime;
extern int *fileMap;

// Input data array
extern double *volume_fraction;

// Output data
extern double *moving_avg;
extern double *stationary_avg;

// dom_struct doms
extern dom_struct dom;

/**** FUNCTIONS ****/
// set up directory structure
void directory_init(int argc, char *argv[]);
 
// set up batch run
void batch_setup(int argc, char *argv[]);

// read input file
void main_read_input(void);

// read and sort files
void init_input_files(void);

// Merge sort
void merge_sort(double *A, int n, int *A2);
void merge(double *A, int n, int m, int *A2);

// Create directory for output data, init output files
void create_output(void);

// initialize dom_struct binDom
void domain_init(void);

// Read flow_struct data
void cgns_fill_input(void);

// show binDom and bc structures
void show_domain(void);

// get sigfigs
void get_sigfigs(void);

// write reconsructed data
void write_field(double curr_location);

// Free parts
void free_vars(void);

#endif
