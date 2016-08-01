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

extern int sigFigPre;
extern int sigFigPost;

// Solution variables
extern int nBins;
extern int nparts;
extern double *volume_fraction;
extern double *part_w;
extern int *histogram_vf;
extern int *histogram_wp;
extern int *bihistogram_vf_wp;

// struct init
extern dom_struct dom;

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

// get sigfigs
void get_sigfigs(void);

// write reconsructed data
void write_field(void);

// Free parts
void free_vars(void);

#endif
