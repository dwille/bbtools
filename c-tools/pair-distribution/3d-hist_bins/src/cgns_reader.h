#ifndef _CGNS_READER_H
#define _CGNS_READER_H

#include "main.h"
#include <cgnslib.h>
#include <dirent.h>

/**** STRUCTURES ****/
// part_struct
typedef struct part_struct {
  double r;
  double x;
  double y;
  double z;
} part_struct;

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

// boundary condition structure
typedef struct BC {
  int pW;
  int pE;
  int pS;
  int pN;
  int pB;
  int pT;
} BC;

/**** VARIABLES ****/
// File Variables
extern int nFiles;
extern char **flowFiles;
extern double *flowFileTime;
extern int *flowFileMap;

// Number of Particles
extern int nparts;
extern double meanR;

// Reconstruction params
extern int nBinsTh;
extern int nBinsR;
extern double cutoff;
extern int gtr;
extern int less;
extern double dr;
extern double dth;
extern double *evalTh;
extern double *evalR;
extern double *gHist;
extern double *volume_fraction;

// part_struct parts;
extern part_struct *parts;

// host and dev dom_struct doms
extern dom_struct dom;

// host and dev bc struct bcs
extern BC bc;

/**** FUNCTIONS ****/
// set up directory structure
void directory_init(int argc, char *argv[]);
 
// read input file
void main_read_input(void);

// read and sort part/flow files
void init_files(void);

// Merge sort
void merge_sort(double *A, int n, int *A2);
void merge(double *A, int n, int m, int *A2);

// Create directory for output data, init output files
void create_output(void);

// initialize parts, domain, and histogram
void parts_init(void);
void domain_init(void);
void histogram_init(void);

// Read part_struct data
void cgns_fill(void);

// show binDom and bc structures
void show_domain(void);

// write
void write_coeffs(int in, int ll, int nn);
void write_reconstruct(void);
void write_avg_reconstruct(void);

// Free parts
void free_vars(void);

#endif
