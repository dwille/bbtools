#ifndef _READER_H
#define _READER_H

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
extern char **partFiles;
extern double *partFileTime;
extern int *partFileMap;

// Number of Particles
extern int nparts;
extern double meanR;


// host and dev part_struct parts;
extern part_struct *parts;

// host and dev dom_struct doms
extern dom_struct dom;

// host and dev bc struct bcs
extern BC bc;

/**** FUNCTIONS ****/
// read input file
void main_read_input(void);

// read and sort part/flow files
void init_part_files(void);

// Merge sort
void merge_sort(double *A, int n, int *A2);
void merge(double *A, int n, int m, int *A2);

// Create directory for output data, init output files
void create_output(void);

// Read nparts
int cgns_read_nparts(void);

// initialize part_struct and flow vars
void parts_init(void);

// initialize dom_struct binDom
void domain_init(void);

// Read part_struct data
void cgns_fill_parts(void);

// show binDom and bc structures
void show_domain(void);

// write
void get_sigfigs(void);
void write_coeffs(int in);
void write_reconstruct(void);

// Free parts
void free_vars(void);

#endif
