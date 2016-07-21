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
extern char **flowFiles;
extern double *flowFileTime;
extern int *flowFileMap;

// Fluid velocities
extern double *uf;
extern double *vf;
extern double *wf;
extern int *phase;

extern fftw_complex *uf_k;
extern fftw_complex *vf_k;
extern fftw_complex *wf_k;

// Reconstruct variables
extern double *evalZ;
extern double *wf_rec_Re;
extern double *wf_rec_Im;

// dev dom_struct doms
extern dom_struct dom;

/**** FUNCTIONS ****/
// set up directory structure
void directory_init(int argc, char *argv[]);
 
// read input file
void main_read_input(void);

// read and sort flow files
void init_flow_files(void);

// Merge sort
void merge_sort(double *A, int n, int *A2);
void merge(double *A, int n, int m, int *A2);

// Create directory for output data, init output files
void create_output(void);

// initialize dom_struct binDom
void domain_init(void);

// Read flow_struct data
void cgns_fill_flow(void);

// Test
void test_fill(void);
void test_out(void);

// show binDom and bc structures
void show_domain(void);

// write
void write_coeffs(int in);
void write_reconstruct(void);

// Free parts
void free_vars(void);

#endif
