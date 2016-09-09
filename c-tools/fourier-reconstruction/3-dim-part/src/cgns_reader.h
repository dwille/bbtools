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
  double u;
  double v;
  double w;
  double ox;
  double oy;
  double oz;
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

/**** VARIABLES ****/
// File Variables
extern int nFiles;
extern char **flowFiles;
extern double *flowFileTime;
extern int *flowFileMap;

// Solution variables
extern int nparts;
// extern double *x;
// extern double *y;
// extern double *z;
// extern double *up;
// extern double *vp;
// extern double *wp;
// extern double *oxp;
// extern double *oyp;
// extern double *ozp;

extern int *phase;
extern double *uf;
extern double *vf;
extern double *wf;

// FFTW variables
extern fftw_complex *chi;
extern fftw_complex *up_field;             // for taking FFT, size(dom.Gcc.s3)
extern fftw_complex *vp_field;
extern fftw_complex *wp_field;
extern fftw_complex *ke_field;
// extern fftw_complex *uop_field;
// extern fftw_complex *vop_field;
// extern fftw_complex *wop_field;

// struct init
extern part_struct *parts;
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

// Read data
// void cgns_fill_part(void);
void cgns_fill_flow(void);

// Test
void test_fill(void);
void test_out(void);

// show binDom and bc structures
void show_domain(void);

// get sigfigs
void get_sigfigs(void);

// write reconsructed data
void cgns_write_field(void);

// Free parts
void free_vars(void);

#endif
