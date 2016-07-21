#ifndef _READER_H
#define _READER_H

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
// output directory flowFiles and flowFileTime
extern int nFiles;
extern char **flowFiles;
extern double *flowFileTime;
extern int *fileMap;
extern int sigFigPre;
extern int sigFigPost;
extern double *uf;
extern double *_uf;
extern double *vf;
extern double *_vf;
extern double *wf;
extern double *_wf;
extern int *phase;
extern int *_phase;
extern double *phaseAvgUf;
extern double *phaseAvgVf;
extern double *phaseAvgWf;

// host and dev dom_struct doms
extern dom_struct dom;
extern dom_struct *_dom;

/**** FUNCTIONS ****/
// set up directory structure
void directory_init(int argc, char *argv[]);
 
// read tetrad.config input file
void main_read_input(void);

// read and sort output directory
void init_input_files(void);

// Merge sort
void merge_sort(double *A, int n, int *A2);
void merge(double *A, int n, int m, int *A2);

// create output dir
void create_output_dir(void);

// allocate phase averaged results
void alloc_result(void);

// read CUDA_VISIBLE_DEVICES
int read_devices(void);

// initialize dom_struct binDom
void domain_init(void);

// Read flow_struct data
void cgns_fill_flow(void);

// show binDom and bc structures
void show_domain(void);

// get sigfigs of last file 
//void get_sigfigs(void);

// write phaseAveraged
void write_averaged(void);

//// write each timestep
//void write_timestep(void);

// Free flow vars
void free_flow(void);

#endif
