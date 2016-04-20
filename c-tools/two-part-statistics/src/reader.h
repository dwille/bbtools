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
  double u;
  double v;
  double w;
  int bin;
  int flipCountX;
  int flipCountY;
  int flipCountZ;
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

// pair structure
typedef struct pair_struct {
  int N1;
  int N2;
  int N1_initFlip_X;
  int N1_initFlip_Y;
  int N1_initFlip_Z;
  int N2_initFlip_X;
  int N2_initFlip_Y;
  int N2_initFlip_Z;
} pair_struct;

/**** VARIABLES ****/
// output directory partFiles and partFileTime
extern int nFiles;
extern char **partFiles;
extern double *partFileTime;
extern double *simTime;
extern int *fileMap;
extern int sigFigPre;
extern int sigFigPost;
extern char runDir[256];

// nparts
extern int nparts;
extern double rmax;

// binLength
extern double binLength;

// nMax -- maximum number of particles in 27 bins
extern int nMax;

// host and dev part_struct parts;
extern part_struct *parts;
extern part_struct *_parts;
extern part_struct *_partsPrev;

// host and dev dom_struct doms
extern dom_struct dom;
extern dom_struct *_dom;

// host and dev dom_struct binDoms
extern dom_struct binDom;
extern dom_struct *_binDom;

// host and dev bc struct bcs
extern BC bc;
extern BC *_bc;

// host and dev pair_structures
extern pair_struct *pairs;
extern pair_struct *_pairs;

// Host and Dev pair data
extern double *rSep;        // pair separation distance
extern double *_rSep;

/**** FUNCTIONS ****/
// read main.config input file
void read_input(void);

// read and sort output directory
void init_input_files(void);

// Merge sort
void merge_sort(double *A, int n, int *A2);
void merge(double *A, int n, int m, int *A2);

// create output dir
void create_output(void);
void init_output(void);

// read CUDA_VISIBLE_DEVICES
int read_devices(void);

// Read nparts
int cgns_read_nparts(void);

// initialize part_struct
void parts_init(void);

// initialize dom_struct binDom
void domain_init(void);

// show binDom and bc structures
void show_domain(void);

// Allocate memory for pair arrays
void alloc_pair_arrays(void);

// Read part_struct data
void cgns_fill_part_struct(void);

// get sigfigs of last file 
void get_sigfigs(void);

// write nodes
void write_info(void);

// write each timestep
void write_timestep(void);

// Free parts
void free_parts(void);

#endif
