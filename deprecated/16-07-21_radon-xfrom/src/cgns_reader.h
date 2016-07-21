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
// host and dev dom_struct doms
extern dom_struct dom;

/**** FUNCTIONS ****/
// read input file
void main_read_input(void);

// Create directory for output data, init output files
void create_output(void);

// initialize dom_struct binDom
void domain_init(void);

// show binDom and bc structures
void show_domain(void);

// write
void write_radon(void);

// Free parts
void free_vars(void);

#endif
