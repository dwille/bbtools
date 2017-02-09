#ifndef _READER_H
#define _READER_H

#include "main.h"
#include <cgnslib.h>
#include <dirent.h>

/**** VARIABLES ****/
// File Variables
extern int nFiles;
extern char **partFiles;
extern double *partFileTime;
extern int *partFileMap;

// Number of Particles
extern int nparts;

/**** FUNCTIONS ****/
// set up directory structure
void directory_init(int argc, char *argv[]);

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

// Read part_struct data
int cgns_fill_parts(int tt);

// show binDom and bc structures
void show_domain(void);

// write
void write_mean(int tt, int sum);

// Free parts
void free_vars(void);

#endif
