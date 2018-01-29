#ifndef _READER_H
#define _READER_H

#include "main.h"
#include <cgnslib.h>
#include <dirent.h>

/*******************/
/**** VARIABLES ****/
/*******************/
extern char *SIM_ROOT_DIR;
extern char *ANALYSIS_DIR;
extern int n_files;
extern char **part_files;
extern double *part_file_time;
extern int *file_map;

/*******************/
/**** FUNCTIONS ****/
/*******************/

// Initialize data directory
void directory_init(int argc, char *argv[]);

// Read input file
void read_input(void);

// Parse sim output directory for cgns files
void init_cgns_files(void);

// Merge sort
void merge_sort(double *A, int n, int *A2);
void merge(double *A, int n, int m, int *A2);

// Create output directory
void create_output(void);

// Initilize part_struct and fill for first timestep
void parts_init(void);

// Init domain
void domain_init(void);

// Alloc host memory for single part stats
void alloc_host_mem(void);

// Fill particle structure at new timestep
void cgns_fill_part_struct(void);

// Write results
void write_output(void);

// Free memory
void free_mem(void);

#endif // _READER_H
