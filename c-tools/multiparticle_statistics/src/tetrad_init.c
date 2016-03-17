#include "tetrad_init.h"

// Define global variables declared in header file
int dev_start;
double tStart;
double tEnd;
double R0_a;
double eps_a;
double varCutLow;
double varCutHigh;
double shapeCutLow;
double shapeCutHigh;
int findTets;
int nRegular;
int tt;

double meanR2;
double meanVar;
double meanShape;
double stdR2;
double stdVar;
double stdShape;

double mean_g1_s1;   // Alignment of principal shape axes with initial strain axes
double mean_g2_s1;
double mean_g3_s1;
double mean_g1_s2;
double mean_g2_s2;
double mean_g3_s2;
double mean_g1_s3;
double mean_g2_s3;
double mean_g3_s3;
double mean_g1_z;    // Alignment of shape, strain, vorticity with gravity
double mean_g2_z;
double mean_g3_z;
double mean_s1_z;
double mean_s2_z;
double mean_s3_z;
double mean_w_z;
double mean_w_g1;    // Alignment of vorticity with initial shape, strain axes
double mean_w_g2;
double mean_w_g3;
double mean_w_s1;
double mean_w_s2;
double mean_w_s3;
double mean_vortMag;

double *_g1_s1;   // Alignment of principal shape axes with initial strain axes
double *_g2_s1;
double *_g3_s1;
double *_g1_s2;
double *_g2_s2;
double *_g3_s2;
double *_g1_s3;
double *_g2_s3;
double *_g3_s3;
double *_g1_z;    // Alignment of shape, strain, vorticity with gravity
double *_g2_z;
double *_g3_z;
double *_s1_z;
double *_s2_z;
double *_s3_z;
double *_w_z;
double *_w_g1;    // Alignment of vorticity with initial shape, strain axes
double *_w_g2;
double *_w_g3;
double *_w_s1;
double *_w_s2;
double *_w_s3;
double *_vortMag;

int main(void) 
{
  // Read tetrad input file
  tetrad_read_input();

  // Read and sort output directory for finding files within our time limits
  // and also initialize time array
  init_input_files();

  // Create output directory and init mean output file
  create_output_dir();
  init_stat_output();

  // Messy hack for taking advantage of CUDA_VISIBLE_DEVICES in SLURM
  dev_start = read_devices();

  // Get number of particles
  nparts = cgns_read_nparts();

  // Initialize part_struct and fill for first timestep
  parts_init();

  // Initialize domain, binDom and fill
  domain_init(); 

  // Allocate device memory
  cuda_dev_malloc();
  cuda_dom_push(); 
  cuda_part_push(); 

  // Bin particles, find tetrads
  if (findTets == 1) {
    cuda_find_tetrads();
  } else {
    printf("  uniqueNodes already exists, reading in file.\n");
    printf("  Except this feature doesn't exist yet.\n");
    printf("  See so-18737117, 12733105\n");
    exit(EXIT_FAILURE);
  }

  // Allocate host memory for tetrad statistics
  alloc_tetrad_arrays();

  // Allocate device memory for tetrads
  cuda_tetrad_malloc();

  // get sig figs for writing
  get_sigfigs();

  // Calculate tetrad stats
  for (tt = 0; tt < nFiles; tt++) {
    // Read in new data and push to device
    cgns_fill_part_struct();
    cuda_part_push();

    // Calculate tetrads
    cuda_tetrad_stats();
    
    //  Write to file
    write_timestep();
  }

  // write uniqueNodes to file
  write_nodes();
   
  // Free and exit
  free_parts();
  cuda_dev_free();
  return EXIT_SUCCESS;
}


