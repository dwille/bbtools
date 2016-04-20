#include "main.h"

// Define global variables declared in header file
int dev_start;
double tStart;
double tEnd;
double R0_a;
double eps_a;
int multRuns;
int nPairs;
int tt;

int main(void) 
{
  // Read main input file
  read_input();

  // Read and sort output directory for finding files within our time limits
  // and also initialize time array
  init_input_files();

  // get sig figs for writing
  get_sigfigs();

  // Create output directory and init mean output file
  create_output();
  init_output();

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

  // Bin particles, find pairs
  cuda_find_pairs();

  // Allocate host memory for pair statistics
  alloc_pair_arrays();

  // Allocate device memory for pairs
  cuda_pair_malloc();

  // Calculate pair stats
//  for (tt = 0; tt < nFiles; tt++) {
//    // Read in new data and push to device
//    cgns_fill_part_struct();
//    cuda_part_push();
//
//    // Take care of periodicity
//    //  -- Tetrads that are intialized over periodic boundaries need an
//    //      extra flip -- this is determined in check tolerances, and
//    //      done in pair_geometry
//    //  -- All particles need to be checked at tt > 0 to see if they
//    //      went over a periodic boundary -- this is done in cuda_periodic_flip
//    if (tt > 0) {
//      cuda_periodic_flip();
//    }
//
//    // Calculate pairs
//    cuda_pair_stats();
//    
//    // Write to file
//    write_timestep();
//
//    // Save current partstruct data to previous timestep
//    cuda_save_parts_prev();
//  }

  // write regular nodes to file, as well as nTetrads and nTsteps
//  write_info();
   
  // Free and exit
  free_parts();
  cuda_dev_free();
  return EXIT_SUCCESS;
}

double m_rSep[4];
double m_g1_s1[4];  // Alignment of principal shape axes with initial 
double *_g1_s1;   // Alignment of principal shape axes with initial strain axes
