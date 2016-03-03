#include "tetrad_init.h"

// Define global variables declared in header file
int dev_start;
double tStart;
double tEnd;
double R0_a;
double eps_a;
double varCut;
double shapeCut;
int findTets;
int nUnique;
int *uniqueNodes;
int *_uniqueNodes;

int main(void) 
{
  // Read tetrad input file
  tetrad_read_input();

  // Read and sort output directory for finding files within our time limits
  init_input_files();

  // Create output directory
  create_output_dir();

  // Messy hack for taking advantage of CUDA_VISIBLE_DEVICES in SLURM
  dev_start = read_devices();

  // Get number of particles
  nparts = cgns_read_nparts();

  // Initialize part_struct and fill for first timestep
  parts_init(nparts);

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

  // Allocate memory for tetrads
  cuda_tetrad_malloc();
  cuda_tetrad_init();

  // Calculate tetrad stats
  for (int tt = 0; tt < nFiles; tt++) {
    // Read in new data and push to device
    cgns_fill_part_struct(nparts, tt);
    cuda_part_push();

    // Calculate tetrads
    cuda_tetrad_stats();

    //  Write to file
    write_timestep(tt);
  }

  // write uniqueNodes to file
  write_nodes();
   
  // Free and exit
  free_parts();
  cuda_dev_free();
  return EXIT_SUCCESS;
}


