#include "main.h"

// Define global variables declared in header file
int dev_start;
double tStart;
double tEnd;
int tt;

int main(void) 
{
  // Read input file
  main_read_input();

  // Read and sort output directory for finding files within our time limits
  init_part_files();
  init_flow_files();

  // Create output directory
  create_output_dir();

  // Get number of particles
  nparts = cgns_read_nparts();

  // Initialize domain and flow arrays
  domain_init(); 

  // Initialize partstruct and flow vars
  parts_init();
  flow_init();

  // Messy hack for taking advantage of CUDA_VISIBLE_DEVICES in SLURM
  dev_start = read_devices();

  // Allocate device memory
  cuda_dev_malloc();
  cuda_dom_push(); 
  //cuda_part_push();
  //cuda_part_pull();

  // Calculate tetrad stats
  for (tt = 0; tt < nFiles; tt++) {
    // Read in new data and push to device
    cgns_fill_flow();
    cuda_flow_push();

    // Calculate phase average
    cuda_phase_averaged_vel();
  }

  // write phaseAvereaged to file
  write_part_data();
   
  // Free and exit
  free_vars();
  cuda_dev_free();
  return EXIT_SUCCESS;
}


