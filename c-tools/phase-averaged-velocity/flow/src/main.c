#include "main.h"

// Define global variables declared in header file
char *SIM_ROOT_DIR;
char *ANALYSIS_DIR;
int dev_start;
double tStart;
double tEnd;
int tt;

int main(int argc, char *argv[]) 
{
  // Set directory structure
  directory_init(argc, argv);


  // Read main input file
  main_read_input();

  // Read and sort output directory for finding files within our time limits
  init_input_files(); // >> nFiles

  // Create output directory
  create_output_dir();

  // Allocate result
  alloc_result();

  // Messy hack for taking advantage of CUDA_VISIBLE_DEVICES in SLURM
  dev_start = read_devices();

  // Initialize domain and flow arrays
  domain_init(); 

  // Allocate device memory
  cuda_dev_malloc();
  cuda_dom_push(); 

  // Calculate tetrad stats
  for (tt = 0; tt < nFiles; tt++) {
    // Read in new data and push to device
    cgns_fill_flow();
    cuda_flow_push();

    // Calculate phase average
    cuda_phase_averaged_vel();
  }

  // write phaseAvereaged to file
  write_averaged();
   
  // Free and exit
  free_flow();
  cuda_dev_free();
  return EXIT_SUCCESS;
  printf("Done!!\n");
}


