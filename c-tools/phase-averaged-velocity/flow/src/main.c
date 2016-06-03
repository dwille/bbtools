#include "main.h"

// Define global variables declared in header file
#ifdef BATCH
  char *ROOT_DIR;
  char *SIM_ROOT_DIR;
#endif
int dev_start;
double tStart;
double tEnd;
int tt;

int main(int argc, char *argv[]) 
{
  #ifdef BATCH
    SIM_ROOT_DIR = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
    ROOT_DIR = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));

    // arg[0] = program name
    // arg[1] = SIM_ROOT_DIR
    if (argc == 2) {          // should be 2 for batch execution
      sprintf(SIM_ROOT_DIR, "%s", argv[1]);
      sprintf(ROOT_DIR, "%s/flow_vel", SIM_ROOT_DIR);
    } else if (argc != 2) {   // print batch usage
      printf("usage: %s SIM_ROOT_DIR\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  #else   // To prevent compiler warnings
    argc = argc;
    argv = argv;
  #endif

  printf("\n SIM_ROOT_DIR = %s\n", SIM_ROOT_DIR);
  printf(" ROOT_DIR = %s\n\n", ROOT_DIR);

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
}


