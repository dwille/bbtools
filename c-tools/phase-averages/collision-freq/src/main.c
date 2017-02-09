#include "main.h"

// Define global variables declared in header file
char *SIM_ROOT_DIR;
char *ANALYSIS_DIR;
double tStart;
double tEnd;

int main(int argc, char *argv[]) 
{
  // Set directory structure
  directory_init(argc, argv);

  // Read input file
  main_read_input();

  // Read and sort output directory for finding files within our time limits
  init_part_files();

  // Get number of particles
  nparts = cgns_read_nparts();

  // Initialize domain and flow arrays
  #ifdef DEBUG
    show_domain();
  #endif

  // Create output directory
  create_output();

  // Loop over time
  //double inparts = 1./nparts;
  //double inFiles = 1./nFiles;
  for (int tt = 0; tt < nFiles; tt++) {
    int sum = cgns_fill_parts(tt);
    printf("t = %lf, sum = %d\n", partFileTime[tt], sum);

    // write to file
    write_mean(tt, sum);
  }

  // Free and exit
  free_vars();
  return EXIT_SUCCESS;
}
