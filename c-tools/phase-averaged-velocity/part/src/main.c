#include "main.h"

// Define global variables declared in header file
#ifdef BATCH
  char *ROOT_DIR;
  char *SIM_ROOT_DIR;
#endif
double tStart;
double tEnd;
int tt;

double upMean, vpMean, wpMean;
double upSdev, vpSdev, wpSdev;

int main(int argc, char *argv[]) 
{
  // setup batch submission
  #ifdef BATCH
    SIM_ROOT_DIR = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
    ROOT_DIR = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));

    // arg[0] = program name
    // arg[1] = SIM_ROOT_DIR
    if (argc == 2) {
      sprintf(SIM_ROOT_DIR, "%s", argv[1]);
      sprintf(ROOT_DIR, "%s/part_vel", SIM_ROOT_DIR);
    } else if (argc != 2) {
      printf("usage: %s SIM_ROOT_DIR\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  #else           // prevent compiler warning
    argc = argc;
    argv = argv;
  #endif

  printf("\n SIM_ROOT_DIR = %s\n", SIM_ROOT_DIR);
  printf(" ROOT_DIR = %s\n\n", ROOT_DIR);

  // Read input file
  main_read_input();

  // Read and sort output directory for finding files within our time limits
  init_part_files();

  // Get number of particles
  nparts = cgns_read_nparts();

  // Initialize partstruct and flow vars
  parts_init();

  // Initialize domain and flow arrays
  #ifdef DEBUG
    show_domain()
  #endif

  // Create output directory
  create_output();

  // Loop over time, avg
  double inparts = 1./nparts;
  for (tt = 0; tt < nFiles; tt++) {
    // Fill parts
    cgns_fill_parts();
    
    upMean = 0.;
    vpMean = 0.;
    wpMean = 0.;

    // Calculate mean
    for (int pp = 0; pp < nparts; pp++) {
      upMean += up[pp];
      vpMean += vp[pp];
      wpMean += wp[pp];
    }
    upMean *= inparts;
    vpMean *= inparts;
    wpMean *= inparts;

    // Calculate standard dev
    upSdev = 0.;
    vpSdev = 0.;
    wpSdev = 0.;
    for (int pp = 0; pp < nparts; pp++) {
      upSdev += (up[pp] - upMean)*(up[pp] - upMean);
      vpSdev += (vp[pp] - vpMean)*(vp[pp] - vpMean);
      wpSdev += (wp[pp] - wpMean)*(wp[pp] - wpMean);
    }
    upSdev = sqrt(upSdev * inparts);
    vpSdev = sqrt(vpSdev * inparts);
    wpSdev = sqrt(wpSdev * inparts);

    // write to file
    write_mean();
  }

   
  // Free and exit
  free_vars();

  return EXIT_SUCCESS;
}
