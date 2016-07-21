#include "main.h"

// Define global variables declared in header file
char *SIM_ROOT_DIR;
char *ANALYSIS_DIR;
double tStart;
double tEnd;
int tt;

double upMean, vpMean, wpMean;
double upSdev, vpSdev, wpSdev;

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

  printf("Done!!\n");
  return EXIT_SUCCESS;
}
