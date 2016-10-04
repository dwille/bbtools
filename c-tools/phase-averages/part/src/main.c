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

  // Loop over time, avg at each time step
  // Also, keep track of average over all particles and time
  double up_overall_mean = 0.;
  double vp_overall_mean = 0.;
  double wp_overall_mean = 0.;
  double inparts = 1./nparts;
  double inFiles = 1./nFiles;
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

    // Keep track of mean over time
    up_overall_mean += upMean;
    vp_overall_mean += vpMean;
    wp_overall_mean += wpMean;

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

  // Normalize mean
  up_overall_mean *= inFiles;
  vp_overall_mean *= inFiles;
  wp_overall_mean *= inFiles;

  // Loop again, this time calculating variance
  double vertical_var = 0;
  double horizontal_var = 0;
  for (tt = 0; tt < nFiles; tt++) {
    // Fill parts
    cgns_fill_parts();

    // Loop over particles, add to variance
    // vertical = sqrt( < sum [ ( w(n,t) - <w>_t )^2 ] >_t )
    // horizont = sqrt( 0.5 * < sum [ (u(n,t) - <u>_t)^2 + 
    //                                (v(n,t) - <v>_t)^2 ] >_t )
    for (int pp = 0; pp < nparts; pp++) {
      vertical_var += (wp[pp] - wp_overall_mean)*(wp[pp] - wp_overall_mean);
      horizontal_var += (up[pp] - up_overall_mean)*(up[pp] - up_overall_mean) +
                        (vp[pp] - vp_overall_mean)*(vp[pp] - vp_overall_mean);
    }
    // Variance
    // horizontal

  }
  // Normalize mean
  vertical_var *= inFiles * inparts;
  horizontal_var *= 0.5 * inFiles * inparts;

  vertical_var = sqrt(vertical_var);
  horizontal_var = sqrt(horizontal_var);

  // Print to stdout
  printf("  Starting time = %lf\n", tStart);
  if (floor(tStart) == 0) {
    printf("  WARNING: Starting time set to %lf. If the simulation has a \n"
           "    start up time, the variance will reflect this. Choose a \n"
           "    steady-state time to avoid this issue\n", tStart);
  }
  printf("  Vertical variance = %lf\n", vertical_var);
  printf("  Horizontal variance = %lf\n", horizontal_var);
   
  // Free and exit
  free_vars();

  printf("Done!!\n");
  return EXIT_SUCCESS;
}
