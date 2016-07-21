#include "main.h"

// take care of batch job submission
char *SIM_ROOT_DIR;
char *ANALYSIS_DIR;

// Define global variables declared in header file
double tStart;
double tEnd;
int tt;

double initial_location;
double averaging_distance;
double radius;
double wavespeed;

// Main loop
int main(int argc, char *argv[]) 
{
  // Set directory structure
  directory_init(argc, argv);

  printf("Initializing...\n"); 
  fflush(stdout);
  // Read input file
  main_read_input();

  // Read and sort output directory for finding files within our time limits
  init_input_files();

  // Initialize domain and arrays
  domain_init(); 

  // Create output directory
  create_output();
  get_sigfigs();

  // Initialize variables to be used in the loop
  double dt = 0.;                             // Time step
  double dz = 0.;                             // change in z
  int dist_index = averaging_distance/dom.dz; // Averaging distance, index
  double norm = 1./(2.*dist_index);           // for normalization of averagin
  int ks, ke;                                 // indices to loop over

  double curr_location = initial_location;    // Current Location
  int init_index = (curr_location - dom.zs)/dom.dz;   // initial array index
  int curr_index = init_index;                        // Current array index

  // Loop over time
  printf("Looping...\n"); 
  fflush(stdout);
  for (tt = 0; tt < nFiles; tt++) {
    printf("  Timestep = %d of %d, Current Location (index) = %lf (%d)\n", 
      tt+1, nFiles, curr_location, curr_index);
    fflush(stdout);

    // Fill flow
    cgns_fill_input();

    // Moving average over width of 2*averaging_distance
    for (int i = 0; i < dom.xn; i++) {
      for (int j = 0; j < dom.yn; j++) {
        int pp2 = j + i*dom.yn;

        // Average over moving plane
        ks = curr_index - dist_index;
        ke = curr_index + dist_index;
        for (int k = ks; k <= ke; k++) {
          // Correct starting / ending indices for periodicity
          if (k < 0) {
            k += dom.zn;
          } else if (k > dom.zn) {
            k -= dom.zn;
          }
          int pp3 = k + j*dom.Gcc.s1 + i*dom.Gcc.s2;
          moving_avg[pp2] += norm*volume_fraction[pp3];
        }

        // Average over stationary plane
        ks = init_index - dist_index;
        ke = init_index + dist_index;
        for (int k = ks; k <= ke; k++) {
          // Correct starting / ending indices for periodicity
          if (k < 0) {
            k += dom.zn;
          } else if (k > dom.zn) {
            k -= dom.zn;
          }
          int pp3 = k + j*dom.Gcc.s1 + i*dom.Gcc.s2;
          stationary_avg[pp2] += norm*volume_fraction[pp3];
        }
      }
    }

    // write to file
    write_field(curr_location);

    // find next current location for every except last tstep
    if (tt < (nFiles - 1)) {
      // Find vertical shift based on wavespeed
      dt = fileTime[tt+1] - fileTime[tt];
      dz = wavespeed * dt; 

      // Set to next timestep -- keep location increasing relative to start
      curr_location += dz;
      curr_index = (curr_location - dom.zs)/dom.dz;   // Current array index
      
      // Correct index for periodicity
      if (curr_index > dom.zn) {
        curr_index -= dom.zn;
      } else if (curr_index < 0) {    // This should never trigger unless c < 0
        curr_index += dom.zn;
      }
    }
  }

  // Free and exit
  printf("Remember to copy grid.cgns to %s/%s!\n", SIM_ROOT_DIR, DATA_DIR);
  free_vars();
  printf("Done!\n");
  return EXIT_SUCCESS;
}
