#include "main.h"

// Set up batch job submission
char *SIM_ROOT_DIR;
char *ANALYSIS_DIR;

// Define global variables declared in header file
double tStart;
double tEnd;
double R0;
int nPointsR;
int nPointsTh;
int tt;

double r_ab;
double th_ab;
double mean_vf = 0.;
double sdev_vf = 0.;

int main(int argc, char *argv[]) 
{
  // Set directory structure
  directory_init(argc, argv);

  // Read input file
  main_read_input();

  // Read and sort output directory for finding files within our time limits
  init_files();

  // Initialize parts, domain, and hist
  parts_init();
  domain_init(); 
  histogram_init();

  // Create output directory
  create_output();

  // Pull first timestep, find mean and sdev of volume fraction
  vfrac_stats();

  // Loop over time, calc coeffs, and reconstruct
  int rInd, thInd;
  double idr = 1./dr;
  double idth = 1./dth;
  int filecount = 0;
  int partcount = 0;
  for (tt = 0; tt < nFiles; tt++) {
    filecount++;
    // Fill data from cgns file
    cgns_fill();

    // Loop over part-pair combinations
    for (int alpha = 0; alpha < nparts; alpha++) {
      for (int beta = 0; beta < nparts; beta++) {
        // don't include beta == alpha
        if (beta == alpha) {
          continue;
        } else {
          // find part-pair geometry
          part_pair_geometry(alpha,beta);

          // if we're outside the R0 range, don't do work
          if (r_ab < R0) {

            /* Find correct indices for r
             * r = 2*a + i*dr 
             * if ind < 0 (i.e. during contact) set ind = 0
             * this is equivalent to r <= 2r
             */
            rInd = floor((r_ab - 2.*meanR)*idr)*(r_ab >= 2.*meanR);

            /* Correct indices for the
             * th = dth*j
             */
            thInd = floor(th_ab*idth);

            // Increment the bincount
            gHist[thInd + nBinsTh*rInd]++;
            partcount++;
          }
        }
      }
    }
  }
  printf("  Found %d pairs\n", partcount);

  // Prefactor and divide by nfiles
  double inFiles = 1./filecount;
  double prefactor = dom.xl*dom.yl*dom.zl/(nparts*(nparts-1.));
  for (int rr = 0; rr < nBinsR; rr++) {
    for (int th = 0; th < nBinsTh; th++) {
      gHist[th + rr*nBinsTh] *= inFiles * prefactor;
    }
  }

  // Write the hist to a file
  write_avg_reconstruct();
   
  // Free and exit
  free_vars();
  printf("  ... Done!\n");
  return EXIT_SUCCESS;
}
