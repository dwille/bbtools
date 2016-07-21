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
  domain_init(); 

  // Create output directory
  create_output();

  // Loop over time, calc coeffs, and reconstruct
  int rInd, thInd;
  double idr = 1./dr;
  double idth = 1./dth;
  int tcount = 0;
  for (tt = 0; tt < nFiles; tt+=50) {
    tcount++;
    // Fill parts
    cgns_fill_parts();

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
            //th_ab = 0.;
            //r_ab = 0.25*R0;

            // Find correct index
            // r = 2*a + i*dr --if < 0 (contact) put ind = 0
            // th = dth*j
            rInd = floor((r_ab - 2.*meanR)*idr)*(r_ab >= 2.*meanR);
            // rInd = rInd*(rInd >= 0);
            thInd = floor(th_ab*idth);

            gHist[thInd + nBinsTh*rInd]++;
          }
        }
      }
    }
  }
  printf("tcount = %d\n", tcount);

  // Prefactor and divide by nfiles
  //double inFiles = 1./nFiles;
  double inFiles = 1./tcount;
  double prefactor = dom.xl*dom.yl*dom.zl/(nparts*(nparts-1.));
  for (int rr = 0; rr < nBinsR; rr++) {
    for (int th = 0; th < nBinsTh; th++) {
      gHist[th + rr*nBinsTh] *= inFiles * prefactor;
      //gHist[th + rr*nBinsTh] *= prefactor;
    }
  }

  write_avg_reconstruct();
   
  // Free and exit
  free_vars();
  printf("... Done!\n");
  return EXIT_SUCCESS;
}
