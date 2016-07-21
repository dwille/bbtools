#include "main.h"

// Set up batch job submission
char *SIM_ROOT_DIR;
char *ANALYSIS_DIR;

// Define global variables declared in header file
double tStart;
double tEnd;
double R0;
int legendreOrder;
int laguerreOrder;
int printAvgFlag;    // 1 prints avg, 0 prints each tstep
int nPointsR;
int nPointsTh;
int tt;

double r_ab;
double mu_ab;
double P_ell;
double L_enn;

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

  // Allocate arrays
  alloc_arrays();

  // Create output directory
  create_output();
  get_sigfigs();

  // Loop over time, calc coeffs, and reconstruct
  double inFiles = 1./((double) nFiles);
  double prefactor = dom.xl*dom.yl*dom.zl/(4.*PI*R0*nparts*(nparts-1.));
  for (tt = 0; tt < nFiles; tt++) {
    // Fill parts
    cgns_fill_parts();

    // Reset summation variables
    reset_sums();

    /* Calculate coefficients */
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
          if (r_ab <= R0) {
            // loop over legendre order -- evens
            for (int ll = 0; ll <= legendreOrder; ll += 2) {

              // calculate P_l(mu_ab)
              eval_legendre(mu_ab, ll);

              // loop over laguerre order
              for (int nn = 0; nn <= laguerreOrder; nn++) {
                // calculate L_n(r_ab)
                eval_laguerre(r_ab, nn);

                // Calculate g_ln
                calc_coeffs(ll,nn);
              }
            }
          }
        }
      }
    }

    // Loop over legendre, fourier, and reconstruct
    for (int ll = 0; ll <= legendreOrder; ll += 2) {
      for (int nn = 0; nn <= laguerreOrder; nn++) {
        eval_ln(ll, nn);
      }
    }

    // Prefactor and add to averaged
    if (printAvgFlag == 1) {           // add to averaged
      for (int rr = 0; rr < nPointsR; rr++) {
        for (int th = 0; th < nPointsTh; th++) {
          g_rec_avg[th + rr*nPointsTh] += inFiles * prefactor * g_rec[th + rr*nPointsTh];
        }
      }

    } else {                        // keep discrete and output each step
      for (int rr = 0; rr < nPointsR; rr++) {
        for (int th = 0; th < nPointsTh; th++) {
          g_rec[th + rr*nPointsTh] *= prefactor;
        }
      }
      // write to file
      write_reconstruct();
    }

  }
  write_avg_reconstruct();
   
  // Free and exit
  free_vars();
  printf("... Done!\n");
  return EXIT_SUCCESS;
}
