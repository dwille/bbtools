#include "main.h"

// Set up batch job submission
char *SIM_ROOT_DIR;
char *ANALYSIS_DIR;

// Define global variables declared in header file
double tStart;
double tEnd;
double R0;
int legendreOrder;
int fourierOrder;
int printAvgFlag;    // 1 prints avg, 0 prints each tstep
int printCoeffsFlag; // 1 prints coeffs
int nPointsR;
int nPointsTh;
int tt;

double r_ab;
double mu_ab;
double p_ell;

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
  if (printCoeffsFlag == 1) write_coeffs(-1, -1, -1);

  // Loop over time, calc coeffs, and reconstruct
  double inFiles = 1./((double) nFiles);
  double prefactor = dom.xl*dom.yl*dom.zl/(4.*PI*(R0 - 2*meanR)*nparts*(nparts-1.));
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
          if (r_ab > R0) {
            continue;
          } else {

            // loop over legendre order -- evens
            for (int ll = 0; ll <= legendreOrder; ll += 2) {

              // calculate P_l(mu_ab)
              eval_legendre(mu_ab, ll);

              // Calculate g_l0
              const_coeffs(ll);

              // loop over fourier order
              for (int nn = 1; nn <= fourierOrder; nn++) {
                 // Calculate even and odd coefficients
                 calc_coeffs(ll, nn);
              }
            }
          }
        }
      }
    }

    // Loop over legendre, fourier, and reconstruct
    for (int ll = 0; ll <= legendreOrder; ll += 2) {
      eval_l0(ll);
      if (printCoeffsFlag == 1) write_coeffs(0, ll, -1);
      for (int nn = 1; nn <= fourierOrder; nn++) {
        eval_ln(ll,nn);
        if (printCoeffsFlag == 1) write_coeffs(1, ll, nn);
      }
    }
    if (printCoeffsFlag == 1) write_coeffs(2, -1, -1); // endline

    // Prefactor and add to averaged
    double r_norm, th_norm;
    if (printAvgFlag == 1) {           // add to averaged
      for (int rr = 0; rr < nPointsR; rr++) {
        r_norm = evalR[rr] + 0.5*dr;
        for (int th = 0; th < nPointsTh; th++) {
          th_norm = evalTh[th] + 0.5*dth;
          iNorm = 1./(2.*PI*r_norm*r_norm*dr*sin(th_norm)*dth);
          g_ces_avg[th + rr*nPointsTh] += inFiles * prefactor * iNorm *
            g_ces[th + rr*nPointsTh];
        }
      }

    } else {                        // keep discrete and output each step
      for (int rr = 0; rr < nPointsR; rr++) {
        r_norm = evalR[rr] + 0.5*dr;
        for (int th = 0; th < nPointsTh; th++) {
          th_norm = evalTh[th] + 0.5*dth;
          iNorm = 1./(2.*PI*r_norm*r_norm*dr*sin(th_norm)*dth);
          g_ces[th + rr*nPointsTh] *= prefactor * iNorm;
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
