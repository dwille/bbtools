#include "main.h"

// Define global variables declared in header file
double tStart;
double tEnd;
double R0;
int legendreOrder;
int fourierOrder;
int nPointsR;
int nPointsTh;
int tt;

double r_ab;
double mu_ab;
double p_ell;

int main(void) 
{
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
  for (tt = 0; tt < nFiles; tt++) {
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
          if (r_ab > R0) {
            continue;
          } else {

            // loop over legendre order
            for (int ll = 0; ll <= legendreOrder; ll++) {

              // calculate P_l(mu_ab)
              eval_legendre(mu_ab, ll);

              // Calculate g_l0
              const_coeffs();

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
    for (int ll = 0; ll <= legendreOrder; ll++) {
      eval_l0(ll);
      for (int nn = 1; nn <= fourierOrder; nn++) {
        eval_ln(ll,nn);
      }
    }

    // Prefactor
    for (int rr = 0; rr < nPointsR; rr++) {
      for (int th = 0; th < nPointsTh; th++) {
        g_ces[th + rr*nPointsTh] *= dom.xl*dom.yl*dom.zl/(4.*PI*nparts*(nparts-1.));
      }
    }


    // write to file
    write_reconstruct();
    // write the rest of the coefficients
    // TODO: write once
    //write_coeffs(-1);
  }


   
  // Free and exit
  free_vars();
  return EXIT_SUCCESS;
}
