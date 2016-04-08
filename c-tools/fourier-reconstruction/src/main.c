#include "main.h"

// Define global variables declared in header file
double tStart;
double tEnd;
int order;
int npoints;
int tt;

double pref;

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

  // Prefactor
  pref = 2./(dom.xl*dom.yl*dom.zl);

  // Loop over time, calc coeffs, and reconstruct
  for (tt = 0; tt < nFiles; tt++) {
    // Fill parts
    cgns_fill_parts();

    // Calculate nq0; init cesaro sum to this value
    const_coeffs(ones, n_ces);
    const_coeffs(up, nu_ces);
    const_coeffs(vp, nv_ces);
    const_coeffs(wp, nw_ces);

    // Calclate nql, add to cesaro sum
    for (int ll = 1; ll <= order; ll++) {
      double ell = (double) ll;

      // Calculate wavenumber
      double k_ell = 2.*PI*ell/dom.zl;

      // Calculate even and odd coefficients
      calc_coeffs(&nl_even[ll], &nl_odd[ll], ones, parts, k_ell);
      calc_coeffs(&nul_even[ll], &nul_odd[ll], up, parts, k_ell);
      calc_coeffs(&nvl_even[ll], &nvl_odd[ll], vp, parts, k_ell);
      calc_coeffs(&nwl_even[ll], &nwl_odd[ll], wp, parts, k_ell);

      // Calculate volume fraction coeffs from number density
      eval_vfrac(vFrac_ces, nl_even[ll], nl_odd[ll], evalZ, ell, k_ell);

      // Evaluate series at z
      eval_series(n_ces, nl_even[ll], nl_odd[ll], evalZ, ell, k_ell);
      eval_series(nu_ces, nul_even[ll], nul_odd[ll], evalZ, ell, k_ell);
      eval_series(nv_ces, nvl_even[ll], nvl_odd[ll], evalZ, ell, k_ell);
      eval_series(nw_ces, nwl_even[ll], nwl_odd[ll], evalZ, ell, k_ell);
    }
  }

  // write to file
  write_reconstruct();
   
  // Free and exit
  free_vars();
  return EXIT_SUCCESS;
}
