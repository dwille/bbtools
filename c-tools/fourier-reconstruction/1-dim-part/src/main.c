#include "main.h"

// take care of batch job submission
char *SIM_ROOT_DIR;
char *ANALYSIS_DIR;

// Define global variables declared in header file
double tStart;
double tEnd;
int order_s;
int order_e;
int n_order;
int coeffsOut;
int npoints;
int tt;

double pref;

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

  // Prefactor
  pref = 2./(dom.xl*dom.yl*dom.zl);

  // Loop over time, calc coeffs, and reconstruct
  for (tt = 0; tt < nFiles; tt++) {
    // Fill parts
    cgns_fill_parts();

    // Calclate nql, add to cesaro sum
    for (int ll = order_s; ll <= order_e; ll++) {
      // Calculate nq0; init cesaro sum to this value
      if (ll == 0) {
        const_coeffs(ones, mean_vol, vfrac_ces);
        const_coeffs(wp, mean_vol, vfrac_wp_ces);
        const_coeffs(ones, 1., n_ces);
        const_coeffs(up, 1., nu_ces);
        const_coeffs(vp, 1., nv_ces);
        const_coeffs(wp, 1., nw_ces);
        const_coeffs(ke, 1., nke_ces);

        // write constant coeffs to file
        if (coeffsOut == 1) write_coeffs(0);
      } else {
        // Calculate wavenumber
        double ell = (double) ll;
        double k_ell = 2.*PI*ell/dom.zl;

        // Calculate even and odd coefficients
        calc_coeffs(&nl_even[ll], &nl_odd[ll], ones, parts, k_ell);
        calc_coeffs(&nul_even[ll], &nul_odd[ll], up, parts, k_ell);
        calc_coeffs(&nvl_even[ll], &nvl_odd[ll], vp, parts, k_ell);
        calc_coeffs(&nwl_even[ll], &nwl_odd[ll], wp, parts, k_ell);
        calc_coeffs(&nkel_even[ll], &nkel_odd[ll], ke, parts, k_ell);

        // Calculate volume fraction coeffs from number density
        eval_phase_avg(vfrac_ces, nl_even[ll], nl_odd[ll], evalZ, ell, k_ell);
        eval_phase_avg(vfrac_wp_ces, nwl_even[ll], nwl_odd[ll], evalZ, ell, k_ell);

        // Evaluate series at z
        eval_series(n_ces, nl_even[ll], nl_odd[ll], evalZ, ell, k_ell);
        eval_series(nu_ces, nul_even[ll], nul_odd[ll], evalZ, ell, k_ell);
        eval_series(nv_ces, nvl_even[ll], nvl_odd[ll], evalZ, ell, k_ell);
        eval_series(nw_ces, nwl_even[ll], nwl_odd[ll], evalZ, ell, k_ell);
        eval_series(nke_ces, nkel_even[ll], nkel_odd[ll], evalZ, ell, k_ell);
      }
    }
    // Normalize nq by n to find q
    normalize(nu_ces, n_ces);
    normalize(nv_ces, n_ces);
    normalize(nw_ces, n_ces);
    normalize(nke_ces, n_ces);

    // write the rest of the coefficients
    if (coeffsOut == 1) write_coeffs(-1);
  }

  // write to file
  write_reconstruct();
   
  // Free and exit
  free_vars();
  printf("Done!!\n");
  return EXIT_SUCCESS;
}
