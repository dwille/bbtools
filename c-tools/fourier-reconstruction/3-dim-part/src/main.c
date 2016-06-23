#include "main.h"

// Define global variables declared in header file
double tStart;
double tEnd;
int orderL;
int orderM;
int orderN;
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
  get_sigfigs();
  create_output();

  /* MAIN LOOP */
  // Loop over time
//  for (tt = 0; tt < nFiles; tt++) {
//    // Fill parts with new info
    tt = 0;
    cgns_fill_parts();

    // Loop over all coefficients l,m,n (x,y,z)
    for (int ll = 0; ll <= orderL; ll++) {
      double kl = 2.*PI*((double) ll)/dom.xl;

      for (int mm = 0; mm <= orderM; mm++) {
        double km = 2.*PI*((double) mm)/dom.yl;

        for (int nn = 0; nn <= orderN; nn++) {
          double kn = 2.*PI*((double) nn)/dom.zl;

          int corder = nn + (orderN + 1)*mm + (orderN + 1)*(orderM + 1)*ll;

          // Calculate coefficients n_lmn
          // TODO: need diff for vfrac
          calc_coeffs(n_lmn, ones, parts, kl, km, kn, corder);
          printf("(ll,mm,nn) = (%d,%d,%d)\n", ll, mm , nn);
          //printf("n_lmn[%d] = %f + %fi\n", corder, creal(n_lmn[corder]), cimag(n_lmn[corder]));
          // does it make sense to not even store n_lmn?

          // evaluate series for n_lmn at x,y,z
          // TODO: need diff for vfrac
          // TODO: parallelize? probably not, is not TOO slow
          // TODO: is n_rec always real?
          eval_series(n_rec, n_lmn, kl, km, kn, corder);

        }
      }
    }
    // Normalize nq by n to find q
    // make sure to divide by volume...
    //normalize(nu_ces, n_ces);

//  // write to file -- TODO take specific nq_rec as input
    cgns_write_field();
//  }
   
  // Free and exit
  free_vars();
  return EXIT_SUCCESS;
}
