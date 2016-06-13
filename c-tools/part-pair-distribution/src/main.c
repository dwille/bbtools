#include "main.h"

// Set up batch job submission
#ifdef BATCH
  char *ROOT_DIR;
  char *SIM_ROOT_DIR;
#endif

// Define global variables declared in header file
double tStart;
double tEnd;
double R0;
int legendreOrder;
int fourierOrder;
int printAvgFlag;    // 1 prints avg, 0 prints each tstep
int nPointsR;
int nPointsTh;
int tt;

double r_ab;
double mu_ab;
double p_ell;

int main(int argc, char *argv[]) 
{
  // Batch submission
  #ifdef BATCH
    SIM_ROOT_DIR = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
    ROOT_DIR = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));

    // arg[0] = program name
    // arg[1] = SIM_ROOT_DIR
    if (argc == 2) {
      sprintf(SIM_ROOT_DIR, "%s", argv[1]);
      sprintf(ROOT_DIR, "%s/part-pair", SIM_ROOT_DIR);
    } else if (argc != 2) {
      printf("usage: %s SIM_ROOT_DIR\n", argv[0]);
      exit(EXIT_FAILURE);
    }
    printf("\n SIM_ROOT_DIR = %s\n", SIM_ROOT_DIR);
    printf(" ROOT_DIR = %s\n\n", ROOT_DIR);
  #else           // prevent compiler warning
    argc = argc;
    argv = argv;
  #endif

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
          if (r_ab > R0) {
            continue;
          } else {

            // loop over legendre order
            // TODO: only even order, because symmetric over th=pi/2?
            for (int ll = 0; ll <= legendreOrder; ll++) {

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
    for (int ll = 0; ll <= legendreOrder; ll++) {
      eval_l0(ll);
      for (int nn = 1; nn <= fourierOrder; nn++) {
        eval_ln(ll,nn);
      }
    }

    // Prefactor and add to averaged
    if (printAvgFlag == 1) {           // add to averaged
      for (int rr = 0; rr < nPointsR; rr++) {
        for (int th = 0; th < nPointsTh; th++) {
          g_ces_avg[th + rr*nPointsTh] += inFiles * prefactor * g_ces[th + rr*nPointsTh];
        }
      }

    } else {                        // keep discrete and output each step
      for (int rr = 0; rr < nPointsR; rr++) {
        for (int th = 0; th < nPointsTh; th++) {
          g_ces[th + rr*nPointsTh] *= prefactor;
        }
      }
      // write to file
      write_reconstruct();
    }

  }
  write_avg_reconstruct();


   
  // Free and exit
  free_vars();
  return EXIT_SUCCESS;
}
