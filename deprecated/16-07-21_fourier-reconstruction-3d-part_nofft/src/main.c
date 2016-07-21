#include "main.h"

// take care of batch job submission
#ifdef BATCH
  char *ROOT_DIR;
  char *SIM_ROOT_DIR;
#endif

// Define global variables declared in header file
double tStart;
double tEnd;
int orderL;
int orderM;
int orderN;
int tt;

double pref;

int main(int argc, char *argv[]) 
{
  // Set up batch submission
  #ifdef BATCH
    SIM_ROOT_DIR = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
    ROOT_DIR = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));

    // arg[0] = program name
    // arg[1] = SIM_ROOT_DIR
    if (argc == 2) {
      sprintf(SIM_ROOT_DIR, "%s", argv[1]);
      sprintf(ROOT_DIR, "%s/f-rec-part-3D", SIM_ROOT_DIR);
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
  get_sigfigs();
  create_output();

  /* MAIN LOOP */
  double kl, km, kn;
  int ll, mm, nn, corder;
  // Loop over time
  for (tt = 0; tt < nFiles; tt++) {
    // Fill parts with new info
    cgns_fill_parts();

    // Loop over all coefficients l,m,n (x,y,z)
    for (ll = 0; ll <= orderL; ll++) {
      kl = 2.*PI*ll/dom.xl;
      for (mm = 0; mm <= orderM; mm++) {
        km = 2.*PI*mm/dom.yl;
        for (nn = 0; nn <= orderN; nn++) {
          kn = 2.*PI*nn/dom.zl;
          corder = nn + (orderN + 1)*mm + (orderN + 1)*(orderM + 1)*ll;

          // Calculate coefficients n_lmn
          // TODO: need diff for vfrac
          calc_coeffs(n_lmn, ones, parts, kl, km, kn, corder);
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
  }
   
  // Free and exit
  free_vars();
  printf("... Done!\n");
  return EXIT_SUCCESS;
}
