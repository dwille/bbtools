#include "main.h"

// take care of batch job submission
char *SIM_ROOT_DIR;
char *ANALYSIS_DIR;

// Define global variables declared in header file
double tStart;
double tEnd;
int orderX;
int orderY;
int orderZ;
int npoints;
int tt;

fftw_plan chi2phi_k;
fftw_plan phi_k2phi;

int main(int argc, char *argv[]) 
{
  // Set directory structure
  directory_init(argc, argv);

  printf("Initializing...\n"); 
  fflush(stdout);
  // Read input file
  main_read_input();

  // Read and sort output directory for finding files within our time limits
  init_flow_files();

  // Initialize domain and flow arrays
  domain_init(); 

  // Create output directory
  create_output();
  get_sigfigs();

  // Init FFT
  printf("Initializing fftw plans...\n"); 
  fflush(stdout);
  /* Set up fftw plan
    -- n0, n1, n2   -- n2 varies the fastest
    -- uf, uf       -- in, out arrays -- overwrites!!
    -- FFTW_MEASURE -- find optimal plan
    -- these are UNNORMALIZED
   */
  chi2phi_k = fftw_plan_dft_3d(dom.zn, dom.yn, dom.xn, chi, phi_k, 1, FFTW_MEASURE);
  phi_k2phi = fftw_plan_dft_3d(dom.zn, dom.yn, dom.xn, phi_k, phi, -1, FFTW_MEASURE);

  // Loop over time, calc coeffs, and reconstruct
  printf("Looping...\n"); 
  fflush(stdout);

  // Prefactor
  double iN3 = 1./(dom.Gcc.s3);
  int cc, xCheck, yCheck, zCheck, zeroCheck;

  for (tt = 0; tt < nFiles; tt++) {
    printf("  Timestep = %d of %d\n", tt+1, nFiles);
    fflush(stdout);
    // Fill flow
    cgns_fill_flow();

    // Find coefficients
    fftw_execute(chi2phi_k);

    // Debug tests
    //test_fill();
    //test_out();

    /* Output format
      -- Array of nk,nj,ni is returned as nk,nj,(1+ni/2)
      -- k_(nn,mm,ll) is retrieved by
      --  cc = ll + (1+ni/2)*mm + (1+ni/2)*nj*nn
      -- Can loop from 0 <= nn <= (desired order) for all nn,mm,ll
      -- proper normalization in N3 
        -- take k == 0, then phi_k should be sum(phase > -1) and 
            phi = sum(phase > -1)/sizeof(phase)
     */

    // Zero coefficients higher than we'd like, also divide by N3
    for (int nn = 0; nn < dom.zn; nn++) {
      zCheck = (nn <= orderZ) || (nn >= (dom.zn - orderZ));
      for (int mm = 0; mm < dom.yn; mm++) {
        yCheck = (mm <= orderY) || (mm >= (dom.yn - orderY));
        for (int ll = 0; ll < dom.xn; ll++) {
          xCheck = (ll <= orderX) || (ll >= (dom.xn - orderX));
          cc = ll +  dom.Gcc.s1*mm + dom.Gcc.s2*nn;
          
          zeroCheck = xCheck*yCheck*zCheck;
          phi_k[cc][0] *= iN3*zeroCheck;
          phi_k[cc][1] *= iN3*zeroCheck;
        }
      }
    }
    // TODO: perhaps output k=0 and check volume fraction?

    // Back transform
    fftw_execute(phi_k2phi);

    // find max and print 
   /*double max = -DBL_MAX;
    double min = DBL_MAX;
    for (int i = 0; i < dom.Gcc.s3; i++) {
      if (phi[i][0] > max) {    // Max
        max = phi[i][0];
      }
      if (phi[i][0] < min) {    // Min
        min = phi[i][0];
      }
      if (fabs(phi[i][0] - chi[i][0]) >= 1e-16) { // check reconstruction
        printf("phi[%d] = %lf, chi[%d] = %lf\n", i, phi[i][0], i, chi[i][0]);
      }
    }
    printf("range = [%lf, %lf]\n", min, max); */

    // write to file
    cgns_write_field();
  }

  // Free and exit
  printf("Remember to copy grid.cgns to %s/%s!\n", SIM_ROOT_DIR, DATA_DIR);
  fftw_destroy_plan(chi2phi_k);
  fftw_destroy_plan(phi_k2phi);
  free_vars();
  fftw_cleanup();
  printf("Done!\n");
  return EXIT_SUCCESS;
}
