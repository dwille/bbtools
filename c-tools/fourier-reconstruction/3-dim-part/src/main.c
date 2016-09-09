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

fftw_plan to_spectral;
fftw_plan to_cartesian;

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
  // get_sigfigs();

  // Init FFT
  printf("Initializing fftw plans...\n"); 
  fflush(stdout);
  /* Set up fftw plan
    -- n0, n1, n2   -- n2 varies the fastest
    -- uf, uf       -- in, out arrays -- overwrites!!
    -- FFTW_MEASURE -- find optimal plan
    -- these are UNNORMALIZED
   */
  to_spectral = fftw_plan_dft_3d(dom.zn, dom.yn, dom.xn, chi, chi, 1, 
    FFTW_MEASURE);
  to_cartesian = fftw_plan_dft_3d(dom.zn, dom.yn, dom.xn, chi, chi, -1, 
    FFTW_MEASURE);

  // Loop over time, calc coeffs, and reconstruct
  printf("Looping...\n"); 
  fflush(stdout);


  for (tt = 0; tt < nFiles; tt++) {
    printf("  Timestep = %d of %d\n", tt+1, nFiles);
    fflush(stdout);
    // Fill flow
    // cgns_fill_part();
    cgns_fill_flow();

    // Find coefficients
    fftw_execute(to_spectral);
    fftw_execute_dft(to_spectral, up_field, up_field);
    fftw_execute_dft(to_spectral, vp_field, vp_field);
    fftw_execute_dft(to_spectral, wp_field, wp_field);
    fftw_execute_dft(to_spectral, ke_field, ke_field);

    // Debug tests
    //test_fill();
    //test_out();

    // Zero out coeffs
    zero_coefficients();

    // Back transform
    fftw_execute(to_cartesian);
    fftw_execute_dft(to_cartesian, up_field, up_field);
    fftw_execute_dft(to_cartesian, vp_field, vp_field);
    fftw_execute_dft(to_cartesian, wp_field, wp_field);
    fftw_execute_dft(to_cartesian, ke_field, ke_field);

    // Normalize
    normalize();

    // write to file
    cgns_write_field();
  }

  // Free and exit
  printf("Remember to copy grid.cgns to %s/%s!\n", ANALYSIS_DIR, DATA_DIR);
  fftw_destroy_plan(to_spectral);
  fftw_destroy_plan(to_cartesian);
  free_vars();
  fftw_cleanup();
  printf("Done!\n");
  return EXIT_SUCCESS;
}

// Zero Coefficients
void zero_coefficients()
{
  /* Output format
    * Array of nk,nj,ni is returned as nk,nj,ni
    * k_(nn,mm,ll) is retrieved by cc = ll + ni*mm + ni*nj*nn
    * frequencies:
        ind:  0 1 2 3 4 ... N-3 N-2 N-1
        freq: 0 1 2 3 4 .... -3 -2 -1
    * proper normalization in N3 
      ** take k == 0, then phi_k should be sum(phase > -1) and 
         phi = sum(phase > -1)/sizeof(phase)
  */

  // Prefactor
  double iN3 = 1./(dom.Gcc.s3);
  int cc; 
  double cesaroX, cesaroY, cesaroZ;

  /* TODO:
   * Create zeroed array, get rid of triple for loop by only looping
   * over indices we know need to be filled -- O((2n)^3) rather than
   * order N^3
   */

  // Zero coefficients higher than we'd like, also divide by N3
  for (int nn = 0; nn < dom.zn; nn++) {
    cesaroZ = (1. - nn/(orderZ + 1.))*(nn <= orderZ)
            + (1. - (dom.zn - nn)/(orderZ + 1.))*(nn >= (dom.zn - orderZ));

    for (int mm = 0; mm < dom.yn; mm++) {
      cesaroY = (1. - mm/(orderY + 1.))*(mm <= orderY)
              + (1. - (dom.yn - mm)/(orderY + 1.))*(mm >= (dom.yn - orderY));


      for (int ll = 0; ll < dom.xn; ll++) {
        cesaroX = (1. - ll/(orderX + 1.))*(ll <= orderX)
                + (1. - (dom.xn - ll)/(orderX + 1.))*(ll >= (dom.xn - orderX));


        cc = ll +  dom.Gcc.s1*mm + dom.Gcc.s2*nn;
        
        chi[cc][0] *= iN3*cesaroX*cesaroY*cesaroZ;
        up_field[cc][0] *= iN3*cesaroX*cesaroY*cesaroZ;
        vp_field[cc][0] *= iN3*cesaroX*cesaroY*cesaroZ;
        wp_field[cc][0] *= iN3*cesaroX*cesaroY*cesaroZ;
        ke_field[cc][0] *= iN3*cesaroX*cesaroY*cesaroZ;

        chi[cc][1] *= iN3*cesaroX*cesaroY*cesaroZ;
        up_field[cc][1] *= iN3*cesaroX*cesaroY*cesaroZ;
        vp_field[cc][1] *= iN3*cesaroX*cesaroY*cesaroZ;
        wp_field[cc][1] *= iN3*cesaroX*cesaroY*cesaroZ;
        ke_field[cc][1] *= iN3*cesaroX*cesaroY*cesaroZ;
      }
    }
  }
}

void normalize()
{
  // Normalize q = (n*q)/n, or in this case q = (phi * q)/phi
  for (int k = 0; k < dom.zn; k++) {
    for (int j = 0; j < dom.yn; j++) {
      for (int i = 0; i < dom.xn; i++) {
        int cc = i + dom.Gcc.s1*j + dom.Gcc.s2*k;
        
        up_field[cc][0] /= chi[cc][0];
        vp_field[cc][0] /= chi[cc][0];
        wp_field[cc][0] /= chi[cc][0];
        ke_field[cc][0] /= chi[cc][0];

        #ifdef DEBUG
          up_field[cc][1] /= chi[cc][1];
          vp_field[cc][1] /= chi[cc][1];
          wp_field[cc][1] /= chi[cc][1];
          ke_field[cc][1] /= chi[cc][1];
        #endif
      }
    }
  }
}


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

