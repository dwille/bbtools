#include "main.h"

// take care of batch job submission
#ifdef BATCH
  char *ROOT_DIR;
  char *SIM_ROOT_DIR;
#endif

// Define global variables declared in header file
double tStart;
double tEnd;
int order;
int coeffsOut;
int npoints;
int tt;

fftw_plan pU;
fftw_plan pV;
fftw_plan pW;

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
      sprintf(ROOT_DIR, "%s/f-rec-flow-1D", SIM_ROOT_DIR);
    } else if (argc != 2) {
      printf("usage: %s SIM_ROOT_DIR\n", argv[0]);
      exit(EXIT_FAILURE);
    }
    printf("\n SIM_ROOT_DIR = %s\n", SIM_ROOT_DIR);
    printf(" ROOT_DIR = %s\n\n", ROOT_DIR);
    fflush(stdout);
  #else           // prevent compiler warning
    argc = argc;
    argv = argv;
  #endif

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

  // Init FFT
  printf("Initializing fftw plans...\n"); 
  fflush(stdout);
  /* Set up fftw plan
    -- n0, n1, n2   -- n2 varies the fastest
    -- uf, uf       -- in, out arrays -- overwrites!!
    -- FFTW_MEASURE -- find optimal plan
    -- these are UNNORMALIZED so need to divide by 1/sqrt(nx*ny*nz)
   */
  pU = fftw_plan_dft_r2c_3d(dom.zn, dom.yn, dom.xn, uf, uf_k, FFTW_MEASURE);
  pV = fftw_plan_dft_r2c_3d(dom.zn, dom.yn, dom.xn, vf, vf_k, FFTW_MEASURE);
  pW = fftw_plan_dft_r2c_3d(dom.zn, dom.yn, dom.xn, wf, wf_k, FFTW_MEASURE);

  // Loop over time, calc coeffs, and reconstruct
  printf("Looping...\n"); 
  fflush(stdout);

  // Prefactor
  double iV = 1./(dom.xl*dom.yl*dom.zl);
  double k_enn, wf_Re, wf_Im;
  int kx = 0; kx = kx;
  int ky = 0; ky = ky;
  int cc;
  int halfIn = 0.5*dom.xn + 1;

  for (tt = 0; tt < nFiles; tt++) {
    // Fill flow
    cgns_fill_flow();

    // Find coefficients
    //fftw_execute(pU);
    //fftw_execute(pV);
    fftw_execute(pW);

    // Debug tests
    // test_fill();
    // test_out();

    /* Output format
      -- Array of nk,nj,ni is returned as nk,nj,(1+ni/2)
      -- k_(nn,mm,ll) is retrieved by
      --  cc = ll + (1+ni/2)*mm + (1+ni/2)*nj*nn
      -- Can loop from 0 <= nn <= (desired order) for all nn,mm,ll
      -- for 1-D reconstuction, take kx=ky=0
     */
    for (int nn = 0; nn <= order; nn++) {
      cc = halfIn*dom.Gcc.jn*nn;
      k_enn = 2.*PI*nn/dom.zl;

      wf_Re = wf_k[cc][0];
      wf_Im = wf_k[cc][1];
      
      for (int zz = 0; zz < npoints; zz++) {
        wf_rec_Re[zz + tt*npoints] += iV*(wf_Re*cos(k_enn*evalZ[zz]) 
                                            - wf_Im*sin(k_enn*evalZ[zz]));
        wf_rec_Im[zz + tt*npoints] += iV*(wf_Re*sin(k_enn*evalZ[zz]) 
                                            + wf_Im*cos(k_enn*evalZ[zz]));
      }
    }
  }

  // write to file
  write_reconstruct();
   
  // Free and exit
  fftw_destroy_plan(pU);
  fftw_destroy_plan(pV);
  fftw_destroy_plan(pW);
  free_vars();
  fftw_cleanup();
  printf("Done!\n");
  return EXIT_SUCCESS;
}
