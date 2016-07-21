#include "main.h"

// take care of batch job submission
char *SIM_ROOT_DIR;
char *ANALYSIS_DIR;

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

  // Init FFT
  printf("Initializing fftw plans...\n"); 
  fflush(stdout);
  /* Set up fftw plan
    -- n0, n1, n2   -- n2 varies the fastest
    -- uf, uf       -- in, out arrays -- overwrites!!
    -- FFTW_MEASURE -- find optimal plan
    -- these are UNNORMALIZED
   */
  pU = fftw_plan_dft_r2c_3d(dom.zn, dom.yn, dom.xn, uf, uf_k, FFTW_MEASURE);
  pV = fftw_plan_dft_r2c_3d(dom.zn, dom.yn, dom.xn, vf, vf_k, FFTW_MEASURE);
  pW = fftw_plan_dft_r2c_3d(dom.zn, dom.yn, dom.xn, wf, wf_k, FFTW_MEASURE);

  // Loop over time, calc coeffs, and reconstruct
  printf("Looping...\n"); 
  fflush(stdout);

  // Prefactor
  double iN3 = 1./dom.Gcc.s3;
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
    // Need 0 <= k <= kz
    for (int nn = 0; nn <= order; nn++) {
      cc = halfIn*dom.Gcc.jn*nn;
      k_enn = 2.*PI*nn/dom.zl;

      wf_Re = wf_k[cc][0];
      wf_Im = wf_k[cc][1];
      
      for (int zz = 0; zz < npoints; zz++) {
        wf_rec_Re[zz + tt*npoints] += iN3*(wf_Re*cos(k_enn*evalZ[zz]) 
                                            - wf_Im*sin(k_enn*evalZ[zz]));
        wf_rec_Im[zz + tt*npoints] += iN3*(wf_Re*sin(k_enn*evalZ[zz]) 
                                            + wf_Im*cos(k_enn*evalZ[zz]));
      }
    }
    // Need zn - kz <= k < zn because of wrap around
    for (int nn = dom.zn - order; nn < dom.zn; nn++) {
      cc = halfIn*dom.Gcc.jn*nn;
      k_enn = -2.*PI*(dom.zn - nn)/dom.zl;

      wf_Re = wf_k[cc][0];
      wf_Im = wf_k[cc][1];
      
      for (int zz = 0; zz < npoints; zz++) {
        wf_rec_Re[zz + tt*npoints] += iN3*(wf_Re*cos(k_enn*evalZ[zz]) 
                                            - wf_Im*sin(k_enn*evalZ[zz]));
        wf_rec_Im[zz + tt*npoints] += iN3*(wf_Re*sin(k_enn*evalZ[zz]) 
                                            + wf_Im*cos(k_enn*evalZ[zz]));
      }
    }
  }

  // write to file
  for (int zz = 0; zz < npoints; zz++) {
    printf("zz = %lf: wf = %lf + %lf i\n", dom.zs+zz*dom.dz, 
      wf_rec_Re[zz], wf_rec_Im[zz]);
  }
  //write_reconstruct();
   
  // Free and exit
  fftw_destroy_plan(pU);
  fftw_destroy_plan(pV);
  fftw_destroy_plan(pW);
  free_vars();
  fftw_cleanup();
  printf("Done!\n");
  return EXIT_SUCCESS;
}
