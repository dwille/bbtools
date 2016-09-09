#include "main.h"

// take care of batch job submission
char *SIM_ROOT_DIR;
char *ANALYSIS_DIR;

// Define global variables declared in header file
double tStart;
double tEnd;
int tt;

/* Histogram arrays defined at bottom of the file */

int main(int argc, char *argv[]) 
{
  // Set directory structure
  directory_init(argc, argv);

  printf("Initializing...\n"); 
  fflush(stdout);
  // Read input file
  main_read_input();

  // Read and sort output directory for finding files within our time limits
  init_files();

  // Initialize domain and flow arrays
  domain_init(); 

  // Create output directory
  create_output();

  // Pull initial time step
  tt = 0;
  cgns_fill();

  // Find min and set up dbin
  minmax();
  bin_init();

  // normalization for means
  double norm = 1./(nparts * nFiles);

  // Loop over time and bin
  printf("Looping...\n"); 
  fflush(stdout);

  int ind;
  for (tt = 0; tt < nFiles; tt++) {
    //printf("  Timestep = %d of %d\n", tt+1, nFiles);
    fflush(stdout);
    // Fill data from cgns file
    cgns_fill();

    // Loop over space
    for (int cc = 0; cc < nparts; cc++) {
      /* U */
      bins.U.mean += norm*U[cc];
      // calculate index
      ind = floor(U[cc]/bins.U.dBin);
      // if ind < 0, make it zero
      // if ind > nBins, make it nBins
      ind = ind*(ind < nBins) + (nBins - 1)*(ind >= nBins);
      hist_U[ind]++;

      /* Kinetic Energy */
      bins.ke.mean += norm*ke[cc];
      ind = floor(ke[cc]/bins.ke.dBin);
      ind = ind*(ind < nBins) + (nBins - 1)*(ind >= nBins);
      hist_ke[ind]++;

      /* Kinetic Energy (rotational) */
      bins.ke_rot.mean += norm*ke_rot[cc];
      ind = floor(ke_rot[cc]/bins.ke_rot.dBin);
      ind = ind*(ind < nBins) + (nBins - 1)*(ind >= nBins);
      hist_ke_rot[ind]++;

      /* Kinetic Energy (translational) */
      bins.ke_trans.mean += norm*ke_trans[cc];
      ind = floor(ke_trans[cc]/bins.ke_trans.dBin);
      ind = ind*(ind < nBins) + (nBins - 1)*(ind >= nBins);
      hist_ke_trans[ind]++;

      /* Temperature (vertical) */
      bins.T_z.mean += norm*T_z[cc];
      ind = floor(T_z[cc]/bins.T_z.dBin);
      ind = ind*(ind < nBins) + (nBins - 1)*(ind >= nBins);
      hist_T_z[ind]++;

      /* Temperature (perpendicular) */
      bins.T_perp.mean += norm*T_perp[cc];
      ind = floor(T_perp[cc]/bins.T_perp.dBin);
      ind = ind*(ind < nBins) + (nBins - 1)*(ind >= nBins);
      hist_T_perp[ind]++;

      /* Temperature (total) */
      bins.T.mean += norm*T[cc];
      ind = floor(T[cc]/bins.T.dBin);
      ind = ind*(ind < nBins) + (nBins - 1)*(ind >= nBins);
      hist_T[ind]++;

      /* ux */
      bins.ux.mean += norm*up[cc];
      ind = floor((up[cc] - bins.ux.min)/bins.ux.dBin);
      ind = (ind*(ind < nBins) + (nBins - 1)*(ind >= nBins))*(ind >= 0);
      hist_ux[ind]++;

      /* vy */
      bins.vy.mean += norm*vp[cc];
      ind = floor((vp[cc] - bins.vy.min)/bins.vy.dBin);
      ind = (ind*(ind < nBins) + (nBins - 1)*(ind >= nBins))*(ind >= 0);
      hist_vy[ind]++;

      /* wz */
      bins.wz.mean += norm*wp[cc];
      ind = floor((wp[cc] - bins.wz.min)/bins.wz.dBin);
      ind = (ind*(ind < nBins) + (nBins - 1)*(ind >= nBins))*(ind >= 0);
      hist_wz[ind]++;
    }
  }
  // write to file
  write_field();

  // Free and exit
  printf("Done!\n");
  free_vars();
  return EXIT_SUCCESS;
}

/**** FUNCTIONS ****/
void minmax(void)
{
  bins.U.min = 0.;
  bins.U.max = -DBL_MAX;
  bins.ke_trans.min = 0.;
  bins.ke_trans.max = -DBL_MAX;
  bins.ke_rot.min = 0.;
  bins.ke_rot.max = -DBL_MAX;
  bins.ke.min = 0.;
  bins.ke.max = -DBL_MAX;
  bins.T_z.min = 0.;
  bins.T_z.max = -DBL_MAX;
  bins.T_perp.min = 0.;
  bins.T_perp.max = -DBL_MAX;
  bins.T.min = 0.;
  bins.T.max = -DBL_MAX;
  
  bins.ux.min = DBL_MAX;
  bins.ux.max = -DBL_MAX;
  bins.vy.min = DBL_MAX;
  bins.vy.max = -DBL_MAX;
  bins.wz.min = DBL_MAX;
  bins.wz.max = -DBL_MAX;

  for (tt = 0; tt < 100; tt++) {
    for (int i = 0; i < nparts; i++) {
      // Find max
      bins.U.max = (U[i] >= bins.U.max)*U[i] + 
                   (U[i] <  bins.U.max)*bins.U.max; 
      bins.ke_trans.max = (ke_trans[i] >= bins.ke_trans.max)*ke_trans[i] + 
                          (ke_trans[i] <  bins.ke_trans.max)*bins.ke_trans.max; 
      bins.ke_rot.max = (ke_rot[i] >= bins.ke_rot.max)*ke_rot[i] + 
                        (ke_rot[i] <  bins.ke_rot.max)*bins.ke_rot.max; 
      bins.ke.max = (ke[i] >= bins.ke.max)*ke[i] + 
                    (ke[i] <  bins.ke.max)*bins.ke.max; 
      bins.T_z.max = (T_z[i] >= bins.T_z.max)*T_z[i] + 
                     (T_z[i] <  bins.T_z.max)*bins.T_z.max; 
      bins.T_perp.max = (T_perp[i] >= bins.T_perp.max)*T_perp[i] + 
                        (T_perp[i] <  bins.T_perp.max)*bins.T_perp.max; 
      bins.T.max = (T[i] >= bins.T.max)*T[i] + 
                   (T[i] <  bins.T.max)*bins.T.max; 
      bins.ux.max = (up[i] >= bins.ux.max)*up[i] + 
                    (up[i] <  bins.ux.max)*bins.ux.max; 
      bins.vy.max = (vp[i] >= bins.vy.max)*vp[i] + 
                    (vp[i] <  bins.vy.max)*bins.vy.max; 
      bins.wz.max = (wp[i] >= bins.wz.max)*wp[i] + 
                    (wp[i] <  bins.wz.max)*bins.wz.max; 

      // Find min
      bins.ux.min = (up[i] <= bins.ux.min)*up[i] + 
                    (up[i] >  bins.ux.min)*bins.ux.min; 
      bins.vy.min = (vp[i] <= bins.vy.min)*vp[i] + 
                    (vp[i] >  bins.vy.min)*bins.vy.min; 
      bins.wz.min = (wp[i] <= bins.wz.min)*wp[i] + 
                    (wp[i] >  bins.wz.min)*bins.wz.min; 
    }
  }
}

void bin_init(void)
{
  /* bins go from min:dBin:max and length is nBins
   * the outer two bins will hold data that does is <min or >max
   * bins are created so that:
   ** edge_j <= bin j < edge_{j+1}
   * TODO max mean std for ALL files?
  */
  bins.U.mean = 0.;
  bins.U.dBin = 0;
  bins.ke_trans.mean = 0.;
  bins.ke_trans.dBin = 0;
  bins.ke_rot.mean = 0.;
  bins.ke_rot.dBin = 0;
  bins.ke.mean = 0.;
  bins.ke.dBin = 0;
  bins.T_z.mean = 0.;
  bins.T_z.dBin = 0;
  bins.T_perp.mean = 0.;
  bins.T_perp.dBin = 0;
  bins.T.mean = 0.;
  bins.T.dBin = 0;
  bins.ux.mean = 0.;
  bins.ux.dBin = 0.;
  bins.vy.mean = 0.;
  bins.vy.dBin = 0.;
  bins.wz.mean = 0.;
  bins.wz.dBin = 0.;


  bins.U.dBin = bins.U.max / nBins;
  bins.ke_trans.dBin = bins.ke_trans.max / nBins;
  bins.ke_rot.dBin = bins.ke_rot.max / nBins;
  bins.ke.dBin = bins.ke.max / nBins;
  bins.T_z.dBin = bins.T_z.max / nBins;
  bins.T_perp.dBin = bins.T_perp.max / nBins;
  bins.T.dBin = bins.T.max / nBins;
  bins.ux.dBin = (bins.ux.max - bins.ux.min)/nBins;
  bins.vy.dBin = (bins.vy.max - bins.vy.min)/nBins;
  bins.wz.dBin = (bins.wz.max - bins.wz.min)/nBins;
}
