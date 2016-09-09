#include "main.h"

// take care of batch job submission
char *SIM_ROOT_DIR;
char *ANALYSIS_DIR;

// Define global variables declared in header file
double tStart;
double tEnd;
int npoints;
int tt;

/* Volume Fraction */
double min_vf = DBL_MAX;
double mean_vf = 0.;
double max_vf = -DBL_MAX;
double dBin_vf;
double binStart_vf;
double binEnd_vf;

/* Vertical Velocity */
double min_wp = DBL_MAX;
double mean_wp = 0.;
double max_wp = -DBL_MAX;
double dBin_wp;
double binStart_wp;
double binEnd_wp;

/* Kinetic Energy */
double min_ke = 0.;
double mean_ke = 0.;
double max_ke = -DBL_MAX;
double dBin_ke;
double binStart_ke;
double binEnd_ke;

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
  get_sigfigs();

  // Pull initial time step
  tt = 0;
  cgns_fill();

  // Find min and max and set up bins
  minmax(&min_vf, &max_vf, volume_fraction);
  minmax(&min_wp, &max_wp, part_w);
  minmax(&min_ke, &max_ke, part_ke);
  bin_init(min_vf, max_vf, &dBin_vf, &binStart_vf, &binEnd_vf);
  bin_init(min_wp, max_wp, &dBin_wp, &binStart_wp, &binEnd_wp);
  bin_init(min_ke, max_ke, &dBin_ke, &binStart_ke, &binEnd_ke);

  #ifdef DEBUG
    printf("  Bin Start (vf) = %lf\n", binStart_vf);
    printf("  Min (vf) = %lf\n", min_vf);
    printf("  dBin (vf) = %lf\n", dBin_vf);
    printf("  Max (vf) = %lf\n", max_vf);
    printf("  Bin End (vf) = %lf\n", binEnd_vf);
  #endif
    printf("  Bin Start (ke) = %lf\n", binStart_ke);
    printf("  Min (ke) = %lf\n", min_ke);
    printf("  dBin (ke) = %lf\n", dBin_ke);
    printf("  Max (ke) = %lf\n", max_ke);
    printf("  Bin End (ke) = %lf\n", binEnd_ke);

  // normalization for means
  double norm = 1./(dom.Gcc.s3 * nFiles);

  // Loop over time and bin
  printf("Looping...\n"); 
  fflush(stdout);
  int ind, ind_vf, ind_wp;

  for (tt = 0; tt < nFiles; tt++) {
    printf("  Timestep = %d of %d\n", tt+1, nFiles);
    fflush(stdout);
    // Fill data from cgns file
    cgns_fill();

    // Loop over space
    for (int cc = 0; cc < dom.Gcc.s3; cc++) {
      /* Volume Fraction */
      mean_vf += norm*volume_fraction[cc];
      // calculate index
      ind = (volume_fraction[cc] - binStart_vf)/dBin_vf;
      // if ind < 0, make it zero
      // if ind > nBins + 1, make it nBins +1
      ind = ind*(ind >= 0);         
      ind = ind*(ind <= (nBins + 1)) + (nBins + 1)*(ind > (nBins + 1));
      histogram_vf[ind]++;

      ind_vf = ind;

      /* Part Wp */
      mean_wp += norm*part_w[cc];
      ind = (part_w[cc] - binStart_wp)/dBin_wp;
      ind = ind*(ind >= 0);         
      ind = ind*(ind <= (nBins + 1)) + (nBins + 1)*(ind > (nBins + 1));
      histogram_wp[ind]++;

      ind_wp = ind;

      /* Bivariate */
      bihistogram_vf_wp[ind_wp + (nBins + 2)*ind_vf]++;

      /* Kinetic Energy */
      mean_ke += norm*part_ke[cc];
      ind = (part_ke[cc] - binStart_ke)/dBin_ke;
      ind = ind*(ind >= 0);         
      ind = ind*(ind <= (nBins + 1)) + (nBins + 1)*(ind > (nBins + 1));
      histogram_ke[ind]++;
    }

    // keep track of min max mean std for ALL files?
    // have bin padding as an input parameter?

  }
  // write to file
  write_field();

  // Free and exit
  printf("Done!\n");
  free_vars();
  return EXIT_SUCCESS;
}

/**** FUNCTIONS ****/
void minmax(double *min, double *max, double *array)
{
  // temporary min/max
  double tMin = *min;
  double tMax = *max;
  for (int i = 0; i < dom.Gcc.s3; i++) {
    tMin = (array[i] <= tMin)*array[i] +
          (array[i] > tMin)*tMin;
    tMax = (array[i] >= tMax)*array[i] +
          (array[i] < tMax)*tMax;
  }
  *min = tMin;
  *max = tMax;
}

void bin_init(double min, double max, double *dBin, double *binStart,
  double *binEnd)
{
  /* bins go from (min - dBin):dBin:(max + dBin) and length is nBins+2
   * the outer two bins will hold data that does is <min or >max
   * bins are created so that:
    * edge_j <= bin j < edge_{j+1}
  */
  *dBin = (max - min)/nBins;
  *binStart = min - *dBin;
  *binEnd = max + *dBin;
}

// if error checking on the indexing is needed
//      if (ind < 0) {
//        printf("Index is less than zero!\n");
//        printf("ind = %d\n", ind);
//        printf("volume_fraction[%d] = %lf\n", cc, volume_fraction[cc]); 
//        exit(EXIT_FAILURE);
//      } else if (ind > nBins + 1) {
//        printf("Index is > %d!\n", nBins + 1);
//        printf("ind = %d\n", ind);
//        printf("volume_fraction[%d] = %lf\n", cc, volume_fraction[cc]); 
//        exit(EXIT_FAILURE);
//      }
//
