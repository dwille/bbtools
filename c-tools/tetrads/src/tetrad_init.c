#include "tetrad_init.h"

// Define global variables declared in header file
int dev_start;
double tStart;
double tEnd;
double R0_a;
double eps_a;
double EVarCutLow;
double EVarCutHigh;
double shapeCutLow;
double shapeCutHigh;
int findTets;
int multRuns;
int nRegular;
int tt;

int main(void) 
{
  // Read tetrad input file
  tetrad_read_input();

  // Read and sort output directory for finding files within our time limits
  // and also initialize time array
  init_input_files();

  // get sig figs for writing
  get_sigfigs();

  // Create output directory and init mean output file
  create_output_dir();
  init_stat_output();

  // Messy hack for taking advantage of CUDA_VISIBLE_DEVICES in SLURM
  dev_start = read_devices();

  // Get number of particles
  nparts = cgns_read_nparts();

  // Initialize part_struct and fill for first timestep
  parts_init();

  // Initialize domain, binDom and fill
  domain_init(); 

  // Allocate device memory
  cuda_dev_malloc();
  cuda_dom_push(); 
  cuda_part_push(); 

  // Bin particles, find tetrads
  if (findTets == 1) {
    cuda_find_tetrads();
  } else {
    printf("  regularNodes already exists, reading in file.\n");
    printf("  Except this feature doesn't exist yet.\n");
    printf("  See so-18737117, 12733105\n");
    exit(EXIT_FAILURE);
  }

  // Allocate host memory for tetrad statistics
  alloc_tetrad_arrays();

  // Allocate device memory for tetrads
  cuda_tetrad_malloc();

  // Calculate tetrad stats
  for (tt = 0; tt < nFiles; tt++) {
    // Read in new data and push to device
    cgns_fill_part_struct();
    cuda_part_push();

    // Take care of periodicity
    //  -- Tetrads that are intialized over periodic boundaries need an
    //      extra flip -- this is determined in check tolerances, and
    //      done in tetrad_geometry
    //  -- All particles need to be checked at tt > 0 to see if they
    //      went over a periodic boundary -- this is done in cuda_periodic_flip
    if (tt > 0) {
      cuda_periodic_flip();
    }

    // Calculate tetrads
    cuda_tetrad_stats();
    
    // Write to file
    write_timestep();

    // Save current partstruct data to previous timestep
    cuda_save_parts_prev();
  }

  // write regular nodes to file, as well as nTetrads and nTsteps
  write_info();
   
  // Free and exit
  free_parts();
  cuda_dev_free();
  return EXIT_SUCCESS;
}

double m_RoG[4];
double m_EVar[4];
double m_Shape[4];
double m_I1[4];
double m_I2[4];
double m_I3[4];
double m_S11[4];
double m_S22[4];
double m_S33[4];

double m_g1_s1[4];  // Alignment of principal shape axes with initial 
double m_g2_s1[4];  // strain axes
double m_g3_s1[4];
double m_g1_s2[4];
double m_g2_s2[4];
double m_g3_s2[4];
double m_g1_s3[4];
double m_g2_s3[4];
double m_g3_s3[4];
double m_g1_z[4];    // Alignment of shape, strain, vorticity with gravity
double m_g2_z[4];
double m_g3_z[4];
double m_s1_z[4];
double m_s2_z[4];
double m_s3_z[4];
double m_w_z[4];
double m_w_g1[4];    // Alignment of vorticity with initial shape, strain axes
double m_w_g2[4];
double m_w_g3[4];
double m_w_s1[4];
double m_w_s2[4];
double m_w_s3[4];
double m_vortMag[4];

double *_g1_s1;   // Alignment of principal shape axes with initial strain axes
double *_g2_s1;
double *_g3_s1;
double *_g1_s2;
double *_g2_s2;
double *_g3_s2;
double *_g1_s3;
double *_g2_s3;
double *_g3_s3;
double *_g1_z;    // Alignment of shape, strain, vorticity with gravity
double *_g2_z;
double *_g3_z;
double *_s1_z;
double *_s2_z;
double *_s3_z;
double *_w_z;
double *_w_g1;    // Alignment of vorticity with initial shape, strain axes
double *_w_g2;
double *_w_g3;
double *_w_s1;
double *_w_s2;
double *_w_s3;
double *_vortMag;
