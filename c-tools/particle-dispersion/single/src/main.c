#include "main.h"

// Global variables defined at bottom of this file


int main(int argc, char *argv[])
{
  // Set directory structure
  directory_init(argc, argv);

  // Read input config file
  read_input();

  // Parse sim output directory for cgns files
  init_cgns_files();

  // Create output directory
  create_output();

  // Initialize part struct and fill for first time step
  parts_init();

  // Initialize domain
  domain_init();

  // Assign GPU with CUDA_VISIBLE_DEVICES
  cuda_set_device();

  // Allocate dev memory and push
  cuda_dev_malloc();
  cuda_dom_push();
  cuda_part_push();

  // Allocate host memory for single particle statistics
  alloc_host_mem();

  // Fill initial particle positions
  cuda_fill_initial_positions(); 

  // Loop over time series and calculate stats
  for (tt = 0; tt < n_files; tt++) {
    // Read in new data and push to device
    cgns_fill_part_struct();
    cuda_part_push();

    // Deal with periodicity when particle move across periodic boundaries
    cuda_periodic_flip();

    // Calculate single-particle separations and average over particles
    cuda_find_separations();
    printf("time step = %d of %d\r", tt, n_files);
  }

  // Write results
  // TODO set this up to average over different runs (e.g. init times)
  write_output();

  // Free memory
  cuda_free();
  free_mem();
}

// Declarations
double t_start;
double t_end;
double *sim_time;
int nparts;
int tt;

part_struct *parts;
part_struct *_parts;

dom_struct dom;
dom_struct *_dom;

double *r2_total;
double *r2_horiz;
double *r2_verti;
double *_x0;      
double *_y0;      
double *_z0;      
double *_x;      
double *_y;      
double *_z;      
double *_r2_total;
double *_r2_horiz;
double *_r2_verti;
