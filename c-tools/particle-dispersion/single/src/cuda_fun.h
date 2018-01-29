#ifndef _CUDA_FUN_H
#define _CUDA_FUN_H

extern "C"
{
#include "main.h"
}

// Fill initial particle positions
__global__ void fill_init_positions(part_struct *parts, double *x0, double *y0,
  double *z0, int nparts);

// Fix periodicity
__global__ void flip_kernel(part_struct *parts, double *x, double *y, double *z,
  dom_struct  *dom, int nparts);

// find separations r_total, r_horiz, r_vert
__global__ void find_separation(double *x, double *y, double *z, double *x0, 
  double *y0, double *z0, double *r2_total, double *r2_horiz, double *r2_verti,
  int nparts);

#endif
