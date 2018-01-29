#include "cuda_fun.h"

// Fill init particle positions
__global__ void fill_init_positions(part_struct *parts, double *x0, double *y0,
  double *z0, int nparts)
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;

  if (pp < nparts) {
    x0[pp] = parts[pp].x;
    y0[pp] = parts[pp].y;
    z0[pp] = parts[pp].z;
  }
}

// Fix periodicity
__global__ void flip_kernel(part_struct *parts, double *x, double *y, double *z,
  dom_struct *dom, int nparts)
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;

  double dx; double dy; double dz;
  int cmp_i; int cmp_j; int cmp_k;
  int sgn_i; int sgn_j; int sgn_k;
  
  if (pp < nparts) {
    // Before checking, correct particle position with old flip count
    parts[pp].x += dom->xl*parts[pp].flip_count_i;
    parts[pp].y += dom->yl*parts[pp].flip_count_j;
    parts[pp].z += dom->zl*parts[pp].flip_count_k;

    // Position change since last timestep
    dx = x[pp] - parts[pp].x; 
    dy = y[pp] - parts[pp].y; 
    dz = z[pp] - parts[pp].z;

    // See if particle has crossed domain since last timestep
    //  1 iff it did
    //  0 iff it did not
    cmp_i = (fabs(dx) >= 0.5*dom->xl);
    cmp_j = (fabs(dy) >= 0.5*dom->yl);
    cmp_k = (fabs(dz) >= 0.5*dom->zl);

    // Get sign of dx,dy,dz
    //  -- if (prev - curr) > 0, went R->L, (+) a domain length
    //  -- if (prev - curr) < 0, went L->R, (-) a domain length
    sgn_i = (dx > 0.) - (dx < 0.);
    sgn_j = (dy > 0.) - (dy < 0.);
    sgn_k = (dz > 0.) - (dz < 0.);

    // increment or decrement the appropriate flipcount if so
    parts[pp].flip_count_i += cmp_i*sgn_i;
    parts[pp].flip_count_j += cmp_j*sgn_j;
    parts[pp].flip_count_k += cmp_k*sgn_k;

    // perform one more flip if necessary
    parts[pp].x += dom->xl*cmp_i*sgn_i;
    parts[pp].y += dom->yl*cmp_j*sgn_j;
    parts[pp].z += dom->zl*cmp_k*sgn_k;

    // save into x,y,z
    x[pp] = parts[pp].x;
    y[pp] = parts[pp].y;
    z[pp] = parts[pp].z;
  }
}

__global__ void find_separation(double *x, double *y, double *z, 
                                double *x0, double *y0, double *z0, 
                                double *r2_total, double *r2_horiz,
                                  double *r2_verti,
                                int nparts)
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;

  if (pp < nparts) {
    // Find separation in each component
    double dx = x[pp] - x0[pp];
    double dy = y[pp] - y0[pp];
    double dz = z[pp] - z0[pp];

    // Find total mean squared separation
    r2_total[pp] = dx*dx + dy*dy + dz*dz;
    r2_verti[pp] = dz*dz;
    r2_horiz[pp] = dx*dx + dy*dy;
  }
}
