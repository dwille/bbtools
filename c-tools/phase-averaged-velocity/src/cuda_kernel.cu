#include "myCuda.h"

// create phase mask, multiply by vel
__global__ void phase_mask(double *uf, double *vf, double *wf, int *phase, 
  int N3)
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;;

  if (pp < N3) {
    phase[pp] = (phase[pp] == -1);
    uf[pp] *= phase[pp];
    vf[pp] *= phase[pp];
    wf[pp] *= phase[pp];
  }
}
