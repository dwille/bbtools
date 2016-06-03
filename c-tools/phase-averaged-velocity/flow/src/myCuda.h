#ifndef _myCuda_H
#define _myCuda_H

extern "C"
{
#include "main.h"
#include "reader.h"
}

__global__ void phase_mask(double *uf, double *vf, double *wf, int *phase,
  int N3);

#endif
