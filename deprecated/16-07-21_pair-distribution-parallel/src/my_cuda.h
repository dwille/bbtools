#ifndef _MY_CUDA_H
#define _MY_CUDA_H

extern "C"
{
#include "main.h"
#include "reader.h"
}

// Fill the particle bins arrays -- partBin and partInd
__global__ void bin_fill(int *partInd, int *partBin, int nparts,
  part_struct *parts, dom_struct *binDom, BC bc);

// generic gpu array initialization
__global__ void init(int *array, int length, int val);

// find start and end indices of each in sorted array
__global__ void bin_start(int *binStart, int *binEnd, int *partBin, int nparts);

// find potential tetrad nodes for each particle
__global__ void find_nodes(part_struct *parts, int nparts, dom_struct *dom, 
  BC bc, int *binStart, int *binEnd, int *partBin, int *partInd, 
  dom_struct *binDom, int *neighborList, int *neighborCount,
  int orderL, int orderN);

// calculate the P_ln(r_ij, mu_ij)
__device__ double eval_legendre_poly(double mu_ij, int ell);

#endif
