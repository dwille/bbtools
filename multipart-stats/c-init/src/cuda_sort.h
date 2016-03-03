#ifndef _CUDA_SORT_H
#define _CUDA_SORT_H

extern "C"
{
#include "tetrad_init.h"
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
  dom_struct *binDom, int *neighborList, int *neighborCount, int nMax);

// find nChoose3 for each particle based on neighborCount
__global__ void choose3(int *neighborCount, int *nChoose3, int nparts);

// find combinations for each particle
__global__ void combine_nodes(int *_neighborList, int *_neighborCount,
  int *nodes, int *strides, int nparts, int nMax);

// sort combinations
__global__ void sort_combos(int *nodes, int nPerms);

// find unique combinations
__global__ void find_unique(int *nodes, int base, int remainder,
  int *isUnique);

// find unique combinations
__global__ void find_unique2(int *nodes, int *isUnique, int nPerms);

// pull unique combinations
__global__ void pull_unique(int *uniqueNodes, int *nodes, int *isUnique, 
  int nPerms, int *uniquePrefix, int nUnique);

// fill nodes
__global__ void fill_nodes(tetrad_struct *tetrads, int *uniqueNodes, 
  int nUnique);

// tetrad geometry
__global__ void tetrad_geometry(part_struct *parts, tetrad_struct *tetrads, 
  dom_struct *dom, int nUnique);

// matrix tests if necessary
__global__ void matrixTests(void);

// fix periodicity
__device__ void periodic_flip(double *r1, double *r2, double *r3, double *r4,
  tetrad_struct tetrads, part_struct *parts, double xl, double yl, double zl);

// matrix determinat
__device__ double matrixDet3(double *A);

// matrix trace
__device__ double matrixTrace3(double *A);

// trace(A^2)
__device__ double matrixSquaredTrace3(double *A);

// matrix inverse on device (3x3)
__device__ void matrixInverse3(double *A, double *invA);

// matrix multiplication on device
__device__ void matrixMult3(double *A, double *B, double *R);

// jacobi eigenvalue method
__device__ void jacobiEig3(double *a, double *d, double *v, int *nrot);
/*  --  a     -- Matrix to find eigenvalues and eigenvectors of
    --  d     -- Vector of eigenvalues in descending order
    --  v     -- Matrix of eigenvectors, columns contain corresponding vectors
    --  nrot  -- number of iterations to convergence
*/

// rotation kernel for jacobi eigenvalue method
__device__ void rot(double *a, double s, double tau, int i, int j, int k, 
  int l);

// eigenvalue sort for jacobi eigenvalue method
__device__ void eigsrt(double *d, double *v);

#endif
