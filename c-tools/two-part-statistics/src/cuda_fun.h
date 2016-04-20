#ifndef _CUDA_FUN_H
#define _CUDA_FUN_H

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

// find potential pair nodes for each particle
__global__ void find_nodes(part_struct *parts, int nparts, dom_struct *dom, 
  BC bc, int *binStart, int *binEnd, int *partBin, int *partInd, 
  dom_struct *binDom, int *partI, int *partJ, int *keepIJ, int *iMatches,
  int *initFlipCount, int nMax, double R0_a, double eps_a, double rmax);

// find strides
__global__ void find_strides(int *keepIJ, int *strides, int nparts, int nMax);

// find nChoose3 for each particle based on neighborCount
__global__ void choose2(int *neighborCount, int *nChoose3, int nparts);

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
__global__ void fill_nodes(pair_struct *pairs, int *uniqueNodes, int nPairs);

// check periodicity at ini
__global__ void flip_check(part_struct *parts, pair_struct *pairs,
  dom_struct *dom, int nPairs);

// Fix periodicity
__global__ void flip_kernel(part_struct *parts, part_struct *partsPrev,
  dom_struct *dom, int nparts);

// pair geometry
__global__ void pair_statistics(part_struct *parts, pair_struct *pairs, 
  dom_struct *dom, double *rSep, int nPairs, int tt);

// fix periodicity
__device__ void periodic_flip(double *r1, double *r2, 
  pair_struct *pairs, part_struct *parts, double xl, double yl, double zl);

// calculate statistical moments of scalar arrays
__global__ void higher_moments_kernel(double *array, double mean, int length,
  double *sum2, double *skew, double *kurt); //double *diff, double *diff2,

// Calculate alignment of vectors
__global__ void  align_vectors(double *gEigVec, double *sEigVec,
  double *vorticity, double *gEigVecInit, double *sEigVecInit, 
  int nRegular,
    double *g1_s1, double *g1_s2, double *g1_s3,
    double *g2_s1, double *g2_s2, double *g2_s3,
    double *g3_s1, double *g3_s2, double *g3_s3,
    double *g1_z, double *g2_z, double *g3_z,
    double *s1_z, double *s2_z, double *s3_z,
    double *w_z,
    double *w_g1, double *w_g2, double *w_g3,
    double *w_s1, double *w_s2, double *w_s3);

#endif
