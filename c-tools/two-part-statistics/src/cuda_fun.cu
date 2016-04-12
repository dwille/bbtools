#include "cuda_fun.h"
#include "time.h"

#include <cuda.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/sort.h>
#include <thrust/scan.h>


extern "C"
void cuda_dev_pull(void)
{
//  cudaMemcpy(neighborList, _neighborList, nMax*nparts*sizeof(int),
//    cudaMemcpyDeviceToHost);
}

extern "C"
void cuda_dom_push(void)
{
  cudaSetDevice(dev_start);
  // copy host data to device
  cudaMemcpy(_dom, &dom, sizeof(dom_struct), cudaMemcpyHostToDevice);
  cudaMemcpy(_binDom, &binDom, sizeof(dom_struct), cudaMemcpyHostToDevice);
}

extern "C"
void cuda_part_push(void)
{
  cudaSetDevice(dev_start);
  cudaMemcpy(_parts, parts, sizeof(part_struct) * nparts, 
    cudaMemcpyHostToDevice);
}

extern "C"
void cuda_dev_malloc(void)
{
  // allocate device memory on device
  cudaSetDevice(dev_start);

  cudaMalloc((void**) &(_parts), sizeof(part_struct) * nparts);
  cudaMalloc((void**) &(_partsPrev), sizeof(part_struct) * nparts);
  cudaMalloc((void**) &(_dom), sizeof(dom_struct));
  cudaMalloc((void**) &(_binDom), sizeof(dom_struct));
}

void cuda_find_pairs()
{
  // set up cuda threads and blocks
  int threads = MAX_THREADS_1D;
  int blocks = (int) ceil((double) nparts / (double) threads);
  if (threads > nparts) {
    threads = nparts;
    blocks = 1;
  }
  dim3 dimBlocks(threads);
  dim3 numBlocks(blocks);

  // set up bins and search for pairs
  if (nparts < NPOINTS) {
    printf("nparts = %d, no pairs to find.\n", nparts);
    exit(EXIT_FAILURE);
  } else if (nparts >= NPOINTS) {

    int nBins = binDom.Gcc.s3;

    // initialize threads for nBin size
    int threads_nb = MAX_THREADS_1D;
    int blocks_nb = (int) ceil((double) nBins / (double) threads_nb);
    if (threads_nb > nBins) {
      threads_nb = nBins;
      blocks_nb = 1;
    }
    dim3 dimBlocks_nb(threads_nb);
    dim3 numBlocks_nb(blocks_nb);

    // Go to each particle and find its bin
    int *_partInd;
    int *_partBin;

    cudaMalloc((void**) &_partInd, nparts*sizeof(int));
    cudaMalloc((void**) &_partBin, nparts*sizeof(int));
    bin_fill<<<numBlocks, dimBlocks>>>(_partInd, _partBin, nparts,
      _parts, _binDom, bc);
      
    /* sort by bin */
    thrust::device_ptr<int> ptr_partBin(_partBin);
    thrust::device_ptr<int> ptr_partInd(_partInd);
    thrust::sort_by_key(ptr_partBin, ptr_partBin + nparts, ptr_partInd);
    _partBin = thrust::raw_pointer_cast(ptr_partBin);
    _partInd = thrust::raw_pointer_cast(ptr_partInd);

    /* calculate start and end index of each bin */
    int *_binStart;
    int *_binEnd;
    cudaMalloc((void**) &_binStart, nBins*sizeof(int));
    cudaMalloc((void**) &_binEnd, nBins*sizeof(int));
    init<<<numBlocks_nb, dimBlocks_nb>>>(_binStart, nBins, -1);
    init<<<numBlocks_nb, dimBlocks_nb>>>(_binEnd, nBins, -1);

    int smemSize = sizeof(int)*(threads + 1);
    bin_start<<<blocks, threads, smemSize>>>(_binStart, _binEnd, _partBin,
      nparts);

    /*  FIND_NODES */
    // Find all neighbors in adjacent bins for each particle; count them
    int *_neighborList;
    int *_neighborCount;
    cudaMalloc((void**) &_neighborList, nMax*nparts*sizeof(int));
    cudaMalloc((void**) &_neighborCount, nparts*sizeof(int));
    init<<<numBlocks, dimBlocks>>>(_neighborCount, nparts, 0);

    printf("\tFinding possible pair permutations... ");
    find_nodes<<<numBlocks, dimBlocks>>>(_parts, nparts, _dom, bc, _binStart,
      _binEnd, _partBin, _partInd, _binDom, _neighborList, _neighborCount, 
      nMax, R0_a, eps_a, rmax);

    /*  CHOOSE2 of the neighboring nodes -- for each particle */
    int *_nChoose2;
    cudaMalloc((void **) &_nChoose2, nparts*sizeof(int));
    choose2<<<numBlocks, dimBlocks>>>(_neighborCount, _nChoose2, nparts);

    // Find total number of permutations -- sum _nChoose2
    thrust::device_ptr<int> ptr_nChoose2(_nChoose2);
    int nPermutations = thrust::reduce(ptr_nChoose2, ptr_nChoose2 + nparts);
    printf("Found %d.\n", nPermutations);

    // Find stride for each pairs
    int *_strides;
    cudaMalloc((void **) &_strides, nparts*sizeof(int));
    thrust::device_ptr<int> ptr_strides(_strides);
    thrust::exclusive_scan(ptr_nChoose2, ptr_nChoose2 + nparts, ptr_strides);

    // create array to hold particle index of pair nodes
    int *_nodes;
    cudaMalloc((void **) &_nodes, 2*nPermutations*sizeof(int));

    int threads_nodes = MAX_THREADS_1D;
    int blocks_nodes = (int) ceil((double) 2*nPermutations / (double) threads_nodes);
    if (threads_nodes > 2*nPermutations) {
      threads_nodes = 2*nPermutations;
      blocks_nodes = 1;
    }
    dim3 dimBlocks_nodes(threads_nodes);
    dim3 numBlocks_nodes(blocks_nodes);

    init<<<numBlocks_nodes, dimBlocks_nodes>>>(_nodes, 2*nPermutations, -1);

    /*  COMBINE_NODES */
    // parallelizing over all particles, find all combitions for each particle
    // also, sort in place
    printf("\tCombining nodes... ");
    combine_nodes<<<numBlocks, dimBlocks>>>(_neighborList, _neighborCount, 
      _nodes, _strides, nparts, nMax);

//    /*  SORT_COMBOS */
//    // Parallelizing over all permutations, sort each
//
//    printf("Done!\n\tSorting permutations... ");
//    sort_combos<<<numBlocks_perms, dimBlocks_perms>>>(_nodes, nPermutations);

    /*  FIND_UNIQUE */
    // compare and find unique combinations
    int threads_perms = MAX_THREADS_1D;
    int blocks_perms = (int) ceil((double) nPermutations / (double) threads_perms);
    if (threads_perms > nPermutations) {
      threads_perms = nPermutations;
      blocks_perms = 1;
    }
    dim3 dimBlocks_perms(threads_perms);
    dim3 numBlocks_perms(blocks_perms);

    int *_isUnique;
    cudaMalloc((void **) &_isUnique, nPermutations*sizeof(int));
    init<<<numBlocks_perms, dimBlocks_perms>>>(_isUnique, nPermutations, 1);

    // Loop over each permutation, then parallelize over the remaining
    printf("Done!\n\tLooping over permutations and finding unique sets...");
    for (int base = 0; base < (nPermutations - 1); base++) {

      // set up threads and blocks
      int remainder = nPermutations - base - 1;
      int remT = MAX_THREADS_1D;
      int remB = (int) ceil((double) remainder / (double) remT);
      if (remT > remainder) {
        remT = remainder;
        remB = 1;
      }
      dim3 dimBlocks_rem(remT);
      dim3 numBlocks_rem(remB);

      // determine whether target node is a duplicate, mark if so
      find_unique<<<numBlocks_rem, dimBlocks_rem>>>(_nodes, base, remainder, 
        _isUnique);  
    }
//    find_unique2<<<numBlocks_perms, dimBlocks_perms>>>(_nodes, _isUnique, 
//      nPermutations);

    // sum to find number of unique combinations
    thrust::device_ptr<int> ptr_isUnique(_isUnique);
    nPairs = thrust::reduce(ptr_isUnique, ptr_isUnique + nPermutations);

    printf("Found %d\n", nPairs);

    /*  PULL UNIQUE NODES */
    // Last entry is trash for finding indices and redirecting
    int *_uniqueNodes;
    cudaMalloc((void**) &_uniqueNodes, NPOINTS*(nPairs + 1)*sizeof(int));

    int threadsU = MAX_THREADS_1D;
    int blocksU = (int) ceil((double) NPOINTS*(nPairs + 1) / (double) threadsU);
    if (threadsU > NPOINTS*(nPairs + 1)) {
      threadsU = NPOINTS*(nPairs + 1);
      blocksU = 1;
    }
    dim3 dimBlocks_U(threadsU);
    dim3 numBlocks_U(blocksU);
    init<<<numBlocks_U, dimBlocks_U>>>(_uniqueNodes, NPOINTS*(nPairs + 1), -1);

    // Prefix sum of _isUnique
    int *_uniquePrefix;
    cudaMalloc((void **) &_uniquePrefix, nPermutations*sizeof(int));
    thrust::device_ptr<int> ptr_uPref(_uniquePrefix);
    thrust::device_ptr<int> ptr_isUn(_isUnique);
    thrust::inclusive_scan(ptr_isUn, ptr_isUn + nPermutations, ptr_uPref);

    printf("\tPulling unique nodes... ");
    pull_unique<<<numBlocks_perms, dimBlocks_perms>>>(_uniqueNodes, _nodes, 
      _isUnique, nPermutations, _uniquePrefix, nPairs);
    printf("Done!\n");
  
    /* FIND REGULAR */
    // Initialize pair struct for all unique pairs
    pair_struct *_pairs;
    cudaMalloc((void**) &(_pairs), sizeof(pair_struct) * nPairs);

    // Set up threads, blocks for each pair
    int threads_pairs = MAX_THREADS_1D;
    int blocks_pairs = (int) ceil((double) nPairs /(double) threads_pairs);
    if (threads_pairs > nPairs) {
      threads_pairs = nPairs;
      blocks_pairs = 1;
    }
    dim3 dimBlocks_pairs(threads_pairs);
    dim3 numBlocks_pairs(blocks_pairs);

    // Fill _pairs with the correct nodes
    fill_nodes<<<numBlocks_pairs, dimBlocks_pairs>>>(_pairs, _uniqueNodes, 
      nPairs);

    // check if initial
    flip_check<<<numBlocks_pairs, dimBlocks_pairs>>>(_parts, _pairs, _dom, 
      nPairs);

    printf("Done.\n");

    // Free variables
    cudaFree(_partInd);
    cudaFree(_partBin);
    cudaFree(_binStart);
    cudaFree(_binEnd);
    cudaFree(_neighborCount);
    cudaFree(_neighborList);
    cudaFree(_nChoose2);
    cudaFree(_strides);
    cudaFree(_nodes);
    cudaFree(_uniquePrefix);
    cudaFree(_isUnique);
    cudaFree(_uniqueNodes);
  }
}      

extern "C"
void cuda_pair_malloc(void)
{
  // Allocate pair struct on host and pull from device
  pairs = (pair_struct*) malloc(nPairs * sizeof(pair_struct));
  // Pull pairs back to host
  cudaMemcpy(pairs, _pairs, nPairs * sizeof(pair_struct), 
      cudaMemcpyDeviceToHost);

  cudaSetDevice(dev_start);

  cudaMalloc((void**) &(_RoG), sizeof(double) * nPairs);
  cudaMalloc((void**) &(_g1_s1), sizeof(double) * nPairs);
}

void cuda_periodic_flip(void)
{
  // Parallize over pairs
  int threads = MAX_THREADS_1D;
  int blocks = (int) ceil((double) nparts / (double) threads);
  if (threads > nparts) {
    threads = nparts;
    blocks = 1;
  }
  dim3 numBlocks(blocks);
  dim3 dimBlocks(threads);

  // Fix periodicity
  flip_kernel<<<numBlocks, dimBlocks>>>(_parts, _partsPrev, _dom, nparts);
}

void cuda_save_parts_prev(void)
{
  cudaMemcpy(_partsPrev, _parts, sizeof(part_struct) * nparts,
    cudaMemcpyDeviceToDevice);
}

void cuda_pair_stats(void)
{ 
  // Matrix tests
  //#ifdef DEBUG
  //  if (tt == 0) {
  //    matrixTests<<<1,1>>>();
  //  }
  //#endif

  // Parallelize over each pair
  int threads_pairs = MAX_THREADS_1D;
  int blocks_pairs = (int) ceil((double) nPairs / (double) threads_pairs);
  if (threads_pairs > nPairs) {
    threads_pairs = nPairs;
    blocks_pairs = 1;
  }
  dim3 dimBlocks_pairs(threads_pairs);
  dim3 numBlocks_pairs(blocks_pairs);

  // Calculate pair geometry and velocity measures
  pair_statistics<<<numBlocks_pairs, dimBlocks_pairs>>>(_parts, _pairs, _dom, _RoG, 
    nPairs, tt);

//  // If first timestep, save vectors for later comparison
//  if (tt == 0) {
//    cudaMemcpy(_gEigVecInit, _gEigVec, 9*sizeof(double)*nPairs,
//      cudaMemcpyDeviceToDevice);
//    cudaMemcpy(_sEigVecInit, _sEigVec, 9*sizeof(double)*nPairs,
//      cudaMemcpyDeviceToDevice);
//  }

  // Copy back raw data to host for writing to file
  cudaMemcpy(RoG, _RoG, sizeof(double) * nPairs, cudaMemcpyDeviceToHost);

  // Calculate alignment of vectors
//  align_vectors<<<numBlocks_pairs, dimBlocks_pairs>>>(_gEigVec, _sEigVec,
//    _vorticity, _gEigVecInit, _sEigVecInit, nPairs,
//    _g1_s1, _g1_s2, _g1_s3, _g2_s1, _g2_s2, _g2_s3, _g3_s1, _g3_s2, _g3_s3,
//    _g1_z, _g2_z, _g3_z, _s1_z, _s2_z, _s3_z, _w_z,
//    _w_g1, _w_g2, _w_g3, _w_s1, _w_s2, _w_s3);

  /* Higher order statistics */
  // Calculate higher stat moments of RoG, EVar, Shape, and Sii
  cuda_higher_moments(_RoG, nPairs, m_RoG);

  // Calculate higher stat moments of alignments
  cuda_higher_moments(_g1_s1, nPairs, m_g1_s1); 
}

void cuda_higher_moments(double *_array, int length, double *moments)
{
  // Adapted from Numerical Recipes, 14.1
  double mean, ep, var, sdev, skew, kurt;
  double N = (double) length;
  double iN = 1./N;

  // Calculate Mean
  thrust::device_ptr<double> ptr_array(_array);
  mean = thrust::reduce(ptr_array, ptr_array + length) * iN;

  // Allocate arrays for higher order stats calcs
  double *_diff;
  double *_diff2;
  double *_skew;
  double *_kurt;
  cudaMalloc((void**) &(_diff), sizeof(double) * nPairs);  // xi - x_bar
  cudaMalloc((void**) &(_diff2), sizeof(double) * nPairs); // (xi - x_bar)^2
  cudaMalloc((void**) &(_skew), sizeof(double) * nPairs);  // (xi - x_bar)^3
  cudaMalloc((void**) &(_kurt), sizeof(double) * nPairs);  // (xi - xbar)^4

  // Parallelize over each pair
  int threads = MAX_THREADS_1D;
  int blocks = (int) ceil((double) length / (double) threads);
  if (threads > length) {
    threads = length;
    blocks = 1;
  }
  dim3 numBlocks(blocks);
  dim3 dimBlocks(threads);

  // Parallelize calculation
  higher_moments_kernel<<<numBlocks, dimBlocks>>>(_array, mean, length,
    _diff, _diff2, _skew, _kurt);
  
  // cuda ptr to arrays
  thrust::device_ptr<double> ptr_diff(_diff);
  thrust::device_ptr<double> ptr_diff2(_diff2);
  thrust::device_ptr<double> ptr_skew(_skew);
  thrust::device_ptr<double> ptr_kurt(_kurt);

  // Parallel reduce arrays
  ep = thrust::reduce(ptr_diff, ptr_diff + length);
  var = thrust::reduce(ptr_diff2, ptr_diff2 + length);
  skew = thrust::reduce(ptr_skew, ptr_skew + length);
  kurt = thrust::reduce(ptr_kurt, ptr_kurt + length);

  // Correct them
  var = (var - ep*ep*iN)/(N - 1.);
  sdev = sqrt(var);

  double iVAR = 1./var;
  double iSDEV = 1./sdev;
  if (var != 0.) {
    skew *= iN*iVAR*iSDEV;
    kurt = kurt*iN*iVAR*iVAR - 3.;
  } else {
    skew = DBL_MAX;  
    kurt = DBL_MAX;  
  }

  // Assign to array
  moments[0] = mean;
  moments[1] = sdev;
  moments[2] = skew;
  moments[3] = kurt;

  // Free arrays
  cudaFree(_diff);
  cudaFree(_diff2);
  cudaFree(_skew);
  cudaFree(_kurt);
}

extern "C"
double cuda_sum(double *_array, int N)
{
  thrust::device_ptr<double> ptr_array(_array);
  double sum = thrust::reduce(ptr_array, ptr_array + N);
  return sum;
}

extern "C"
void cuda_dev_free(void)
{
  cudaFree(_parts);
  cudaFree(_partsPrev);
  cudaFree(_dom);
  cudaFree(_binDom);
  cudaFree(_pairs);

  cudaFree(_RoG);
  cudaFree(_g1_s1);

  cudaDeviceReset();
}
