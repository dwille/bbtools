#include "cuda_sort.h"
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

void cuda_find_tetrads()
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

  // set up bins and search for tetrads
  if (nparts < 4) {
    printf("nparts = %d, no tetrads to find.\n", nparts);
    exit(EXIT_FAILURE);
  } else if (nparts >= 4) {

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

    printf("\tFinding possible tetrad permutations... ");
    find_nodes<<<numBlocks, dimBlocks>>>(_parts, nparts, _dom, bc, _binStart,
      _binEnd, _partBin, _partInd, _binDom, _neighborList, _neighborCount, 
      nMax);

    /*  CHOOSE3 */
    int *_nChoose3;
    cudaMalloc((void **) &_nChoose3, nparts*sizeof(int));
    choose3<<<numBlocks, dimBlocks>>>(_neighborCount, _nChoose3, nparts);

    // Find total number of permutations -- sum _nChoose3
    thrust::device_ptr<int> ptr_nChoose3(_nChoose3);
    int nPerms = thrust::reduce(ptr_nChoose3, 
      ptr_nChoose3 + nparts);
    int totalNodes = 4*nPerms;
    printf("Found %d.\n", nPerms);

    // Find stride for each particle
    int *_strides;
    cudaMalloc((void **) &_strides, nparts*sizeof(int));
    thrust::device_ptr<int> ptr_strides(_strides);
    thrust::exclusive_scan(ptr_nChoose3, ptr_nChoose3 + nparts, ptr_strides);

    // create array to hold particle index of tetrad nodes
    int *_nodes;
    cudaMalloc((void **) &_nodes, totalNodes*sizeof(int));

    int threads_nodes = MAX_THREADS_1D;
    int blocks_nodes = (int) ceil((double) totalNodes / (double) threads_nodes);
    if (threads_nodes > totalNodes) {
      threads_nodes = totalNodes;
      blocks_nodes = 1;
    }
    dim3 dimBlocks_nodes(threads_nodes);
    dim3 numBlocks_nodes(blocks_nodes);

    init<<<numBlocks_nodes, dimBlocks_nodes>>>(_nodes, totalNodes, -1);

    /*  COMBINE_NODES */
    // parallelizing over all particles, find all combitions for each particle
    printf("\tCombining nodes... ");
    combine_nodes<<<numBlocks, dimBlocks>>>(_neighborList, _neighborCount, 
      _nodes, _strides, nparts, nMax);

    /*  SORT_COMBOS */
    // Parallelizing over all permutations, sort each
    int threads_perms = MAX_THREADS_1D;
    int blocks_perms = (int) ceil((double) nPerms / (double) threads_perms);
    if (threads_perms > nPerms) {
      threads_perms = nPerms;
      blocks_perms = 1;
    }
    dim3 dimBlocks_perms(threads_perms);
    dim3 numBlocks_perms(blocks_perms);

    printf("Done!\n\tSorting permutations... ");
    sort_combos<<<numBlocks_perms, dimBlocks_perms>>>(_nodes, nPerms);

    /*  FIND_UNIQUE */
    // compare and find unique ones
    int *_isUnique;
    cudaMalloc((void **) &_isUnique, nPerms*sizeof(int));
    init<<<numBlocks_perms, dimBlocks_perms>>>(_isUnique, nPerms, 1);

    // Loop over each permutations, then parallelize over the remaining
    printf("Done!\n\tLooping over permutations and finding unique sets...");
    for (int base = 0; base < (nPerms - 1); base++) {

      // set up threads and blocks
      int remainder = nPerms - base - 1;
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
//      nPerms);

    // sum to find number of unique combinations
    thrust::device_ptr<int> ptr_isUnique(_isUnique);
    int nUnique = thrust::reduce(ptr_isUnique, ptr_isUnique + nPerms);

    printf("Found %d\n", nUnique);

    /*  PULL UNIQUE NODES */
    // Last entry is trash for finding indices and redirecting
    int *_uniqueNodes;
    cudaMalloc((void**) &_uniqueNodes, 4*(nUnique + 1)*sizeof(int));

    int threadsU = MAX_THREADS_1D;
    int blocksU = (int) ceil((double) 4*(nUnique + 1) / (double) threadsU);
    if (threadsU > 4*(nUnique + 1)) {
      threadsU = 4*(nUnique + 1);
      blocksU = 1;
    }
    dim3 dimBlocks_U(threadsU);
    dim3 numBlocks_U(blocksU);
    init<<<numBlocks_U, dimBlocks_U>>>(_uniqueNodes, 4*(nUnique + 1), -1);

    // Prefix sum of _isUnique
    int *_uniquePrefix;
    cudaMalloc((void **) &_uniquePrefix, nPerms*sizeof(int));
    thrust::device_ptr<int> ptr_uPref(_uniquePrefix);
    thrust::device_ptr<int> ptr_isUn(_isUnique);
    thrust::inclusive_scan(ptr_isUn, ptr_isUn + nPerms, ptr_uPref);

    printf("\tPulling unique nodes... ");
    pull_unique<<<numBlocks_perms, dimBlocks_perms>>>(_uniqueNodes, _nodes, 
      _isUnique, nPerms, _uniquePrefix, nUnique);
    printf("Done!\n");
  
    /* FIND REGULAR */
    // Initialize tetrad struct for all unique tetrads
    tetrad_struct *_allTetrads;
    cudaMalloc((void**) &(_allTetrads), sizeof(tetrad_struct) * nUnique);

    // Set up threads, blocks for each tetrad
    int threads_tetrads = MAX_THREADS_1D;
    int blocks_tetrads = (int) ceil((double) nUnique /(double) threads_tetrads);
    if (threads_tetrads > nUnique) {
      threads_tetrads = nUnique;
      blocks_tetrads = 1;
    }
    dim3 dimBlocks_tetrads(threads_tetrads);
    dim3 numBlocks_tetrads(blocks_tetrads);

    // Init isRegular array
    printf("\tFinding regular tetrads... ");
    int *_isRegular;
    cudaMalloc((void**) &(_isRegular), nUnique * sizeof(int));

    // Fill _allTetrads with the correct nodes and init isRegular
    fill_nodes<<<numBlocks_tetrads, dimBlocks_tetrads>>>(_allTetrads,
      _uniqueNodes, _isRegular, nUnique);

    // Tolerance check on all tetrads
    check_tolerances<<<numBlocks_tetrads, dimBlocks_tetrads>>>(_parts, 
      _allTetrads, _dom, _isRegular, nUnique, EVarCutLow, EVarCutHigh,
      shapeCutLow, shapeCutHigh);

    // Find number of tetrads that meet the regularity tolerance
    thrust::device_ptr<int> ptr_isReg(_isRegular);
    nRegular = thrust::reduce(ptr_isReg, ptr_isReg + nUnique);
    printf("Found %d\n", nRegular);

    printf("\tIntializing regular tetrads... ");
    
    // Prefix sum on _isRegular -- will give indices for smaller array
    int *_regularPrefix;
    cudaMalloc((void **) &(_regularPrefix), nUnique * sizeof(int));
    thrust::device_ptr<int> ptr_rPref(_regularPrefix);
    thrust::inclusive_scan(ptr_isReg, ptr_isReg + nUnique, ptr_rPref);

    // Initialize array to hold indices of regular tetrads
    // -- last index is trash for redirecting output
    int *_regularTetrads;
    cudaMalloc((void**) &(_regularTetrads), (nRegular + 1) * sizeof(int));

    // Pull regular tetrads
    pull_regular<<<numBlocks_tetrads, dimBlocks_tetrads>>>(_regularTetrads,
      _isRegular, _regularPrefix, nUnique, nRegular);

    // Set up threads, blocks for each regular tetrad
    int threads_regular = MAX_THREADS_1D;
    int blocks_regular = (int) ceil((double) nRegular/(double) threads_regular);
    if (threads_regular > nRegular) {
      threads_regular = nRegular;
      blocks_regular = 1;
    }
    dim3 dimBlocks_regular(threads_regular);
    dim3 numBlocks_regular(blocks_regular);

    // Alloc new tetrad struct, and pull indices / nodes
    cudaMalloc((void**) &_tetrads, sizeof(tetrad_struct) * nRegular);
    copy_regular<<<numBlocks_regular, dimBlocks_regular>>>(_tetrads, 
      _allTetrads, _regularTetrads, nRegular, _isRegular);

    printf("Done.\n");

    // Free variables
    cudaFree(_partInd);
    cudaFree(_partBin);
    cudaFree(_binStart);
    cudaFree(_binEnd);
    cudaFree(_neighborCount);
    cudaFree(_neighborList);
    cudaFree(_nChoose3);
    cudaFree(_strides);
    cudaFree(_nodes);
    cudaFree(_uniquePrefix);
    cudaFree(_isUnique);
    cudaFree(_uniqueNodes);

    cudaFree(_isRegular);
    cudaFree(_allTetrads);
    cudaFree(_regularPrefix);
    cudaFree(_regularTetrads);
  }
}      

extern "C"
void cuda_tetrad_malloc(void)
{
  // Allocate tetrad struct on host and pull from device
  tetrads = (tetrad_struct*) malloc(nRegular * sizeof(tetrad_struct));
  // Pull tetrads back to host
  cudaMemcpy(tetrads, _tetrads, nRegular * sizeof(tetrad_struct), 
      cudaMemcpyDeviceToHost);

  cudaSetDevice(dev_start);

  cudaMalloc((void**) &(_RoG), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_EVar), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_shape), sizeof(double) * nRegular);

  cudaMalloc((void**) &(_I1), 3 * sizeof(double) * nRegular);
  cudaMalloc((void**) &(_I2), 3 * sizeof(double) * nRegular);
  cudaMalloc((void**) &(_I3), 3 * sizeof(double) * nRegular);
  cudaMalloc((void**) &(_gEigVec), 9 * sizeof(double) * nRegular);
  cudaMalloc((void**) &(_sEigVal), 3 * sizeof(double) * nRegular);
  cudaMalloc((void**) &(_sEigVec), 9 * sizeof(double) * nRegular);
  cudaMalloc((void**) &(_vorticity), 3 * sizeof(double) * nRegular);
  cudaMalloc((void**) &(_S11), 3 * sizeof(double) * nRegular);
  cudaMalloc((void**) &(_S22), 3 * sizeof(double) * nRegular);
  cudaMalloc((void**) &(_S33), 3 * sizeof(double) * nRegular);
  cudaMalloc((void**) &(_vortMag), sizeof(double) * nRegular);
  
  cudaMalloc((void**) &(_gEigVecInit), 9 * sizeof(double) * nRegular);
  cudaMalloc((void**) &(_sEigVecInit), 9 * sizeof(double) * nRegular);


  cudaMalloc((void**) &(_g1_s1), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_g1_s2), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_g1_s3), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_g2_s1), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_g2_s2), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_g2_s3), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_g3_s1), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_g3_s2), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_g3_s3), sizeof(double) * nRegular);

  cudaMalloc((void**) &(_g1_z), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_g2_z), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_g3_z), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_s1_z), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_s2_z), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_s3_z), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_w_z), sizeof(double) * nRegular);

  cudaMalloc((void**) &(_w_g1), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_w_g2), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_w_g3), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_w_s1), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_w_s2), sizeof(double) * nRegular);
  cudaMalloc((void**) &(_w_s3), sizeof(double) * nRegular);
}

void cuda_periodic_flip(void)
{
  // Parallize over tetrads
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

void cuda_tetrad_stats(void)
{ 
  // Matrix tests
  //#ifdef DEBUG
  //  if (tt == 0) {
  //    matrixTests<<<1,1>>>();
  //  }
  //#endif

  // Parallelize over each tetrad
  int threads_tetrads = MAX_THREADS_1D;
  int blocks_tetrads = (int) ceil((double) nRegular / (double) threads_tetrads);
  if (threads_tetrads > nRegular) {
    threads_tetrads = nRegular;
    blocks_tetrads = 1;
  }
  dim3 dimBlocks_tetrads(threads_tetrads);
  dim3 numBlocks_tetrads(blocks_tetrads);

  // Calculate tetrad geometry and velocity measures
  tetrad_geometry<<<numBlocks_tetrads, dimBlocks_tetrads>>>(_parts, _tetrads,
    _dom, _RoG, _EVar, _shape, _I1, _I2, _I3, _gEigVec, _sEigVal, _sEigVec, 
    _vorticity, _S11, _S22, _S33, _vortMag, nRegular, tt);

  // If first timestep, save vectors for later comparison
  if (tt == 0) {
    cudaMemcpy(_gEigVecInit, _gEigVec, 9*sizeof(double)*nRegular,
      cudaMemcpyDeviceToDevice);
    cudaMemcpy(_sEigVecInit, _sEigVec, 9*sizeof(double)*nRegular,
      cudaMemcpyDeviceToDevice);
  }

  // Copy back raw data to host for writing to file
  cudaMemcpy(RoG, _RoG, sizeof(double) * nRegular, cudaMemcpyDeviceToHost);
  cudaMemcpy(EVar, _EVar, sizeof(double) * nRegular, cudaMemcpyDeviceToHost);
  cudaMemcpy(shape, _shape, sizeof(double) * nRegular, cudaMemcpyDeviceToHost);

  cudaMemcpy(I1, _I1, sizeof(double)*nRegular, cudaMemcpyDeviceToHost);
  cudaMemcpy(I2, _I2, sizeof(double)*nRegular, cudaMemcpyDeviceToHost);
  cudaMemcpy(I3, _I3, sizeof(double)*nRegular, cudaMemcpyDeviceToHost);
  cudaMemcpy(gEigVec,_gEigVec,9*sizeof(double)*nRegular,cudaMemcpyDeviceToHost);
  cudaMemcpy(sEigVal,_sEigVal,3*sizeof(double)*nRegular,cudaMemcpyDeviceToHost);
  cudaMemcpy(sEigVec,_sEigVec,9*sizeof(double)*nRegular,cudaMemcpyDeviceToHost);
  cudaMemcpy(vorticity, _vorticity, 3*sizeof(double)*nRegular, 
    cudaMemcpyDeviceToHost);
  cudaMemcpy(S11, _S11, sizeof(double)*nRegular, cudaMemcpyDeviceToHost);
  cudaMemcpy(S22, _S22, sizeof(double)*nRegular, cudaMemcpyDeviceToHost);
  cudaMemcpy(S33, _S33, sizeof(double)*nRegular, cudaMemcpyDeviceToHost);

  // Calculate alignment of vectors
  align_vectors<<<numBlocks_tetrads, dimBlocks_tetrads>>>(_gEigVec, _sEigVec,
    _vorticity, _gEigVecInit, _sEigVecInit, nRegular,
    _g1_s1, _g1_s2, _g1_s3, _g2_s1, _g2_s2, _g2_s3, _g3_s1, _g3_s2, _g3_s3,
    _g1_z, _g2_z, _g3_z, _s1_z, _s2_z, _s3_z, _w_z,
    _w_g1, _w_g2, _w_g3, _w_s1, _w_s2, _w_s3);

  /* Higher order statistics */
  // Calculate higher stat moments of RoG, EVar, Shape, and Sii
  cuda_higher_moments(_RoG, nRegular, m_RoG);
  cuda_higher_moments(_EVar, nRegular, m_EVar);
  cuda_higher_moments(_shape, nRegular, m_Shape);
  cuda_higher_moments(_I1, nRegular, m_I1);
  cuda_higher_moments(_I2, nRegular, m_I2);
  cuda_higher_moments(_I3, nRegular, m_I3);
  cuda_higher_moments(_S11, nRegular, m_S11);
  cuda_higher_moments(_S22, nRegular, m_S22);
  cuda_higher_moments(_S33, nRegular, m_S33);

  // Calculate higher stat moments of alignments
  cuda_higher_moments(_g1_s1, nRegular, m_g1_s1); 
  cuda_higher_moments(_g1_s2, nRegular, m_g1_s2); 
  cuda_higher_moments(_g1_s3, nRegular, m_g1_s3); 
  cuda_higher_moments(_g2_s1, nRegular, m_g2_s1); 
  cuda_higher_moments(_g2_s2, nRegular, m_g2_s2); 
  cuda_higher_moments(_g2_s3, nRegular, m_g2_s3); 
  cuda_higher_moments(_g3_s1, nRegular, m_g3_s1); 
  cuda_higher_moments(_g3_s2, nRegular, m_g3_s2); 
  cuda_higher_moments(_g3_s3, nRegular, m_g3_s3); 

  cuda_higher_moments(_g1_z, nRegular, m_g1_z); 
  cuda_higher_moments(_g2_z, nRegular, m_g2_z); 
  cuda_higher_moments(_g3_z, nRegular, m_g3_z); 
  cuda_higher_moments(_s1_z, nRegular, m_s1_z); 
  cuda_higher_moments(_s2_z, nRegular, m_s2_z); 
  cuda_higher_moments(_s3_z, nRegular, m_s3_z); 
  cuda_higher_moments(_w_z, nRegular, m_w_z);

  cuda_higher_moments(_w_g1, nRegular, m_w_g1); 
  cuda_higher_moments(_w_g2, nRegular, m_w_g2); 
  cuda_higher_moments(_w_g3, nRegular, m_w_g3); 
  cuda_higher_moments(_w_s1, nRegular, m_w_s1); 
  cuda_higher_moments(_w_s2, nRegular, m_w_s2); 
  cuda_higher_moments(_w_s3, nRegular, m_w_s3); 

  cuda_higher_moments(_vortMag, nRegular, m_vortMag);
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
  cudaMalloc((void**) &(_diff), sizeof(double) * nRegular);  // xi - x_bar
  cudaMalloc((void**) &(_diff2), sizeof(double) * nRegular); // (xi - x_bar)^2
  cudaMalloc((void**) &(_skew), sizeof(double) * nRegular);  // (xi - x_bar)^3
  cudaMalloc((void**) &(_kurt), sizeof(double) * nRegular);  // (xi - xbar)^4

  // Parallelize over each tetrad
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
  cudaFree(_tetrads);

  cudaFree(_RoG);
  cudaFree(_EVar);
  cudaFree(_shape);
  cudaFree(_I1);
  cudaFree(_I2);
  cudaFree(_I3);
  cudaFree(_gEigVec);
  cudaFree(_sEigVal);
  cudaFree(_sEigVec);
  cudaFree(_vorticity);
  cudaFree(_vortMag);
  cudaFree(_S11);
  cudaFree(_S22);
  cudaFree(_S33);

  cudaFree(_gEigVecInit);
  cudaFree(_sEigVecInit);


  cudaFree(_g1_s1);
  cudaFree(_g1_s2);
  cudaFree(_g1_s3);
  cudaFree(_g2_s1);
  cudaFree(_g2_s2);
  cudaFree(_g2_s3);
  cudaFree(_g3_s1);
  cudaFree(_g3_s2);
  cudaFree(_g3_s3);

  cudaFree(_g1_z);
  cudaFree(_g2_z);
  cudaFree(_g3_z);
  cudaFree(_s1_z);
  cudaFree(_s2_z);
  cudaFree(_s3_z);
  cudaFree(_w_z);

  cudaFree(_w_g1);
  cudaFree(_w_g2);
  cudaFree(_w_g3);
  cudaFree(_w_s1);
  cudaFree(_w_s2);
  cudaFree(_w_s3);

  cudaDeviceReset();
}
