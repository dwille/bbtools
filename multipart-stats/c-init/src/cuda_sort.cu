#include "cuda_sort.h"
#include "time.h"

#include <cuda.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/sort.h>
#include <thrust/scan.h>
//#include <thrust/set_operations.h>


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
    printf("\tCombining nodes...\n");
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

    printf("\tSorting permutations...\n");
    sort_combos<<<numBlocks_perms, dimBlocks_perms>>>(_nodes, nPerms);

    /*  FIND_UNIQUE */
    // compare and find unique ones
    int *_isUnique;
    cudaMalloc((void **) &_isUnique, nPerms*sizeof(int));
    init<<<numBlocks_perms, dimBlocks_perms>>>(_isUnique, nPerms, 1);

    // Loop over each permutations, then parallelize over the remaining
    printf("\tLooping over permutations and finding unique sets...");
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
    nUnique = thrust::reduce(ptr_isUnique, ptr_isUnique + nPerms);

    printf("Found %d\n", nUnique);

    /*  PULL_UNIQUE */
    // Last entry is trash for finding indices and redirecting
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

    printf("\tPulling unique nodes...\n");
    pull_unique<<<numBlocks_perms, dimBlocks_perms>>>(_uniqueNodes, _nodes, 
      _isUnique, nPerms, _uniquePrefix, nUnique);
  
    uniqueNodes = (int*) malloc(4*nUnique*sizeof(int));
    cudaMemcpy(uniqueNodes, _uniqueNodes, 4*nUnique*sizeof(int),
      cudaMemcpyDeviceToHost);

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
  }
}      

extern "C"
void cuda_tetrad_malloc(void)
{
  // Allocate tetrad struct on host, then device
  tetrads = (tetrad_struct*) malloc(nUnique * sizeof(tetrad_struct));


  cudaSetDevice(dev_start);
  cudaMalloc((void**) &(_tetrads), sizeof(tetrad_struct) * nUnique);
}

void cuda_tetrad_init(void)
{
  // Set up threads, blocks for each tetrad
  int threads_tetrads = MAX_THREADS_1D;
  int blocks_tetrads = (int) ceil((double) nUnique / (double) threads_tetrads);
  if (threads_tetrads > nUnique) {
    threads_tetrads = nUnique;
    blocks_tetrads = 1;
  }
  dim3 dimBlocks_tetrads(threads_tetrads);
  dim3 numBlocks_tetrads(blocks_tetrads);

  // Fill _tetrads with the correct nodes and init the rest
  fill_nodes<<<numBlocks_tetrads, dimBlocks_tetrads>>>(_tetrads, _uniqueNodes,
    nUnique);

  // Pull tetrads back to host
  cudaMemcpy(tetrads, _tetrads, nUnique * sizeof(tetrad_struct), 
    cudaMemcpyDeviceToHost);
}

void cuda_tetrad_stats(void)
{ 
  // Matrix tests
  #ifdef DEBUG
    matrixTests<<<1,1>>>();
  #endif
  // Parallelize over each tetrad
  int threads_tetrads = MAX_THREADS_1D;
  int blocks_tetrads = (int) ceil((double) nUnique / (double) threads_tetrads);
  if (threads_tetrads > nUnique) {
    threads_tetrads = nUnique;
    blocks_tetrads = 1;
  }
  dim3 dimBlocks_tetrads(threads_tetrads);
  dim3 numBlocks_tetrads(blocks_tetrads);

  tetrad_geometry<<<numBlocks_tetrads, dimBlocks_tetrads>>>(_parts, _tetrads,
    _dom, nUnique);

  // Copy back to host for writing to file
  cudaMemcpy(tetrads, _tetrads, nUnique * sizeof(tetrad_struct), 
    cudaMemcpyDeviceToHost);

}

extern "C"
void cuda_dev_free(void)
{
  cudaFree(_parts);
  cudaFree(_dom);
  cudaFree(_binDom);
  cudaFree(_uniqueNodes);
  cudaFree(_tetrads);

  cudaDeviceReset();
}
