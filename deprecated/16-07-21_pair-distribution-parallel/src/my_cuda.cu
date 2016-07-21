#include "my_cuda.h"
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

  // set up bins and search for tetrads
  if (nparts < 2) {
    printf("nparts = %d, no tetrads to find.\n", nparts);
    exit(EXIT_FAILURE);
  } else if (nparts >= 2) {

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

    find_nodes<<<numBlocks, dimBlocks>>>(_parts, nparts, _dom, bc, _binStart,
      _binEnd, _partBin, _partInd, _binDom, _neighborList, _neighborCount,
      orderN, orderL);
    

    // Free variables
    cudaFree(_partInd);
    cudaFree(_partBin);
    cudaFree(_binStart);
    cudaFree(_binEnd);
    cudaFree(_neighborCount);
    cudaFree(_neighborList);
  }
}      

extern "C"
void cuda_dev_free(void)
{
  cudaFree(_parts);
  cudaFree(_dom);
  cudaFree(_binDom);

  cudaDeviceReset();
}
