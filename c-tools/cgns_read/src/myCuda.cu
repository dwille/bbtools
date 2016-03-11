#include "myCuda.h"
//#include "time.h"

#include <cuda.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>

extern "C"
void cuda_dev_malloc(void)
{
  // allocate device memory on device
  cudaSetDevice(dev_start);

  cudaMalloc((void**) &_parts, sizeof(part_struct) * nparts);
  cudaMalloc((void**) &_dom, sizeof(dom_struct));
  cudaMalloc((void**) &_uf, sizeof(double) * dom.Gcc.s3);
  cudaMalloc((void**) &_vf, sizeof(double) * dom.Gcc.s3);
  cudaMalloc((void**) &_wf, sizeof(double) * dom.Gcc.s3);
  cudaMalloc((void**) &_phase, sizeof(int) * dom.Gcc.s3);

}

extern "C"
void cuda_dom_push(void)
{
  cudaSetDevice(dev_start);
  // copy host data to device
  cudaMemcpy(_dom, &dom, sizeof(dom_struct), cudaMemcpyHostToDevice);
}

extern "C"
void cuda_flow_push(void)
{
  cudaSetDevice(dev_start);

  cudaMemcpy(_uf, uf, sizeof(double) * dom.Gcc.s3, cudaMemcpyHostToDevice);
  cudaMemcpy(_vf, vf, sizeof(double) * dom.Gcc.s3, cudaMemcpyHostToDevice);
  cudaMemcpy(_wf, wf, sizeof(double) * dom.Gcc.s3, cudaMemcpyHostToDevice);
  cudaMemcpy(_phase, phase, sizeof(int) *dom.Gcc.s3, cudaMemcpyHostToDevice);
}

extern "C"
void cuda_flow_pull(void)
{
  cudaMemcpy(uf, _uf, sizeof(double) * dom.Gcc.s3, cudaMemcpyDeviceToHost);
  cudaMemcpy(vf, _vf, sizeof(double) * dom.Gcc.s3, cudaMemcpyDeviceToHost);
  cudaMemcpy(wf, _wf, sizeof(double) * dom.Gcc.s3, cudaMemcpyDeviceToHost);
  cudaMemcpy(phase, _phase, sizeof(int) *dom.Gcc.s3, cudaMemcpyDeviceToHost);
}

extern "C"
void cuda_part_push(void)
{
  cudaSetDevice(dev_start);
  cudaMemcpy(_parts, parts, sizeof(part_struct) * nparts, 
    cudaMemcpyHostToDevice);
}

extern "C"
void cuda_part_pull(void)
{
  cudaMemcpy(parts, _parts, sizeof(part_struct) * nparts, 
    cudaMemcpyDeviceToHost);
}

void cuda_phase_averaged_vel(void)
{
  // Paralleize over flow nodes
  int threads = MAX_THREADS_1D;
  int blocks = (int) ceil((double) dom.Gcc.s3 / (double) threads);
  if (threads > dom.Gcc.s3) {
    threads = dom.Gcc.s3;
    blocks = 1;
  }
  dim3 numBlocks(blocks);
  dim3 dimBlocks(threads);

  // create phase mask
  phase_mask<<<numBlocks, dimBlocks>>>(_uf, _vf, _wf, _phase, dom.Gcc.s3);

}

extern "C"
void cuda_dev_free(void)
{
  cudaFree(_dom);
  cudaFree(_uf);
  cudaFree(_vf);
  cudaFree(_wf);
  cudaFree(_phase);
  cudaFree(_parts);

  cudaDeviceReset();
}
