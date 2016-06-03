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

void cuda_phase_averaged_vel(void)
{
  int threads = MAX_THREADS_1D;
  int blocks = (int) ceil((double) dom.Gcc.s3 / (double) threads);
  if (threads > dom.Gcc.s3) {
    threads = dom.Gcc.s3;
    blocks = 1;
  }
  dim3 numBlocks(blocks);
  dim3 dimBlocks(threads);

  // create phase mask, multiply it by velocity
  phase_mask<<<numBlocks, dimBlocks>>>(_uf, _vf, _wf, _phase, dom.Gcc.s3);

  // sum phase to find number of flow nodes
  thrust::device_ptr<int> ptr_phase(_phase);
  int nNodes = thrust::reduce(ptr_phase, ptr_phase + dom.Gcc.s3);

  // sum _u, _v, _w
  thrust::device_ptr<double> ptr_uf(_uf);
  thrust::device_ptr<double> ptr_vf(_vf);
  thrust::device_ptr<double> ptr_wf(_wf);
  double uSum = thrust::reduce(ptr_uf, ptr_uf + dom.Gcc.s3);
  double vSum = thrust::reduce(ptr_vf, ptr_vf + dom.Gcc.s3);
  double wSum = thrust::reduce(ptr_wf, ptr_wf + dom.Gcc.s3);

  cudaMemcpy(phase, _phase, sizeof(int) *dom.Gcc.s3, cudaMemcpyDeviceToHost);

  double invnNodes = 0;
  if (nNodes > 0) {
    invnNodes = 1./nNodes;
  } else if (nNodes <= 0) {
    printf("nNodes <= 0 (%d) -- probable cause is no particles in domain ", 
      nNodes);
    printf("or missing grid.cgns file.\n");
  }

  phaseAvgUf[tt] = uSum*invnNodes;
  phaseAvgVf[tt] = vSum*invnNodes;
  phaseAvgWf[tt] = wSum*invnNodes;
}

extern "C"
void cuda_dev_free(void)
{
  cudaFree(_dom);
  cudaFree(_uf);
  cudaFree(_vf);
  cudaFree(_wf);
  cudaFree(_phase);

  cudaDeviceReset();
}
