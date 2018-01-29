#include <cuda.h>
#include <helper_cuda.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>

#include "cuda_fun.h"


extern "C"
void cuda_set_device(void)
{
  // Should only be one device, so dev_count = 0
  // Should only be one rank, so rank = 0
  // CUDA_VISIBLE_DEVICES should be len(1), so cudaSetDevice[0] will map to the
  //  only device
  checkCudaErrors(cudaSetDevice(0));
}

extern "C"
void cuda_dev_malloc(void)
{
  checkCudaErrors(cudaMalloc(&_parts, nparts * sizeof(part_struct)));
  checkCudaErrors(cudaMalloc(&_dom, sizeof(dom_struct)));

  checkCudaErrors(cudaMalloc(&_x0, nparts * sizeof(double)));
  checkCudaErrors(cudaMalloc(&_y0, nparts * sizeof(double)));
  checkCudaErrors(cudaMalloc(&_z0, nparts * sizeof(double)));
  checkCudaErrors(cudaMalloc(&_x, nparts * sizeof(double)));
  checkCudaErrors(cudaMalloc(&_y, nparts * sizeof(double)));
  checkCudaErrors(cudaMalloc(&_z, nparts * sizeof(double)));
  checkCudaErrors(cudaMalloc(&_r2_total, nparts * sizeof(double)));
  checkCudaErrors(cudaMalloc(&_r2_horiz, nparts * sizeof(double)));
  checkCudaErrors(cudaMalloc(&_r2_verti, nparts * sizeof(double)));
}

extern "C"
void cuda_dom_push(void)
{
  cudaMemcpy(_dom, &dom, sizeof(dom_struct), cudaMemcpyHostToDevice);
}

extern "C"
void cuda_part_push(void)
{
  checkCudaErrors(cudaMemcpy(_parts, parts, nparts * sizeof(part_struct), cudaMemcpyHostToDevice));
}

extern "C"
void cuda_fill_initial_positions(void)
{
  // exec config
  int threads = MAX_THREADS_1D;
  int blocks = (int) ceil((double) nparts / (double) threads);
  if (threads > nparts) {
    threads = nparts;
    blocks = 1;
  }
  dim3 dimBlocks(threads);
  dim3 numBlocks(blocks);

  // Fill device arrays _x0, _y0, _z0 with particle positions at t0
  fill_init_positions<<<numBlocks, dimBlocks>>>(_parts, _x0, _y0, _z0, nparts);
}

extern "C"
void cuda_periodic_flip(void)
{
  // exec config
  int threads = MAX_THREADS_1D;
  int blocks = (int) ceil((double) nparts / (double) threads);
  if (threads > nparts) {
    threads = nparts;
    blocks = 1;
  }
  dim3 dimBlocks(threads);
  dim3 numBlocks(blocks);

  if (tt > 0) {
    // Fix periodicity
    flip_kernel<<<numBlocks, dimBlocks>>>(_parts, _x, _y, _z, _dom, nparts);
  } else if (tt == 0) { 
    // still need to fill _x, _y, _z
    cudaMemcpy(_x, _x0, nparts * sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(_y, _y0, nparts * sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(_z, _z0, nparts * sizeof(double), cudaMemcpyDeviceToDevice);
  }
}

extern "C"
void cuda_find_separations(void)
{
  // exec config
  int threads = MAX_THREADS_1D;
  int blocks = (int) ceil((double) nparts / (double) threads);
  if (threads > nparts) {
    threads = nparts;
    blocks = 1;
  }
  dim3 dimBlocks(threads);
  dim3 numBlocks(blocks);

  // find separations r_total, r_horiz, r_vert
  find_separation<<<numBlocks, dimBlocks>>>(_x, _y, _z, _x0, _y0, _z0, 
    _r2_total, _r2_horiz, _r2_verti, nparts);

  cudaDeviceSynchronize();

  // get pointers to dev arrays
  thrust::device_ptr<double> t_r2_total(_r2_total);
  thrust::device_ptr<double> t_r2_verti(_r2_verti);
  thrust::device_ptr<double> t_r2_horiz(_r2_horiz);

  // average separations over nparts
  r2_total[tt] = thrust::reduce(t_r2_total, t_r2_total + nparts, 0., 
                                  thrust::plus<double>());
  r2_verti[tt] = thrust::reduce(t_r2_verti, t_r2_verti + nparts, 0., 
                                  thrust::plus<double>());
  r2_horiz[tt] = thrust::reduce(t_r2_horiz, t_r2_horiz + nparts, 0., 
                                  thrust::plus<double>());

  // normalize
  double inparts = 1./nparts;
  r2_total[tt] *= inparts;
  r2_verti[tt] *= inparts;
  r2_horiz[tt] *= inparts;

}

extern "C"
void cuda_free(void)
{
  cudaFree(_parts);
  cudaFree(_dom);
  cudaFree(_x0);
  cudaFree(_y0);
  cudaFree(_z0);
  cudaFree(_r2_total);
  cudaFree(_r2_horiz);
  cudaFree(_r2_verti);
}
