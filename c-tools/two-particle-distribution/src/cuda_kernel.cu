#include "my_cuda.h"

// Fill the particle bins arrays -- partBin and partInd
__global__ void bin_fill(int *partInd, int *partBin, int nparts,
  part_struct *parts, dom_struct *binDom, BC bc) 
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;

  int c;
  int ibin, jbin, kbin;

  // Find the correct bin index for each part and store it
  if (pp < nparts) {
    ibin = floor((parts[pp].x - binDom->xs)/binDom->dx);
    jbin = floor((parts[pp].y - binDom->ys)/binDom->dy);
    kbin = floor((parts[pp].z - binDom->zs)/binDom->dz);
    c = ibin + jbin*binDom->Gcc.s1 + kbin*binDom->Gcc.s2;

    partInd[pp] = pp;         // index of particle
    partBin[pp] = c;          // bin index
    parts[pp].bin = c;        // bin index (stored in particle)
  }
}

__global__ void init(int *array, int length, int val)
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;;
  if (pp < length)  {
    array[pp] =  val;
  }
}

__global__ void bin_start(int *binStart, int *binEnd, int *partBin, int nparts)
{
  // This kernel function was adapted from NVIDIA CUDA 5.5 Examples
  // This software contains source code provided by NVIDIA Corporation
  extern __shared__ int sharedBin[];    //blockSize + 1
  int index = threadIdx.x + blockIdx.x*blockDim.x;
  int bin;

  // for a given bin index, the previous bins's index is stored in sharedBin
  if (index < nparts) {
    bin = partBin[index]; 

    // Load bin data into shared memory so that we can look
    // at neighboring particle's hash value without loading
    // two bin values per thread
    sharedBin[threadIdx.x + 1] = bin;

    if (index > 0 && threadIdx.x == 0) {
      // first thread in block must load neighbor particle bin
      sharedBin[0] = partBin[index - 1];
    }
  }
  __syncthreads();

  if (index < nparts) {
    // If this particle has a different cell index to the previous
    // particle then it must be the first particle in the cell,
    // so store the index of this particle in the cell.
    // As it isn't the first particle, it must also be the cell end of
    // the previous particle's cell
    bin = partBin[index]; 

    if (index == 0 || bin != sharedBin[threadIdx.x]) {
    binStart[bin] = index;

        if (index > 0)
            binEnd[sharedBin[threadIdx.x]] = index;
    }

    if (index == nparts - 1)
    {
        binEnd[bin] = index + 1;
    }
  }
}

__global__ void find_nodes(part_struct *parts, int nparts, dom_struct *dom, 
  BC bc, int *binStart, int *binEnd, int *partBin, int *partInd, 
  dom_struct *binDom, int *neighborList, int *neighborCount, int orderN,
  int orderL)
{
  int index = threadIdx.x + blockIdx.x*blockDim.x;

  if (index < nparts) {
    int i = partInd[index];
    int bin = partBin[index];

    int kbin = floorf(bin/binDom->Gcc.s2);
    int jbin = floorf((bin - kbin*binDom->Gcc.s2)/binDom->Gcc.s1);
    int ibin = bin - kbin*binDom->Gcc.s2 - jbin*binDom->Gcc.s1;

    int l, m, n;                          // adjacent bin iterators
    int target, j;                        // target indices
    int adjBin, adjStart, adjEnd;         // adjacent bin stuff
    int iStride, kStride, jStride;        // how to get to Sesame Street

    // predefine face locations 
    // -1, -2 due to local vs global indexing and defiinition of dom_struct
    int fW = binDom->Gcc.is - 1;
    int fE = binDom->Gcc.ie - 2;
    int fS = binDom->Gcc.js - 1;
    int fN = binDom->Gcc.je - 2;
    int fB = binDom->Gcc.ks - 1;
    int fT = binDom->Gcc.ke - 2;

    // size checks
    int xnBin = (binDom->xn > 2);
    int ynBin = (binDom->yn > 2);
    int znBin = (binDom->zn > 2);

    // particle pair variables
    double xi, yi, zi;
    double xj, yj, zj;
    double rx1, ry1, rz1;
    double rx2, ry2, rz2;
    double rx, ry, rz;
    double r_ij, mu_ij;

    // loop over adjacent bins and take care of periodic conditions 
    for (n = -1; n <= 1; n++) {
      // if on a face and not periodic, continue
      // if on a face and periodic but only 2 bins, continue
      if ((n == -1 && kbin == fB && bc.pB != PERIODIC) || 
          (n ==  1 && kbin == fT && bc.pT != PERIODIC) ||
          (n == -1 && kbin == fB && bc.pB == PERIODIC && znBin == 0) ||
          (n ==  1 && kbin == fT && bc.pT == PERIODIC && znBin == 0)) {
        continue;
      // if on a face and periodic, flip to other side
      } else if (n == -1 && kbin == fB && bc.pB == PERIODIC) {
        kStride = fT*binDom->Gcc.s2;
      } else if (n ==  1 && kbin == fT && bc.pT == PERIODIC) {
        kStride = fB*binDom->Gcc.s2;
      // else, we are in the middle, do nothing special
      } else {
        kStride = (kbin + n)*binDom->Gcc.s2;
      }

      for (m = -1; m <= 1; m++) {
        if ((m == -1 && jbin == fS && bc.pS != PERIODIC) ||
            (m ==  1 && jbin == fN && bc.pN != PERIODIC) ||
            (m == -1 && jbin == fS && bc.pS == PERIODIC && ynBin == 0) ||
            (m ==  1 && jbin == fN && bc.pN == PERIODIC && ynBin == 0)) {
          continue;
        } else if (m == -1 && jbin == fS && bc.pS == PERIODIC) {
          jStride = fN*binDom->Gcc.s1;  
        } else if (m ==  1 && jbin == fN && bc.pN == PERIODIC) {
          jStride = fS*binDom->Gcc.s1;
        } else {
          jStride = (jbin + m)*binDom->Gcc.s1;
        }

        for (l = -1; l <= 1; l++) {
          if ((l == -1 && ibin == fW && bc.pW != PERIODIC) ||
              (l ==  1 && ibin == fE && bc.pE != PERIODIC) ||
              (l == -1 && ibin == fW && bc.pW == PERIODIC && xnBin == 0) ||
              (l ==  1 && ibin == fE && bc.pE == PERIODIC && xnBin == 0)) {
            continue;
          } else if (l == -1 && ibin == fW && bc.pW == PERIODIC) {
            iStride = fE;
          } else if (l ==  1 && ibin == fE && bc.pE == PERIODIC) {
            iStride = fW;
          } else {
            iStride = ibin + l;
          }

          adjBin = iStride + jStride + kStride; 
          adjStart = binStart[adjBin];        // find start and end of bins
          adjEnd = binEnd[adjBin];
          if (adjStart != -1) {               // if bin is not empty
            for (target = adjStart; target < adjEnd; target++) {
              j = partInd[target];
              if (j != i) {                   // if its not original part

                /* Find part separation, check for periodic neighbors */
                // X
                xi = parts[i].x;
                xj = parts[j].x;
                rx = xi - xj;
                // check and correct for separation
                rx1 = xi - (xj + dom->xl);
                rx1 = xi - (xj - dom->xl);
                rx = rx1*(rx1*rx1 < rx*rx) + rx2*(rx2*rx2 < rx*rx);
                rx = (bc.pW == PERIODIC)*rx + (bc.pW != PERIODIC)*(xi - xj);
                
                // Y
                yi = parts[i].y;
                yj = parts[j].y;

                // check and correct for separation
                ry1 = yi - (yj + dom->yl);
                ry1 = yi - (yj - dom->yl);
                ry = ry1*(ry1*ry1 < ry*ry) + ry2*(ry2*ry2 < ry*ry);
                ry = (bc.pS == PERIODIC)*ry + (bc.pS != PERIODIC)*(yi - yj);

                // Z
                zi = parts[i].z;
                zj = parts[j].z;
                rz = zi - zj;
                // check and correct for separation
                rz1 = zi - (zj + dom->zl);
                rz1 = zi - (zj - dom->zl);
                rz = rz1*(rz1*rz1 < rz*rz) + rz2*(rz2*rz2 < rz*rz);
                rz = (bc.pB == PERIODIC)*rz + (bc.pB != PERIODIC)*(zi - zj);

                // corrected separation
                r_ij = sqrt(rx*rx + ry*ry + rz*rz);

                // angle mu = cos(th) = z/r -- TODO: symmetric over pi/2?
                mu_ij = rz/r_ij;

                // Loop over L, N
                for (int enn = 0; enn < orderN; enn++) {
                    kn = enn*PI/L
                  for (int ell = 0; ell < orderL; ell++) {
                    // Calculate coefficients for g_ln
                    // Calculate P_l(mu)
                    P_ell = eval_legendre_poly(mu_ij, ell)
                    // Calculate kernel, p_l(mu)/r^2*sin(kn
                  }
                }

              }
            }
          }
        }
      }
    }
  }
}

__device__ double eval_legendre_poly(double mu_ij, int ell)
{
  // from numerical recipes, page 184
  double p1 = 1.;
  double p2 = 0.;
  for (int j = 0; j < ell; j++) {
    p3 = p2;
    p2 = p1;
    p1 = ((2.*j + 1.)*mu*p2 - j*p3)/(j + 1.);
  }
}
