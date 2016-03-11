#include "cuda_sort.h"

// Fill the particle bins arrays -- partBin and partInd
__global__ void bin_fill(int *partInd, int *partBin, int nparts,
  part_struct *parts, dom_struct *binDom, BC bc) 
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;;

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
  dom_struct *binDom, int *neighborList, int *neighborCount, int nMax)
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

    int cc = 0;                           // counter for inner neighborlist

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

                // Add particle j to neighborList
                int nStride = i*nMax + cc;    // increment a counter
                neighborList[nStride] = j;
                neighborCount[i]++;
                cc++;
              }
            }
          }
        }
      }
    }
  }
}

__global__ void choose3(int *neighborCount, int *nChoose3, int nparts)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;

  if (i < nparts) {
    int c = neighborCount[i];
    nChoose3[i] = c*(c-1)*(c-2)/6;
  }
}

__global__ void combine_nodes(int *neighborList, int *neighborCount,
  int *nodes, int *strides, int nparts, int nMax)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int N1, N2, N3, N4;
  int countPerm = 0;
  int quartStride;


  if (i < nparts) {
    N1 = i;
    for (int j = 0; j < neighborCount[i]; j++) {
      N2 = neighborList[i*nMax + j];
      for (int k = j + 1; k < neighborCount[i]; k++) {
        N3 = neighborList[i*nMax + k];
        for (int l = k + 1; l < neighborCount[i]; l++) {
          N4 = neighborList[i*nMax + l];

          quartStride = 4*strides[i] + 4*countPerm;
          countPerm++;

          nodes[quartStride] = N1;
          nodes[quartStride + 1] = N2;
          nodes[quartStride + 2] = N3;
          nodes[quartStride + 3] = N4;
        }
      }
    }
  }
}

__global__ void sort_combos(int *nodes, int nPerms)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int startInd;
  int vals[4];

  if (i < nPerms) {
    startInd = 4*i; 
    for (int v = 0; v < 4; v++) {
      vals[v] = nodes[startInd + v];
    }
    // Adapated from stackoverflow.com 2786899, 
    //  "Fastest sort of fixed length 6 int array"
    // and "pages.ripco.net/~jgamble/nw.html"
    // and perhaps "Bit Twiddling Hacks" by Sean Eron Anderson @ Stanford
    #define min(x, y) (x*(x < y) + y*(y < x) + x*(x == y))
    #define max(x, y) (x*(x > y) + y*(y > x) + y*(x == y))
    #define SWAP(x,y) { const int a = min(vals[x], vals[y]); \
                        const int b = max(vals[x], vals[y]); \
                        vals[x] = a; vals[y] = b; }
    SWAP(0,1);
    SWAP(2,3);
    SWAP(0,2);
    SWAP(1,3);
    SWAP(1,2);

    #undef SWAP
    #undef min
    #undef max

    for (int v = 0; v < 4; v++) {
      nodes[startInd + v] = vals[v];
    }
  }
}

__global__ void find_unique(int *nodes, int base, int remainder, int *isUnique)
{
  int threadInd = threadIdx.x + blockIdx.x*blockDim.x;

  // only want to look at OTHER sets of nodes, so add (n+1)
  int target = threadInd + (base + 1);// Target set
  int targetNode = 4*target;          // First node of target set

  int baseNode = 4*base;              // First node of base set
  int nSame = 0;

  if (isUnique[base] == 0) {    // if base is not unique, skip
      return;
  }
  __syncthreads();

  if (threadInd < remainder) {
    nSame = (nodes[baseNode] == nodes[targetNode]) 
          + (nodes[baseNode + 1] == nodes[targetNode + 1])
          + (nodes[baseNode + 2] == nodes[targetNode + 2]) 
          + (nodes[baseNode + 3] == nodes[targetNode + 3]);

    // set isUnique to 1 if nSame != 4, else set to zero
    // if isUnique == 1
    //  if nSame != 4 --> isUnique = 1*1 = 1
    //  if nSame == 4 --> isUnique = 0*1 = 0
    // if isUnique == 0
    //  if nSame != 4 --> isUnique = 1*0 = 0
    //  if nSame == 4 --> isUnique = 0*0 = 0
    isUnique[target] = (nSame != 4)*isUnique[target];
  }
}

__global__ void find_unique2(int *nodes, int *isUnique, int nPerms) {
  int base = threadIdx.x + blockIdx.x*blockDim.x;

  if (base < nPerms) {
    int baseNode = 4*base;
    int target, targetNode;
    int nSame = 0;
    for (target = 0; target < nPerms; target++) {
      targetNode = 4*target;

      nSame = (nodes[baseNode] == nodes[targetNode]) 
            + (nodes[baseNode + 1] == nodes[targetNode + 1])
            + (nodes[baseNode + 2] == nodes[targetNode + 2]) 
            + (nodes[baseNode + 3] == nodes[targetNode + 3]);
      // if target == base, set nSame to 1
      nSame = (target == base) + nSame*(target != base);

      // set isUnique to 1 if nSame != 4, else set to zero
      // if isUnique == 1
      //  if nSame != 4 --> isUnique = 1*1 = 1
      //  if nSame == 4 --> isUnique = 0*1 = 0
      // if isUnique == 0
      //  if nSame != 4 --> isUnique = 1*0 = 0
      //  if nSame == 4 --> isUnique = 0*0 = 0
      isUnique[base] *= (nSame != 4);
       
    }
  }
}

__global__ void pull_unique(int *uniqueNodes, int *nodes, int *isUnique, 
  int nPerms, int *uniquePrefix, int nUnique)
{
  int TID = threadIdx.x + blockIdx.x*blockDim.x;
  int ind;

  if (TID < nPerms) {
    // becomes uniquePrefix[TID] - 1 iff is unique
    // becomes nUnique iff is !unique
    ind = (uniquePrefix[TID] - 1) * isUnique[TID] 
        + nUnique*(1 - isUnique[TID]);

    uniqueNodes[4*ind] = nodes[4*TID];
    uniqueNodes[4*ind + 1] = nodes[4*TID + 1];
    uniqueNodes[4*ind + 2] = nodes[4*TID + 2];
    uniqueNodes[4*ind + 3] = nodes[4*TID + 3];

  }
}

__global__ void fill_nodes(tetrad_struct *tetrads, int *uniqueNodes, 
  int nUnique)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;

  if (i < nUnique) {
    tetrads[i].N1 = uniqueNodes[4*i];    
    tetrads[i].N2 = uniqueNodes[4*i + 1];    
    tetrads[i].N3 = uniqueNodes[4*i + 2];    
    tetrads[i].N4 = uniqueNodes[4*i + 3];    
    tetrads[i].tolCheck = 0;
  }
}

__global__ void tetrad_geometry(part_struct *parts, tetrad_struct *tetrads,
  dom_struct *dom, int nUnique)
{
  int tet = threadIdx.x + blockIdx.x*blockDim.x;

  /* Tetrahedron Geometry Variables */
  double XCM = 0;       // Tetrad center of mass -- x
  double YCM = 0;       // Tetrad center of mass -- y
  double ZCM = 0;       // Tetrad center of mass -- x
  double r1[nDim];      // Node 1 relative coordinates
  double r2[nDim];      // Node 2 relative coordinates
  double r3[nDim];      // Node 3 relative coordinates
  double r4[nDim];      // Node 4 relative coordinates

  /* Shape Tensor Variables */
  double g[nDim2];      // Gyration tensor
  double avgLambda;     // Average eigenvalue of g = trace(g)/3 
  double g_hat[nDim2];  // Deviatoric part of g

  /* Velocity variables */

  double gInv[nDim2];   // Gyration tensor inverse
  int nrot;

  double UCM = 0;       // Tetrad center of vel -- u
  double VCM = 0;       // Tetrad center of vel -- v
  double WCM = 0;       // Tetrad center of vel -- w
  double u1[nDim];      // Node 1 relative vel
  double u2[nDim];      // Node 2 relative vel
  double u3[nDim];      // Node 3 relative vel
  double u4[nDim];      // Node 4 relative vel
  double W[nDim2];      // Velocity tensor

  double M[nDim2];      // Coarse-grained vel grad tensor
  double S[nDim2];      // Symmetric part of M
  double O[nDim2];      // Anti-Symmetric part of M

  if (tet < nUnique) {
    /*  POSITION  */
    // Fix periodicity issues
    periodic_flip(r1, r2, r3, r4, tetrads[tet], parts, dom->xl, dom->yl, 
      dom->zl);

    // Calculate tetrad center of mass
    // reference all of them to N1, if > dom.size, flip it
    XCM = 0.25*(r1[0] + r2[0] + r3[0] + r4[0]);
    YCM = 0.25*(r1[1] + r2[1] + r3[1] + r4[1]);
    ZCM = 0.25*(r1[2] + r2[2] + r3[2] + r4[2]);

    // Relate nodes to center of mass
    r1[0] -= XCM;
    r1[1] -= YCM;
    r1[2] -= ZCM;
            
    r2[0] -= XCM;
    r2[1] -= YCM;
    r2[2] -= ZCM;
            
    r3[0] -= XCM;
    r3[1] -= YCM;
    r3[2] -= ZCM;
            
    r4[0] -= XCM;
    r4[1] -= YCM;
    r4[2] -= ZCM;

    // Gyration tensor
    for (int i = 0; i < nDim; i++) {
      for (int j = 0; j < nDim; j++) {
        g[nDim*i + j] = 0.25*(r1[i]*r1[j] + r2[i]*r2[j] 
                            + r3[i]*r3[j] + r4[i]*r4[j]);
      }
    }
    tetrads[tet].R2 = matrixTrace3(g);

    // Calculate average eigenvalue of g
    avgLambda = tetrads[tet].R2/3.;

    // Calculate deviatoric part of g, g_hat
    for (int i = 0; i < nDim; i++) {
      for (int j = 0; j < nDim; j++) {
        int c = nDim*i + j;
        g_hat[c] = g[c] - avgLambda*(i == j);
      }
    }

    // Calculate variance of g's eigenvalues, var = trace(g_hat^2)
    tetrads[tet].var = matrixSquaredTrace3(g_hat);

    // Calculate shape of g using det(g_hat)
    tetrads[tet].det = matrixDet3(g_hat);

    // Calculate I1, I2, I3 and principal directions of shape tensor
    jacobiEig3(g, tetrads[tet].gEigVal, tetrads[tet].gEigVec, &nrot);

    /* Velocity */
    // Reinit g since it was overwritten in last step
    for (int i = 0; i < nDim; i++) {
      for (int j = 0; j < nDim; j++) {
        g[nDim*i + j] = (r1[i]*r1[j] + r2[i]*r2[j] 
                       + r3[i]*r3[j] + r4[i]*r4[j]);
      }
    }

    // find gInv
    matrixInverse3(g, gInv);

    // Calculate tetrad center of vel
    int N1 = tetrads[tet].N1;
    int N2 = tetrads[tet].N2;
    int N3 = tetrads[tet].N3;
    int N4 = tetrads[tet].N4;

    UCM = 0.25*(parts[N1].u + parts[N2].u + parts[N3].u + parts[N4].u);
    VCM = 0.25*(parts[N1].v + parts[N2].v + parts[N3].v + parts[N4].v);
    WCM = 0.25*(parts[N1].w + parts[N2].w + parts[N3].w + parts[N4].w);

    // Relate nodes to center of vel
    u1[0] = parts[N1].u - UCM;
    u1[1] = parts[N1].v - VCM;
    u1[2] = parts[N1].w - WCM;

    u2[0] = parts[N2].u - UCM;
    u2[1] = parts[N2].v - VCM;
    u2[2] = parts[N2].w - WCM;

    u3[0] = parts[N3].u - UCM;
    u3[1] = parts[N3].v - VCM;
    u3[2] = parts[N3].w - WCM;

    u4[0] = parts[N4].u - UCM;
    u4[1] = parts[N4].v - VCM;
    u4[2] = parts[N4].w - WCM;

    // Vel tensor
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        W[3*i + j] = r1[i]*u1[j] + r2[i]*u2[j] + r3[i]*u3[j] + r4[i]*u4[j];
      }
    }

    // Velocity Gradient Tensor
    matrixMult3(gInv, W, M);

    // Decompose M
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        S[3*i + j] = 0.5*(M[3*i + j] + M[3*j + i]);
        O[3*i + j] = 0.5*(M[3*i + j] - M[3*j + i]);
      }
    }
    // Find princpal directions and values of strain tensor
    jacobiEig3(S, tetrads[tet].sEigVal, tetrads[tet].sEigVec, &nrot);

    // pull vorticity vector:
    // See AP, Fluid Dynamics, pages 38 + 40 -- O is -(7.53) bc (8.9)
    // TODO: need to check
    tetrads[tet].vorticity[0] = 2*O[3];
    tetrads[tet].vorticity[1] = 2*O[2];
    tetrads[tet].vorticity[2] = 2*O[7];

  }
}

__global__ void check_tolerances(tetrad_struct *tetrads, double varCutLow,
  double varCutHigh, double shapeCutLow, double shapeCutHigh, int nUnique)
{
  int tet = threadIdx.x + blockIdx.x*blockDim.x;
  if (tet < nUnique) {
    double eigVar = 1.5*tetrads[tet].var/(tetrads[tet].R2*tetrads[tet].R2);
    double shape = 27.*tetrads[tet].det/(tetrads[tet].R2*tetrads[tet].R2*
                                          tetrads[tet].R2);
    tetrads[tet].tolCheck = (eigVar >= varCutLow &&
                             eigVar <= varCutHigh &&
                             shape >= shapeCutLow &&
                             shape <= shapeCutHigh);
  }
}

__global__ void matrixTests(void)
{
  // TEST MATRIX FUNCTIONS
    double A[9];
    double a_in[9];
    double B[9];
    double R[9];
    double invA[9];
    double d[3];
    double v[9];
    int nrot = 0;

  // INITIALIZE MATRICES
    A[0] = 3.; A[1] = 3.; A[2] = 5.; 
    A[3] = 3.; A[4] = 4.; A[5] = 6.;
    A[6] = 5.; A[7] = 6.; A[8] = -6.;

    B[0] = 2; B[1] = 6; B[2] = 1; 
    B[3] = 2; B[4] = 6; B[5] = 1;
    B[6] = 2; B[7] = 6; B[8] = 1;

    R[0] = 0.; R[1] = 0.; R[2] = 0.; 
    R[3] = 0.; R[4] = 0.; R[5] = 0.;
    R[7] = 0.; R[7] = 0.; R[8] = 0.;

    invA[0] = 0.; invA[3] = 0.; invA[6] = 0.; 
    invA[1] = 0.; invA[4] = 0.; invA[7] = 0.;
    invA[2] = 0.; invA[5] = 0.; invA[8] = 0.;

    d[0] = 0.; d[1] = 0.; d[2] = 0.; 

    v[0] = 0.; v[1] = 0.; v[2] = 0.; 
    v[1] = 0.; v[4] = 0.; v[5] = 0.;
    v[2] = 0.; v[7] = 0.; v[8] = 0.;
    

    // Copy input matrix
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        a_in[3*i + j] = A[3*i + j];
      }
    }

    matrixInverse3(A, invA);
    matrixMult3(A, B, R);

    jacobiEig3(a_in, d, v, &nrot);

    printf("A = %lf, %lf, %lf\n\
    %lf, %lf, %lf,\n\
    %lf, %lf, %lf\n", A[0], A[1], A[2],
         A[3], A[4], A[5], A[6], A[7], A[8]);

    printf("invA = %lf, %lf, %lf\n\
       %lf, %lf, %lf,\n\
       %lf, %lf, %lf\n", invA[0], invA[1], invA[2],
         invA[3], invA[4], invA[5], invA[6], invA[7], invA[8]);

    printf("R = %lf, %lf, %lf\n\
    %lf, %lf, %lf,\n\
    %lf, %lf, %lf\n", R[0], R[1], R[2],
         R[3], R[4], R[5], R[6], R[7], R[8]);

    printf("d = %lf, %lf, %lf; nrot = %d\n", d[0], d[1], d[2], nrot);

    printf("v = %lf, %lf, %lf\n\
    %lf, %lf, %lf,\n\
    %lf, %lf, %lf\n", v[0], v[1], v[2],
                      v[3], v[4], v[5], 
                      v[6], v[7], v[8]);
}

__device__ void periodic_flip(double *r1, double *r2, double *r3, double *r4, 
  tetrad_struct tetrads, part_struct *parts, double xl, double yl, double zl)
{
    int N1 = tetrads.N1;
    int N2 = tetrads.N2;
    int N3 = tetrads.N3;
    int N4 = tetrads.N4;
    // branchless min
    #define flip(s1,s2,l,i) \
      {const double standard = s1[i] - s2[i]; \
       const double standard2 = standard*standard; \
       const double var1 = s1[i] - (s2[i] + l); \
       const double var2 = s1[i] - (s2[i] - l); \
       s2[i] += l*((var1*var1 < standard2) - (var2*var2 < standard2)); }

    // Set up position vectors to each particle from origin
    r1[0] = parts[N1].x;
    r1[1] = parts[N1].y;
    r1[2] = parts[N1].z;

    r2[0] = parts[N2].x;
    r2[1] = parts[N2].y;
    r2[2] = parts[N2].z;

    r3[0] = parts[N3].x;
    r3[1] = parts[N3].y;
    r3[2] = parts[N3].z;

    r4[0] = parts[N4].x;
    r4[1] = parts[N4].y;
    r4[2] = parts[N4].z;

    /* X direction */
    flip(r1, r2, xl, 0);
    flip(r1, r3, xl, 0);
    flip(r1, r4, xl, 0);

    /* Y direction */
    flip(r1, r2, yl, 1);
    flip(r1, r3, yl, 1);
    flip(r1, r4, yl, 1);

    /* Z direction */
    flip(r1, r2, zl, 2);
    flip(r1, r3, zl, 2);
    flip(r1, r4, zl, 2);

  #undef flip
}

__device__ double matrixDet3(double *A)
{
  double detA;
  detA = A[0]*(A[8]*A[4] - A[5]*A[7])
       - A[3]*(A[8]*A[1] - A[2]*A[7])
       + A[6]*(A[5]*A[1] - A[2]*A[4]);

  return detA;
}

__device__ double matrixTrace3(double *A)
{
  double trace = 0;
  for (int i = 0; i < 3; i++) {
    trace += A[nDim*i + i];
  }
  return trace;
}

__device__ double matrixSquaredTrace3(double *A)
{
  // computes trace(A^2) by first calculating A^2
  double A2[9] = {0,0,0,0,0,0,0,0,0};

  matrixMult3(A,A,A2);

  double trace = matrixTrace3(A2);

  return trace;
}

__device__ void matrixInverse3(double *A, double *invA)
{
  //     | A11 A12 A13 |   | A[0] A[1] A[2] |
  // A = | A21 A22 A23 | = | A[3] A[4] A[5] |
  //     | A31 A32 A33 |   | A[6] A[7] A[8] |
  double detA = matrixDet3(A);

  double invDETA = 1/detA;             

  // If it's bigger than this, the matrix is probably singular
  // 1 if okay, 0 if not
  int detCheck = (invDETA < 1e10);

  invA[0] = A[8]*A[4] - A[5]*A[7];
  invA[1] = -(A[8]*A[3] - A[6]*A[5]);
  invA[2] = A[7]*A[3] - A[6]*A[4];
  invA[3] = -(A[8]*A[1] - A[7]*A[2]);
  invA[4] = A[8]*A[0] - A[2]*A[6];
  invA[5] = -(A[7]*A[0] - A[1]*A[6]);
  invA[6] = A[5]*A[1] - A[4]*A[2];
  invA[7] = -(A[5]*A[0] - A[3]*A[2]);
  invA[8] = A[4]*A[0] - A[1]*A[3];

  for (int i = 0; i < 9; i++) {
    invA[i] *= detCheck*invDETA;
  }
}

__device__ void matrixMult3(double *A, double *B, double *R)
{
  // R = A*B
  // Rij = Aim*Bmj
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      R[3*i + j] = 0;
      for (int m = 0; m < 3; m++) {
        R[3*i + j] += A[3*i + m]*B[3*m + j];
      }
    }
  }
}

__device__ void jacobiEig3(double *a, double *d, double *v, int *nrot)
{
  // Modified from Numerical Recipes
  int i,j,ip,iq;
  double thresh, theta, tau, t, sm, s, h, g, c;
  int n = 3;
  double b[3], z[3];

  // Initialize v to identity
  for (ip = 0; ip < n; ip++) {
    for (iq = 0; iq < n; iq++) {
      v[3*ip + iq] = (ip == iq);
    }
  }

  // Intialize b,d to diagonal of a; z to 0
  for (ip = 0; ip < n; ip++) {
    b[ip] = a[3*ip + ip];
    d[ip] = a[3*ip + ip];
    z[ip] = 0.;
  }

  // Main loop
  for (i = 1; i < 50; i++) {
    sm = 0.0;
    // Sum magnitude of off-diagonal elements of a
    for (ip = 0; ip < n-1; ip++) {
      for (iq = ip + 1; iq < n; iq++) {
        sm += abs(a[3*ip + iq]);
      }
    }

    // Normal return, relies on quadratic convergence to machine underflow
    if (sm == 0.0) {
      eigsrt(d,v);
      return;
    }

    if (i < 4) {
      thresh = 0.2*sm/(n*n);    // On first 3 sweeps...
    } else {
      thresh = 0.;               // ...thereafter
    }

    for (ip = 0; ip < n - 1; ip++) {
      for (iq = ip + 1; iq < n; iq++) {
        g = 100.0*abs(a[3*ip + iq]);

        // After 4 sweeps, skip the rotation if the off-diagonal element is sm.
        if (i > 4 && g <= DBL_EPSILON*abs(d[ip]) && g <= DBL_EPSILON*abs(d[iq])) {
          a[3*ip + iq] = 0; 
        } else if (abs(a[3*ip + iq]) > thresh) {
          h = d[iq] - d[ip];
          if (g <= DBL_EPSILON*abs(h)) {
            t = a[3*ip + iq]/h;   // t = 1/(2theta)
          } else {
            theta = 0.5*h/a[3*ip + iq];
            t = 1.0/(abs(theta) + sqrt(1.0 + theta*theta));
            if (theta < 0.0) {
              t = -t;
            }
          }
          c = 1.0/sqrt(1.0 + t*t);
          s = t*c;
          tau = s/(1.0 + c);
          h = t*a[3*ip + iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[3*ip + iq] = 0.0;
          for (j = 0; j < ip; j++) {  // Case 0 <= j < p
            rot(a,s,tau,j,ip,j,iq);
          }
          for (j = ip + 1; j < iq; j++) {  // Case p < j < q
            rot(a,s,tau,ip,j,j,iq);
          }
          for (j = iq + 1; j < n; j++) {  // Case q < j < n
            rot(a,s,tau,ip,j,iq,j);
          }
          for (j = 0; j < n; j++) {  // Case 0 <= j < p
            rot(v,s,tau,j,ip,j,iq);
          }
          (*nrot)++;
        }
      }
    }

    // Update d with sum of t*a_pq and reinit z
    for (ip = 0; ip < n; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
//  //printf("Too many iterations in routine jacobi\n");
//  //return EXIT_FAILURE;
}

// rotation kernel
__device__ void rot(double *a, double s, double tau, int i, int j, int k, 
  int l)
{
  double g = a[3*i + j];
  double h = a[3*k + l];
  a[3*i + j] = g - s*(h + g*tau);
  a[3*k + l] = h + s*(g - h*tau);
}

// sort eigen values and eigenvectors into descending order
__device__ void eigsrt(double *d, double *v)
{
  int k;
  int n = 3;
  double p;

  // Loop over all eigenvalues (find max)
  for (int i = 0; i < n - 1; i++) {
    k = i;      // Index of currEig
    p = d[k];   // Value of currEig

    // Loop over all other eigenvalues
    for (int j = i; j < n; j++) {

      // If targetEig > currEig
      if (d[j] >= p) {
        k = j;      // Index of targetEig
        p = d[k];   // Value of targetEig
      }
    }

    // If we've found a targetEig > currEig
    if (k != i) {
      d[k] = d[i];  // Set value of targetEig to currEig
      d[i] = p;     // Set value of currEig to targetEig

      // Arrange eigenvectors
      for (int j = 0; j < n; j++) {
        p = v[3*j + i];
        v[3*j + i] = v[3*j + k];
        v[3*j + k] = p;
      }
    }
  }
}
