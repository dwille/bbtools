#include "cuda_fun.h"

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
  dom_struct *binDom, int *partI, int *partJ, int *keepIJ, int *iMatches,
  int *initFlipCount, int nMax, double R0_a, double eps_a, double rmax)
{
  int index = threadIdx.x + blockIdx.x*blockDim.x;

  if (index < nparts) {
    int i = partInd[index];
    /* NBODY variables */
    int bin = partBin[index];

    int kbin = floorf(bin/binDom->Gcc.s2);
    int jbin = floorf((bin - kbin*binDom->Gcc.s2)/binDom->Gcc.s1);
    int ibin = bin - kbin*binDom->Gcc.s2 - jbin*binDom->Gcc.s1;

    int l, m, n;                          // adjacent bin iterators
    int target, j;                        // target indices
    int adjBin, adjStart, adjEnd;         // adjacent bin stuff
    int iStride, kStride, jStride;        // how to get to Sesame Street

    /* Particle pair variables */
    // Initialize all
    for (int aa = 0; aa < nMax; aa++) {
      partI[aa + nMax*i] = i;             // fill partI with i
      partJ[aa + nMax*i] = -1;
      keepIJ[aa + nMax*i] = 0;
      initFlipCount[0 + 3*aa + nMax*i] = 0;
      initFlipCount[1 + 3*aa + nMax*i] = 0;
      initFlipCount[2 + 3*aa + nMax*i] = 0;
    }
    iMatches[i] = 0;
    int cc = 0;                              // counter for filling arrays

    double ri[3], rj[3];                     // i, j positions
    double rx, ry, rz, r2;                   // i, j separations
    double minR0 = rmax*rmax*(R0_a - eps_a)*(R0_a - eps_a); // tolerance
    double maxR0 = rmax*rmax*(R0_a + eps_a)*(R0_a + eps_a); // tolerance
    int inRange;                             // meets tolerance
    
    // tmp variables for init flip
    int flipX = 0., flipY = 0., flipZ = 0.;  // tmp holder for init flip
    double dx, dx2, flipL, flipR;
    int flipFlag;

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
              
              // check for neighbors when using periodic boundaries
              #define flip(s1,s2,l,i,count) \
                {dx = s1[i] - s2[i]; \
                 dx2 = dx*dx; \
                 flipL = s1[i] - (s2[i] + l); \
                 flipR = s1[i] - (s2[i] - l); \
                 flipFlag = (flipL*flipL < dx2) - (flipR*flipR < dx2); \
                 s2[i] += l*(flipFlag); \
                 count += flipFlag; }

              // Pull positions
              ri[0] = parts[i].x; ri[1] = parts[i].y; ri[2] = parts[i].z;
              rj[0] = parts[j].x; rj[1] = parts[j].y; rj[2] = parts[j].z;

              // Determine if particles are over periodic boundary
              // increment appropriate tmp counter if so
              flip(ri, rj, dom->xl, 0, flipX);
              flip(ri, rj, dom->yl, 1, flipY);
              flip(ri, rj, dom->zl, 2, flipZ);

              // separation in x,y,z and r2
              rx = ri[0] - rj[0];
              ry = ri[1] - rj[1];
              rz = ri[2] - rj[2];
              r2 = rx*rx + ry*ry + rz*rz;

              // If particles are different and separation is in range
              inRange = (r2 < maxR0) && (r2 > minR0);
              if (j != i && inRange == 1) {
                int stride = i*nMax + cc;    // counter for partI,J arrays

                // Add particle j to partI/partJ -- sort
                partI[stride] = i*(i < j) + j*(j < i);
                partJ[stride] = i*(i > j) + j*(j > i);
                keepIJ[stride] = 1;
                iMatches[i] += 1;

                // account for periodicity
                int sX = 0 + 3*cc + i*nMax;
                int sY = 1 + 3*cc + i*nMax;
                int sZ = 2 + 3*cc + i*nMax;
                initFlipCount[sX] += flipX;
                initFlipCount[sY] += flipY;
                initFlipCount[sZ] += flipZ;

                cc++;
              }
              #undef flip
            }
          }
        }
      }
    }
  }
}

__global__ void find_strides(int *keepIJ, int *strides, int nparts, int nMax);

__global__ void choose2(int *neighborCount, int *nChoose2, int nparts)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;

  if (i < nparts) {
    int c = neighborCount[i];
    nChoose2[i] = c*(c-1)/2;
  }
}

__global__ void combine_nodes(int *neighborList, int *neighborCount,
  int *nodes, int *strides, int nparts, int nMax)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int N1, N2;
  int countPair = 0;
  int c;


  if (i < nparts) {
    N1 = i;
    for (int j = 0; j < neighborCount[i]; j++) {
      N2 = neighborList[i*nMax + j];

          c = NPOINTS*strides[i] + NPOINTS*countPair;
          countPair++;

          // sort while adding to nodes
          nodes[c] = N1*(N1 < N2) + N2*(N2 < N1) + INT_MIN*(N1 == N2);
          nodes[c + 1] = N2*(N1 < N2) + N1*(N2 < N1) + INT_MIN*(N1 == N2);
    }
  }
}

__global__ void sort_combos(int *nodes, int nPairs)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int startInd;
  int vals[NPOINTS];

  if (i < nPairs) {
    startInd = NPOINTS*i; 
    for (int v = 0; v < NPOINTS; v++) {
      vals[v] = nodes[startInd + v];
    }
    int a, b;
    // Adapated from stackoverflow.com 2786899, 
    //  "Fastest sort of fixed length 6 int array"
    // and "pages.ripco.net/~jgamble/nw.html"
    // and perhaps "Bit Twiddling Hacks" by Sean Eron Anderson @ Stanford
    #define min(x, y) (x*(x < y) + y*(y < x) + x*(x == y))
    #define max(x, y) (x*(x > y) + y*(y > x) + y*(x == y))
    #define SWAP(x,y) { a = min(vals[x], vals[y]); \
                        b = max(vals[x], vals[y]); \
                        vals[x] = a; vals[y] = b; }
    SWAP(0,1);

    #undef SWAP
    #undef min
    #undef max

    for (int v = 0; v < NPOINTS; v++) {
      nodes[startInd + v] = vals[v];
    }
  }
}

__global__ void find_unique(int *nodes, int base, int remainder, int *isUnique)
{
  int threadInd = threadIdx.x + blockIdx.x*blockDim.x;

  // only want to look at OTHER sets of nodes, so add (n+1)
  int target = threadInd + (base + 1);// Target set
  int targetNode = NPOINTS*target;          // First node of target set

  int baseNode = NPOINTS*base;              // First node of base set
  int nSame = 0;

  if (isUnique[base] == 0) {    // if base is not unique, skip
      return;
  }
  __syncthreads();

  if (threadInd < remainder) {
    nSame = (nodes[baseNode] == nodes[targetNode]) 
          + (nodes[baseNode + 1] == nodes[targetNode + 1]);

    // set isUnique to 1 if nSame != NPOINTS, else set to zero
    // if isUnique == 1
    //  if nSame != NPOINTS --> isUnique = 1*1 = 1
    //  if nSame == NPOINTS --> isUnique = 0*1 = 0
    // if isUnique == 0
    //  if nSame != NPOINTS --> isUnique = 1*0 = 0
    //  if nSame == NPOINTS --> isUnique = 0*0 = 0
    isUnique[target] = (nSame != NPOINTS)*isUnique[target];
  }
}

//__global__ void find_unique2(int *nodes, int *isUnique, int nPerms) {
//  int base = threadIdx.x + blockIdx.x*blockDim.x;
//
//  if (base < nPerms) {
//    int baseNode = NPOINTS*base;;
//    int target, targetNode;
//    int nSame = 0;
//    for (target = 0; target < nPerms; target++) {
//      targetNode = NPOINTS*target;
//
//      nSame = (nodes[baseNode] == nodes[targetNode]) 
//            + (nodes[baseNode + 1] == nodes[targetNode + 1]);
//      // if target == base, set nSame to 1
//      nSame = (target == base) + nSame*(target != base);
//
//      // set isUnique to 1 if nSame != NPOINTS, else set to zero
//      // if isUnique == 1
//      //  if nSame != NPOINTS --> isUnique = 1*1 = 1
//      //  if nSame == NPOINTS --> isUnique = 0*1 = 0
//      // if isUnique == 0
//      //  if nSame != NPOINTS --> isUnique = 1*0 = 0
//      //  if nSame == NPOINTS --> isUnique = 0*0 = 0
//      isUnique[base] *= (nSame != NPOINTS);
//       
//    }
//  }
//}

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

    uniqueNodes[NPOINTS*ind] = nodes[NPOINTS*TID];
    uniqueNodes[NPOINTS*ind + 1] = nodes[NPOINTS*TID + 1];

  }
}

__global__ void fill_nodes(pair_struct *pairs, int *uniqueNodes, int nPairs)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;

  if (i < nPairs) {
    pairs[i].N1 = uniqueNodes[NPOINTS*i];    
    pairs[i].N2 = uniqueNodes[NPOINTS*i + 1];    

    pairs[i].N1_initFlip_X = 0;
    pairs[i].N1_initFlip_Y = 0;
    pairs[i].N1_initFlip_Z = 0;

    pairs[i].N2_initFlip_X = 0;
    pairs[i].N2_initFlip_Y = 0;
    pairs[i].N2_initFlip_Z = 0;
  }
}

__global__ void flip_check(part_struct *parts, pair_struct *pairs,
  dom_struct *dom, int nPairs)
{
  int p = threadIdx.x + blockIdx.x*blockDim.x;

  double r1[nDim]; // TODO: what? maybe only need one
  double r2[nDim];

  if (p < nPairs) {
    /*  POSITION  */
    // Fix periodicity issues for first timestep 
    // Save an initial flip if one is needed
    periodic_flip(r1, r2, &pairs[p], parts, dom->xl, dom->yl, dom->zl);
  }
}

__global__ void flip_kernel(part_struct *parts, part_struct *partsPrev,
  dom_struct *dom, int nparts)
{
  // I almost wonder if the old version of this would be 
  //  fine as long as the initialization flip is taken
  //  care of on a pair-by-pair case -- do timestep 
  //  flip here, and then correct for pairs if need be?

  int n = threadIdx.x + blockIdx.x*blockDim.x;
  double dx; double dy; double dz;
  int cmpX; int cmpY; int cmpZ;
  int sgnDx; int sgnDy; int sgnDz;

  if (n < nparts) {
    // Before comparison, correct particle position with old flip count
    parts[n].x += dom->xl*parts[n].flipCountX;
    parts[n].y += dom->yl*parts[n].flipCountY;
    parts[n].z += dom->zl*parts[n].flipCountZ;

    // Set up dx,dy,dz
    dx = partsPrev[n].x - parts[n].x;
    dy = partsPrev[n].y - parts[n].y;
    dz = partsPrev[n].z - parts[n].z;

    // See if particle has crossed domain since last timestep
    //  1 iff it did
    //  0 iff it did not
    cmpX = (fabs(dx) >= 0.5*dom->xl);
    cmpY = (fabs(dy) >= 0.5*dom->yl);
    cmpZ = (fabs(dz) >= 0.5*dom->zl);

    // Get sign of dx,dy,dz
    //  -- if (prev - curr) > 0, went R->L, (+) a domain length
    //  -- if (prev - curr) < 0, went L->R, (-) a domain length
    sgnDx = (dx > 0.) - (dx < 0.);
    sgnDy = (dy > 0.) - (dy < 0.);
    sgnDz = (dz > 0.) - (dz < 0.);

    // increment or decrement the appropriate flipcount if so
    parts[n].flipCountX += cmpX*sgnDx;
    parts[n].flipCountY += cmpY*sgnDy;
    parts[n].flipCountZ += cmpZ*sgnDz;

    // Perform one more flip if necssary
    parts[n].x += dom->xl*cmpX*sgnDx;
    parts[n].y += dom->yl*cmpY*sgnDy;
    parts[n].z += dom->zl*cmpZ*sgnDz;
  }
}

__global__ void pair_statistics(part_struct *parts, pair_struct *pairs,
  dom_struct *dom, double *rSep, int nPairs, int tt)
{
  int p = threadIdx.x + blockIdx.x*blockDim.x;

  /* Tetrahedron Geometry Variables */
  double r1[3], r2[3];        // position of particles post flip
  double rx, ry, rz;    // separation of particles

  /* Velocity variables */
 // double u1[nDim];      // Node 1 relative vel
 // double u2[nDim];      // Node 2 relative vel

  /* Misc */
  int N1, N2;           // pair nodes
//int nrot;             // number of jacobi rotations

  if (p < nPairs) {
    /*  POSITION  */
    N1 = pairs[p].N1;
    N2 = pairs[p].N2;
    
    // Pull particle positions
    //  -- These have been corrected for periodicity between tsteps, but not
    //      yet for the initialization
    r1[0] = parts[N1].x; r1[1] = parts[N1].y; r1[2] = parts[N1].z;
    r2[0] = parts[N2].x; r2[1] = parts[N2].y; r2[2] = parts[N2].z;

    // Fix periodicity issues for pairs that were initialized over a boundary
    r1[0] += dom->xl * pairs[p].N1_initFlip_X;
    r1[1] += dom->yl * pairs[p].N1_initFlip_Y;
    r1[2] += dom->zl * pairs[p].N1_initFlip_Z;

    r2[0] += dom->xl * pairs[p].N2_initFlip_X;
    r2[1] += dom->yl * pairs[p].N2_initFlip_Y;
    r2[2] += dom->zl * pairs[p].N2_initFlip_Z;

    // separation distance
    rx = r1[0] - r2[0];
    ry = r1[1] - r2[1];
    rz = r1[2] - r2[2];
    rSep[p] = sqrt(rx*rx + ry*ry + rz*rz);

    // Relate nodes to center of vel
  //  u1[0] = parts[N1].u;
  //  u1[1] = parts[N1].v;
  //  u1[2] = parts[N1].w;

  //  u2[0] = parts[N2].u;
  //  u2[1] = parts[N2].v;
  //  u2[2] = parts[N2].w;
  }
}

__device__ void periodic_flip(double *r1, double *r2, pair_struct *pairs, 
  part_struct *parts, double xl, double yl, double zl)
{
  // For first timestep, correct pairs that have crossed periodic boundaries
  //  but don't know it.
  // Also, increment initial flip counter if a specific node needs to be
  //  permanently flipped down the line

  int N1 = pairs->N1;
  int N2 = pairs->N2;

  double dx, dx2;
  double flipL, flipR;
  int flipFlag;

  // branchless min
  #define flip(s1,s2,l,i,count) \
    {dx = s1[i] - s2[i]; \
     dx2 = dx*dx; \
     flipL = s1[i] - (s2[i] + l); \
     flipR = s1[i] - (s2[i] - l); \
     flipFlag = (flipL*flipL < dx2) - (flipR*flipR < dx2); \
     s2[i] += l*(flipFlag); \
     count += flipFlag; }

  // Set up position vectors to each particle from origin
  r1[0] = parts[N1].x;
  r1[1] = parts[N1].y;
  r1[2] = parts[N1].z;

  r2[0] = parts[N2].x;
  r2[1] = parts[N2].y;
  r2[2] = parts[N2].z;

  /* X,Y,Z directions */
  flip(r1, r2, xl, 0, pairs->N2_initFlip_X);
  flip(r1, r2, yl, 1, pairs->N2_initFlip_Y);
  flip(r1, r2, zl, 2, pairs->N2_initFlip_Z);

  #undef flip
}

__global__ void higher_moments_kernel(double *array, double mean, int length,
  double *sum2, double *skew, double *kurt) //double *diff, double *diff2)
{
  int n = threadIdx.x + blockIdx.x*blockDim.x;
  double diff, diff2;

  if (n < length) {
    sum2[n] = array[n]*array[n];  // xi*xi
    diff = array[n] - mean;    // xi - x_bar
    diff2 = diff*diff;         // (xi - x_bar)^2
    skew[n] = diff2*diff;      // (xi - x_bar)^3
    kurt[n] = diff2*diff2;     // (xi - xbar)^4
  } 
}

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
    double *w_s1, double *w_s2, double *w_s3)
{
  // cos(theta) = (a,b)
  #define dot3(a,b) a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

  int p = threadIdx.x + blockIdx.x*blockDim.x;

  double g1[3];
  double g2[3];
  double g3[3];
  double s1[3];
  double s2[3];
  double s3[3];
  double g1_0[3];
  double g2_0[3];
  double g3_0[3];
  double s1_0[3];
  double s2_0[3];
  double s3_0[3];
  double vort[3];
  double z[3];

  if (p < nRegular) {
    for (int i = 0; i < nDim; i++) {
      // Current vectors
      g1[i] = gEigVec[nDim2*p + i*nDim];
      g2[i] = gEigVec[nDim2*p + i*nDim + 1];
      g3[i] = gEigVec[nDim2*p + i*nDim + 2];

      s1[i] = sEigVec[nDim2*p + i*nDim];
      s2[i] = sEigVec[nDim2*p + i*nDim + 1];
      s3[i] = sEigVec[nDim2*p + i*nDim + 2];

      vort[i] = vorticity[nDim*p + i];

      // Initial vectors
      g1_0[i] = gEigVecInit[nDim2*p + i*nDim];
      g2_0[i] = gEigVecInit[nDim2*p + i*nDim + 1];
      g3_0[i] = gEigVecInit[nDim2*p + i*nDim + 2];
                                                     
      s1_0[i] = sEigVecInit[nDim2*p + i*nDim];
      s2_0[i] = sEigVecInit[nDim2*p + i*nDim + 1];
      s3_0[i] = sEigVecInit[nDim2*p + i*nDim + 2];
    }
    z[0] = 0.; z[1] = 0.; z[2] = 1.;

  
    // Alignment of shape axes with initial strain axes
    g1_s1[p] = fabs(dot3(g1,s1_0));
    g1_s2[p] = fabs(dot3(g1,s2_0));
    g1_s3[p] = fabs(dot3(g1,s3_0));

    g2_s1[p] = fabs(dot3(g2, s1_0));
    g2_s2[p] = fabs(dot3(g2, s2_0));
    g2_s3[p] = fabs(dot3(g2, s3_0));

    g3_s1[p] = fabs(dot3(g3, s1_0));
    g3_s2[p] = fabs(dot3(g3, s2_0));
    g3_s3[p] = fabs(dot3(g3, s3_0));

    // Alignment of shape,strain,vorticity with z axis
    g1_z[p] = fabs(dot3(g1, z));
    g2_z[p] = fabs(dot3(g2, z));
    g3_z[p] = fabs(dot3(g3, z));

    s1_z[p] = fabs(dot3(s1, z));
    s2_z[p] = fabs(dot3(s2, z));
    s3_z[p] = fabs(dot3(s3, z));

    // actual direction of vorticity DOES matter -- so don't use fabs
    w_z[p] = dot3(vort, z);

    // Alignment of vorticity with initial shape, strain
    w_g1[p] = dot3(vort, g1_0);
    w_g2[p] = dot3(vort, g2_0);
    w_g3[p] = dot3(vort, g3_0);

    w_s1[p] = dot3(vort, s1_0);
    w_s2[p] = dot3(vort, s2_0);
    w_s3[p] = dot3(vort, s3_0);
  }
  #undef dot3
}
