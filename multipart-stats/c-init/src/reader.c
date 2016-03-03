#include "tetrad_init.h"
#include "reader.h"

int nFiles;
char **partFiles;
double *partFileTime;
int *fileMap;
int nparts;
double binLength;
int nMax;
part_struct *parts;
part_struct *_parts;
dom_struct dom;
dom_struct *_dom;
dom_struct binDom;
dom_struct *_binDom;
BC bc;
BC *_bc;
tetrad_struct *tetrads;
tetrad_struct *_tetrads;

// Read tetrad.config input file
void tetrad_read_input(void)
{
  int fret = 0;
  fret = fret;

  // open config file for reading
  char fname[CHAR_BUF_SIZE] = "";
  sprintf(fname, "%s/%s/tetrad.config", ROOT_DIR, DATA_OUT_DIR);
  FILE *infile = fopen(fname, "r");
  
  // read input
  fret = fscanf(infile, "tStart %lf\n", &tStart);
  fret = fscanf(infile, "tEnd %lf\n", &tEnd);
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "R0/a %lf\n", &R0_a);
  fret = fscanf(infile, "eps/a %lf\n", &eps_a);
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "Variance Cutoff %lf\n", &varCut);
  fret = fscanf(infile, "Shape Cutoff %lf\n", &shapeCut);
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "Find Tetrads %d\n", &findTets);
  fclose(infile);
}

// read and sort output directory
void init_input_files(void) {
  DIR *dir;
  struct dirent *ent;
  char output_path[FILE_NAME_SIZE] = "";
  int fret = 0; fret=fret;
  double time;

  sprintf(output_path, "%s/%s", ROOT_DIR, OUTPUT_DIR);

  int isPart;
  int inRange;
  nFiles = 0;

  // count number of part files in directory that fall in time range
  if ((dir = opendir (output_path)) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      // check if part file (0 if match)
      isPart = (strncmp(ent->d_name, "part", 4) == 0);

      if (isPart == 1) {
        // check if in time range
        fret = sscanf(ent->d_name, "part-%lf.cgns", &time);
        inRange = ((time >= tStart) & (time <= tEnd));
        nFiles += isPart*inRange;
      } else {
        continue;
      }
    }
    closedir (dir);
  } else {
    printf("Output directory does not exist!\n");
    exit(EXIT_FAILURE);
  }

  // store cgns filenames and times within range
  // char **partFiles;
  partFiles = (char**) malloc(nFiles * sizeof(char*));
  for (int i = 0; i < nFiles; i ++) {
    partFiles[i] = (char*) malloc(FILE_NAME_SIZE*sizeof(char));
  }
  partFileTime = (double*) malloc(nFiles * sizeof(double));

  int cc = 0;

  if ((dir = opendir (output_path)) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      isPart = (strncmp(ent->d_name, "part", 4) == 0);

      if (isPart == 1) {
        // check if in time range
        fret = sscanf(ent->d_name, "part-%lf.cgns", &time);
        inRange = ((time >= tStart) & (time <= tEnd));
        
        if (inRange == 1) {
          fret = sscanf(ent->d_name, "%s", partFiles[cc]);
          partFileTime[cc] = time;
          cc++;
        }
      } else {
        continue;
      }
    }
    closedir (dir);
  } else {
    printf("Output directory does not exist!\n");
    exit(EXIT_FAILURE);
  }

  // Sort the resulting array by time
  // create temporary array to sort partFiles by
  fileMap = malloc(nFiles * sizeof(double));
  for (int i = 0; i < nFiles; i++) {
    fileMap[i] = i;
  }

  merge_sort(partFileTime, nFiles, fileMap);
  printf("Found %d files in range [%lf, %lf]\n", nFiles, tStart, tEnd);
}

// entry point for mergesort
void merge_sort(double *A, int n, int *A2) 
{
  if (n < 2)                          // if there is only one element: done
    return;
  int m = 0.5 * n;                    // cut array in half
  merge_sort(A, m, A2);               // recurse on first half
  merge_sort(A + m, n - m, A2 + m);   // recurse on second half
  merge(A, n, m, A2);                 // merge the two halves
}

// merge two arrays
void merge(double *A, int n, int m, int *A2) 
{
  // indices
  int i, j, k;
  // temporary array
  double *B  = malloc(n * sizeof(double));
  int *B2 = malloc(n * sizeof(int));
  i = 0;                      // first-half index
  j = m;                      // second-half index
  // proceed through entire array
  for(k = 0; k < n; k++) {
    if(j == n) {              // if j has reached the end
      B[k] = A[i];            // take the remaining elements of first-half of A
      B2[k] = A2[i];          // take A2 along for the ride
      i++;                    // increment i
    } else if(i == m) {       // else if i has reached half-way
      B[k] = A[j];            // take the remaining elements of second-half of A
      B2[k] = A2[j];          // take A2 along for the ride
      j++;                    // increment j
    } else if(A[j] < A[i]) {  // else compare two halves of A
      B[k] = A[j];            // take the second half if smaller
      B2[k] = A2[j];          // take A2 along for the ride
      j++;                    // increment j
    } else {                  // else
      B[k] = A[i];            // take first half if smaller
      B2[k] = A2[i];          // take A2 along for the ride
      i++;                    // increment i
    }
  }
  // overwrite unsorted A with sorted B
  for(i = 0; i < n; i++) {
    A[i] = B[i];
    A2[i] = B2[i];   // take A2 along for the ride
  }
  free(B);                    // clean up
  free(B2);
}

// read CUDA_VISIBLE_DEVICES from slurm and set dev_start
int read_devices(void)
{
  char *cvdin;
  cvdin = getenv("CUDA_VISIBLE_DEVICES");
  if (cvdin != NULL) {
    // number of devices
    int n_CUDA_VISIBLE_DEVICES = 0.5*(strlen(cvdin) + 1.);

    // list of devices
    int *CUDA_VISIBLE_DEVICES = malloc(n_CUDA_VISIBLE_DEVICES * sizeof(int));

    // fill list of devices assuming single-character separation
    int j = 0;
    for (int i = 0; i < 2*n_CUDA_VISIBLE_DEVICES - 1; i += 2) {
      CUDA_VISIBLE_DEVICES[j] = cvdin[i] - '0';
      j++;
    }

    // use the second device available (devices reindexed by CUDA)
    if (n_CUDA_VISIBLE_DEVICES > 1 ) {
      dev_start = 1;
    // if only one is available try to use it
    } else if (n_CUDA_VISIBLE_DEVICES == 1) {
      dev_start = 0;
    } else { // exit
      printf("Environment variable CUDA_VISIBLE_DEVICES is empty:\n");
      exit(EXIT_FAILURE);
    }

    free(CUDA_VISIBLE_DEVICES);
  } else {
    printf("CUDA_VISIBLE_DEVICES is NULL\n");
    printf("  If running on guinevere, use ");
    printf("export CUDA_VISIBLE_DEVICES=1\n");
    dev_start = -1;
    exit(EXIT_FAILURE);
  }

  return dev_start;
}

// initialize part_struct
void parts_init(int nparts)
{
  parts = (part_struct*) malloc(nparts * sizeof(part_struct));

  for(int p = 0; p < nparts; p++) {
    parts[p].x = -1;
    parts[p].y = -1;
    parts[p].z = -1;
    parts[p].r = -1;
    parts[p].bin = -1;
    parts[p].u = -1;
    parts[p].v = -1;
    parts[p].w = -1;

    // Allocate space for tetrad statistic arrays
    parts[p].g_lambda1 = malloc(nFiles * sizeof(double));
    parts[p].g_lambda2 = malloc(nFiles * sizeof(double));
    parts[p].g_lambda3 = malloc(nFiles * sizeof(double));
    for (int n = 0; n < nFiles; n++) {
      parts[p].g_lambda1[n] = -1; 
      parts[p].g_lambda2[n] = -1; 
      parts[p].g_lambda3[n] = -1; 
    }
  }


  // Open cgns file and get cgns file index number fn
  char buf[FILE_NAME_SIZE];
  sprintf(buf, "%s/%s/%s", ROOT_DIR, OUTPUT_DIR, partFiles[0]);
  int fn;
  cg_open(buf, CG_MODE_READ, &fn);
  
  // Set base, zone, and solutions index numbers
  int bn = 1;
  int zn = 1;
  int sn = 1;

  // Read part coords
  double *x = malloc(nparts * sizeof(double));
  double *y = malloc(nparts * sizeof(double));
  double *z = malloc(nparts * sizeof(double));
  for (int p = 0; p < nparts; p++) {
    x[p] = 0;
    y[p] = 0;
    z[p] = 0;
  }

  cgsize_t range_min = 1;
  cgsize_t range_max = nparts;

  cg_coord_read(fn,bn,zn, "CoordinateX", RealDouble, &range_min, &range_max, x);
  cg_coord_read(fn,bn,zn, "CoordinateY", RealDouble, &range_min, &range_max, y);
  cg_coord_read(fn,bn,zn, "CoordinateZ", RealDouble, &range_min, &range_max, z);

  for (int p = 0; p < nparts; p++) {
    parts[p].x = x[p];
    parts[p].y = y[p];
    parts[p].z = z[p];
  }

  // Read part radius (TODO: only once?)
  double *r = malloc(nparts * sizeof(double));
  for (int p = 0; p < nparts; p++) {
    r[p] = 0;
  }

  cg_field_read(fn,bn,zn,sn, "Radius", RealDouble, &range_min, &range_max, r);

  for (int p = 0; p < nparts; p++) {
    parts[p].r = r[p];
  }

  #ifdef DEBUG
  for (int p = 0; p < nparts; p++) {
    printf("  p = %d, (x,y,z); (%lf, %lf, %lf); %lf\n", p, x[p], y[p], z[p], r[p]);
  }
  #endif

  cg_close(fn);
  free(x);
  free(y);
  free(z);
  free(r);
}

void domain_init(void)
{
  int fret = 0;
  fret = fret;  // prevent compiler warning

  // open config file for reading
  char fname[FILE_NAME_SIZE] = "";
  sprintf(fname, "%s/input/flow.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if (infile == NULL) {
    printf("Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }
  
  // read buffers
  int ibuf = 0;
  char cbuf[CHAR_BUF_SIZE] = "";
  double dbuf = 0;

  // read domain
  fret = fscanf(infile, "DOMAIN\n");
  fret = fscanf(infile, "(Xs, Xe, Xn) %lf %lf %d\n", &dom.xs, &dom.xe, &dom.xn);
  fret = fscanf(infile, "(Ys, Ye, Yn) %lf %lf %d\n", &dom.ys, &dom.ye, &dom.yn);
  fret = fscanf(infile, "(Zs, Ze, Zn) %lf %lf %d\n", &dom.zs, &dom.ze, &dom.zn);
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "GPU DOMAIN DECOMPOSITION\n");
  fret = fscanf(infile, "DEV RANGE %d %d\n", &ibuf, &ibuf);
  fret = fscanf(infile, "n %d\n", &ibuf);
  fret = fscanf(infile, "(Xs, Xe, Xn) %lf %lf %d\n", &dbuf, &dbuf, &ibuf);
  fret = fscanf(infile, "(Ys, Ye, Yn) %lf %lf %d\n", &dbuf, &dbuf, &ibuf);
  fret = fscanf(infile, "(Zs, Ze, Zn) %lf %lf %d\n", &dbuf, &dbuf, &ibuf);
  fret = fscanf(infile, "E %d W %d N %d S %d T %d B %d\n", &ibuf, &ibuf, &ibuf,
    &ibuf, &ibuf, &ibuf);
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "PHYSICAL PARAMETERS\n");
  fret = fscanf(infile, "rho_f %lf\n", &dbuf);
  fret = fscanf(infile, "nu %lf\n", &dbuf);
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "SIMULATION PARAMETERS\n");
  fret = fscanf(infile, "duration %lf\n", &dbuf);
  fret = fscanf(infile, "CFL %lf\n", &dbuf);
  fret = fscanf(infile, "pp_max_iter %d\n", &ibuf);
  fret = fscanf(infile, "pp_residual %lf\n", &dbuf);
  fret = fscanf(infile, "lamb_max_iter %d\n", &ibuf);
  fret = fscanf(infile, "lamb_residual %lf\n", &dbuf);
  fret = fscanf(infile, "lamb_relax %lf\n", &dbuf);
  fret = fscanf(infile, "lamb_cut %lf\n", &dbuf);
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "BOUNDARY CONDITIONS\n");
  fret = fscanf(infile, "vel_tDelay %lf\n", &dbuf);

  fret = fscanf(infile, "PRESSURE\n");
  fret = fscanf(infile, "bc.pW %s", cbuf);
  if(strcmp(cbuf, "PERIODIC") == 0) {
    bc.pW = PERIODIC;
  } else if(strcmp(cbuf, "NEUMANN") == 0) {
    bc.pW = NEUMANN;
    fret = fscanf(infile, "%lf", &dbuf);
  } else {
    printf("flow.config read error -- read %s\n",cbuf);
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "bc.pE %s", cbuf);
  if(strcmp(cbuf, "PERIODIC") == 0) {
    bc.pE = PERIODIC;
  } else if(strcmp(cbuf, "NEUMANN") == 0) {
    bc.pE = NEUMANN;
    fret = fscanf(infile, "%lf", &dbuf);
  } else {
    printf("flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "bc.pS %s", cbuf);
  if(strcmp(cbuf, "PERIODIC") == 0) {
    bc.pS = PERIODIC;
  } else if(strcmp(cbuf, "NEUMANN") == 0) {
    bc.pS = NEUMANN;
    fret = fscanf(infile, "%lf", &dbuf);
  } else {
    printf("flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "bc.pN %s", cbuf);
  if(strcmp(cbuf, "PERIODIC") == 0) {
    bc.pN = PERIODIC;
  } else if(strcmp(cbuf, "NEUMANN") == 0) {
    bc.pN = NEUMANN;
    fret = fscanf(infile, "%lf", &dbuf);
  } else {
    printf("flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "bc.pB %s", cbuf);
  if(strcmp(cbuf, "PERIODIC") == 0) {
    bc.pB = PERIODIC;
  } else if(strcmp(cbuf, "NEUMANN") == 0) {
    bc.pB = NEUMANN;
    fret = fscanf(infile, "%lf", &dbuf);
  } else {
    printf("flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "bc.pT %s", cbuf);
  if(strcmp(cbuf, "PERIODIC") == 0) {
    bc.pT = PERIODIC;
  } else if(strcmp(cbuf, "NEUMANN") == 0) {
    bc.pT = NEUMANN;
    fret = fscanf(infile, "%lf", &dbuf);
  } else {
    printf("flow.config read error.\n");
    exit(EXIT_FAILURE);
  }

  /**** dom ****/
  // Calculate domain sizes
  dom.xl = dom.xe - dom.xs;
  dom.yl = dom.ye - dom.ys;
  dom.zl = dom.ze - dom.zs;

  // Calculate cell size
  dom.dx = dom.xl / dom.xn;
  dom.dy = dom.yl / dom.yn;
  dom.dz = dom.zl / dom.zn;

  // Set up grids
  dom.Gcc.is = 1;
  dom.Gcc.js = 1;
  dom.Gcc.ks = 1;

  dom.Gcc.in = dom.xn;
  dom.Gcc.jn = dom.yn;
  dom.Gcc.kn = dom.zn;

  dom.Gcc.ie = dom.Gcc.is + dom.xn;
  dom.Gcc.je = dom.Gcc.js + dom.yn;
  dom.Gcc.ke = dom.Gcc.ks + dom.zn;

  dom.Gcc.s1 = dom.Gcc.in;
  dom.Gcc.s2 = dom.Gcc.jn * dom.Gcc.s1;
  dom.Gcc.s3 = dom.Gcc.kn * dom.Gcc.s2;

  /**** binDom ****/
  // Calculate domain size
  binDom.xs = dom.xs;
  binDom.ys = dom.ys;
  binDom.zs = dom.zs;

  binDom.xe = dom.xe;
  binDom.ye = dom.ye;
  binDom.ze = dom.ze;

  binDom.xl = dom.xl;
  binDom.yl = dom.yl;
  binDom.zl = dom.zl;

  // Calculate cell sizes
  double rmax = 0;
  for (int i = 0; i < nparts; i++) {
    rmax = rmax + (parts[i].r > rmax)*(parts[i].r - rmax);
  }
  binLength = 2.*(R0_a*rmax + eps_a*rmax)/3.;

  binDom.xn = floor(binDom.xl/(binLength));
  binDom.yn = floor(binDom.yl/(binLength));
  binDom.zn = floor(binDom.zl/(binLength));

  // Avoid /0
  if (binDom.xn == 0) {
    binDom.xn = 1;
    binDom.dx = binDom.xl;
  } else {
    binDom.dx = binDom.xl / binDom.xn; 
  }
  if (binDom.yn == 0) {
    binDom.yn = 1;
    binDom.dy = binDom.yl;
  } else {
    binDom.dy = binDom.yl / binDom.yn; 
  }
  if (binDom.zn == 0) {
    binDom.zn = 1;
    binDom.dz = binDom.zl;
  } else {
    binDom.dz = binDom.zl / binDom.zn; 
  }

  // Set up grids
  binDom.Gcc.is = 1;
  binDom.Gcc.js = 1;
  binDom.Gcc.ks = 1;

  binDom.Gcc.in = binDom.xn;
  binDom.Gcc.jn = binDom.yn;
  binDom.Gcc.kn = binDom.zn;
  
  binDom.Gcc.ie = binDom.Gcc.is + binDom.Gcc.in;
  binDom.Gcc.je = binDom.Gcc.js + binDom.Gcc.jn;
  binDom.Gcc.ke = binDom.Gcc.ks + binDom.Gcc.kn;

  binDom.Gcc.s1 = binDom.Gcc.in;
  binDom.Gcc.s2 = binDom.Gcc.s1 * binDom.Gcc.jn;
  binDom.Gcc.s3 = binDom.Gcc.s2 * binDom.Gcc.kn;

  // alpha = n*Vp/Vd
  nMax = ceil((ALPHA_MAX*8.*(R0_a + eps_a)*(R0_a + eps_a)*(R0_a + eps_a))/
    (4./3.*PI));

  fclose(infile);

  #ifdef DEBUG
    show_domain();
  #endif
}

// read number of particles from cgns file
int cgns_read_nparts(void)
{
  // Open cgns file and get cgns file index number fn
  char buf[FILE_NAME_SIZE];
  sprintf(buf, "%s/%s/%s", ROOT_DIR, OUTPUT_DIR, partFiles[0]);
  int fn;
  cg_open(buf, CG_MODE_READ, &fn);

  // Set base index nuber and zone index number (only one, so is 1)
  int bn = 1;
  int zn = 1;

  // Read zone to find cgsize_t *size, or nparts
  char zonename[FILE_NAME_SIZE] = "";
  cgsize_t nparts = 0;
  cg_zone_read(fn, bn, zn, zonename, &nparts);

  cg_close(fn);

  return nparts;
}

// Read part_struct data
void cgns_fill_part_struct(int nparts, int tt)
{
  // Open cgns file and get cgns file index number fn
  char buf[FILE_NAME_SIZE];
  sprintf(buf, "%s/%s/%s", ROOT_DIR, OUTPUT_DIR, partFiles[tt]);
  int fn;
  cg_open(buf, CG_MODE_READ, &fn);
  
  // Set base, zone, and solutions index numbers
  int bn = 1;
  int zn = 1;
  int sn = 1;

  // Read part coords
  double *x = malloc(nparts * sizeof(double));
  double *y = malloc(nparts * sizeof(double));
  double *z = malloc(nparts * sizeof(double));
  for (int p = 0; p < nparts; p++) {
    x[p] = 0;
    y[p] = 0;
    z[p] = 0;
  }

  cgsize_t range_min = 1;
  cgsize_t range_max = nparts;

  cg_coord_read(fn,bn,zn, "CoordinateX", RealDouble, &range_min, &range_max, x);
  cg_coord_read(fn,bn,zn, "CoordinateY", RealDouble, &range_min, &range_max, y);
  cg_coord_read(fn,bn,zn, "CoordinateZ", RealDouble, &range_min, &range_max, z);

  for (int p = 0; p < nparts; p++) {
    parts[p].x = x[p];
    parts[p].y = y[p];
    parts[p].z = z[p];
  }

  // Read part vel
  double *u = malloc(nparts * sizeof(double));
  double *v = malloc(nparts * sizeof(double));
  double *w = malloc(nparts * sizeof(double));
  for (int p = 0; p < nparts; p++) {
    w[p] = 0;
    v[p] = 0;
    w[p] = 0;
  }
  cg_field_read(fn,bn,zn,sn, "VelocityX", RealDouble, &range_min, &range_max, u);
  cg_field_read(fn,bn,zn,sn, "VelocityY", RealDouble, &range_min, &range_max, v);
  cg_field_read(fn,bn,zn,sn, "VelocityZ", RealDouble, &range_min, &range_max, w);

  for (int p = 0; p < nparts; p++) {
    parts[p].u = u[p];
    parts[p].v = v[p];
    parts[p].w = w[p];
  }

  free(x);
  free(y);
  free(z);
  free(u);
  free(v);
  free(w);
  
  cg_close(fn);
}

// Show binDom
void show_domain(void)
{
  printf("Domain:\n");
  printf("  X: (%f, %f), dX = %f\n", dom.xs, dom.xe, dom.dx);
  printf("  Y: (%f, %f), dY = %f\n", dom.ys, dom.ye, dom.dy);
  printf("  Z: (%f, %f), dZ = %f\n", dom.zs, dom.ze, dom.dz);
  printf("  Xn = %d, Yn = %d, Zn = %d\n", dom.xn, dom.yn, dom.zn);
  printf("Domain Grids:\n");
  printf(" dom.Gcc:\n");
  printf("    is = %d, ie = %d, in = %d\n", dom.Gcc.is, dom.Gcc.ie, 
     dom.Gcc.in);
  printf("    js = %d, je = %d, jn = %d\n", dom.Gcc.js, dom.Gcc.je, 
     dom.Gcc.jn);
  printf("    ks = %d, ke = %d, kn = %d\n", dom.Gcc.ks, dom.Gcc.ke, 
     dom.Gcc.kn);
  printf("    s1 = %d, s2 = %d, s3 = %d\n", dom.Gcc.s1, dom.Gcc.s2,
     dom.Gcc.s3);
  printf("\n");

  printf("Bin Domain:\n");
  printf("  X: (%f, %f), dX = %f\n", binDom.xs, binDom.xe, binDom.dx);
  printf("  Y: (%f, %f), dY = %f\n", binDom.ys, binDom.ye, binDom.dy);
  printf("  Z: (%f, %f), dZ = %f\n", binDom.zs, binDom.ze, binDom.dz);
  printf("  Xn = %d, Yn = %d, Zn = %d\n", binDom.xn, binDom.yn, binDom.zn);
  printf("binDomain Grids:\n");
  printf("  binDom.Gcc:\n");
  printf("    is = %d, ie = %d, in = %d\n", binDom.Gcc.is, binDom.Gcc.ie, 
    binDom.Gcc.in);
  printf("    js = %d, je = %d, jn = %d\n", binDom.Gcc.js, binDom.Gcc.je, 
    binDom.Gcc.jn);
  printf("    ks = %d, ke = %d, kn = %d\n", binDom.Gcc.ks, binDom.Gcc.ke, 
    binDom.Gcc.kn);
  printf("    s1 = %d, s2 = %d, s3 = %d\n", binDom.Gcc.s1, binDom.Gcc.s2,
    binDom.Gcc.s3);
  printf("\n");

  printf("Boundary Condition Structure\n");
  printf("  PERIODIC = 0\n");
  printf("  DIRICHLET 1\n");
  printf("  NEUMANN = 2\n");
  printf("    bc.pW = %d\n", bc.pW);
  printf("    bc.pE = %d\n", bc.pE);
  printf("    bc.pS = %d\n", bc.pS);
  printf("    bc.pN = %d\n", bc.pN);
  printf("    bc.pB = %d\n", bc.pB);
  printf("    bc.pT = %d\n", bc.pT);
  printf("\n");


  printf("Input Parameters\n");
  printf("  tStart %lf\n", tStart);
  printf("  tEnd %lf\n", tEnd);
  printf("  R0/a %lf\n", R0_a);
  printf("  eps/a %lf\n", eps_a);
  printf("  Variance Cutoff %lf\n", varCut);
  printf("  Shape Cutoff %lf\n", shapeCut);
  printf("  Find Tetrads %d\n", findTets);
  printf("  Number of parts = %d\n", nparts);
  printf("  Max parts per shell = %d\n", nMax);
}

// Write nodes
void write_nodes(void)
{
  int count = 0;
  printf("\tWriting to file ""uniqueNodes""... ");

  // Set up output file name
  char fname[CHAR_BUF_SIZE] = "";
  sprintf(fname, "%s/%s/uniqueNodes", ROOT_DIR, DATA_OUT_DIR);
  FILE *fnode = fopen(fname, "w");
  if (fnode == NULL) {
    printf("Error opening file!\n");
    exit(EXIT_FAILURE);
  }

  fprintf(fnode, "n1,n2,n3,n4\n");
  for (int i = 0; i < nUnique; i++) {
  // TODO: add input parameters on these limits
    if (fabs(tetrads[i].eigVar) <= varCut) {
      if (fabs(tetrads[i].shape) <= shapeCut) {
        fprintf(fnode, "%d,%d,%d,%d\n",
          tetrads[i].N1, tetrads[i].N2, tetrads[i].N3, tetrads[i].N4);
        count++;
      }
    }
  }
  fclose(fnode);

  printf("Wrote %d tetrads\n", count);
}

void create_output_dir (void) {
  // Create output directory if it doesn't exist
  // From stackoverflow-7430248
  struct stat st = {0};
  char buf[CHAR_BUF_SIZE];
  sprintf(buf, "%s/%s", ROOT_DIR, DATA_OUT_DIR);
  if (stat(buf, &st) == -1) {
    mkdir(buf, 0700);
  }
}

// Write data at each timestep
void write_timestep(int tt)
{
  // Set up output file name
  char format[CHAR_BUF_SIZE] = "";
  char fnameall[CHAR_BUF_SIZE] = "";
  char fnameall2[CHAR_BUF_SIZE] = "";

  double tout = partFileTime[tt];
  int sigfigs = ceil(log10(1. / tout));
  if (sigfigs < 1) sigfigs = 1;

  sprintf(format, "%%.%df", sigfigs);
  sprintf(fnameall2, "%s/%s/tetrad-data-%s", ROOT_DIR, DATA_OUT_DIR, format);
  sprintf(fnameall, fnameall2, tout);

  FILE *fdat = fopen(fnameall, "w");
  if (fdat == NULL) {
    printf("Error opening file %s!\n", fnameall);
    exit(EXIT_FAILURE);
  }

  fprintf(fdat, "shape,eigVar,R2\n");
  for (int i = 0; i < nUnique; i++) {
    if (fabs(tetrads[i].eigVar) <= varCut) {
      if (fabs(tetrads[i].shape) <= shapeCut) {
        fprintf(fdat, "%lf,%lf,%lf\n",
          tetrads[i].shape, tetrads[i].eigVar, tetrads[i].R2);
      }
    }
  }
  fclose(fdat);
}

// Free parts
void free_parts(void)
{
  for (int p = 0; p < nparts; p++) {
    free(parts[p].g_lambda1);
    free(parts[p].g_lambda2);
    free(parts[p].g_lambda3);
  }
  free(parts);
  free(uniqueNodes);
  free(tetrads);
  for (int i = 0; i < nFiles; i++) {
    free(partFiles[i]);
  }
  free(partFiles);
  free(fileMap);
  free(partFileTime);
}

