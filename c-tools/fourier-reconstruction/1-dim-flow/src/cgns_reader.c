#include "main.h"
#include "cgns_reader.h"

int nFiles;
char **flowFiles;
double *flowFileTime;
int *flowFileMap;
int sigFigPre;
int sigFigPost;
dom_struct dom;

// fftw variables
double *uf;
double *vf;
double *wf;
int *phase;

fftw_complex *uf_k;
fftw_complex *vf_k;
fftw_complex *wf_k;

double *evalZ;
double *wf_rec_Re;
double *wf_rec_Im;

// Set up directory structure
void directory_init(int argc, char *argv[])
{
  SIM_ROOT_DIR = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
  ANALYSIS_DIR = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));

  // arg[0] = program name
  // arg[1] = SIM_ROOT_DIR
  if (argc == 2) {
  sprintf(SIM_ROOT_DIR, "%s", argv[1]);
  sprintf(ANALYSIS_DIR, "%s/analysis/%s/%s", SIM_ROOT_DIR, FREC_DIR, ANALYSIS);
  } else if (argc != 2) {
  printf("usage: %s SIM_ROOT_DIR\n", argv[0]);
  exit(EXIT_FAILURE);
  }
  printf("\n SIM_ROOT_DIR = %s\n", SIM_ROOT_DIR);
  printf(" ANALYSIS_DIR = %s\n\n", ANALYSIS_DIR);
  fflush(stdout);
}

// Read main.config input file
void main_read_input(void)
{
  int fret = 0;
  fret = fret;

  // open config file for reading
  char fname[CHAR_BUF_SIZE] = "";
  sprintf(fname, "%s/%s", ANALYSIS_DIR, CONFIG_FILE);
  FILE *infile = fopen(fname, "r");
  if (infile == NULL) {
    printf("Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  } else {
    printf("Reading config from %s...\n", fname);
  }
  
  // read input
  fret = fscanf(infile, "Starting Time %lf\n", &tStart);
  fret = fscanf(infile, "Ending Time %lf\n", &tEnd);
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "Fourier Order %d\n", &order);
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "Output Coefficients %d\n", &coeffsOut);
  fclose(infile);
}

// read and sort flow files directory
void init_flow_files(void) 
{
  DIR *dir;
  struct dirent *ent;
  char output_path[FILE_NAME_SIZE] = "";
  int fret = 0; fret=fret;
  double time;

  sprintf(output_path, "%s/%s", SIM_ROOT_DIR, OUTPUT_DIR);

  int isFlow;
  int inRange;
  nFiles = 0;

  // count number of flow files in directory that fall in time range
  if ((dir = opendir (output_path)) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      // check if flow file (0 if match)
      isFlow = (strncmp(ent->d_name, "flow", 4) == 0);

      if (isFlow == 1) {
        // check if in time range
        fret = sscanf(ent->d_name, "flow-%lf.cgns", &time);
        inRange = ((time >= tStart) & (time <= tEnd));
        nFiles += isFlow*inRange;
      } else {
        continue;
      }
    }
    closedir (dir);
  } else {
    printf("Output directory %s does not exist!\n", output_path);
    exit(EXIT_FAILURE);
  }

  // store cgns filenames and times within range
  flowFiles = (char**) malloc(nFiles * sizeof(char*));
  for (int i = 0; i < nFiles; i ++) {
    flowFiles[i] = (char*) malloc(FILE_NAME_SIZE*sizeof(char));
  }
  flowFileTime = (double*) malloc(nFiles * sizeof(double));

  int cc = 0;

  if ((dir = opendir (output_path)) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      isFlow = (strncmp(ent->d_name, "flow", 4) == 0);

      if (isFlow == 1) {
        // check if in time range
        fret = sscanf(ent->d_name, "flow-%lf.cgns", &time);
        inRange = ((time >= tStart) & (time <= tEnd));
        
        if (inRange == 1) {
          fret = sscanf(ent->d_name, "%s", flowFiles[cc]);
          flowFileTime[cc] = time;
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
  // create temporary array to sort flowFiles by
  flowFileMap = malloc(nFiles * sizeof(double));
  for (int i = 0; i < nFiles; i++) {
    flowFileMap[i] = i;
  }

  merge_sort(flowFileTime, nFiles, flowFileMap);
  printf("Found %d flow files in range [%lf, %lf]\n", nFiles, tStart, tEnd);
  if (nFiles == 0) {
    printf("Check something...\n");
    exit(EXIT_FAILURE);
  }
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

// Create direcotry for output data
void create_output(void) 
{
  // Create data directory if it doesn't exist
  // From stackoverflow-7430248
  struct stat st = {0};
  char buf[CHAR_BUF_SIZE];
  sprintf(buf, "%s/%s", ANALYSIS_DIR, DATA_DIR);
  if (stat(buf, &st) == -1) {
    mkdir(buf, 0700);
  }

  // Create output files
  char path2file[FILE_NAME_SIZE] = "";

  /* Create eval/time file */
  sprintf(path2file, "%s/%s/info", ANALYSIS_DIR, DATA_DIR);
  FILE *file = fopen(path2file, "w");
  if (file == NULL) {
    printf("Could not open file %s\n", path2file);
    exit(EXIT_FAILURE);
  }

  // Print time on first row
  fprintf(file, "time ");
  for (int t = 0; t < nFiles; t++) {
    fprintf(file, "%lf ", flowFileTime[t]);
  }
  fprintf(file, "\n");

  // Print evalZ on second row
  fprintf(file, "evalZ ");
  for (int i = 0; i < npoints; i++) {
    fprintf(file, "%lf ", evalZ[i]);
  }

  fclose(file);

//  /* volume fraction */
//  sprintf(path2file, "%s/%s/volume-fraction", ANALYSIS_DIR, DATA_DIR);
//  file = fopen(path2file, "w");
//  if (file == NULL) {
//    printf("Could not open file %s\n", path2file);
//    exit(EXIT_FAILURE);
//  }

//  /* only output coefficients if desired */
//  if (coeffsOut == 1) {
//
//    // flow-w coefficients even and odd
//    sprintf(path2file, "%s/%s/flow-w-coeffs-even", ANALYSIS_DIR, DATA_DIR);
//    file = fopen(path2file, "w");
//    if (file == NULL) {
//      printf("Could not open file %s\n", path2file);
//      exit(EXIT_FAILURE);
//    }
//    fclose(file);
//    sprintf(path2file, "%s/%s/flow-w-coeffs-odd", ANALYSIS_DIR, DATA_DIR);
//    file = fopen(path2file, "w");
//    if (file == NULL) {
//      printf("Could not open file %s\n", path2file);
//      exit(EXIT_FAILURE);
//    }
//    fclose(file);
//  }

}

// Read domain
void domain_init(void)
{
  int fret = 0;
  fret = fret;  // prevent compiler warning

  // open config file for reading
  char fname[FILE_NAME_SIZE] = "";
  sprintf(fname, "%s/input/flow.config", SIM_ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if (infile == NULL) {
    printf("Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }
  
  // read buffers
  int ibuf = 0;
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

  // Set up flow
  uf = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  vf = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  wf = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  phase = (int*) malloc(dom.Gcc.s3 * sizeof(int));

  for (int ii = 0; ii < dom.Gcc.s3; ii++) {
    uf[ii] = 0.;
    vf[ii] = 0.;
    wf[ii] = 0.;
    phase[ii] = 0;
  }

  //int outN = dom.Gcc.in*dom.Gcc.jn*(0.5*dom.Gcc.kn + 1);
  int outN = (0.5*dom.Gcc.in + 1)*dom.Gcc.jn*dom.Gcc.kn;
  uf_k = (fftw_complex*) fftw_malloc(outN * sizeof(fftw_complex));
  vf_k = (fftw_complex*) fftw_malloc(outN * sizeof(fftw_complex));
  wf_k = (fftw_complex*) fftw_malloc(outN * sizeof(fftw_complex));

  for (int ii = 0; ii < outN; ii++) {
    uf_k[ii][0] = 0.;
    uf_k[ii][1] = 0.;
    vf_k[ii][0] = 0.;
    vf_k[ii][1] = 0.;
    wf_k[ii][0] = 0.;
    wf_k[ii][1] = 0.;
  }

  // Calculate order using sampling theorem -- 2 points / wave
  if (order == -1) {
    order = (int) floor(0.5*dom.zn);
  } // else if order > ...
  printf("Order = %d\n", order);

  // Number of points to reconstruct at
  // Just because we can only resolve xn/2 doesnt mean we should underresolve 
  //  the reconstruction...
  // npoints = 4*dom.zn;
  npoints = dom.zn;
  double dz = dom.zl/npoints;

  evalZ = (double*) malloc(npoints * sizeof(double));
  for (int zz = 0; zz < npoints; zz++) {
    evalZ[zz] = dom.zs + dz*zz;
  }

  // Reconstruct variables
  wf_rec_Re = (double*) malloc(npoints*nFiles * sizeof(double));
  wf_rec_Im = (double*) malloc(npoints*nFiles * sizeof(double));
  for (int tt = 0; tt < nFiles; tt++) {
    for (int zz = 0; zz < npoints; zz++) {
      wf_rec_Re[zz + tt*npoints] = 0.;
      wf_rec_Im[zz + tt*npoints] = 0.;
    }
  }

  #ifdef DEBUG
    show_domain();
  #endif
  fclose(infile);
}

// Read part_struct data
void cgns_fill_flow(void)
{
  // Open cgns file and get cgns file index number fn
  char buf[FILE_NAME_SIZE];
  sprintf(buf, "%s/%s/%s", SIM_ROOT_DIR, OUTPUT_DIR, flowFiles[flowFileMap[tt]]);
  printf("file = %s\n", buf);
  fflush(stdout);
  int fn;
  int ier = cg_open(buf, CG_MODE_READ, &fn);
  if (ier != 0 ) {
    printf("CGNS Error - double check grid.cgns exists in output\n");
    cg_error_exit();
  }
  fflush(stdout);
  
  // Set base, zone, and solutions index numbers
  int bn = 1;
  int zn = 1;
  int sn = 1;

  // Size of array to pull
  cgsize_t range_min[3];
  range_min[0] = 1;
  range_min[1] = 1;
  range_min[2] = 1;
  cgsize_t range_max[3];
  range_max[0] = dom.xn;
  range_max[1] = dom.yn;
  range_max[2] = dom.zn;

  // Read and fill
  //cg_field_read(fn,bn,zn,sn, "VelocityX", RealDouble, range_min, range_max, uf);
  //cg_field_read(fn,bn,zn,sn, "VelocityY", RealDouble, range_min, range_max, vf);
  cg_field_read(fn,bn,zn,sn, "VelocityZ", RealDouble, range_min, range_max, wf);
  cg_field_read(fn,bn,zn,sn, "Phase", Integer, range_min, range_max, phase);

  // Apply phase mask
  for (int kk = 0; kk < dom.zn; kk++) {
    for (int jj = 0; jj < dom.yn; jj++) {
      for (int ii = 0; ii < dom.xn; ii++) {
        int stride = ii + dom.Gcc.s1*jj + dom.Gcc.s2*kk;
        //uf[stride] *= phase[stride];
        //vf[stride] *= phase[stride];
        // wf[stride] *= (phase[stride] == -1);
        wf[stride] = (double) (phase[stride] == -1);
      }
    }
  }
  cg_close(fn);
}

// Test fill
void test_fill(void)
{
  double x = 0.;
  double y = 0.;
  double z = 0.;
  double xstar = 0.;
  double ystar = 0.;
  double zstar = 0.;
  x = x; xstar = xstar;
  y = y; ystar = ystar;
  z = z; zstar = zstar;
  int cc = 0;
  for (int kk = 0; kk < dom.zn; kk++) {
    for (int jj = 0; jj < dom.yn; jj++) {
      for (int ii = 0; ii < dom.xn; ii++) {
        x = dom.xs + ii*dom.dx;
        y = dom.ys + jj*dom.dy;
        z = dom.zs + kk*dom.dz;
        
        cc = ii + dom.Gcc.s1*jj + dom.Gcc.s2*kk;

        xstar = 2.*PI*x/dom.xl;
        ystar = 2.*PI*y/dom.yl;
        zstar = 2.*PI*z/dom.zl;

        //wf[cc] = sin(3.*zstar) + sin(2.*zstar);
        //wf[cc] = sin(3.*x) + sin(2.*x);
        wf[cc] = sin(1.*xstar) + sin(2.*xstar) + sin(3.*xstar) +
          sin(1.*ystar) + sin(2.*ystar) + sin(3.*ystar) +
          sin(1.*zstar) + sin(2.*zstar) + sin(3.*zstar) + 2.;
      }
    }
    printf("z = %lf, sin(kzz) = %lf\n", z, wf[cc]);
  }
  fftw_execute(pW);
}

void test_out(void)
{
  int halfIN = 0.5*dom.Gcc.in + 1;
  printf("All Wavenumbers\n");
  for (int nn = 0; nn < dom.Gcc.kn; nn++) {   // z
    for (int mm = 0; mm < dom.Gcc.jn; mm++) { // y
      for (int ll = 0; ll < halfIN; ll++) {   // x
        int cc = ll + halfIN*mm + halfIN*dom.Gcc.jn*nn;

        if ((fabs(wf_k[cc][0]) > 1e-8) || (fabs(wf_k[cc][1]) > 1e-8)) {
          printf("wf_k[%d + %d*%d + %d*%d] = %lf + i%lf\n", 
            ll, halfIN, mm, halfIN*dom.Gcc.jn, nn, wf_k[cc][0], wf_k[cc][1]);
          fflush(stdout);
        }

      }
    }
  }

  printf("\n\n");

  int kx = 1;
  int ky = 1;
  int kz = 3;
  printf("Choose (kx,ky,kz) = (%d,%d,%d)\n", kx, ky, kz);
  for (int nn = 0; nn <= kz; nn++) {
    for (int mm = 0; mm <= ky; mm++) {
      for (int ll = 0; ll <= kx; ll++) {
        int cc = ll + halfIN*mm + halfIN*dom.Gcc.jn*nn;

        printf("wf_k[%d + %d*%d + %d*%d] = %lf + i%lf\n", 
          ll, halfIN, mm, halfIN*dom.Gcc.jn, nn, wf_k[cc][0], wf_k[cc][1]);
        fflush(stdout);
      }
    }
  }
}

// Show domain
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

  printf("Input Parameters\n");
  printf("  tStart %lf\n", tStart);
  printf("  tEnd %lf\n", tEnd);
  printf("  order %d\n", order);
  printf("  coeffsOut %d\n", coeffsOut);
  printf("  nPoints %d\n", npoints);
}

// Write reconstructed data
void write_reconstruct(void)
{
  char fname[CHAR_BUF_SIZE] = "";

  /* number desnity */
  sprintf(fname, "%s/%s/rec-wf-re", ANALYSIS_DIR, DATA_DIR);
  FILE *fileRe = fopen(fname, "w");
  if (fileRe == NULL) {
    printf("Error opening file %s!\n", fname);
    exit(EXIT_FAILURE);
  }

  sprintf(fname, "%s/%s/rec-wf-im", ANALYSIS_DIR, DATA_DIR);
  FILE *fileIm = fopen(fname, "w");
  if (fileIm == NULL) {
    printf("Error opening file %s!\n", fname);
    exit(EXIT_FAILURE);
  }

  // each reconstructed timestep is a row
  for (int t = 0; t < nFiles; t++) {
    for (int zz = 0; zz < npoints; zz++) {
      int cc = zz + t*npoints;
      fprintf(fileRe, "%lf ", wf_rec_Re[cc]);
      fprintf(fileIm, "%lf ", wf_rec_Im[cc]);
    }
    fprintf(fileRe, "\n");
    fprintf(fileIm, "\n");
  }

  fclose(fileRe);
  fclose(fileIm);
}

// Free parts
void free_vars(void)
{
  for (int i = 0; i < nFiles; i++) {
    free(flowFiles[i]);
  }
  free(flowFiles);
  free(flowFileMap);
  free(flowFileTime);

  free(uf);
  free(vf);
  free(wf);
  free(phase);

  fftw_free(uf_k);
  fftw_free(vf_k);
  fftw_free(wf_k);

  free(evalZ);
  free(wf_rec_Re);
  free(wf_rec_Im);

  free(SIM_ROOT_DIR);
  free(ANALYSIS_DIR);
}

