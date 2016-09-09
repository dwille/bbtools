#include "main.h"
#include "cgns_reader.h"

int nFiles;
char **partFiles;
double *partFileTime;
int *partFileMap;
int sigFigPre;
int sigFigPost;
int nparts;
double meanR;
double *up;
double *vp;
double *wp;
double *ke;
part_struct *parts;
dom_struct dom;

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

// read and sort part files directory
void init_part_files(void) 
{
  DIR *dir;
  struct dirent *ent;
  char output_path[FILE_NAME_SIZE] = "";
  int fret = 0; fret=fret;
  double time;

  sprintf(output_path, "%s/%s", SIM_ROOT_DIR, OUTPUT_DIR);

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
  partFileMap = malloc(nFiles * sizeof(double));
  for (int i = 0; i < nFiles; i++) {
    partFileMap[i] = i;
  }

  merge_sort(partFileTime, nFiles, partFileMap);
  printf("Found %d part files in range [%lf, %lf]\n", nFiles, tStart, tEnd);
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
    fprintf(file, "%lf ", partFileTime[t]);
  }
  fprintf(file, "\n");

  // Print evalZ on second row
  fprintf(file, "evalZ ");
  for (int i = 0; i < npoints; i++) {
    fprintf(file, "%lf ", evalZ[i]);
  }

  fclose(file);

  /* number density */
  sprintf(path2file, "%s/%s/number-density", ANALYSIS_DIR, DATA_DIR);
  file = fopen(path2file, "w");
  if (file == NULL) {
    printf("Could not open file %s\n", path2file);
    exit(EXIT_FAILURE);
  }
  fclose(file);

  /* volume fraction */
  sprintf(path2file, "%s/%s/volume-fraction", ANALYSIS_DIR, DATA_DIR);
  file = fopen(path2file, "w");
  if (file == NULL) {
    printf("Could not open file %s\n", path2file);
    exit(EXIT_FAILURE);
  }
  fclose(file);

  /* part-u */
  sprintf(path2file, "%s/%s/part-u", ANALYSIS_DIR, DATA_DIR);
  file = fopen(path2file, "w");
  if (file == NULL) {
    printf("Could not open file %s\n", path2file);
    exit(EXIT_FAILURE);
  }
  fclose(file);

  /* part-v */
  sprintf(path2file, "%s/%s/part-v", ANALYSIS_DIR, DATA_DIR);
  file = fopen(path2file, "w");
  if (file == NULL) {
    printf("Could not open file %s\n", path2file);
    exit(EXIT_FAILURE);
  }
  fclose(file);

  /* part-w */
  sprintf(path2file, "%s/%s/part-w", ANALYSIS_DIR, DATA_DIR);
  file = fopen(path2file, "w");
  if (file == NULL) {
    printf("Could not open file %s\n", path2file);
    exit(EXIT_FAILURE);
  }
  fclose(file);
  
  /* particle kinetic energy */
  sprintf(path2file, "%s/%s/part-kinetic-energy", ANALYSIS_DIR, DATA_DIR);
  file = fopen(path2file, "w");
  if (file == NULL) {
    printf("Could not open file %s\n", path2file);
    exit(EXIT_FAILURE);
  }
  fclose(file);

  /* only output coefficients if desired */
  if (coeffsOut == 1) {
    // part-w coefficients even and odd
    sprintf(path2file, "%s/%s/number-dens-coeffs-even", ANALYSIS_DIR, DATA_DIR);
    file = fopen(path2file, "w");
    if (file == NULL) {
      printf("Could not open file %s\n", path2file);
      exit(EXIT_FAILURE);
    }
    fclose(file);
    sprintf(path2file, "%s/%s/number-dens-coeffs-odd", ANALYSIS_DIR, DATA_DIR);
    file = fopen(path2file, "w");
    if (file == NULL) {
      printf("Could not open file %s\n", path2file);
      exit(EXIT_FAILURE);
    }
    fclose(file);

    // part-w coefficients even and odd
    sprintf(path2file, "%s/%s/part-w-coeffs-even", ANALYSIS_DIR, DATA_DIR);
    file = fopen(path2file, "w");
    if (file == NULL) {
      printf("Could not open file %s\n", path2file);
      exit(EXIT_FAILURE);
    }
    fclose(file);
    sprintf(path2file, "%s/%s/part-w-coeffs-odd", ANALYSIS_DIR, DATA_DIR);
    file = fopen(path2file, "w");
    if (file == NULL) {
      printf("Could not open file %s\n", path2file);
      exit(EXIT_FAILURE);
    }
    fclose(file);
  }

}

// read number of particles from cgns file
int cgns_read_nparts(void)
{
  // Open cgns file and get cgns file index number fn
  char buf[FILE_NAME_SIZE];
  sprintf(buf, "%s/%s/%s", SIM_ROOT_DIR, OUTPUT_DIR, partFiles[partFileMap[0]]);
  int fn;
  int ier = cg_open(buf, CG_MODE_READ, &fn);
  if (ier != 0 ) {
    printf("\nCGNS Error on file %s\n", buf);
    cg_error_exit();
  }
  fflush(stdout);

  // Set base index nuber and zone index number (only one, so is 1)
  int bn = 1;
  int zn = 1;

  // Read zone to find cgsize_t *size, or nparts
  char zonename[FILE_NAME_SIZE] = "";
  cgsize_t nparts = 0;
  ier = cg_zone_read(fn, bn, zn, zonename, &nparts);
  if (ier != 0 ) {
    printf("\nCGNS Error on file %s\n", buf);
    cg_error_exit();
  }
  fflush(stdout);

  cg_close(fn);

  return nparts;
}

// initialize part_struct
void parts_init(void)
{
  parts = (part_struct*) malloc(nparts * sizeof(part_struct));
  up = (double*) malloc(nparts * sizeof(double));
  vp = (double*) malloc(nparts * sizeof(double));
  wp = (double*) malloc(nparts * sizeof(double));
  ke = (double*) malloc(nparts * sizeof(double));

  for(int p = 0; p < nparts; p++) {
    parts[p].x = -1;
    parts[p].y = -1;
    parts[p].z = -1;
    parts[p].r = -1;
    up[p] = 0.;
    vp[p] = 0.;
    wp[p] = 0.;
    ke[p] = 0.;
  }

  // Open cgns file and get cgns file index number fn
  char buf[FILE_NAME_SIZE];
  sprintf(buf, "%s/%s/%s", SIM_ROOT_DIR, OUTPUT_DIR, partFiles[partFileMap[0]]);
  int fn;
  int ier = cg_open(buf, CG_MODE_READ, &fn);
  if (ier != 0 ) {
    printf("\nCGNS Error on file %s\n", buf);
    cg_error_exit();
  }
  fflush(stdout);
  
  // Set base, zone, and solutions index numbers
  int bn = 1;
  int zn = 1;
  int sn = 1;

  cgsize_t range_min = 1;
  cgsize_t range_max = nparts;

  // Read part radius
  double *r = malloc(nparts * sizeof(double));
  for (int p = 0; p < nparts; p++) {
    r[p] = 0;
  }

  ier = cg_field_read(fn,bn,zn,sn, "Radius", RealDouble, &range_min, &range_max, r);
  if (ier != 0 ) {
    printf("\nCGNS Error on file %s\n", buf);
    cg_error_exit();
  }
  fflush(stdout);

  meanR = 0.;
  for (int p = 0; p < nparts; p++) {
    parts[p].r = r[p];
    meanR += parts[p].r;
  }
  meanR /= nparts;

  cg_close(fn);
  free(r);
}

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

  // Calculate order using sampling theorem -- 1 wave / 1 diameters
  if (order == -1) {
    order = (int) floor(dom.zn / (2. * meanR) );
  }
  printf("Order = %d\n", order);
  // Number of points to reconstruct at
  // Just because we can only resolve xn/2 doesnt mean we should underresolve 
  //  the reconstruction...
  // TODO: how does this affect?
  npoints = dom.zn;
  printf("npoints = %d\n", npoints);

  #ifdef DEBUG
    show_domain();
  #endif
  fclose(infile);
}

// Read part_struct data
void cgns_fill_parts(void)
{
  // Open cgns file and get cgns file index number fn
  char buf[FILE_NAME_SIZE];
  sprintf(buf, "%s/%s/%s", SIM_ROOT_DIR, OUTPUT_DIR, partFiles[partFileMap[tt]]);
  int fn;
  int ier = cg_open(buf, CG_MODE_READ, &fn);
  if (ier != 0 ) {
    printf("\nCGNS Error on file %s\n", buf);
    free_vars();
    cg_close(fn);
    cg_error_exit();
  }
  fflush(stdout);
  
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
  cg_field_read(fn,bn,zn,sn, "VelocityX", RealDouble, &range_min, &range_max, up);
  cg_field_read(fn,bn,zn,sn, "VelocityY", RealDouble, &range_min, &range_max, vp);
  cg_field_read(fn,bn,zn,sn, "VelocityZ", RealDouble, &range_min, &range_max, wp);

  // Calculate kinetic energy
  for (int p = 0; p < nparts; p++) {
    ke[p] = 0.5*(up[p]*up[p] + vp[p]*vp[p] + wp[p]*wp[p]);
    //ke[p] = 0.5*wp[p]*wp[p];
  }

  free(x);
  free(y);
  free(z);
  
  cg_close(fn);
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

void write_coeffs(int in)
{
  char fnameEven[CHAR_BUF_SIZE] = "";
  char fnameOdd[CHAR_BUF_SIZE] = "";

  // number density
  sprintf(fnameEven, "%s/%s/number-dens-coeffs-even", ANALYSIS_DIR, DATA_DIR);
  sprintf(fnameOdd, "%s/%s/number-dens-coeffs-odd", ANALYSIS_DIR, DATA_DIR);
  FILE *fileEven = fopen(fnameEven, "a");
  FILE *fileOdd = fopen(fnameOdd, "a");
  if (fileEven == NULL) {
    printf("Error opening file %s!\n", fnameEven);
    exit(EXIT_FAILURE);
  } else if (fileOdd == NULL) {
    printf("Error opening file %s!\n", fnameOdd);
    exit(EXIT_FAILURE);
  }

  if (in == 0) {            // write constant coeffs
    int cc = 0 + npoints*tt;
    fprintf(fileEven, "%lf", 0.5*n_ces[cc]);
    fprintf(fileOdd, "%lf", 0.5*n_ces[cc]);
  } else if (in == -1) {    // write the rest
    for (int i = 1; i <= order; i++) {
      fprintf(fileEven, " %lf", nl_even[i]);
      fprintf(fileOdd, " %lf", nl_odd[i]);
    }
    fprintf(fileEven, "\n");
    fprintf(fileOdd, "\n");
  }

  fclose(fileEven);
  fclose(fileOdd);

  // part-w coeffs
  sprintf(fnameEven, "%s/%s/part-w-coeffs-even", ANALYSIS_DIR, DATA_DIR);
  sprintf(fnameOdd, "%s/%s/part-w-coeffs-odd", ANALYSIS_DIR, DATA_DIR);
  fileEven = fopen(fnameEven, "a");
  fileOdd = fopen(fnameOdd, "a");
  if (fileEven == NULL) {
    printf("Error opening file %s!\n", fnameEven);
    exit(EXIT_FAILURE);
  } else if (fileOdd == NULL) {
    printf("Error opening file %s!\n", fnameOdd);
    exit(EXIT_FAILURE);
  }

  if (in == 0) {            // write constant coeffs
    int cc = 0 + npoints*tt;
    fprintf(fileEven, "%lf", 0.5*nw_ces[cc]);
    fprintf(fileOdd, "%lf", 0.5*nw_ces[cc]);
  } else if (in == -1) {    // write the rest
    for (int i = 1; i <= order; i++) {
      fprintf(fileEven, " %lf", nwl_even[i]);
      fprintf(fileOdd, " %lf", nwl_odd[i]);
    }
    fprintf(fileEven, "\n");
    fprintf(fileOdd, "\n");
  }

  fclose(fileEven);
  fclose(fileOdd);
}

// Write reconstructed data
void write_reconstruct(void)
{
  char fname[CHAR_BUF_SIZE] = "";

  /* number desnity */
  sprintf(fname, "%s/%s/number-density", ANALYSIS_DIR, DATA_DIR);
  FILE *file = fopen(fname, "w");
  if (file == NULL) {
    printf("Error opening file %s!\n", fname);
    exit(EXIT_FAILURE);
  }

  // each reconstructed timestep is a row
  for (int t = 0; t < nFiles; t++) {
    for (int zz = 0; zz < npoints; zz++) {
      int cc = zz + t*npoints;
      fprintf(file, "%lf ", n_ces[cc]);
    }
    fprintf(file, "\n");
  }
  fclose(file);

  /* volume fraction */
  sprintf(fname, "%s/%s/volume-fraction", ANALYSIS_DIR, DATA_DIR);
  file = fopen(fname, "w");
  if (file == NULL) {
    printf("Error opening file %s!\n", fname);
    exit(EXIT_FAILURE);
  }

  // each reconstructed timestep is a row
  for (int t = 0; t < nFiles; t++) {
    for (int zz = 0; zz < npoints; zz++) {
      int cc = zz + t*npoints;
      fprintf(file, "%lf ", vFrac_ces[cc]);
    }
    fprintf(file, "\n");
  }
  fclose(file);

  /* part-u */
  sprintf(fname, "%s/%s/part-u", ANALYSIS_DIR, DATA_DIR);
  file = fopen(fname, "w");
  if (file == NULL) {
    printf("Error opening file %s!\n", fname);
    exit(EXIT_FAILURE);
  }

  // each reconstructed timestep is a row
  for (int t = 0; t < nFiles; t++) {
    for (int zz = 0; zz < npoints; zz++) {
      int cc = zz + t*npoints;
      fprintf(file, "%lf ", nu_ces[cc]);
    }
    fprintf(file, "\n");
  }
  fclose(file);
  
  /* part-v */
  sprintf(fname, "%s/%s/part-v", ANALYSIS_DIR, DATA_DIR);
  file = fopen(fname, "w");
  if (file == NULL) {
    printf("Error opening file %s!\n", fname);
    exit(EXIT_FAILURE);
  }

  // each reconstructed timestep is a row
  for (int t = 0; t < nFiles; t++) {
    for (int zz = 0; zz < npoints; zz++) {
      int cc = zz + t*npoints;
      fprintf(file, "%lf ", nv_ces[cc]);
    }
    fprintf(file, "\n");
  }
  fclose(file);

  /* part-w */
  sprintf(fname, "%s/%s/part-w", ANALYSIS_DIR, DATA_DIR);
  file = fopen(fname, "w");
  if (file == NULL) {
    printf("Error opening file %s!\n", fname);
    exit(EXIT_FAILURE);
  }

  // each reconstructed timestep is a row
  for (int t = 0; t < nFiles; t++) {
    for (int zz = 0; zz < npoints; zz++) {
      int cc = zz + t*npoints;
      fprintf(file, "%lf ", nw_ces[cc]);
    }
    fprintf(file, "\n");
  }
  fclose(file);

  /* kinetic energy */
  sprintf(fname, "%s/%s/part-kinetic-energy", ANALYSIS_DIR, DATA_DIR);
  file = fopen(fname, "w");
  if (file == NULL) {
    printf("Error opening file %s!\n", fname);
    exit(EXIT_FAILURE);
  }

  // each reconstructed timestep is a row
  for (int t = 0; t < nFiles; t++) {
    for (int zz = 0; zz < npoints; zz++) {
      int cc = zz + t*npoints;
      fprintf(file, "%lf ", nke_ces[cc]);
    }
    fprintf(file, "\n");
  }
  fclose(file);
}

// Free parts
void free_vars(void)
{
  for (int i = 0; i < nFiles; i++) {
    free(partFiles[i]);
  }
  free(partFiles);
  free(partFileMap);
  free(partFileTime);
  free(SIM_ROOT_DIR);
  free(ANALYSIS_DIR);

  free(parts);
  free(up);
  free(vp);
  free(wp);
  free(ke);

  free(evalZ);

  free(nl_even);
  free(nl_odd);
  free(nul_even);
  free(nul_odd);
  free(nvl_even);
  free(nvl_odd);
  free(nwl_even);
  free(nwl_odd);
  free(nkel_even);
  free(nkel_odd);

  free(nwl_avg_even);
  free(nwl_avg_odd);

  free(n_ces);
  free(vFrac_ces);
  free(nu_ces);
  free(nv_ces);
  free(nw_ces);
  free(nke_ces);

  free(nw_avg_ces);

  free(ones);
}

