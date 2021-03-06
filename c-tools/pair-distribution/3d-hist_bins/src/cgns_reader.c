#include "main.h"
#include "cgns_reader.h"

int nFiles;
char **flowFiles;
double *flowFileTime;
int *flowFileMap;
int nparts;
double meanR;

part_struct *parts;
dom_struct dom;
BC bc;

int nBinsTh;
int nBinsR;
double cutoff;
int gtr;
int less;
double dr;
double dth;
double *evalTh;
double *evalR;
double *gHist;
double *volume_fraction;

// set up directory structure
void directory_init(int argc, char *argv[])
{
  SIM_ROOT_DIR = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
  ANALYSIS_DIR = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));

  // arg[0] = program name
  // arg[1] = sim_root_dir
  if (argc == 2) {
    sprintf(SIM_ROOT_DIR, "%s", argv[1]);
    sprintf(ANALYSIS_DIR, "%s/analysis/%s/%s", SIM_ROOT_DIR, PAIR_DIR, ANALYSIS);
  } else if (argc != 2) {
    printf("usage: %s sim_root_dir\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  #ifdef DEBUG
    printf("\n SIM_ROOT_DIR = %s\n", SIM_ROOT_DIR);
    printf(" ANALYSIS_DIR = %s\n\n", ANALYSIS_DIR);
    fflush(stdout);
  #endif
}

// read main.config input file
void main_read_input(void)
{
  int fret = 0;
  fret = fret;

  // open config file for reading
  char fname[CHAR_BUF_SIZE] = "";
  sprintf(fname, "%s/%s", ANALYSIS_DIR, CONFIG_FILE);
  FILE *infile = fopen(fname, "r");
  if (infile == NULL) {
    printf("could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  } else {
    #ifdef DEBUG
      printf("  Reading config from %s...\n", fname);
      fflush(stdout);
    #endif
  }
  
  // read input
  fret = fscanf(infile, "tStart %lf\n", &tStart);
  fret = fscanf(infile, "tEnd %lf\n", &tEnd);
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "R0/a %lf\n", &R0);
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "ntheta bins %d\n", &nBinsTh);
  fret = fscanf(infile, "nr bins %d\n", &nBinsR);
  fret = fscanf(infile, "Greater Than %d\n", &gtr);
  fret = fscanf(infile, "Less Than %d\n", &less);
  fret = fscanf(infile, "Cutoff (sdev) %lf\n", &cutoff);
  fclose(infile);

  // Checks on gtr, less flags
  if (gtr > 1 || gtr < 0 || less < 0 || less > 1) {
    printf("Greater and Less Than flags must be 0 or 1\n");
    exit(EXIT_FAILURE);
  }
  if (gtr ^ less) {
  } else {
    printf("Greater Than is %d and Less Than is %d\n", gtr, less);
    printf("Only one of the two may be 1, the other must be 0\n");
    exit(EXIT_FAILURE);
  }
  if (gtr) {
    printf("  Running for greater than %.2lf standard deviations from the mean\n",
      cutoff);
  } else {
    printf("  Running for less than %.2lf standard deviations from the mean\n",
      cutoff);
  }
  printf("  %d bins in theta, %d bins in r\n", nBinsTh, nBinsR);
}

// read and sort flow files directory
void init_files(void) {
  DIR *dir;
  struct dirent *ent;
  char cgns_file[FILE_NAME_SIZE] = "";
  char file_name[FILE_NAME_SIZE] = "";
  char output_path[FILE_NAME_SIZE] = "";
  int fret = 0; fret=fret;
  double time;

  sprintf(cgns_file, "f-rec-3D-");
  sprintf(file_name, "%s%%lf.cgns", cgns_file);
  sprintf(output_path, "%s/analysis/fourier-reconstruction/3-dim-part/data", 
    SIM_ROOT_DIR);

  int isFlow;
  int inRange;
  nFiles = 0;

  // count number of flow files in directory that fall in time range
  if ((dir = opendir (output_path)) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      // check if flow file (0 if match)
      isFlow = (strncmp(ent->d_name, cgns_file, strlen(cgns_file)) == 0);

      if (isFlow == 1) {
        // check if in time range
        fret = sscanf(ent->d_name, file_name, &time);
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
      isFlow = (strncmp(ent->d_name, cgns_file, strlen(cgns_file)) == 0);

      if (isFlow == 1) {
        // check if in time range
        fret = sscanf(ent->d_name, file_name, &time);
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
  if (nFiles == 0) {
    printf("  Found %d flow files in range [%lf, %lf]\n", nFiles, tStart, tEnd);
    printf("    cgns_file = %s\n", cgns_file);
    printf("    file_name = %s\n", file_name);
    printf("    output_path = %s\n", output_path);
    printf("  Quitting...\n");
    exit(EXIT_FAILURE);
  } else {
    printf("  Found %d flow files in range [%lf, %lf]\n", nFiles, tStart, tEnd);
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

// Create directory for output data
void create_output(void) 
{
  // Create output directory if it doesn't exist
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
  sprintf(path2file, "%s/%s/bininfo-r%d-th%d", ANALYSIS_DIR, DATA_DIR, nBinsR, 
    nBinsTh);
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
  fprintf(file, "evalR ");
  for (int i = 0; i < nBinsR; i++) {
    fprintf(file, "%lf ", evalR[i]);
  }
  fprintf(file, "\n");

  // Print evalTh on second row
  fprintf(file, "evalTh ");
  for (int i = 0; i < nBinsTh; i++) {
    fprintf(file, "%lf ", evalTh[i]);
  }

  fclose(file);
}

// initialize part_struct
void parts_init(void)
{
  // REad nparts from part.config file
  int fret = 0;
  fret = fret;  // prevent compiler warning

  // read nparts from part.config
  char fname[FILE_NAME_SIZE] = "";
  sprintf(fname, "%s/input/part.config", SIM_ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if (infile == NULL) {
    printf("Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "n %d\n", &nparts);
  fclose(infile);

  printf("  Found %d particles\n", nparts);
  fflush(stdout);

  // Init part_struct
  parts = (part_struct*) malloc(nparts * sizeof(part_struct));

  for(int p = 0; p < nparts; p++) {
    parts[p].x = -1;
    parts[p].y = -1;
    parts[p].z = -1;
    parts[p].r = -1;
  }

  // Create part-%lf.cgns from flow-%lf.cgns
  char timeString1[CHAR_BUF_SIZE] = "";
  sprintf(timeString1, flowFiles[flowFileMap[0]]);

  // We expect f-rec-3D-%lf.cgns, so split along "-"
  char *timeString2;
  timeString2 = strtok(timeString1, "-");   // f
  timeString2 = strtok(NULL, "-");          // rec
  timeString2 = strtok(NULL, "-");          // 3D
  timeString2 = strtok(NULL, "-");          // %lf.cgns

  // Create part file with directory path
  char partFile[FILE_NAME_SIZE];
  sprintf(partFile, "%s/%s/part-%s", SIM_ROOT_DIR, OUTPUT_DIR, timeString2);

  // Open cgns file and get cgns file index number fn
  int fn;
  int ier = cg_open(partFile, CG_MODE_READ, &fn);
  if (ier != 0) {
    printf("CGNS error opening %s in %s at %d\n", partFile, __FILE__, __LINE__);
    printf("  Check existance of grid.cgns\n");
    cg_error_exit();
  }

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
  if (ier != 0) {
    printf("CGNS error opening %s in %s at %d\n", partFile, __FILE__, __LINE__);
    printf("  Check existance of grid.cgns\n");
    cg_error_exit();
  }

  // Calculate mean radius
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
  
  fclose(infile);
}

void histogram_init(void)
{
  /* R_0/a input -- need to convert to distance
    * R_0 <= 0,                       use 0.5*(minimum domain extent)
    * 0 < R_0*meanR <= (min extent),  use R_0*meanR
    * (min extent) < R_0*mean,        use (min extent)
  */

  double minExtent = (dom.xl <= dom.yl) ? dom.xl : dom.yl;
  minExtent = (minExtent <= dom.zl) ? minExtent : dom.zl;
  if (R0 > 0.) {
    R0 *= meanR;
    if (R0 > 0.5*minExtent) {
      R0 = 0.5*minExtent;
    }
  }
  else if (R0 <= 0.) {
    R0 = 0.5*minExtent;
  }

  // Init bins
  evalR = (double*) malloc(nBinsR * sizeof(double));
  evalTh = (double*) malloc(nBinsTh * sizeof(double));
  gHist = (double*) malloc(nBinsR*nBinsTh * sizeof(double));

  // Calculate bin sizes given nBinsTh, nBinsR
  dr = (R0 - 2.*meanR)/nBinsR;
  dth = 0.5*PI/nBinsTh;

  // init locations for eval reconstruct
  for (int ii = 0; ii < nBinsR; ii++) {
    for (int jj = 0; jj < nBinsTh; jj++) {
      gHist[jj + ii*nBinsTh] = 0.;
      evalTh[jj] = jj*dth;
    }
    evalR[ii] = 2.*meanR + ii*dr;  
  }

  #ifdef DEBUG
    show_domain();
  #endif

  // Init vol fraction
  volume_fraction = (double*) malloc(dom.Gcc.s3 * sizeof(double));
}

// Read part_struct data
void cgns_fill(void)
{
  /* Read part-%lf.cgns */
  // Create part-%lf.cgns from flow-%lf.cgns
  char timeString1[CHAR_BUF_SIZE] = "";
  sprintf(timeString1, flowFiles[flowFileMap[tt]]);

  // We expect f-rec-3D-%lf.cgns, so split along "-"
  char *timeString2;
  timeString2 = strtok(timeString1, "-");   // f
  timeString2 = strtok(NULL, "-");          // rec
  timeString2 = strtok(NULL, "-");          // 3D
  timeString2 = strtok(NULL, "-");          // %lf.cgns

  // Create part file with directory path
  char partFile[FILE_NAME_SIZE];
  sprintf(partFile, "%s/%s/part-%s", SIM_ROOT_DIR, OUTPUT_DIR, timeString2);

  // Open cgns file and get cgns file index number fn
  int fn;
  int ier = cg_open(partFile, CG_MODE_READ, &fn);
  if (ier != 0) {
    printf("CGNS error opening %s in %s at %d\n", partFile, __FILE__, __LINE__);
    printf("  Check existance of grid.cgns\n");
    cg_error_exit();
  }
  
  // Set base, zone, and solutions index numbers
  int bn = 1;
  int zn = 1;

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

  ier = cg_coord_read(fn,bn,zn, "CoordinateX", RealDouble, &range_min, &range_max, x);
  ier = cg_coord_read(fn,bn,zn, "CoordinateY", RealDouble, &range_min, &range_max, y);
  ier = cg_coord_read(fn,bn,zn, "CoordinateZ", RealDouble, &range_min, &range_max, z);
  if (ier != 0) {     // only error check one... lazy bum
    printf("CGNS error opening %s in %s at %d\n", partFile, __FILE__, __LINE__);
    printf("  Check existance of grid.cgns\n");
    cg_error_exit();
  }

  for (int p = 0; p < nparts; p++) {
    parts[p].x = x[p];
    parts[p].y = y[p];
    parts[p].z = z[p];
  }

  free(x);
  free(y);
  free(z);
  
  cg_close(fn);

  /* Read f-rec-3D-%lf.cgns */
  char buf[FILE_NAME_SIZE] = "";
  sprintf(buf, "%s/analysis/fourier-reconstruction/3-dim-part", SIM_ROOT_DIR);
  sprintf(buf, "%s/data/f-rec-3D-%s", buf, timeString2);
  ier = cg_open(buf, CG_MODE_READ, &fn);
  if (ier != 0) {
    printf("CGNS error opening %s in %s at %d\n", buf, __FILE__, __LINE__);
    printf("  Check existance of grid.cgns\n");
    cg_error_exit();
  }

  // Set base, zone, and solution index numbers
  bn = 1;
  zn = 1;
  int sn = 1;

  // Size of array to pull
  cgsize_t range_min_3d[3];
  range_min_3d[0] = 1;
  range_min_3d[1] = 1;
  range_min_3d[2] = 1;
  cgsize_t range_max_3d[3];
  range_max_3d[0] = dom.xn;
  range_max_3d[1] = dom.yn;
  range_max_3d[2] = dom.zn;

  // Read and fill
  ier = cg_field_read(fn,bn,zn,sn,"Volume Fraction Real", RealDouble, range_min_3d, 
    range_max_3d, volume_fraction);
  if (ier != 0) {
    printf("CGNS error reading %s in %s at %d\n", buf, __FILE__, __LINE__);
    printf("  Check existance of grid.cgns\n");
    cg_error_exit();
  }

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

  printf("Parameters\n");
  printf("  nparts %d\n", nparts);
  printf("  tStart %lf\n", tStart);
  printf("  tEnd %lf\n", tEnd);
  printf("  R0 %lf\n", R0);
  printf("\n");

  printf("  nBinsR = %d\n", nBinsR);
  printf("  nBinsTh = %d\n", nBinsTh);
  printf("  dr = %lf\n", dr);
  printf("  dth = %lf\n", dth);
  printf("\n");

  printf("  Cutoff = %lf\n", cutoff);
  printf("  Gtr Flag = %d\n", gtr);
  printf("  Less Flag = %d\n", less);
  printf("\n");
}

// Write averaged data
void write_avg_reconstruct(void)
{
  // Set up the filename
  char fname[CHAR_BUF_SIZE] = "";
  if (gtr == 1) {
    sprintf(fname, "%s/%s/bins-r%d-th%d_gtr-%.1fsdev", ANALYSIS_DIR, DATA_DIR, 
      nBinsR, nBinsTh, cutoff);
  } else if (less == 1) {
    sprintf(fname, "%s/%s/bins-r%d-th%d_lss-%.1fsdev", ANALYSIS_DIR, DATA_DIR, 
      nBinsR, nBinsTh, cutoff);
  }

  FILE *fdat = fopen(fname, "w");
  if (fdat == NULL) {
    printf("Error opening file %s!\n", fname);
    exit(EXIT_FAILURE);
  }

  // print reconstructed data -- each r is a row
  for (int rr = 0; rr < nBinsR; rr++) {
    for (int th = 0; th < nBinsTh; th++) {
      int cc = th + rr*nBinsTh;
      fprintf(fdat, "%lf ", gHist[cc]);
    }
    fprintf(fdat, "\n");
  }
  fclose(fdat);

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

  free(parts);

  free(evalR);
  free(evalTh);
  free(gHist);
  free(volume_fraction);

  free(SIM_ROOT_DIR);
  free(ANALYSIS_DIR);
}
