#include "main.h"
#include "cgns_reader.h"

int nFiles;
char **files;
double *fileTime;
int *fileMap;
int sigFigPre;
int sigFigPost;
dom_struct dom;

double *volume_fraction;
double *moving_avg;
double *stationary_avg;

// Set up script for batch run
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
  fret = fscanf(infile, "Initial Location %lf\n", &initial_location);
  fret = fscanf(infile, "Total Averaging Distance (2z/a) %lf\n", &averaging_distance);
  fret = fscanf(infile, "Radius %lf\n", &radius);
  fret = fscanf(infile, "Wavespeed (mm/ms) %lf", &wavespeed);

  // Un-normalize averaging distance -- is equal to distance above or below
  averaging_distance *= radius*0.5;

  fclose(infile);
}

// read and sort flow files directory
void init_input_files(void) 
{
  DIR *dir;
  struct dirent *ent;
  char output_path[FILE_NAME_SIZE] = "";
  int fret = 0; fret=fret;
  double time;

  sprintf(output_path, "%s/data", ANALYSIS_DIR);

  int isFlow;
  int inRange;
  nFiles = 0;

  // count number of files in directory that fall in time range
  if ((dir = opendir (output_path)) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      // check if flow file (0 if match)
      isFlow = (strncmp(ent->d_name, "f-rec", 5) == 0);

      if (isFlow == 1) {
        // check if in time range
        fret = sscanf(ent->d_name, "f-rec-part-phase-3D-%lf.cgns", &time);
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
  files = (char**) malloc(nFiles * sizeof(char*));
  for (int i = 0; i < nFiles; i ++) {
    files[i] = (char*) malloc(FILE_NAME_SIZE*sizeof(char));
  }
  fileTime = (double*) malloc(nFiles * sizeof(double));

  int cc = 0;

  if ((dir = opendir (output_path)) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      isFlow = (strncmp(ent->d_name, "f-rec", 4) == 0);

      if (isFlow == 1) {
        // check if in time range
        fret = sscanf(ent->d_name, "f-rec-part-phase-3D-%lf.cgns", &time);
        inRange = ((time >= tStart) & (time <= tEnd));
        
        if (inRange == 1) {
          fret = sscanf(ent->d_name, "%s", files[cc]);
          fileTime[cc] = time;
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
  // create temporary array to sort files by
  fileMap = malloc(nFiles * sizeof(double));
  for (int i = 0; i < nFiles; i++) {
    fileMap[i] = i;
  }

  merge_sort(fileTime, nFiles, fileMap);
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

  fclose(infile);

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

  // init data arrays
  volume_fraction = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  moving_avg = (double*) malloc(dom.Gcc.s2 * sizeof(double));
  stationary_avg = (double*) malloc(dom.Gcc.s2 * sizeof(double));
  for (int i = 0; i < dom.Gcc.s3; i++) {
    volume_fraction[i] = 0.;
  }

  #ifdef DEBUG
    show_domain();
  #endif
}

// Read data
void cgns_fill_input(void)
{
  // Open cgns file and get cgns file index number fn
  char buf[FILE_NAME_SIZE];
  sprintf(buf, "%s/data/%s", ANALYSIS_DIR, files[fileMap[tt]]);
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
  cg_field_read(fn,bn,zn,sn, "Volume Fraction Real", RealDouble, range_min, 
    range_max, volume_fraction);

  cg_close(fn);

  // Zero the averages
  for (int i = 0; i < dom.Gcc.s2; i++) {
    moving_avg[i] = 0.;
    stationary_avg[i] = 0.;
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
  printf("  Initial Location %lf\n", initial_location);
  printf("  Total Averaging Distance (2z/a) %lf\n", averaging_distance);
}

// Get sigfigs
void get_sigfigs(void) 
{
  // print the last file name (longest time) to temp variable
  char tempLastFile[CHAR_BUF_SIZE] = "";
  sprintf(tempLastFile, files[fileMap[nFiles - 1]]);

  // We expect flow-*.*.cgns, so split along "-"
  // We expect f-rec-part-phase-3D-*.*.cgns, so split along "-"
  char *tempDash;
  tempDash = strtok(tempLastFile, "-"); // "f"
  tempDash = strtok(NULL, "-");         // "rec"
  tempDash = strtok(NULL, "-");         // "part"
  tempDash = strtok(NULL, "-");         // "phase"
  tempDash = strtok(NULL, "-");         // "3D"
  tempDash = strtok(NULL, "-");         // "*.*.cgns"

  // Now expect *.*.cgns, split along "."
  char *tempDot;
  tempDot = strtok(tempDash, ".");
  sigFigPre = strlen(tempDot);
  tempDot = strtok(NULL, ".");
  sigFigPost = strlen(tempDot);
}

// Write averaged data
void write_field(double curr_location)
{
  // Sig figs stuff
  char format[CHAR_BUF_SIZE] = "";
  char fnameall[CHAR_BUF_SIZE] = "";
  char fnameall2[CHAR_BUF_SIZE] = "";
  sprintf(format, "%%0%d.%df", sigFigPre + sigFigPost + 1, sigFigPost);

  /* Set up the filename for MOVING AVERAGE */
  sprintf(fnameall2, "%s/%s/vfrac-moving-avg-%s", ANALYSIS_DIR, DATA_DIR, 
    format);
  sprintf(fnameall, fnameall2, fileTime[tt]);

  // Open file error checking
  FILE *file = fopen(fnameall, "w");
  if (file == NULL) {
    printf("Error opening file %s!\n", fnameall);
    exit(EXIT_FAILURE);
  }

  // Write array
  for (int i = 0; i < dom.xn; i++) {
    for (int j = 0; j < dom.yn; j++) {
      // j is increasing across a row, i is increasing down columns 
      fprintf(file, "%lf ", moving_avg[j + dom.yn*i]);
    }
    fprintf(file, "\n");
  }
  fclose(file);

  /* Set up the filename for STATIONARY AVERAGE */
  sprintf(fnameall2, "%s/%s/vfrac-stationary-avg-%s", ANALYSIS_DIR, DATA_DIR, 
    format);
  sprintf(fnameall, fnameall2, fileTime[tt]);

  // Open file error checking
  file = fopen(fnameall, "w");
  if (file == NULL) {
    printf("Error opening file %s!\n", fnameall);
    exit(EXIT_FAILURE);
  }

  // Write array
  for (int i = 0; i < dom.xn; i++) {
    for (int j = 0; j < dom.yn; j++) {
      // j is increasing across a row, i is increasing down columns 
      fprintf(file, "%lf ", stationary_avg[j + dom.yn*i]);
    }
    fprintf(file, "\n");
  }
  fclose(file);

  // Fill up the info file with time and location
  char fname[CHAR_BUF_SIZE] = "";
  sprintf(fname, "%s/%s/info", ANALYSIS_DIR, DATA_DIR);

  if (tt == 0) {
    file = fopen(fname, "w");
    fprintf(file, "time z_distance\n");
    fclose(file);
  }
  file = fopen(fname, "a");
  fprintf(file, "%lf %lf\n", fileTime[tt], curr_location);

  fclose(file);
}

// Free vars
void free_vars(void)
{
  for (int i = 0; i < nFiles; i++) {
    free(files[i]);
  }
  free(files);
  free(fileMap);
  free(fileTime);

  free(SIM_ROOT_DIR);
  free(ANALYSIS_DIR);

  free(volume_fraction);
  free(moving_avg);
  free(stationary_avg);
}

