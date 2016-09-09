#include "main.h"
#include "cgns_reader.h"

// File variables
int nFiles;
char **files;
double *fileTime;
int *fileMap;

bin_struct bins;

// part stuff
int nparts;
double mass;
double meanR;
double meanI;
double *up;
double *vp;
double *wp;
double *ox;
double *oy;
double *oz;
double *U;
double *ke;
double *ke_rot;
double *ke_trans;
double *T;
double *T_perp;
double *T_z;

// hists
int nBins;
int *hist_U;
int *hist_ke;
int *hist_ke_rot;
int *hist_ke_trans;
int *hist_T;
int *hist_T_perp;
int *hist_T_z;
int *hist_ux;
int *hist_vy;
int *hist_wz;

// Set up directory structure
void directory_init(int argc, char *argv[])
{
  SIM_ROOT_DIR = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
  ANALYSIS_DIR = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));

  // arg[0] = program name
  // arg[1] = SIM_ROOT_DIR
  if (argc == 2) {
  sprintf(SIM_ROOT_DIR, "%s", argv[1]);
  sprintf(ANALYSIS_DIR, "%s/analysis/%s/%s", SIM_ROOT_DIR, PHASE_DIR, ANALYSIS);
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
  fret = fscanf(infile, "nBins %d\n", &nBins);
  fclose(infile);
}

// read and sort files directory
void init_files(void) 
{
  DIR *dir;
  struct dirent *ent;
  char output_path[FILE_NAME_SIZE] = "";  // where the cgns files are stored
  char file_name[FILE_NAME_SIZE] = "";    // CGNS_FILE%lf.cgns
  int fret = 0; fret=fret;
  double time;

  sprintf(output_path, "%s/output", SIM_ROOT_DIR);
  sprintf(file_name, "%s%%lf.cgns", CGNS_FILE);

  int isFile;
  int inRange;
  nFiles = 0;

  // count number of files in directory that fall in time range
  if ((dir = opendir (output_path)) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      // check if file (0 if match)
      isFile = (strncmp(ent->d_name, CGNS_FILE, strlen(CGNS_FILE)) == 0);

      if (isFile == 1) {
        // check if in time range
        fret = sscanf(ent->d_name, file_name, &time);
        inRange = ((time >= tStart) & (time <= tEnd));
        nFiles += isFile*inRange;
      } else {
        continue;
      }
    }
    closedir (dir);
  } else {
    printf("Output directory %s does not exist! (Line %d in %s)\n", 
      output_path, __LINE__, __FILE__);
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
      isFile = (strncmp(ent->d_name, CGNS_FILE, strlen(CGNS_FILE)) == 0);

      if (isFile == 1) {
        // check if in time range
        fret = sscanf(ent->d_name, file_name, &time);
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
    printf("Output directory %s does not exist! (Line %d in %s)\n", 
      output_path, __LINE__, __FILE__);
    exit(EXIT_FAILURE);
  }

  // Sort the resulting array by time
  // create temporary array to sort files by
  fileMap = malloc(nFiles * sizeof(double));
  for (int i = 0; i < nFiles; i++) {
    fileMap[i] = i;
  }

  merge_sort(fileTime, nFiles, fileMap);
  printf("Found %d files in range [%lf, %lf]\n", nFiles, tStart, tEnd);
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

  // Set up variables
  up = (double*) malloc(nparts * sizeof(double));
  vp = (double*) malloc(nparts * sizeof(double));
  wp = (double*) malloc(nparts * sizeof(double));
  ox = (double*) malloc(nparts * sizeof(double));
  oy = (double*) malloc(nparts * sizeof(double));
  oz = (double*) malloc(nparts * sizeof(double));
  U = (double*) malloc(nparts * sizeof(double));
  ke = (double*) malloc(nparts * sizeof(double));
  ke_rot = (double*) malloc(nparts * sizeof(double));
  ke_trans = (double*) malloc(nparts * sizeof(double));
  T = (double*) malloc(nparts * sizeof(double));
  T_perp = (double*) malloc(nparts * sizeof(double));
  T_z = (double*) malloc(nparts * sizeof(double));
  for (int i = 0; i < nparts; i++) {
    up[i] = 0.; 
    vp[i] = 0.;
    wp[i] = 0.;
    ox[i] = 0.;
    oy[i] = 0.;
    oz[i] = 0.;
    U[i] = 0.;
    ke[i] = 0.;
    ke_rot[i] = 0.;
    ke_trans[i] = 0.;
    T[i] = 0.;
    T_perp[i] = 0.;
    T_z[i] = 0.;
  }
  hist_U = (int*) malloc(nBins * sizeof(int));
  hist_ke = (int*) malloc(nBins * sizeof(int));
  hist_ke_rot = (int*) malloc(nBins * sizeof(int));
  hist_ke_trans = (int*) malloc(nBins * sizeof(int));
  hist_T = (int*) malloc(nBins * sizeof(int));
  hist_T_perp = (int*) malloc(nBins * sizeof(int));
  hist_T_z = (int*) malloc(nBins * sizeof(int));
  hist_ux = (int*) malloc(nBins * sizeof(int));
  hist_vy = (int*) malloc(nBins * sizeof(int));
  hist_wz = (int*) malloc(nBins * sizeof(int));
  for (int i = 0; i < nBins; i++) {
    hist_U[i] = 0;
    hist_ke[i] = 0;
    hist_ke_rot[i] = 0;
    hist_ke_trans[i] = 0;
    hist_T[i] = 0;
    hist_T_perp[i] = 0;
    hist_T_z[i] = 0;
    hist_ux[i] = 0;
    hist_vy[i] = 0;
    hist_wz[i] = 0;
  }

  // Init parts
  char buf[FILE_NAME_SIZE];
  sprintf(buf, "%s/%s/%s", SIM_ROOT_DIR,OUTPUT_DIR, files[fileMap[0]]);
  int fn;
  cg_open(buf, CG_MODE_READ, &fn);

  // SEt base, zone, soln
  int bn = 1;
  int zn = 1;
  int sn = 1;

  cgsize_t range_min = 1;
  cgsize_t range_max = nparts;

  // Read part radius
  double *r = malloc(nparts * sizeof(double));
  double *rho = malloc(nparts * sizeof(double));
  cg_field_read(fn,bn,zn,sn, "Radius", RealDouble, &range_min, &range_max, r);
  cg_field_read(fn,bn,zn,sn, "Density", RealDouble, &range_min, &range_max, rho);

  // Calculate meanR, rho, meanI, mass
  meanR = 0.;
  double meanRho = 0.;
  for (int p = 0; p < nparts; p++) {
    meanR += r[p];
    meanRho += rho[p];
  }
  meanR /= nparts;
  meanRho /= nparts;
  mass = 4./3.*PI*meanR*meanR*meanR * meanRho;
  meanI = 2.*mass*meanR*meanR/5.;

  cg_close(fn);
  free(r);
  free(rho);
}

// Read data
void cgns_fill(void)
{
  // Open cgns file and get cgns file index number fn
  char buf[FILE_NAME_SIZE];
  sprintf(buf, "%s/output/%s", SIM_ROOT_DIR, files[fileMap[tt]]);
  int fn;
  int ier = cg_open(buf, CG_MODE_READ, &fn);
  if (ier != 0 ) {
    printf("CGNS Error - double check grid.cgns exists");
    printf(" in same directory as %s\n", buf);
    cg_error_exit();
  }
  fflush(stdout);
  
  // Set base, zone, and solutions index numbers
  int bn = 1;
  int zn = 1;
  int sn = 1;

  // Size of array to pull
  cgsize_t range_min = 1;
  cgsize_t range_max = nparts;;

  // Read and fill
  cg_field_read(fn,bn,zn,sn, "VelocityX", RealDouble, &range_min, 
    &range_max, up);
  cg_field_read(fn,bn,zn,sn, "VelocityY", RealDouble, &range_min, 
    &range_max, vp);
  cg_field_read(fn,bn,zn,zn, "VelocityZ", RealDouble, &range_min, 
    &range_max, wp);
  cg_field_read(fn,bn,zn,sn, "AngularVelocityX", RealDouble, &range_min, 
    &range_max, ox);
  cg_field_read(fn,bn,zn,sn, "AngularVelocityY", RealDouble, &range_min, 
    &range_max, oy);
  cg_field_read(fn,bn,zn,zn, "AngularVelocityZ", RealDouble, &range_min, 
    &range_max, oz);

  // Calculate mean
  double up_mean = 0.;
  double vp_mean = 0.;
  double wp_mean = 0.;
  double inparts = 1./nparts;
  for (int p = 0; p < nparts; p++) {
    up_mean += inparts*up[p];
    vp_mean += inparts*vp[p];
    wp_mean += inparts*wp[p];
  }

  for (int p = 0; p < nparts; p++) {
    // subtract mean
    up[p] -= up_mean;
    vp[p] -= vp_mean;
    wp[p] -= wp_mean;

    // Calculate temperatures
    T_perp[p] = mass*(up[p]*up[p] + vp[p]*vp[p]);
    T_z[p]    = mass*(wp[p]*wp[p]);
    T[p]      = (T_perp[p] + T_z[p])/3.;

    // Speed
    U[p] = sqrt(up[p]*up[p] + vp[p]*vp[p] + wp[p]*wp[p]);

    // Ke
    ke_rot[p] = 0.5*meanI*(ox[p]*ox[p] + oy[p]*oy[p] + oz[p]*oz[p]);
    ke_trans[p] = 1.5*T[p];
    ke[p] = ke_rot[p] + ke_trans[p];
  }

  cg_close(fn);
}

// Show domain
void show_domain(void)
{
  printf("Input Parameters\n");
  printf("  tStart %lf\n", tStart);
  printf("  tEnd %lf\n", tEnd);
  printf("  nBins %d\n", nBins);
}

// Write reconstructed data
void write_field(void)
{
  // set up file
  char fname[FILE_NAME_SIZE] = "";

  // INFO
  sprintf(fname, "%s/%s/%s", ANALYSIS_DIR, DATA_DIR, INFO_FILE);
  
  FILE *file = fopen(fname, "w");
  if (file == NULL) {
    printf("Error opening file %s! (%s at %d)\n", fname, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  fprintf(file, "array min dbin max mean\n");
  fprintf(file, "U %e %e %e %e\n", 
    bins.U.min, bins.U.dBin, bins.U.max, bins.U.mean);
  fprintf(file, "ke %e %e %e %e\n", 
    bins.ke.min, bins.ke.dBin, bins.ke.max, bins.ke.mean);
  fprintf(file, "ke_rot %e %e %e %e\n", 
    bins.ke_rot.min, bins.ke_rot.dBin, bins.ke_rot.max, bins.ke_rot.mean);
  fprintf(file, "ke_trans %e %e %e %e\n", 
    bins.ke_trans.min, bins.ke_trans.dBin, bins.ke_trans.max, bins.ke_trans.mean);
  fprintf(file, "T_z %e %e %e %e\n",
    bins.T_z.min, bins.T_z.dBin, bins.T_z.max, bins.T_z.mean);
  fprintf(file, "T_perp %e %e %e %e\n",
    bins.T_perp.min, bins.T_perp.dBin, bins.T_perp.max, bins.T_perp.mean);
  fprintf(file, "T %e %e %e %e\n",
    bins.T.min, bins.T.dBin, bins.T.max, bins.T.mean);
  fprintf(file, "ux %e %e %e %e\n",
    bins.ux.min, bins.ux.dBin, bins.ux.max, bins.ux.mean);
  fprintf(file, "vy %e %e %e %e\n",
    bins.vy.min, bins.vy.dBin, bins.vy.max, bins.vy.mean);
  fprintf(file, "wz %e %e %e %e\n",
    bins.wz.min, bins.wz.dBin, bins.wz.max, bins.wz.mean);

  fclose(file);

  // velocity counts
  sprintf(fname, "%s/%s/%s", ANALYSIS_DIR, DATA_DIR, "U_hist");
  file = fopen(fname, "w");
  for (int c = 0; c < nBins; c++) {
    fprintf(file, "%d ", hist_U[c]);
  }
  fclose(file);

  // ke counts
  sprintf(fname, "%s/%s/%s", ANALYSIS_DIR, DATA_DIR, "ke_hist");
  file = fopen(fname, "w");
  for (int c = 0; c < nBins; c++) {
    fprintf(file, "%d ", hist_ke[c]);
  }
  fclose(file);

  // ke_rot counts
  sprintf(fname, "%s/%s/%s", ANALYSIS_DIR, DATA_DIR, "ke_hist_rot");
  file = fopen(fname, "w");
  for (int c = 0; c < nBins; c++) {
    fprintf(file, "%d ", hist_ke_rot[c]);
  }
  fclose(file);

  // ke_trans counts
  sprintf(fname, "%s/%s/%s", ANALYSIS_DIR, DATA_DIR, "ke_hist_trans");
  file = fopen(fname, "w");
  for (int c = 0; c < nBins; c++) {
    fprintf(file, "%d ", hist_ke_trans[c]);
  }
  fclose(file);

  // T_z counts
  sprintf(fname, "%s/%s/%s", ANALYSIS_DIR, DATA_DIR, "T_hist_z");
  file = fopen(fname, "w");
  for (int c = 0; c < nBins; c++) {
    fprintf(file, "%d ", hist_T_z[c]);
  }
  fclose(file);

  // T_perp counts
  sprintf(fname, "%s/%s/%s", ANALYSIS_DIR, DATA_DIR, "T_hist_perp");
  file = fopen(fname, "w");
  for (int c = 0; c < nBins; c++) {
    fprintf(file, "%d ", hist_T_perp[c]);
  }
  fclose(file);

  // T counts
  sprintf(fname, "%s/%s/%s", ANALYSIS_DIR, DATA_DIR, "T_hist");
  file = fopen(fname, "w");
  for (int c = 0; c < nBins; c++) {
    fprintf(file, "%d ", hist_T[c]);
  }
  fclose(file);

  // ux counts
  sprintf(fname, "%s/%s/%s", ANALYSIS_DIR, DATA_DIR, "ux_hist");
  file = fopen(fname, "w");
  for (int c = 0; c < nBins; c++) {
    fprintf(file, "%d ", hist_ux[c]);
  }
  fclose(file);

  // vy counts
  sprintf(fname, "%s/%s/%s", ANALYSIS_DIR, DATA_DIR, "vy_hist");
  file = fopen(fname, "w");
  for (int c = 0; c < nBins; c++) {
    fprintf(file, "%d ", hist_vy[c]);
  }
  fclose(file);

  // wz counts
  sprintf(fname, "%s/%s/%s", ANALYSIS_DIR, DATA_DIR, "wz_hist");
  file = fopen(fname, "w");
  for (int c = 0; c < nBins; c++) {
    fprintf(file, "%d ", hist_wz[c]);
  }
  fclose(file);
}

// Free vars
void free_vars(void)
{
  free(ANALYSIS_DIR);
  free(SIM_ROOT_DIR);

  for (int i = 0; i < nFiles; i++) {
    free(files[i]);
  }
  free(files);
  free(fileMap);
  free(fileTime);

  free(up);
  free(vp);
  free(wp);
  free(ox);
  free(oy);
  free(oz);
  free(ke);
  free(ke_rot);
  free(ke_trans);
  free(U);
  free(hist_U);
  free(hist_ke);
  free(hist_ke_rot);
  free(hist_ke_trans);
  free(T_perp);
  free(T_z);
  free(T);
  free(hist_T);
  free(hist_T_z);
  free(hist_T_perp);
  free(hist_ux);
  free(hist_vy);
  free(hist_wz);
}
