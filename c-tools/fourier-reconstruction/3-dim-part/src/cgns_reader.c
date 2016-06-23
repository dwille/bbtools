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
part_struct *parts;
dom_struct dom;

// Read main.config input file
void main_read_input(void)
{
  int fret = 0;
  fret = fret;

  // open config file for reading
  char fname[CHAR_BUF_SIZE] = "";
  sprintf(fname, "%s/f-rec.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  
  // read input
  fret = fscanf(infile, "tStart %lf\n", &tStart);
  fret = fscanf(infile, "tEnd %lf\n", &tEnd);
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "Fourier Order in (X,Y,Z) (%d,%d,%d)\n", &orderL, 
    &orderM, &orderN);
  fclose(infile);
}

// read and sort part files directory
void init_part_files(void) {
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
void create_output(void) {
  // Create DATA_DIR directory if it doesn't exist
  // From stackoverflow-7430248
  struct stat st = {0};
  char buf[CHAR_BUF_SIZE];
  sprintf(buf, "%s/%s", ROOT_DIR, DATA_DIR);
  if (stat(buf, &st) == -1) {
    mkdir(buf, 0700);
  }
  // Create DATA_SUBDIR if it doesn't exist
  sprintf(buf, "%s/%s/%s", ROOT_DIR, DATA_DIR, DATA_SUBDIR);
  if (stat(buf, &st) == -1) {
    mkdir(buf, 0700);
  }

//  // Create output files
//  char path2file[FILE_NAME_SIZE] = "";
//
//  // number density
//  sprintf(path2file, "%s/%s/%s/number-density", ROOT_DIR, DATA_DIR, 
//    DATA_SUBDIR);
//  FILE *file = fopen(path2file, "w");
//  if (file == NULL) {
//    printf("Could not open file %s\n", path2file);
//  }
//
//  fclose(file);
}

// read number of particles from cgns file
int cgns_read_nparts(void)
{
  // Open cgns file and get cgns file index number fn
  char buf[FILE_NAME_SIZE];
  sprintf(buf, "%s/%s/%s", ROOT_DIR, OUTPUT_DIR, partFiles[partFileMap[0]]);
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

// initialize part_struct
void parts_init(void)
{
  parts = (part_struct*) malloc(nparts * sizeof(part_struct));
//  up = (double*) malloc(nparts * sizeof(double));
//  vp = (double*) malloc(nparts * sizeof(double));
//  wp = (double*) malloc(nparts * sizeof(double));

  for(int p = 0; p < nparts; p++) {
    parts[p].x = -1;
    parts[p].y = -1;
    parts[p].z = -1;
    parts[p].r = -1;
    //up[p] = 0.;
    //vp[p] = 0.;
    //wp[p] = 0.;
  }

  // Open cgns file and get cgns file index number fn
  char buf[FILE_NAME_SIZE];
  sprintf(buf, "%s/%s/%s", ROOT_DIR, OUTPUT_DIR, partFiles[partFileMap[0]]);
  int fn;
  cg_open(buf, CG_MODE_READ, &fn);
  
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

  cg_field_read(fn,bn,zn,sn, "Radius", RealDouble, &range_min, &range_max, r);

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
  sprintf(fname, "%s/input/flow.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if (infile == NULL) {
    printf("Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }
  
  // read domain
  fret = fscanf(infile, "DOMAIN\n");
  fret = fscanf(infile, "(Xs, Xe, Xn) %lf %lf %d\n", &dom.xs, &dom.xe, &dom.xn);
  fret = fscanf(infile, "(Ys, Ye, Yn) %lf %lf %d\n", &dom.ys, &dom.ye, &dom.yn);
  fret = fscanf(infile, "(Zs, Ze, Zn) %lf %lf %d\n", &dom.zs, &dom.ze, &dom.zn);
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
  dom.Gcc.is = 0;
  dom.Gcc.js = 0;
  dom.Gcc.ks = 0;

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
  if (orderL == -1) {
    orderL = (int) floor(dom.xn / (2. * meanR) );
  }
  printf("Fourier order in X = %d\n", orderL);

  if (orderM == -1) {
    orderM = (int) floor(dom.yn / (2. * meanR) );
  }
  printf("Fourier order in Y = %d\n", orderM);

  if (orderN == -1) {
    orderN = (int) floor(dom.zn / (2. * meanR) );
  }
  printf("Fourier order in Z = %d\n", orderN);

  // Number of points to reconstruct at
  // Just because we can only resolve xn/2 doesnt mean we should underresolve 
  //  the reconstruction...
  // TODO: how does this affect?
  // TODO: make this = grid size, so that dx=dy=dz??
  // TODO: I should really make a domain struct for this...
  // TODO: for now, is just going to be dom.*n
  //nPointsX = 2*dom.yn + 1;
  //nPointsY = 2*dom.xn + 1;
  //nPointsZ = 2*dom.zn + 1;

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
  sprintf(buf, "%s/%s/%s", ROOT_DIR, OUTPUT_DIR, partFiles[partFileMap[tt]]);
  int fn;
  cg_open(buf, CG_MODE_READ, &fn);
  
  // Set base, zone, and solutions index numbers
  int bn = 1;
  int zn = 1;
  //int sn = 1;

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

//  cg_field_read(fn,bn,zn,sn, "VelocityX", RealDouble, &range_min, &range_max, up);
//  cg_field_read(fn,bn,zn,sn, "VelocityY", RealDouble, &range_min, &range_max, vp);
//  cg_field_read(fn,bn,zn,sn, "VelocityZ", RealDouble, &range_min, &range_max, wp);
//

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
}

// Get sigfigs
void get_sigfigs(void) 
{
  // print the last file name (longest time) to temp variable
  char tempLastFile[CHAR_BUF_SIZE] = "";
  sprintf(tempLastFile, partFiles[partFileMap[nFiles - 1]]);

  // We expect part-*.*.cgns, so split along "-"
  char *tempDash;
  tempDash = strtok(tempLastFile, "-");
  tempDash = strtok(NULL, "-");

  // Now expect *.*.cgns, split along "."
  char *tempDot;
  tempDot = strtok(tempDash, ".");
  sigFigPre = strlen(tempDot);
  tempDot = strtok(NULL, ".");
  sigFigPost = strlen(tempDot);
}

// Write reconstructed data
void cgns_write_field(void)
{
  // Set up the cgns filename
  char format[CHAR_BUF_SIZE] = "";
  char fnameall[CHAR_BUF_SIZE] = "";
  char fnameall2[CHAR_BUF_SIZE] = "";
  sprintf(format, "%%0%d.%df", sigFigPre + sigFigPost + 1, sigFigPost);

  sprintf(fnameall2, "%s/%s/%s/volume-fraction-%s.cgns", ROOT_DIR, DATA_DIR, 
    DATA_SUBDIR, format);
  sprintf(fnameall, fnameall2, partFileTime[tt]);

  // Set up grid filename
  char gname[FILE_NAME_SIZE] = "";
  char gnameall[FILE_NAME_SIZE] = "";
  sprintf(gname, "grid.cgns");
  sprintf(gnameall, "%s/output/%s", ROOT_DIR, "grid.cgns");
  printf("Remember to copy %s to %s/%s/%s!\n", gnameall, ROOT_DIR, DATA_DIR, 
    DATA_SUBDIR);

  // Set up cgns file paths
  char snodename[CHAR_BUF_SIZE] = "";
  char snodenameall[CHAR_BUF_SIZE] = "";
  sprintf(snodename, "Solution-");
  sprintf(snodenameall, "/Base/Zone0/Solution-");
  sprintf(snodename, "%s%s", snodename, format);
  sprintf(snodenameall, "%s%s", snodenameall, format);
  sprintf(snodename, snodename, partFileTime[tt]);
  sprintf(snodenameall, snodenameall, partFileTime[tt]);

  // CGNS vairables
  int fn; // solution file name
  int bn; // solution base name
  int zn; // solution zone name
  int sn; // solution sol name
  int fn_frec;
  cg_open(fnameall, CG_MODE_WRITE, &fn);
  cg_base_write(fn, "Base", 3, 3, &bn);
  cgsize_t size[9];
  size[0] = dom.xn+1; // cells -> vertices
  size[1] = dom.yn+1;
  size[2] = dom.zn+1;
  size[3] = dom.xn;
  size[4] = dom.yn;
  size[5] = dom.zn;
  size[6] = 0;
  size[7] = 0;
  size[8] = 0;
  cg_zone_write(fn, bn, "Zone0", size, Structured, &zn);
  cg_goto(fn, bn, "Zone_t", zn, "end");
  
  cg_link_write("GridCoordinates", gname, "Base/Zone0/GridCoordinates");

  cg_sol_write(fn, bn, zn, "Solution", CellCenter, &sn);
  double *vfout = malloc(dom.Gcc.s3 * sizeof(double));

  cg_field_write(fn, bn, zn, sn, RealDouble, "Volume Fraction", n_rec, &fn_frec);

  cg_user_data_write("Etc");
  cg_goto(fn, bn, "Zone_t", zn, "Etc", 0, "end");
  cgsize_t *N = malloc(sizeof(cgsize_t));
  N[0] = 1;
  cg_array_write("Time", RealDouble, 1, N, &partFileTime[tt]);
  free(N);

  cg_close(fn);
  free(vfout);

}

// Free variables
void free_vars(void)
{
  for (int i = 0; i < nFiles; i++) {
    free(partFiles[i]);
  }
  free(partFiles);
  free(partFileMap);
  free(partFileTime);

  free(parts);
//  free(up);
//  free(vp);
//  free(wp);

  free(n_lmn);

  free(n_rec);

  free(ones);
}

