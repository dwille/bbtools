#include "main.h"
#include "reader.h"

// Declare variables
char *SIM_ROOT_DIR;
char *ANALYSIS_DIR;
int n_files;
char **part_files;
double *part_file_time;
int *file_map;

// Initialize data directory
void directory_init(int argc, char *argv[])
{
  SIM_ROOT_DIR = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
  ANALYSIS_DIR = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));

  // arg[0] = program name
  // arg[1] = SIM_ROOT_DIR
  if (argc == 2) {
    sprintf(SIM_ROOT_DIR, "%s", argv[1]);
    sprintf(ANALYSIS_DIR, "%s/analysis/%s", SIM_ROOT_DIR, 
      "single-particle-dispersion");
  } else if (argc != 2) {
    printf("usage: %s SIM_ROOT_DIR\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  printf("\n SIM_ROOT_DIR = %s\n", SIM_ROOT_DIR);
  printf(" ANALYSIS_DIR = %s\n", ANALYSIS_DIR);
  fflush(stdout);
}

// Read input file
void read_input(void)
{
  int fret = 0;
  fret = fret; // to prevent compiler warning

  // open config file for reading
  char fname[CHAR_BUF_SIZE] = "";
  sprintf(fname, "%s/%s", ANALYSIS_DIR, "input.config");
  FILE *infile = fopen(fname, "r");
  if (infile == NULL) {
    printf("Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  } else {
    printf(" CONFIG_FILE  = %s\n", fname);
  }

  // Read input
  fret = fscanf(infile, "t_start %lf\n", &t_start);
  fret = fscanf(infile, "t_end %lf\n", &t_end);

  fclose(infile);
}

// Parse sim output directory for cgns files
void init_cgns_files(void)
{
  DIR *dir;
  struct dirent *ent;
  char output_path[FILE_NAME_SIZE] = "";
  int fret = 0;
  fret = fret;
  double time;

  sprintf(output_path, "%s/%s", SIM_ROOT_DIR, "output");
  printf(" OUTPUT_DIR   = %s\n", output_path);

  int is_part_file;
  int in_time_range;
  n_files = 0;

  // Count number of part.cgns files that are within time range
  if ((dir = opendir (output_path)) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      // check if part file (0 if match)
      is_part_file = (strncmp(ent->d_name, "part", 4) == 0);

      if (is_part_file == 1) {
        // check if in time range
        fret = sscanf(ent->d_name, "part-%lf.cgns", &time);
        in_time_range = ((time >= t_start) & (time <= t_end));
        n_files += is_part_file*in_time_range;
      } else {
        continue;
      }
    }
    closedir (dir);
  } else {
    printf("Output directory does not exist!\n");
    exit(EXIT_FAILURE);
  }

  // Store cgns filenames and times
  part_files = (char**) malloc(n_files * sizeof(char*));
  for (int i = 0; i < n_files; i ++) {
    part_files[i] = (char*) malloc(FILE_NAME_SIZE*sizeof(char));
  }
  part_file_time = (double*) malloc(n_files * sizeof(double));

  int cc = 0;

  if ((dir = opendir (output_path)) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      is_part_file = (strncmp(ent->d_name, "part", 4) == 0);

      if (is_part_file == 1) {
        // check if in time range
        fret = sscanf(ent->d_name, "part-%lf.cgns", &time);
        in_time_range = ((time >= t_start) & (time <= t_end));
        
        if (in_time_range == 1) {
          fret = sscanf(ent->d_name, "%s", part_files[cc]);
          part_file_time[cc] = time;
          cc++;
        }
      } else {
        continue;
      }
    }
    closedir (dir);
  }

  // Sort the resulting array by time using temp array
  file_map = malloc(n_files * sizeof(double));
  for (int i = 0; i < n_files; i++) {
    file_map[i] = i;
  }

  merge_sort(part_file_time, n_files, file_map);
  printf("\n");
  printf(" Found %d files in range [%lf, %lf]\n", n_files, t_start, t_end);
  if (n_files == 0) {
    printf("Quitting...\n");
    exit(EXIT_FAILURE);
  }

  sim_time = (double*) malloc(n_files * sizeof(double));
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

// Create output directory
void create_output(void)
{
  // Create output directory if it doesn't exist
  // From stackoverflow-7430248
  struct stat st = {0};
  char buf[CHAR_BUF_SIZE];
  sprintf(buf, "%s/%s", ANALYSIS_DIR, "data");
  if (stat(buf, &st) == -1) {
    mkdir(buf, 0700);
  }

  // XXX See tetrads/src/reader.c for how to
  // create diretories if using multiple runs

}

// Initilize part_struct and fill for first timestep
void parts_init(void)
{
  // Read number of particles
  char buf[FILE_NAME_SIZE];
  sprintf(buf, "%s/%s/%s", SIM_ROOT_DIR, "output", part_files[file_map[0]]);
  int fn;
  cg_open(buf, CG_MODE_READ, &fn);

  // Set base index nuber and zone index number (only one, so is 1)
  int bn = 1;
  int zn = 1;

  // Read zone to find cgsize_t *size, or nparts
  char zonename[FILE_NAME_SIZE] = "";
  cgsize_t nparts_tmp = 0;
  cg_zone_read(fn, bn, zn, zonename, &nparts_tmp);

  // Cast nparts as int
  nparts = nparts_tmp;

  cg_close(fn);

  // Allocate new part struct
  parts = (part_struct*) malloc(nparts * sizeof(part_struct));

  for(int p = 0; p < nparts; p++) {
    parts[p].flip_count_i = 0;
    parts[p].flip_count_j = 0;
    parts[p].flip_count_k = 0;
  }

  // Open first cgns file
  sprintf(buf, "%s/%s/%s", SIM_ROOT_DIR, "output", part_files[file_map[0]]);
  cg_open(buf, CG_MODE_READ, &fn);
  
  // Set base, zone, and solutions index numbers
  bn = 1;
  zn = 1;
  int sn = 1;

  // Read part coords
  double *x = malloc(nparts * sizeof(double));
  double *y = malloc(nparts * sizeof(double));
  double *z = malloc(nparts * sizeof(double));
  for (int p = 0; p < nparts; p++) {
    x[p] = 0.;
    y[p] = 0.;
    z[p] = 0.;
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

  // Read part radius
  double *r = malloc(nparts * sizeof(double));
  for (int p = 0; p < nparts; p++) {
    r[p] = 0.;
  }

  cg_field_read(fn,bn,zn,sn, "Radius", RealDouble, &range_min, &range_max, r);

  for (int p = 0; p < nparts; p++) {
    parts[p].r = r[p];
  }

  cg_close(fn);
  free(x);
  free(y);
  free(z);
  free(r);
}

// Init domain
void domain_init(void) {
  int fret = 0;
  fret = fret;  // prevent compiler warning

  // open config file for reading
  char fname[FILE_NAME_SIZE] = "";
  sprintf(fname, "%s/%s/flow.config", SIM_ROOT_DIR, "input");
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
  } else if(strcmp(cbuf, "NEUMANN") == 0) {
    printf("Warning: single particle statistics only implemented for triply-periodic domains");
    fret = fscanf(infile, "%lf", &dbuf);
  } else {
    printf("flow.config read error -- read %s\n",cbuf);
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "bc.pE %s", cbuf);
  if(strcmp(cbuf, "PERIODIC") == 0) {
  } else if(strcmp(cbuf, "NEUMANN") == 0) {
    printf("Warning: single particle statistics only implemented for triply-periodic domains");
    fret = fscanf(infile, "%lf", &dbuf);
  } else {
    printf("flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "bc.pS %s", cbuf);
  if(strcmp(cbuf, "PERIODIC") == 0) {
  } else if(strcmp(cbuf, "NEUMANN") == 0) {
    printf("Warning: single particle statistics only implemented for triply-periodic domains");
    fret = fscanf(infile, "%lf", &dbuf);
  } else {
    printf("flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "bc.pN %s", cbuf);
  if(strcmp(cbuf, "PERIODIC") == 0) {
  } else if(strcmp(cbuf, "NEUMANN") == 0) {
    printf("Warning: single particle statistics only implemented for triply-periodic domains");
    fret = fscanf(infile, "%lf", &dbuf);
  } else {
    printf("flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "bc.pB %s", cbuf);
  if(strcmp(cbuf, "PERIODIC") == 0) {
  } else if(strcmp(cbuf, "NEUMANN") == 0) {
    printf("Warning: single particle statistics only implemented for triply-periodic domains");
    fret = fscanf(infile, "%lf", &dbuf);
  } else {
    printf("flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "bc.pT %s", cbuf);
  if(strcmp(cbuf, "PERIODIC") == 0) {
  } else if(strcmp(cbuf, "NEUMANN") == 0) {
    printf("Warning: single particle statistics only implemented for triply-periodic domains");
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

// Alloc host memory for single part stats
void alloc_host_mem(void)
{
  r2_total = malloc(n_files * sizeof(double));
  r2_horiz = malloc(n_files * sizeof(double));
  r2_verti = malloc(n_files * sizeof(double));
}


// Fill particle structure at new timestep
void cgns_fill_part_struct(void)
{
  // Open cgns file and get cgns file index number fn
  char buf[FILE_NAME_SIZE];
  sprintf(buf, "%s/%s/%s", SIM_ROOT_DIR, "output", part_files[file_map[tt]]);
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

  // Read actual time
  cg_goto(fn, bn, "Zone_t", zn, "Etc", 0, "end");
  cg_array_read(1, &sim_time[tt]);

  free(x);
  free(y);
  free(z);
  free(u);
  free(v);
  free(w);
  
  cg_close(fn);
}

// Write results
void write_output(void)
{
  // TODO set up for multple runs (see tetrads)

  // Output file
  char fname[CHAR_BUF_SIZE] = "";
  sprintf(fname, "%s/%s/data.csv", ANALYSIS_DIR, "data");

  // Open file for reading
  FILE *file = fopen(fname, "w");
  if (file == NULL) {
    printf("Error opening file %s!\n", fname);
    exit(EXIT_FAILURE);
  }

  // Write
  fprintf(file, "time,r_total,r_verti,r_horiz\n");
  for (int i = 0; i < n_files; i++) {
    // We take sqrt here of the mean-squared value r2_*
    fprintf(file, "%lf,%lf,%lf,%lf\n", sim_time[i], sqrt(r2_total[i]), 
      sqrt(r2_verti[i]), sqrt(r2_horiz[i]));
  }

  fclose(file);
}

// Free memory
void free_mem(void)
{
  for (int i = 0; i < n_files; i++) {
    free(part_files[i]);
  }
  free(part_files);
  free(file_map);
  free(part_file_time);
  free(sim_time);
  free(parts);

  free(r2_total);
  free(r2_horiz);
  free(r2_verti);
}
