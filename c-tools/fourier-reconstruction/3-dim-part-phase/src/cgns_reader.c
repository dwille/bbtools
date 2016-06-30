#include "main.h"
#include "cgns_reader.h"

int nFiles;
char **flowFiles;
double *flowFileTime;
int *flowFileMap;
int sigFigPre;
int sigFigPost;
dom_struct dom;
BC bc;

int *phase;
// fftw variables
fftw_complex *chi;    //indicator function
fftw_complex *phi;
fftw_complex *phi_k;

// Read main.config input file
void main_read_input(void)
{
  int fret = 0;
  fret = fret;

  // open config file for reading
  char fname[CHAR_BUF_SIZE] = "";
  sprintf(fname, "%s/%s", ROOT_DIR, CONFIG_FILE);
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
  fret = fscanf(infile, "Fourier Order in (X,Y,Z) (%d,%d,%d)\n", &orderX, 
    &orderY, &orderZ);
  fclose(infile);
  // TODO: add check for order > dom.n/2
  printf("Fourier Order (%d,%d,%d)\n", orderX, orderY, orderZ);
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
  sprintf(buf, "%s/%s", ROOT_DIR, DATA_DIR);
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

  // Set up flow
  phase = (int*) malloc(dom.Gcc.s3 * sizeof(int));
  chi = (fftw_complex*) fftw_malloc(dom.Gcc.s3 * sizeof(fftw_complex));
  phi = (fftw_complex*) fftw_malloc(dom.Gcc.s3 * sizeof(fftw_complex));
  phi_k = (fftw_complex*) fftw_malloc(dom.Gcc.s3 * sizeof(fftw_complex));

  for (int ii = 0; ii < dom.Gcc.s3; ii++) {
    phase[ii] = 0.;
    chi[ii][0] = 0.;
    chi[ii][1] = 0.;
    phi[ii][0] = 0.;
    phi[ii][1] = 0.;
    phi_k[ii][0] = 0.;
    phi_k[ii][1] = 0.;
  }

  #ifdef DEBUG
    show_domain();
  #endif
  fclose(infile);
}

// Read flow_struct data
void cgns_fill_flow(void)
{
  // Open cgns file and get cgns file index number fn
  char buf[FILE_NAME_SIZE];
  sprintf(buf, "%s/%s/%s", SIM_ROOT_DIR, OUTPUT_DIR, flowFiles[flowFileMap[tt]]);
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
  cg_field_read(fn,bn,zn,sn, "Phase", Integer, range_min, range_max, phase);

  // Apply phase mask -- 1 iff particle
  for (int kk = 0; kk < dom.zn; kk++) {
    for (int jj = 0; jj < dom.yn; jj++) {
      for (int ii = 0; ii < dom.xn; ii++) {
        int stride = ii + dom.Gcc.s1*jj + dom.Gcc.s2*kk;
        chi[stride][0] = (double) (phase[stride] > -1);
        chi[stride][1] = 0.;

        //double z = dom.zs + kk*dom.dz;
        //chi[stride][0] = sin(2.*PI*z/dom.zl);
        //chi[stride][0] = (double) (kk > dom.zn*0.5);
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
        chi[cc][0] = sin(1.*xstar) + sin(2.*xstar) /*+ sin(3.*xstar)*/ +
          2.*(sin(1.*ystar) + sin(2.*ystar) /*+ sin(3.*ystar)*/) +
          3.*(sin(1.*zstar) + sin(2.*zstar) /*+ sin(3.*zstar)*/) +
          1.;
      }
    }
  }
  fftw_execute(chi2phi_k);

  double iV = 1./(dom.xl*dom.yl*dom.zl);
  for (int kk = 0; kk < dom.zn; kk++) {
    for (int jj = 0; jj < dom.yn; jj++) {
      for (int ii = 0; ii < dom.xn; ii++) {
        cc = ii + dom.Gcc.s1*jj + dom.Gcc.s2*kk;
        phi_k[cc][0] *= iV;
        phi_k[cc][1] *= iV;
      }
    }
  }
}

void test_out(void)
{
  printf("All Wavenumbers\n");
  for (int nn = 0; nn < dom.Gcc.kn; nn++) {   // z
    for (int mm = 0; mm < dom.Gcc.jn; mm++) { // y
      for (int ll = 0; ll < dom.Gcc.in; ll++) {   // x
        int cc = ll + dom.Gcc.s1*mm + dom.Gcc.s2*nn;

        if ((fabs(phi_k[cc][0]) > 1e-8) || (fabs(phi_k[cc][1]) > 1e-8)) {
          printf("phi_k[%d + %d*%d + %d*%d] = %lf + i%lf\n", 
            ll, dom.Gcc.s1, mm, dom.Gcc.s2, nn, phi_k[cc][0], phi_k[cc][1]);
          fflush(stdout);
        }

      }
    }
  }

  printf("\n\n");

  int kx = 2;
  int ky = 2;
  int kz = 2;
  printf("Choose (kx,ky,kz) = (%d,%d,%d)\n", kx, ky, kz);
  for (int nn = 0; nn <= kz; nn++) {
    for (int mm = 0; mm <= kx; mm++) {
      for (int ll = 0; ll <= ky; ll++) {
        int cc = ll + dom.Gcc.s1*mm + dom.Gcc.s2*nn;
        printf("phi_k[%d + %d*%d + %d*%d] = %lf + i%lf\n", 
          ll, dom.Gcc.s1, mm, dom.Gcc.s2, nn, phi_k[cc][0], phi_k[cc][1]);
        fflush(stdout);
      }
    }
  }
  for (int nn = 0; nn <= kz; nn++) {
    for (int mm = 0; mm <= ky; mm++) {
      for (int ll = 0; ll <= kx; ll++) {
      int xFlip = (dom.xn - ll)*(ll > 0);
      int yFlip = (dom.yn - mm)*(mm > 0);
      int zFlip = (dom.zn - nn)*(nn > 0);
        int cc = xFlip
                +yFlip*dom.Gcc.s1
                +zFlip*dom.Gcc.s2;
        printf("phi_k[%d + %d*%d + %d*%d] = %lf + i%lf\n", 
          xFlip, dom.Gcc.s1, 
          yFlip, dom.Gcc.s2, 
          zFlip, phi_k[cc][0], phi_k[cc][1]);
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
  printf("  orderX %d\n", orderX);
  printf("  orderY %d\n", orderY);
  printf("  orderZ %d\n", orderZ);
}

// Get sigfigs
void get_sigfigs(void) 
{
  // print the last file name (longest time) to temp variable
  char tempLastFile[CHAR_BUF_SIZE] = "";
  sprintf(tempLastFile, flowFiles[flowFileMap[nFiles - 1]]);

  // We expect flow-*.*.cgns, so split along "-"
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

  sprintf(fnameall2, "%s/%s/f-rec-part-phase-3D-%s.cgns", ROOT_DIR, DATA_DIR, 
    format);
  sprintf(fnameall, fnameall2, flowFileTime[tt]);

  // Set up grid filename
  char gname[FILE_NAME_SIZE] = "";
  char gnameall[FILE_NAME_SIZE] = "";
  sprintf(gname, "grid.cgns");
  sprintf(gnameall, "%s/output/%s", SIM_ROOT_DIR, "grid.cgns");

  // Set up cgns file paths
  char snodename[CHAR_BUF_SIZE] = "";
  char snodenameall[CHAR_BUF_SIZE] = "";
  sprintf(snodename, "Solution-");
  sprintf(snodenameall, "/Base/Zone0/Solution-");
  sprintf(snodename, "%s%s", snodename, format);
  sprintf(snodenameall, "%s%s", snodenameall, format);
  sprintf(snodename, snodename, flowFileTime[tt]);
  sprintf(snodenameall, snodenameall, flowFileTime[tt]);

  // CGNS vairables
  int fn; // solution file name
  int bn; // solution base name
  int zn; // solution zone name
  int sn; // solution sol name
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

  // Write real data
  int fn_real;
  double *real_out = malloc(dom.Gcc.s3 * sizeof(double));
  for (int i = 0; i < dom.Gcc.s3; i++) {
    real_out[i] = phi[i][0];
  }
  cg_field_write(fn, bn, zn, sn, RealDouble, "Volume Fraction Real", real_out, &fn_real);
  free(real_out);

  // Write imaginary
  #ifdef DEBUG
    int fn_imag;  //
    double *imag_out = malloc(dom.Gcc.s3 * sizeof(double));
    for (int i = 0; i < dom.Gcc.s3; i++) {
      imag_out[i] = phi[i][1];
    }
    cg_field_write(fn, bn, zn, sn, RealDouble, "Volume Fraction Imag", imag_out, &fn_imag);
    free(imag_out);
  #endif

  cg_user_data_write("Etc");
  cg_goto(fn, bn, "Zone_t", zn, "Etc", 0, "end");
  cgsize_t *N = malloc(sizeof(cgsize_t));
  N[0] = 1;
  cg_array_write("Time", RealDouble, 1, N, &flowFileTime[tt]);
  free(N);

  cg_close(fn);

}

// Free vars
void free_vars(void)
{
  for (int i = 0; i < nFiles; i++) {
    free(flowFiles[i]);
  }
  free(flowFiles);
  free(flowFileMap);
  free(flowFileTime);

  free(phase);

  fftw_free(chi);
  fftw_free(phi);
  fftw_free(phi_k);

  #ifdef BATCH
    free(ROOT_DIR);
    free(SIM_ROOT_DIR);
  #endif
}

