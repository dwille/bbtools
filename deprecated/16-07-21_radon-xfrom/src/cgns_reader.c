#include "main.h"
#include "cgns_reader.h"

// Read main.config input file
void main_read_input(void)
{
  int fret = 0;
  fret = fret;

  // open config file for reading
  char fname[CHAR_BUF_SIZE] = "";
  sprintf(fname, "%s/radon.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  
  // read input
  // fret = fscanf(infile, "tStart %lf\n", &tStart);
  fclose(infile);
}

// Create direcotry for output data
void create_output(void) {
  // Create data directory if it doesn't exist
  // From stackoverflow-7430248
  struct stat st = {0};
  char buf[CHAR_BUF_SIZE];
  sprintf(buf, "%s/%s", ROOT_DIR, DATA_DIR);
  if (stat(buf, &st) == -1) {
    mkdir(buf, 0700);
  }

  // Create output files
  char path2file[FILE_NAME_SIZE] = "";

  // radon transform
  sprintf(path2file, "%s/%s/radon-xform", ROOT_DIR, DATA_DIR);
  FILE *file = fopen(path2file, "w");
  if (file == NULL) {
    printf("Could not open file %s\n", path2file);
  }
  fclose(file);

  /* Create eval/time file */
  sprintf(path2file, "%s/%s/radon-info", ROOT_DIR, DATA_DIR);
  file = fopen(path2file, "w");
  if (file == NULL) {
    printf("Could not open file %s\n", path2file);
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
}

void domain_init(void)
{
  // TODO: perhaps set up some variables here, like dz, dt, total z, total t

  #ifdef DEBUG
    show_domain();
  #endif
}

// Show domain
void show_domain(void)
{
  // print the info from domain_init, or print min(evalZ), max(evalZ),etc
  //printf("Domain:\n");
  //printf("  X: (%f, %f), dX = %f\n", dom.xs, dom.xe, dom.dx);
  //printf("\n");
}

// Write reconstructed data
{
  char fname[CHAR_BUF_SIZE] = "";

  /* radon xform */
  sprintf(fname, "%s/%s/radon-xform", ROOT_DIR, DATA_DIR);
  FILE *file = fopen(fname, "a");
  if (file == NULL) {
    printf("Error opening file %s!\n", fname);
    exit(EXIT_FAILURE);
  }

  /* // each reconstructed timestep is a row
  for (int t = 0; t < nFiles; t++) {
    for (int zz = 0; zz < npoints; zz++) {
      int cc = zz + t*npoints;
      fprintf(file, "%lf ", n_ces[cc]);
    }
    fprintf(file, "\n");
  }
  */
  fclose(file);
}

// Free parts
void free_vars(void)
{

  
  free(evalZ);

}

