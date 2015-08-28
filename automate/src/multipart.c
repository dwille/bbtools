#include "multipart.h"

// Define global variables declared in header file
int nparts;
double ts;
double te;
double dt;
int tetrad_flag;
int triad_flag;
int volume_flag;
int rog_flag;
int shape_factor_flag;
int lambda_flag;
int invari_flag;

double *fileTime;

int main(void) {

  // Read input file and get number of particles
  read_input();
  read_nparts();

  // Read list of files in ./output/part-*.*.cgns
  read_time();
  // strrep part- and .cgns, save remainder to tName array
  
  // pare down tName array to ts < tName < te
  // pare down tName array to multiples of dt, or thereabouts

  // read first file, find nparts
  // construct part_struct[nparts].vec[length(tName)]

  // read part-*.** to part_struct array

  // loop over all r0 TODO: list r0, tol in input file
  //   loop over all particles and find tetrad apirs
  //    create *T struct with first tetrad pair
  //    double n[0], n[1], n[2], n[3]
  //    *prev = NULL
  //    *next = 
  //    create a *T, loop over *prev, check if unique or not
  //    if unique, change *T.next = *currT
  //    *T = *currT
  //    fill out V, R^2, etc. to *T
  // continue loops


  return EXIT_SUCCESS;
}


