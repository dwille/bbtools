#include "template.h"

// Define global variables declared in header file
double ts;
double te;
double *fileTime;

int nparts;

int main(void) {

  // Check if running for all times
  printf("Run for all time? (y/N)\n\t");
  fflush(stdout);
  int tFlag = getchar();
  getchar();
  if (tFlag != 'y' && tFlag != 'Y') {
    char buf[CHAR_BUF_SIZE] = "";
    // Ask for start and ending times
    printf("Starting time?\n\t");
    fflush(stdout);
    fgets(buf, CHAR_BUF_SIZE, stdin);
    sscanf(buf, "%lf\n", &ts);

    printf("Ending time?\n\t");
    fflush(stdout);
    fgets(buf, CHAR_BUF_SIZE, stdin);
    sscanf(buf, "%lf\n", &te);
  } else if (tFlag == 'y' || tFlag == 'Y') {
    ts = 0;
    te = INT_MAX;
  } else {
    printf("Something went wrong.\n");
    return EXIT_FAILURE;
  }

  // TODO: sort *fileTime in read_*fileTime
  // Read part time from output directory
  int partCount = 0;
  find_partCount(&partCount);

  double *partFileTime; 
  partFileTime = malloc(partCount * sizeof(double));
  read_partFileTime(partFileTime);

  // Read flow time from output directory
  int flowCount = 0;
  find_flowCount(&flowCount);

  double *flowFileTime;
  flowFileTime = malloc(flowCount * sizeof(double));
  read_flowFileTime(flowFileTime);

  /* read cgns files
   * wf, phase, wp
   * phaseMask (phase == -1)
   * nf = sum(phaseMask)
   * flowField (Wf.*phaseMask)
   * wf_bar = sum(sum(sum(flowField)))./nf
   * wp_rel = wp - wf_bar
   * wset = mean(wp_rel)
   */

  return EXIT_SUCCESS;
}


