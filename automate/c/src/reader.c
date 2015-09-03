#include "multipart.h"
#include "reader.h"
#include "array.h"

// Read tool input file
void read_input(void)
{
  printf("\nRunning multiparticle statistics analysis tool...\n");
  printf("Reading the tool input file...\n");

  int fret = 0;
  fret = fret; // prevent compiler warning

  // open config file for reading
  char fname[FILE_NAME_SIZE] = "";
  sprintf(fname, "%s/input/input.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if (infile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  //char buf[CHAR_BUF_SIZE] = ""; // character read buffer

  // read input
  fret = fscanf(infile, "Starting Time %lf\n", &ts);
  fret = fscanf(infile, "Ending Time %lf\n", &te);
  fret = fscanf(infile, "Timestep %lf\n", &dt);
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "Tetrads %d\n", &tetrad_flag);
  fret = fscanf(infile, "Triads %d\n", &triad_flag);
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "Tools\n");
  fret = fscanf(infile, "Volume %d\n", &volume_flag);
  fret = fscanf(infile, "Radius of Gyration %d\n", &rog_flag);
  fret = fscanf(infile, "Shape Factors %d\n", &shape_factor_flag);
  fret = fscanf(infile, "Lambda %d\n", &lambda_flag);
  fret = fscanf(infile, "Invariants %d\n", &invari_flag);

  printf("\tStarting time\t%.2lf\n", ts);
  printf("\tEnding time\t%.2lf\n", te);
  printf("\tTimestep\t%.2lf\n", dt);
  printf("\tTetrads\t\t%d\n", tetrad_flag);
  printf("\tTetrads\t\t%d\n", triad_flag);
  printf("\t  Volume\t%d\n", volume_flag);
  printf("\t  Rad of Gyrtn\t%d\n", rog_flag);
  printf("\t  Shape Factors\t%d\n", shape_factor_flag);
  printf("\t  Lambda\t%d\n", lambda_flag);
  printf("\t  Invariants\t%d\n", invari_flag);
}

void read_time(void)
{
  DIR *dir;
  struct dirent *ent;
  char output_path[FILE_NAME_SIZE] = "";
  sprintf(output_path, "%s/output", ROOT_DIR);

  // count number of part files in directory
  int partCount = 0;
  int isPart = -3;
  if ((dir = opendir (output_path)) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      isPart = strncmp(ent->d_name, "part", 4);
      partCount = partCount + !isPart;
    }
    closedir (dir);
  } else {
    perror("");
    exit(EXIT_FAILURE);
  }

  // create vector containing all times
  fileTime = malloc(partCount * sizeof(double));
  int fret = 0;
  fret = fret;
  int ind = 0;
  if ((dir = opendir (output_path)) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      isPart = strncmp(ent->d_name, "part", 4);
      if (!isPart) {
        fret = sscanf(ent->d_name, "part-%lf.cgns", &fileTime[ind]); 
        ind++;
      }
    }
    closedir (dir);
  } else {
    perror("");
    exit(EXIT_FAILURE);
  }

  // sort the array
  qsort(fileTime, partCount, sizeof(double), compare);
  for (int i = 0; i < partCount; i++) {
    printf("%lf\n", fileTime[i]);
  }
}

//void read_nparts(void)
//{
//  // read list of files, parse out times -- save
//}

