#include "template.h"
#include "reader.h"
#include "array.h"

// Find number of part cgns files in output directory
void find_partCount(int *partCount)
{
  DIR *dir;
  struct dirent *ent;
  char output_path[FILE_NAME_SIZE] = "";
  int fret = 0;
  double time;

  fret = fret;
  sprintf(output_path, "%s/%s", ROOT_DIR, OUTPUT_DIR);

  int isPart;
  int inTime;

  // count number of part files in directory that fall in time range
  if ((dir = opendir (output_path)) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      // check if part file (0 if match)
      isPart = (strncmp(ent->d_name, "part", 4) == 0);

      if (isPart == 1) {
        // check if in time range
        fret = sscanf(ent->d_name, "part-%lf.cgns", &time);
        inTime = ((time >= ts) & (time <= te));
        *partCount = *partCount + isPart*inTime;
      } else {
        continue;
      }

    }
    closedir (dir);
  } else {
    printf("Output directory does not exist!\n");
    exit(EXIT_FAILURE);
  }
}

// Read part times from cgns filenames
void read_partFileTime(double *partFileTime)
{
  DIR *dir;
  struct dirent *ent;
  char output_path[FILE_NAME_SIZE] = "";
  int fret = 0;
  double time;

  fret = fret;
  sprintf(output_path, "%s/%s", ROOT_DIR, OUTPUT_DIR);

  int isPart;
  int inTime;

  int ind = 0;
  // create vector of times in file name
  if ((dir = opendir (output_path)) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      // check if part file (0 if match)
      isPart = (strncmp(ent->d_name, "part", 4) == 0);

      if (isPart == 1) {
        // check if in time range
        fret = sscanf(ent->d_name, "part-%lf.cgns", &time);
        inTime = ((time >= ts) & (time <= te));

        // if it's in time range
        if (inTime == 1) {
        //TODO: maybe read this as a char instead of float
          fret = sscanf(ent->d_name, "part-%lf.cgns", &partFileTime[ind]);
          ind++;
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

  /* // TODO: sort the array
  qsort(fileTime, partCount, sizeof(double), compare);
  for (int i = 0; i < partCount; i++) {
    printf("%lf\n", fileTime[i]);
  } */
}

// Find number of flow cgns files in output directory
void find_flowCount(int *flowCount)
{
  DIR *dir;
  struct dirent *ent;
  char output_path[FILE_NAME_SIZE] = "";
  int fret = 0;
  double time;

  fret = fret;
  sprintf(output_path, "%s/%s", ROOT_DIR, OUTPUT_DIR);

  int isFlow;
  int inTime;

  // count number of flow files in directory that fall in time range
  if ((dir = opendir (output_path)) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      // check if flow file (0 if match)
      isFlow = (strncmp(ent->d_name, "flow", 4) == 0);

      if (isFlow == 1) {
        // check if in time range
        fret = sscanf(ent->d_name, "flow-%lf.cgns", &time);
        inTime = ((time >= ts) && (time <= te));
        *flowCount = *flowCount + isFlow*inTime;
      } else {
        continue; 
      }

    }
    closedir (dir);
  } else {
    printf("Output directory does not exist!\n");
    exit(EXIT_FAILURE);
  }
}

// Read flow times from cgns filenames
void read_flowFileTime(double *flowFileTime)
{
  DIR *dir;
  struct dirent *ent;
  char output_path[FILE_NAME_SIZE] = "";
  int fret = 0;
  double time;

  fret = fret;
  sprintf(output_path, "%s/%s", ROOT_DIR, OUTPUT_DIR);

  int isFlow;
  int inTime;

  int ind = 0;
  // create vector of times in file name
  if ((dir = opendir (output_path)) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      // check if part file (0 if match)
      isFlow = (strncmp(ent->d_name, "flow", 4) == 0);

      if (isFlow == 1) {
        // check if in time range
        fret = sscanf(ent->d_name, "flow-%lf.cgns", &time);
        inTime = ((time >= ts) & (time <= te));

        // if partfile and in time range
        if (inTime == 1) {
        //TODO: maybe read this as a char instead of float
          fret = sscanf(ent->d_name, "flow-%lf.cgns", &flowFileTime[ind]); 
          ind++;
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

  /*// TODO: sort the array
  qsort(fileTime, partCount, sizeof(double), compare);
  for (int i = 0; i < partCount; i++) {
    printf("%lf\n", fileTime[i]);
  }
  */
}

void read_cgnsPart(double *partFileTime, int partCount)
{
  char fname[FILE_NAME_SIZE] = "";
  int fn;
  int nbases;
  int B = 1;
  char basename[FILE_NAME_SIZE] = "Base";
  int cell_dim = 3;
  int phys_dim = 3;
  int Z = 1;

  // Loop over each file
  for (int ff = 0; ff < partCount; ff++) {
    printf(fname, "%s/output/part-%d.cgns", partFileTime[ff]);
    
    cg_open(fname, CG_MODE_READ, &fn);      // open file for reading
    //cg_nbases(fn, &nbases);               // get number of bases (should be 1)
    //cg_base_read(fn, B, basename, cell_dim, phys_dim);  // rea
  }
}

void read_cgnsFlow(double *flowFileTime, int flowCount)
{
  char fname[FILE_NAME_SIZE] = "";
  int fn;
  int nbases;
  int B = 1;
  char Bn[FILE_NAME_SIZE] = "Base";
  int cell_dim = 3;
  int phys_dim = 3;
  int Z = 1;
  char Zn[FILE_NAME_SIZE] = "Zone0";
  int size;
  int S = 1;
  char Sn[FILE_NAME_SIZE] = "Solution";

  // Loop over each file
  for (int ff = 0; ff < flowCount; ff++) {
    printf(fname, "%s/output/flow-%d.cgns", flowFileTime[ff]);
    
    cg_open(fname, CG_MODE_READ, &fn);      // open file for reading
    //cg_nbases(fn, &nbases);               // get number of bases (should be 1)
    //cg_base_read(fn, B, basename, cell_dim, phys_dim);  // rea
    cg_zone_read(fn, B, Z, &Zn, &size);
    // size[0] = NVertexI; size[1] = NVertexJ; size[2] = NVertexK
    // size[3] = NCellI; size[4] = NCellJ; size[5] = NCellK
    // size[6]=NBoundVertexI; size[7]=NBoundVertexJ; size[8]=NBoundVertexK = 0
    cg_field_read(fn, B, Z, S, "VelocityX", RealDouble, 0, NCellI?, solnarray);
  }
}

