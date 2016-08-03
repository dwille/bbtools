#ifndef _READER_H
#define _READER_H

#include "template.h"
#include <dirent.h>
#include <cgnslib.h>

// Find number of part cgns files in output directory
void find_partCount(int *partCount);

// Read part times from cgns filenames
void read_partFileTime(double *partFileTime);

// Find number of flow cgns files in output directory
void find_flowCount(int *flowCount);

// Read flow times from cgns filenames
void read_flowFileTime(double *flowFileTime);

// Read cgns part file
void read_cgnsPart(double *partFileTime, int partCount);

// Read cgns flow file
void read_cgnsFlow(double *flowFileTime, int flowCount);


#endif
