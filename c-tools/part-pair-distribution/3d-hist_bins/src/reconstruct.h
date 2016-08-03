#ifndef _RECONSTRUCT_H
#define _RECONSTRUCT_H

#include "main.h"
#include "cgns_reader.h"

/**** variables ****/
/**** FUNCTIONS ****/
// Find mean, sedev
void vfrac_stats(void);

// calculate part-pair geometery
void part_pair_geometry(int alpha, int beta);

#endif
