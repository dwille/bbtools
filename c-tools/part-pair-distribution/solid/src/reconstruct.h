#ifndef _RECONSTRUCT_H
#define _RECONSTRUCT_H

#include "main.h"
#include "cgns_reader.h"

/**** variables ****/
extern double *evalR;        // locations to evaluate series at
extern double *evalTh;       // locations to evaluate series at

// number of even legendre polys
extern int nLegendreEven;   // = (floor(legendreOrder/2) + 1)

// Coefficients
extern double *g_ll;          // coeffs array

// Reconstructions
extern double *g_rec;        // reconstructed pair distr function -- time
extern double *g_rec_avg;    // reconstructed pair distr function -- avg

/**** FUNCTIONS ****/
// Allocate variables
void alloc_arrays();

// reset sums
void reset_sums();

// calculate part-pair geometery
void part_pair_geometry(int alpha, int beta);

// Calculate P_l(mu)
void eval_legendre(double mu, int ll);

// Calculate even and odd coefficents
void calc_coeffs(int ll);

// Evaluate series at desired points
void eval_ll(int ll);

#endif
