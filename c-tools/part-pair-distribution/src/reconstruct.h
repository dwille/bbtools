#ifndef _RECONSTRUCT_H
#define _RECONSTRUCT_H

#include "main.h"
#include "cgns_reader.h"

/**** variables ****/
extern double *evalR;        // locations to evaluate series at
extern double *evalTh;       // locations to evaluate series at

// Coefficients
extern double *g_l0;          // const coeffs
extern double *g_ln_even;     // even coeffs
extern double *g_ln_odd;      // odd coeffs

// Reconstructions
extern double *g_ces;        // particle number density -- cesaro sum

/**** FUNCTIONS ****/
// Allocate variables
void alloc_arrays();

// calculate part-pair geometery
void part_pair_geometry(int alpha, int beta);

// Calculate P_l(mu)
void eval_legendre(double mu, int ll);

// Calculate constant coefficents
void const_coeffs();

// Calculate even and odd coefficents
void calc_coeffs(int ll, int nn);

// evaluate series at constant term, nn = 0
void eval_l0(int ll);

// Evaluate series at desired points
void eval_ln(int ll, int nn);

#endif
