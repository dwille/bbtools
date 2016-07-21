#ifndef _RECONSTRUCT_H
#define _RECONSTRUCT_H

#include "main.h"

/**** variables ****/
// Coefficients
extern complex double *n_lmn;        // particle number density coeffs

// Reconstructions
extern complex double *n_rec;        // particle number density

// Other
extern double *ones;         // array of ones to place hold q for number density

/**** FUNCTIONS ****/
// Allocate variables
void alloc_arrays();

// Calculate even and odd coefficents
void calc_coeffs(complex double *nq_lmn, double *q, part_struct *parts, 
  double kl, double km, double kn, int corder);

// Evaluate series at desired points
void eval_series(complex double *nq_sum, complex double *nq_lmn, 
  double kl, double km, double kn, int corder);

// evaluate volume fraction for current order and add to cesaro sum
//void eval_vfrac(double *vFrac_ces, double nql_even, double nql_odd, 
  //double *evalZ, double ell, double k_ell);

// normalize nq by n
//void normalize(double *cesaro_sum, double *norm);

#endif
