#ifndef _RECONSTRUCT_H
#define _RECONSTRUCT_H

#include "main.h"

/**** variables ****/
extern double *evalZ;        // locations to evaluate series at

// Coefficients
extern double *nl_even;      // particle number density -- even coeffs
extern double *nl_odd;       // particle number density -- odd coeffs
extern double *nul_even;     // particle u-vel
extern double *nul_odd;      // particle u-vel
extern double *nvl_even;     // particle v-vel
extern double *nvl_odd;      // particle v-vel
extern double *nwl_even;     // particle w-vel
extern double *nwl_odd;      // particle w-vel
extern double *nkel_even;     // kinetic energy
extern double *nkel_odd;      // kinetic energy

extern double *nwl_avg_odd;    // phase averaged w-vel
extern double *nwl_avg_even;   // phase averaged w-vel

// Reconstructions
extern double *n_ces;        // particle number density -- cesaro sum
extern double *vfrac_ces;    // volume fraction -- cesaro sum
extern double *vfrac_wp_ces;    // (volume fraction * vertical particle vel) -- cesaro sum
extern double *nu_ces;       // u-vel -- cesaro sum
extern double *nv_ces;       // v-vel -- cesaro sum
extern double *nw_ces;       // w-vel -- cesaro sum
extern double *nke_ces;      // kinetic energy -- cesaro sum

extern double *nw_avg_ces;     // phase averged wvel reconstruct

extern double *ones;         // array of ones to place hold q for number density

/**** FUNCTIONS ****/
// Allocate variables
void alloc_arrays();

// Calculate constant coefficents
void const_coeffs(double *q, double q_const, double *cesaro_sum);

// Calculate even and odd coefficents
void calc_coeffs(double *nql_even, double *nql_odd, double *q, 
  part_struct *parts, double k_ell);

// Evaluate series at desired points
void eval_series(double *cesaro_sum, double nql_even, double nql_odd, 
  double *evalZ, double ell, double k_ell);

// evaluate volume fraction for current order and add to cesaro sum
void eval_phase_avg(double *phase_avg_ces, double nql_even, double nql_odd, 
  double *evalZ, double ell, double k_ell);

// normalize nq by n
void normalize(double *cesaro_sum, double *norm);

#endif
