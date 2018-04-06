#include "main.h"

/* Variables */
double *evalZ;

double *nl_even;
double *nl_odd;
double *nul_even;     // particle u-vel
double *nul_odd;      // particle u-vel
double *nvl_even;     // particle v-vel
double *nvl_odd;      // particle v-vel
double *nwl_even;     // particle w-vel
double *nwl_odd;      // particle w-vel
double *nkel_even;     // kinetic energy
double *nkel_odd;      // kinetic energy

double *nwl_avg_odd;    // phase averaged w-vel
double *nwl_avg_even;   // phase averaged w-vel

double *n_ces;
double *vfrac_ces;
double *vfrac_wp_ces;    // reconstruction of (phi * up)
double *nu_ces;       // u-vel -- cesaro sum
double *nv_ces;       // v-vel -- cesaro sum
double *nw_ces;       // w-vel -- cesaro sum
double *nke_ces;      // kinetic energy -- cesaro sum

double *nw_avg_ces;     // phase averged wvel reconstruct

double *ones;

/* FUNCTION DEFINITIONS */

// Allocate arrays
void alloc_arrays()
{
  /* Init points where will evaluate the reconstruction */
  evalZ = (double*) malloc(npoints * sizeof(double));

  double dz = dom.zl / ((double) (npoints - 1.));
  for (int i = 0; i < npoints; i++) {
    evalZ[i] = dom.zs + ((double) i)*dz;
  }

    #ifdef DEBUG
      printf("  evalZ[0] = %lf\n", evalZ[0]);
      printf("  dz = %lf\n", dz);
      printf("  evalZ[%d] = %lf\n", npoints - 1, evalZ[npoints - 1]);
    #endif
  
  /* Init coefficent arrays */
  // arrays are size (order+1) because includes 0 and order
  //  but zero index will be trash -- this is given by constant coefficient
  // TODO: is allocating much more memory than needed if we're only looking at
  // one wavenumber...
  nl_even = (double*) malloc((order_e + 1)* sizeof(double));
  nl_odd = (double*) malloc((order_e + 1)* sizeof(double));
  nul_even = (double*) malloc((order_e + 1)* sizeof(double));
  nul_odd = (double*) malloc((order_e + 1)* sizeof(double));
  nvl_even = (double*) malloc((order_e + 1)* sizeof(double));
  nvl_odd = (double*) malloc((order_e + 1)* sizeof(double));
  nwl_even = (double*) malloc((order_e + 1)* sizeof(double));
  nwl_odd = (double*) malloc((order_e + 1)* sizeof(double));
  nkel_even = (double*) malloc((order_e + 1)* sizeof(double));
  nkel_odd = (double*) malloc((order_e + 1)* sizeof(double));

  nwl_avg_even = (double*) malloc((order_e + 1)* sizeof(double));
  nwl_avg_odd = (double*) malloc((order_e + 1)* sizeof(double));

  for (int i = 0; i <(order_e + 1); i++) {
    nl_even[i] = 0.;
    nl_odd[i] = 0.;
    nul_even[i] = 0.;
    nul_odd[i] = 0.;
    nvl_even[i] = 0.;
    nvl_odd[i] = 0.;
    nwl_even[i] = 0.;
    nwl_odd[i] = 0.;
    nkel_even[i] = 0.;
    nkel_odd[i] = 0.;

    nwl_avg_even[i] = 0.;
    nwl_avg_odd[i] = 0.;
  }

  /* Init cesaro sum result arrays */
  // size is npoints x nFiles
  n_ces = (double*) malloc(npoints * nFiles * sizeof(double));
  vfrac_ces = (double*) malloc(npoints * nFiles * sizeof(double));
  vfrac_wp_ces = (double*) malloc(npoints * nFiles * sizeof(double));
  nu_ces = (double*) malloc(npoints * nFiles * sizeof(double));
  nv_ces = (double*) malloc(npoints * nFiles * sizeof(double));
  nw_ces = (double*) malloc(npoints * nFiles * sizeof(double));
  nke_ces = (double*) malloc(npoints * nFiles * sizeof(double));
  nw_avg_ces = (double*) malloc(npoints * nFiles * sizeof(double));

  //double vfrac0 = 4./3.*PI*mean_r*mean_r*mean_r * nparts / (dom.xl*dom.yl*dom.zl);
  for (int j = 0; j < nFiles; j++) {
    for (int i = 0; i < npoints; i++) {
      int c = i + j*npoints;

      //if (order_s == 0) {
      //  vfrac_ces[c] = vfrac0;  // set this explicitly
      //} else {
      //  vfrac_ces[c] = 0.;  // set this explicitly
      //}

      vfrac_ces[c] = 0.;
      n_ces[c] = 0.;
      nu_ces[c] = 0.;
      nv_ces[c] = 0.;
      nw_ces[c] = 0.;
      nke_ces[c] = 0.;
      nw_avg_ces[c] = 0.;
    }
  }

  /* q-array for number density is just 1 */
  ones = (double*) malloc(nparts * sizeof(double));

  for (int i = 0; i < nparts; i++) {
    ones[i] = 1.;
  }
}

// Calculate constant coefficients
// q          -- particle quantity to reconstruct
// q_const    -- constant to multiply q by
// cesaro_sum -- sum to add const coeff to
void const_coeffs(double *q, double q_const, double *cesaro_sum)
{
  double sum = 0.;
  for (int nn = 0; nn < nparts; nn++) {
    sum += q[nn]*q_const;
  }
  double nq0 = sum/(dom.xl * dom.yl * dom.zl);

  for (int zz = 0; zz < npoints; zz++) {
    int cc = zz + npoints*tt;
    cesaro_sum[cc] = nq0;
  }
}

// Calculate even and odd coefficents
void calc_coeffs(double *nql_even, double *nql_odd, double *q, 
  part_struct *parts, double k_ell)
{
  double tmp_nql_even = 0.;
  double tmp_nql_odd = 0.;

  // Calculate nql even and odd
  for (int nn = 0; nn < nparts; nn++) {
    tmp_nql_even += q[nn] * cos(k_ell * parts[nn].z); 
    tmp_nql_odd += q[nn] * sin(k_ell * parts[nn].z); 
  }

  // multiply by prefactor constant
  *nql_even = pref*tmp_nql_even;
  *nql_odd = pref*tmp_nql_odd;
}

// evaluate series for current order and add to cesaro sum
void eval_series(double *cesaro_sum, double nql_even, double nql_odd,
  double *evalZ, double ell, double k_ell)
{
  // Cesaro sum weight -- set = 1 for not cesaro summation
  double weight = (1. - ell / (order_e + 1.));
  
  for (int zz = 0; zz < npoints; zz++) {
    int cc = zz + npoints*tt;

    // Add even
    cesaro_sum[cc] += weight*nql_even*cos(k_ell * evalZ[zz]);

    // Add odd
    cesaro_sum[cc] += weight*nql_odd*sin(k_ell * evalZ[zz]);
  }
}

// evaluate volume fraction for current order and add to cesaro sum
void eval_phase_avg(double *phase_avg_ces, double nql_even, double nql_odd, 
  double *evalZ, double ell, double k_ell)
{
  double weight = (1. - ell / (order_e + 1.));
  double correction = 4.*PI/(k_ell*k_ell*k_ell)*(sin(k_ell*mean_r) 
                                                - k_ell*mean_r*cos(k_ell*mean_r));
  for (int zz = 0; zz < npoints; zz++) {
    int cc = zz + npoints*tt;

    // Add even
    phase_avg_ces[cc] += weight*correction*nql_even*cos(k_ell * evalZ[zz]);

    // Add odd
    phase_avg_ces[cc] += weight*correction*nql_odd*sin(k_ell * evalZ[zz]);
  }
}

// normalize nq by n
void normalize(double *cesaro_sum, double *norm)
{
  for (int zz = 0; zz < npoints; zz++) {
    int cc = zz + npoints*tt;

    // Normalize
    cesaro_sum[cc] /= norm[cc];
  }
}
