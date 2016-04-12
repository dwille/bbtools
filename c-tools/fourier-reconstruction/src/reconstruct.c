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

double *n_ces;
double *vFrac_ces;
double *nu_ces;       // u-vel -- cesaro sum
double *nv_ces;       // v-vel -- cesaro sum
double *nw_ces;       // w-vel -- cesaro sum

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
  nl_even = (double*) malloc((order + 1) * sizeof(double));
  nl_odd = (double*) malloc((order + 1) * sizeof(double));
  nul_even = (double*) malloc((order + 1) * sizeof(double));
  nul_odd = (double*) malloc((order + 1) * sizeof(double));
  nvl_even = (double*) malloc((order + 1) * sizeof(double));
  nvl_odd = (double*) malloc((order + 1) * sizeof(double));
  nwl_even = (double*) malloc((order + 1) * sizeof(double));
  nwl_odd = (double*) malloc((order + 1) * sizeof(double));

  for (int i = 0; i < order; i++) {
    nl_even[i] = 0.;
    nl_odd[i] = 0.;
    nul_even[i] = 0.;
    nul_odd[i] = 0.;
    nvl_even[i] = 0.;
    nvl_odd[i] = 0.;
    nwl_even[i] = 0.;
    nwl_odd[i] = 0.;
  }

  /* Init cesaro sum result arrays */
  // size is npoints x nFiles
  n_ces = (double*) malloc(npoints * nFiles * sizeof(double));
  vFrac_ces = (double*) malloc(npoints * nFiles * sizeof(double));
  nu_ces = (double*) malloc(npoints * nFiles * sizeof(double));
  nv_ces = (double*) malloc(npoints * nFiles * sizeof(double));
  nw_ces = (double*) malloc(npoints * nFiles * sizeof(double));

  // initialize volume fraction
  double vFrac0 = 4./3.*PI*meanR*meanR*meanR * nparts / (dom.xl*dom.yl*dom.zl);
  for (int j = 0; j < nFiles; j++) {
    for (int i = 0; i < npoints; i++) {
      int c = i + j*npoints;
      vFrac_ces[c] = vFrac0;
    }
  }

  /* q-array for number density is just 1 */
  ones = (double*) malloc(nparts * sizeof(double));

  for (int i = 0; i < nparts; i++) {
    ones[i] = 1.;
  }
}

// Calculate constant coefficients
void const_coeffs(double *q, double *cesaro_sum)
{
  double sum = 0.;
  for (int nn = 0; nn < nparts; nn++) {
    sum += q[nn];
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
  double weight = (1. - ell / (order + 1.));
  
  for (int zz = 0; zz < npoints; zz++) {
    int cc = zz + npoints*tt;

    // Add even
    cesaro_sum[cc] += weight*nql_even*cos(k_ell * evalZ[zz]);

    // Add odd
    cesaro_sum[cc] += weight*nql_odd*sin(k_ell * evalZ[zz]);
  }
}

// evaluate volume fraction for current order and add to cesaro sum
void eval_vfrac(double *vFrac_ces, double nql_even, double nql_odd, 
  double *evalZ, double ell, double k_ell)
{
  double weight = (1. - ell / (order + 1.));
  double correction = 4.*PI/(k_ell*k_ell*k_ell)*(sin(k_ell*meanR) 
                                                - k_ell*meanR*cos(k_ell*meanR));
  for (int zz = 0; zz < npoints; zz++) {
    int cc = zz + npoints*tt;

    // Add even
    vFrac_ces[cc] += weight*correction*nql_even*cos(k_ell * evalZ[zz]);

    // Add odd
    vFrac_ces[cc] += weight*correction*nql_odd*sin(k_ell * evalZ[zz]);
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
