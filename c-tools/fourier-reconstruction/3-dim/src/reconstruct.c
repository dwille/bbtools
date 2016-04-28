#include "main.h"

/* Variables */
_Complex double *n_lmn;
_Complex double *n_rec;

double *ones;

/* FUNCTION DEFINITIONS */

// Allocate arrays
void alloc_arrays()
{
 
  /* Init coefficent arrays */
  // arrays are size (order+1) because includes 0 and order
  int nOrders = (orderL + 1)*(orderM + 1)*(orderN + 1);
  n_lmn = (_Complex double*) malloc(nOrders * sizeof(_Complex double));

  for (int i = 0; i < nOrders; i++) {
    n_lmn[i] = 0.;
  }

  /* Init resulting arrays */
  n_rec = (_Complex double*) malloc(dom.Gcc.s3 * sizeof(_Complex double));

  for (int i = 0; i < dom.Gcc.s3; i++) {
    n_rec[i] = 0.;
  }

  // initialize volume fraction
//  double vFrac0 = 4./3.*PI*meanR*meanR*meanR * nparts / (dom.xl*dom.yl*dom.zl);
//  for (int j = 0; j < nFiles; j++) {
//    for (int i = 0; i < npoints; i++) {
//      int c = i + j*npoints;
//      vFrac_ces[c] = vFrac0;
//    }
//  }

  /* q-array for number density is just 1 */
  ones = (double*) malloc(nparts * sizeof(double));

  for (int i = 0; i < nparts; i++) {
    ones[i] = 1.;
  }
}

// Calculate even and odd coefficents
void calc_coeffs(_Complex double *nq_lmn, double *q, part_struct *parts, 
  double kl, double km, double kn, int corder)
{
  // Calculate nq_lmn even and odd
  for (int nn = 0; nn < nparts; nn++) {
    // exp(i(k dot x))
    nq_lmn[corder] += q[nn]* 
      cexp(J*(kl*parts[nn].x + km*parts[nn].y + kn*parts[nn].z));
  }
  // Normalize by 1/V 
  nq_lmn[corder] /= (dom.xl * dom.yl * dom.zl);
}

// evaluate series for current order and add to cesaro sum
void eval_series(_Complex double *nq_sum, _Complex double *nq_lmn,
  double kl, double km, double kn, int corder)
{
  for (int zz = 0; zz < dom.zn; zz++) {
    double z = dom.zs + dom.dz*zz;

    for (int yy = 0; yy < dom.yn; yy++) {
      double y = dom.ys + dom.dy*yy;

      for (int xx = 0; xx < dom.xn; xx++) {
        double x = dom.xs + dom.dx*xx;
        int cc = xx + yy*dom.Gcc.s1 + zz*dom.Gcc.s2;

        // Add to sum
        nq_sum[cc] += nq_lmn[corder] * cexp(-J * (kl*x + km*y + kn*z));
      }
    }
  }
}

// evaluate volume fraction for current order and add to cesaro sum
//void eval_vfrac(double *vFrac_ces, double nql_even, double nql_odd, 
//  double *evalZ, double ell, double k_ell)
//{
//  double weight = (1. - ell / (order + 1.));
//  double correction = 4.*PI/(k_ell*k_ell*k_ell)*(sin(k_ell*meanR) 
//                                                - k_ell*meanR*cos(k_ell*meanR));
//  for (int zz = 0; zz < npoints; zz++) {
//    int cc = zz + npoints*tt;
//
//    // Add even
//    vFrac_ces[cc] += weight*correction*nql_even*cos(k_ell * evalZ[zz]);
//
//    // Add odd
//    vFrac_ces[cc] += weight*correction*nql_odd*sin(k_ell * evalZ[zz]);
//  }
//}
//
//// normalize nq by n
//void normalize(double *cesaro_sum, double *norm)
//{
//  for (int zz = 0; zz < npoints; zz++) {
//    int cc = zz + npoints*tt;
//
//    // Normalize
//    cesaro_sum[cc] /= norm[cc];
//  }
//}
