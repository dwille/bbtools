#include "main.h"

/* Variables */
double *evalR;
double *evalTh;

int nLegendreEven;

double *g_l0;
double *g_ln_even;
double *g_ln_odd;

double *g_ces;
double *g_ces_avg;

/* FUNCTION DEFINITIONS */

// Allocate arrays
void alloc_arrays()
{
  /* Init points where will evaluate the reconstruction */
  evalR = (double*) malloc(nPointsR * sizeof(double));
  evalTh = (double*) malloc(nPointsTh * sizeof(double));

  double dr = (R0 - 2.*meanR) / ((double) (nPointsR - 1.));
  for (int i = 0; i < nPointsR; i++) {
    evalR[i] = ((double) i)*dr + 2.*meanR;
  }

  double dth = 0.5*PI / (nPointsTh - 1.);
  for (int i = 0; i < nPointsTh; i++) {
    evalTh[i] = ((double) i)*dth;
  }

  #ifdef DEBUG
    printf("  evalR[0] = %lf\n", evalR[0]);
    printf("  dr = %lf\n", dr);
    printf("  evalR[%d] = %lf\n", nPointsR - 1, evalR[nPointsR - 1]);
    printf("\n");
    printf("  evalTh[0] = %lf\n", evalTh[0]);
    printf("  dth = %lf\n", dth);
    printf("  evalTh[%d] = %lf\n", nPointsTh - 1, evalTh[nPointsTh - 1]);
  #endif
  
  /* Init coefficent arrays */
  nLegendreEven = floor(0.5 * legendreOrder) + 1;

  // nn = 0, ll = 0:2:legendreOrder
  // size = nLegendreEven
  g_l0 = (double*) malloc(nLegendreEven * sizeof(double));

  // nn = 1:fourierOrder, ll = 0:2:legendreOrder
  // size = fourierOrder*nLegendreEven
  g_ln_even = (double*) malloc(fourierOrder * nLegendreEven * sizeof(double));
  g_ln_odd = (double*) malloc(fourierOrder * nLegendreEven * sizeof(double));

  for (int i = 0; i < fourierOrder * nLegendreEven; i++) {
    g_ln_even[i] = 0.;
    g_ln_odd[i] = 0.;
  }
  for (int i = 0; i < nLegendreEven; i++) {
    g_l0[i] = 0.;
  }

  /* Init cesaro sum result arrays */
  // size is npointsR x npointsTh
  g_ces = (double*) malloc(nPointsR * nPointsTh * sizeof(double));
  g_ces_avg = (double*) malloc(nPointsR * nPointsTh * sizeof(double));
  for (int rr = 0; rr < nPointsR; rr++) {
    for (int th = 0; th < nPointsTh; th++) {
      g_ces[th + rr*nPointsTh] = 0.;
      g_ces_avg[th + rr*nPointsTh] = 0.;
    }
  }
}

// set sums to zero
void reset_sums(void)
{
  for (int i = 0; i < fourierOrder * nLegendreEven; i++) {
    g_ln_even[i] = 0.;
    g_ln_odd[i] = 0.;
  }

  for (int i = 0; i < nLegendreEven; i++) {
    g_l0[i] = 0.;
  }

  for (int rr = 0; rr < nPointsR; rr++) {
    for (int th = 0; th < nPointsTh; th++) {
      g_ces[th + rr*nPointsTh] = 0.;
    }
  }
}

// calculate part-pair geometery
void part_pair_geometry(int alpha, int beta)
{
  double ri[3], rj[3];
  double rx, ry, rz;
  double dx, dx2, flipL, flipR;
  int flipFlag;

  // Pull positions
  ri[0] = parts[alpha].x; 
  ri[1] = parts[alpha].y; 
  ri[2] = parts[alpha].z;
  rj[0] = parts[beta].x; 
  rj[1] = parts[beta].y; 
  rj[2] = parts[beta].z;

  #define flip(s1,s2,l,i) \
    {dx = s1[i] - s2[i]; \
     dx2 = dx*dx; \
     flipL = s1[i] - (s2[i] + l); \
     flipR = s1[i] - (s2[i] - l); \
     flipFlag = (flipL*flipL < dx2) - (flipR*flipR < dx2); \
     s2[i] += l*(flipFlag); }

  // Flip particles
  flip(ri, rj, dom.xl, 0);
  flip(ri, rj, dom.yl, 1);
  flip(ri, rj, dom.zl, 2);

  // separation in x,y,z and r2
  rx = ri[0] - rj[0];
  ry = ri[1] - rj[1];
  rz = ri[2] - rj[2];

  // r_ab
  r_ab = sqrt(rx*rx + ry*ry + rz*rz);

  // mu_ab = cos(th_ab) = z/r
  mu_ab = rz/r_ab;

  #undef flip
}

// Calculate P_l(mu_ab)
void eval_legendre(double mu, int ll)
{
  // from numerical recipes, page 184
  double p1 = 1.;
  double p2 = 0.;
  double p3;
  for (int j = 1; j <= ll; j++) {
    p3 = p2;
    p2 = p1;
    p1 = ((2.*j - 1.)*mu*p2 - (j - 1.)*p3)/((double) j);
  }

  // p1 has answer
  p_ell = p1;
}

// Calculate constant coefficients -- each pair has contribution
void const_coeffs(int ll)
{
  int cc = 0.5*ll;
  g_l0[cc] += p_ell / (r_ab*r_ab);
}

// Calculate even and odd coefficents
void calc_coeffs(int ll, int nn)
{
  double k_enn = 2.*PI*nn/(R0- 2*meanR);
  double kernel = p_ell/(r_ab * r_ab);
  /* stride: 
    -- ll goes up by two, need to divide
    -- nn loop doesn't include 0, need to subtract 1 from fOrder
    -- nn starts from 1, need to get to zero at start
   */
  int stride = (0.5*ll)*(fourierOrder - 1) + nn - 1;

  g_ln_even[stride] +=  2.*kernel*cos(k_enn*(r_ab - 2*meanR));
  g_ln_odd[stride] += 2.*kernel*sin(k_enn*(r_ab - 2*meanR));

}

// evaluate series at constant term, nn = 0
void eval_l0(int ll)
{
  double ell_const = 2.*((double) ll + 1.);

  for (int rr = 0; rr < nPointsR; rr++) {
    for (int th = 0; th < nPointsTh; th++) {
      int cc = th + rr*nPointsTh;
      int cl = 0.5*ll;

      double mu = cos(evalTh[th]);
      eval_legendre(mu, ll);

      g_ces[cc] += ell_const*p_ell*g_l0[cl];
    }
  }
}

// evaluate series for current order and add to cesaro sum
void eval_ln(int ll, int nn)
{
  double weight = (1. - (double) nn / ((double)fourierOrder + 1.));
  double k_enn = 2.*PI*nn/(R0 - 2*meanR);
  double ell_const = 2.*((double) ll + 1.);
  
  for (int rr = 0; rr < nPointsR; rr++) {
    double rEval = evalR[rr];
    
    for (int th = 0; th < nPointsTh; th++) {
      double mu = cos(evalTh[th]);
      eval_legendre(mu, ll);

      int cc = th + rr*nPointsTh;
      int stride = (0.5*ll)*(fourierOrder - 1) + nn - 1;

      // Add even
      g_ces[cc] += ell_const*p_ell*weight*g_ln_even[stride]*cos(k_enn*(rEval - 2*meanR));

      // Add odd
      g_ces[cc] += ell_const*p_ell*weight*g_ln_odd[stride]*sin(k_enn*(rEval - 2*meanR));
    }
  }
}

