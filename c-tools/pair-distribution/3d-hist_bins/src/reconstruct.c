#include "main.h"

/* Variables */

/* FUNCTION DEFINITIONS */

// Find mean, sedev
void vfrac_stats(void)
{
  tt = 0;
  cgns_fill();
  double id3 = 1./dom.Gcc.s3;
  double min_vf = DBL_MAX;
  double max_vf = -DBL_MAX;

  for (int i = 0; i < dom.Gcc.s3; i++) {
    mean_vf += volume_fraction[i];

    min_vf = min_vf*(min_vf < volume_fraction[i]) +
             volume_fraction[i]*(volume_fraction[i] < min_vf);
    max_vf = max_vf*(max_vf > volume_fraction[i]) +
             volume_fraction[i]*(volume_fraction[i] > max_vf);
  }

  mean_vf *= id3;

  for (int i = 0; i < dom.Gcc.s3; i++) {
    sdev_vf += (volume_fraction[i] - mean_vf)*
               (volume_fraction[i] - mean_vf);
  }
  sdev_vf *= id3;
  sdev_vf = sqrt(sdev_vf);

  #ifdef DEBUG
    printf("  Volume Fraction Mean: %lf\n", mean_vf);
    printf("  Volume Fraction SDEV: %lf\n", sdev_vf);
    printf("  Volume Fraction Max: %lf\n", max_vf);
    printf("  Volume Fraction Min: %lf\n", min_vf);
  #endif
}


// calculate part-pair geometery
void part_pair_geometry(int alpha, int beta)
{
  double ra[3], rb[3];
  double rx, ry, rz;
  double dx, dx2, flipL, flipR;
  int flipFlag;

  // Pull positions
  ra[0] = parts[alpha].x; 
  ra[1] = parts[alpha].y; 
  ra[2] = parts[alpha].z;
  rb[0] = parts[beta].x; 
  rb[1] = parts[beta].y; 
  rb[2] = parts[beta].z;

  /* Check periodicity:
    -- Compute distance between particles
    -- Compute distance between particles if:
      -- second particle is flipped to the left side of the domain (+domlength)
      -- second particle is flipped to the right side of the domain (+domlength)
   */
  #define flip(s1,s2,l,i) \
    {dx = s1[i] - s2[i]; \
     dx2 = dx*dx; \
     flipL = s1[i] - (s2[i] + l); \
     flipR = s1[i] - (s2[i] - l); \
     flipFlag = (flipL*flipL < dx2) - (flipR*flipR < dx2); \
     s2[i] += l*((double) flipFlag); }

  // Flip particles
  flip(ra, rb, dom.xl, 0);
  flip(ra, rb, dom.yl, 1);
  flip(ra, rb, dom.zl, 2);

  #undef flip

  // separation in x,y,z and r2
  rx = ra[0] - rb[0];
  ry = ra[1] - rb[1];
  rz = fabs(ra[2] - rb[2]);
  
  // r_ab
  r_ab = sqrt(rx*rx + ry*ry + rz*rz);
  // th_ab -- between 0, pi/2
  th_ab = acos(rz / r_ab);

  // Volume fraction check
  int i1,j1,k1;
  int i2,j2,k2;
  int check1 = 0;
  int check2 = 0;
  int cc1, cc2;

  // Determine grid location of particle centers
  i1 = floor((parts[alpha].x - dom.xs)/dom.dx);
  j1 = floor((parts[alpha].y - dom.ys)/dom.dy);
  k1 = floor((parts[alpha].z - dom.zs)/dom.dz);
  cc1 = k1 + dom.Gcc.s1*j1 + dom.Gcc.s2*i1;

  i2 = floor((parts[beta].x - dom.xs)/dom.dx);
  j2 = floor((parts[beta].y - dom.ys)/dom.dy);
  k2 = floor((parts[beta].z - dom.zs)/dom.dz);
  cc2 = k2 + dom.Gcc.s1*j2 + dom.Gcc.s2*i2;

  // Check to see if vfrac is greater than or less than the desired
  // cutoff
  if (gtr == 1) {
    // If volume fraction is greater than number
    check1 = (volume_fraction[cc1] >= (mean_vf + cutoff*sdev_vf));
    check2 = (volume_fraction[cc2] >= (mean_vf + cutoff*sdev_vf));
  } else if (less == 1) {
    // If volume fraction is less than number
    check1 = (volume_fraction[cc1] <= (mean_vf + cutoff*sdev_vf));
    check2 = (volume_fraction[cc2] <= (mean_vf + cutoff*sdev_vf));
  }

  /* if both r_a and r_b are not in the desired volume fraction
   *  regime, make r_ab huge so the  if statement following this function
   *  doesn't trigger
   */
  if ((check1 == 1) && (check2 == 1)) {
    // Do nothing
  } else {
    r_ab = DBL_MAX;
  }


}

