#include "main.h"

/* Variables */

/* FUNCTION DEFINITIONS */


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

  // r_plane 
  double r_plane = sqrt(rx*rx + ry*ry);

  // th_ab = asin(r_plane/r_ab)
  // asin[-1:1] -> -pi/2:pi/2
  // r_plane/r_ab >=0, so 0 <= th <= pi/2
  th_ab = asin(r_plane/r_ab);

  #undef flip
}

