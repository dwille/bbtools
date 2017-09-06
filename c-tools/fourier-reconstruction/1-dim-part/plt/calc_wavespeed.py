#!/usr/bin/env python3

# README
#
# Purpose:
#   Determine the wavespeed of the kinematic/continuity waves from the 1-D
#   Fourier reconstruction. The wavespeed is calculated from cross-correlations
#   of phi at differing locations, e.g.
#     (phi(z_i, t), phi(z_i + dz, t + dt)),
#   where dz ranges over the entire domain. This can also be written as
#     (phi(z_i, t), phi(z_j, t + dt))
#   Then, this is averaged over all z_i to find phi(dz, dt). This is an averaged
#   view of the  kinematic waves, and the slope of the maximum of the waves is
#   the wavespeed.
#
#   This was previously separated into two other scripts,
#     - avg-crosscorr-spacetime.py, which calculates the averaged phi(dz, dt)
#     - avg-plot-spacetime-vfrac.py, which finds the slope based on the maxima
#   However, these only returned the average slope. In order to find the error
#   (e.g., the standard deviations of the velocities before averaging over all
#   z_i), the slope needs to be found within the outside for-loop over z_i.
#
# Usage:
#   ./calc_wavespeed.py triply_per/<nparts>/<density> <sim steady state time>
#
# 2017 Daniel Willen daniel.willen@jhu.edu

from setup import *
from scipy import signal
from matplotlib.ticker import MultipleLocator

# Initialize simulation parameters
(radius, nparts, rho,  _, simdir, ts) = simParams(sys)
(root, simdir, datadir, imgdir) = directoryStructure(simdir)
(time, ind_ts, tn, z, zn) = initData(datadir, ts)
printSimulationData(radius, root, simdir, datadir)
dz = z - z[0]

# Load volume fraction data after steady state time
vfrac_file = datadir + "volume-fraction"
vfrac = np.genfromtxt(vfrac_file).T[:, ind_ts:]

# Initialize data arrays
vfrac_xcorr_mean = np.zeros((zn, tn))
slope = np.zeros(zn)
tmp = np.zeros((zn, tn))

## Main cross correlation loop ##
for zi in np.arange(0, zn):       # Loop over starting locations
  for dz_i in np.arange(0, zn):   # Loop over all dz
    # Correctly loop through domain (correct for periodicity)
    if zi + dz_i >= zn:
      zj = zi + dz_i - zn
    else:
      zj = zi + dz_i

    # Perform cross-correlation (time-reverse z_j for xcorr)
    y1 = vfrac[zi, :] - vfrac[zi, :].mean()
    y2 = vfrac[zj, :] - vfrac[zj, :].mean()

    result = signal.fftconvolve(y2[::-1], y1, mode="full")

    # Flip result (strictly doesn't matter due to wrap around ordering)
    result = result[::-1]

    # Pull positive portion
    tmp[dz_i, :] = result[(len(result) // 2):]

  # end for dz_i...

  # Normalize xcorr and add to mean
  tmp /= tmp[0,0]
  vfrac_xcorr_mean += tmp/zn

  # Find slope for current phi(z_i, tau) (e.g. tmp)
  max_1 = np.zeros((zn, 3)) # [time, z, xcorr]
  max_2 = np.zeros((zn, 3))
  max_3 = np.zeros((zn, 3))
  for zz in np.arange(0, zn):
    # Find maxima locations by searching for 2nd derivative changes
    maxima_loc = (np.diff(np.sign(np.diff(tmp[zz,:]))) < 0).nonzero()[0] + 1
    maxima =  tmp[zz, maxima_loc]

    if np.size(maxima_loc) == 0:
      max_1[zz,0] = np.nan
      max_1[zz,1] = dz[zz]
      max_1[zz,2] = np.nan

      max_2[zz,0] = np.nan
      max_2[zz,1] = dz[zz]
      max_2[zz,2] = np.nan
    elif np.size(maxima_loc) > 0:         # Find first maximum
      max_1[zz,0] = time[maxima_loc[0]]
      max_1[zz,1] = dz[zz]
      max_1[zz,2] = tmp[zz,maxima_loc[0]]

      if np.size(maxima_loc) > 1:         # Second maximum
        max_2[zz,0] = time[maxima_loc[1]]
        max_2[zz,1] = dz[zz]
        max_2[zz,2] = tmp[zz,maxima_loc[1]]

        if np.size(maxima_loc) > 2:       # Third maximum
          max_3[zz,0] = time[maxima_loc[2]]
          max_3[zz,1] = dz[zz]
          max_3[zz,2] = tmp[zz,maxima_loc[2]]

    if tmp[zz,0] > tmp[zz,1]:     # Deal with maxima at t=t0
      if np.size(maxima_loc) > 2:
        max_3[zz,0] = max_2[zz,0]
        max_3[zz,1] = max_2[zz,1]
        max_3[zz,2] = max_2[zz,2]

      if np.size(maxima_loc) > 1:
        max_2[zz,0] = max_1[zz,0]
        max_2[zz,1] = max_1[zz,1]
        max_2[zz,2] = max_1[zz,2]

      max_1[zz,0] = time[0]
      max_1[zz,1] = dz[zz]
      max_1[zz,2] = tmp[zz,0]

  # end for zz..

  # Find a continuous maxima to fit a curve to
  # Assume that the firt maxima is continuous, then loop through
  # and double check. If two adjacent maxima are too far apart,
  # try to use another maxima
  tau = -np.ones(zn)
  zeta = -np.ones(zn)
  tau[0] = max_1[0,0]
  zeta[0] = max_1[0,1]

  tstar = 1 ## XXX if non-dimensionalized, change
  tauDistance = 0.02/tstar  # to make sure maximas are near each other
  for zz in np.arange(0, zn-1):
    tau_next_1 = max_1[zz+1, 0]
    tau_next_2 = max_2[zz+1, 0]
    tau_next_3 = max_3[zz+1, 0]

    # Use first max if it is close
    if np.abs(tau[zz] - tau_next_1) < tauDistance:
      tau[zz+1] = tau_next_1
      zeta[zz+1] = max_1[zz+1,1]
    # if it's too far, try the second max
    elif np.abs(tau[zz] - tau_next_2) < tauDistance:
      tau[zz+1] = tau_next_2
      zeta[zz+1] = max_2[zz+1,1]
    # if it's too far, try the third
    elif np.abs(tau[zz] - tau_next_3) < tauDistance:
      tau[zz+1] = tau_next_3
      zeta[zz+1] = max_3[zz+1,1]
    # if that is not good, quit
    else:
      break;

  tau = tau[tau != -1]
  zeta = zeta[zeta != -1]

  # Fit curve, assume 0 intercept -- p is slope
  p, _, _, _ = np.linalg.lstsq(tau[:,np.newaxis], zeta)
  yFit = p*tau
  slope[zi] += p

print("%d/%.1f: slope mean = %.3lf" % (nparts, rho, np.mean(slope)))
print("%d/%.1f: slope sdev = %.3lf" % (nparts, rho, np.std(slope)))

# XXX save data

# XXX TEMPORARY XXX
nu = 0.0175 # mm^2/ms
a = 2.1     # mm
#tau        # ms
#dz         # mm

dz /= (2. * a)              # mm / mm
time *= 1000. * nu / (2. * a)**2.   # 1000 * s * (mm^2/ms) / mm^2

zeta /= (2. * a)            # mm / mm
yFit /= (2. * a)            # mm / mm
tau *= 1000. * nu / (2. * a)**2.    # 1000 * s * (mm^2/ms) / mm^2

print(np.max(vfrac_xcorr_mean))
print(np.min(vfrac_xcorr_mean))

fig = plt.figure(figsize=(3.25, 1.625))
plt.imshow(vfrac_xcorr_mean, origin='lower', aspect='auto', interpolation='none',
  extent=[time[0], time[-1], dz[0], dz[-1]],
  vmin=-0.5, vmax=0.5, cmap='coolwarm')
plt.colorbar()
plt.plot(tau, zeta, 'k-')
plt.plot(tau, yFit, 'w--')
plt.xlabel(r'$\nu \Delta t / (2a)^2$')
plt.xlim([0, 6])
plt.ylabel(r'$\Delta z / 2a$')
plt.ylim([dz[0], dz[-1]])
plt.gca().yaxis.set_label_coords(-0.15, 0.5)
plt.gca().xaxis.set_major_locator(MultipleLocator(2))
plt.gca().xaxis.set_minor_locator(MultipleLocator(1))
plt.gca().yaxis.set_major_locator(MultipleLocator(5))
plt.gca().yaxis.set_minor_locator(MultipleLocator(2.5))
plt.savefig(imgdir + "wavespeed-slope-fit.png", bbox_inches='tight', format='png')
