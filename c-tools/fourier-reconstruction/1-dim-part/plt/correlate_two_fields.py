#!/usr/bin/env python2
from setup import *
import matplotlib.ticker as tick
from matplotlib.ticker import MultipleLocator
os.system('clear')

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "              Correlation of two Fields"
print ""

# Setup simulation parameters
(partR, simdir, tstart) = simParams(sys)
nparts = int(simdir.partition('/')[0])
rho = float(simdir.partition('/')[2][-4:-1])
vFracMean = nparts*(4./3.)*np.pi*(partR**3.)/(42.*42.*126.)
print "      Mean Volume Fraction = %.4f" % vFracMean

# Setup directory structures
(root, simdir, datadir, imgdir) = directoryStructureMarcc(simdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

# Make new img directory for field correlations
imgdir += "field_correlations/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# Pull data to arrays
vFrac = np.genfromtxt(datadir + "volume-fraction").T[:,tsInd:]
up = np.genfromtxt(datadir + "part-u").T[:,tsInd:]
vp = np.genfromtxt(datadir + "part-v").T[:,tsInd:]
uperp = np.sqrt(up**2 + vp**2)
wp = np.genfromtxt(datadir + "part-w").T[:,tsInd:]

# Correlations
vFrac_up = np.zeros((nz, nt))
vFrac_vp = np.zeros((nz, nt))
vFrac_wp = np.zeros((nz, nt))
vFrac_uperp = np.zeros((nz, nt))
wp_up = np.zeros((nz, nt))
wp_vp = np.zeros((nz, nt))

for zz,_ in enumerate(evalZ):
  # Correlate vFrac with up, vp, wp
  vFrac_up[zz,:] = CrossCorrelationFFT(vFrac[zz,:], up[zz,:])
  vFrac_vp[zz,:] = CrossCorrelationFFT(vFrac[zz,:], vp[zz,:])
  vFrac_wp[zz,:] = CrossCorrelationFFT(vFrac[zz,:], wp[zz,:])
  vFrac_uperp[zz,:] = CrossCorrelationFFT(vFrac[zz,:], uperp[zz,:])

  # Correlati wp with vp, up
  wp_up[zz,:] = CrossCorrelationFFT(wp[zz,:], up[zz,:])
  wp_vp[zz,:] = CrossCorrelationFFT(wp[zz,:], vp[zz,:])

  # Normalize
  #vFrac_up[zz,:] /= vFrac_up[zz,0]
  #vFrac_vp[zz,:] /= vFrac_vp[zz,0]
  #vFrac_wp[zz,:] /= vFrac_wp[zz,0]
  #wp_up[zz,:] /= wp_up[zz,0]
  #wp_vp[zz,:] /= wp_vp[zz,0]

# Plot!
fsize = (6,3)
vscale = 0.15
# Plot of volume fraction and up #
fig1 = plt.figure(figsize=fsize)
plt.imshow(vFrac_up, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=-vscale, vmax=vscale, cmap='seismic')

plt.colorbar()

plt.title(r"Correlation of $\phi$ and $u_p$")
plt.xlabel(r"$Time\ [s]$")
plt.xlim([0, time[-1]])
plt.ylabel(r"$z\ [mm]$")
plt.ylim([evalZ[0], evalZ[-1]])

plt.gca().xaxis.set_major_locator(MultipleLocator(2))
plt.gca().xaxis.set_minor_locator(MultipleLocator(1))

imgname = imgdir + "vfrac-up"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

# Plot of volume fraction and vp #
fig2 = plt.figure(figsize=fsize)
plt.imshow(vFrac_vp, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=-vscale, vmax=vscale, cmap='seismic')

plt.colorbar()

plt.title(r"Correlation of $\phi$ and $v_p$")
plt.xlabel(r"$Time\ [s]$")
plt.xlim([0, time[-1]])
plt.ylabel(r"$z\ [mm]$")
plt.ylim([evalZ[0], evalZ[-1]])

plt.gca().xaxis.set_major_locator(MultipleLocator(2))
plt.gca().xaxis.set_minor_locator(MultipleLocator(1))

imgname = imgdir + "vfrac-vp"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

# Plot of volume fraction and wp #
fig3 = plt.figure(figsize=fsize)
plt.imshow(vFrac_wp, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=-1., vmax=1., cmap='seismic')

plt.colorbar()

plt.title(r"Correlation of $\phi$ and $w_p$")
plt.xlabel(r"$Time\ [s]$")
plt.xlim([0, time[-1]])
plt.ylabel(r"$z\ [mm]$")
plt.ylim([evalZ[0], evalZ[-1]])

plt.gca().xaxis.set_major_locator(MultipleLocator(2))
plt.gca().xaxis.set_minor_locator(MultipleLocator(1))

imgname = imgdir + "vfrac-wp"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

# Plot of volume fraction and uperp #
fig3 = plt.figure(figsize=fsize)
plt.imshow(vFrac_uperp, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=-vscale, vmax=vscale, cmap='seismic')

plt.colorbar()

plt.title(r"Correlation of $\phi$ and $u_\perp$")
plt.xlabel(r"$Time\ [s]$")
plt.xlim([0, time[-1]])
plt.ylabel(r"$z\ [mm]$")
plt.ylim([evalZ[0], evalZ[-1]])

plt.gca().xaxis.set_major_locator(MultipleLocator(2))
plt.gca().xaxis.set_minor_locator(MultipleLocator(1))

imgname = imgdir + "vfrac-uperp"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

# Plot of wp and up #
fig4 = plt.figure(figsize=fsize)
plt.imshow(wp_up, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=-vscale, vmax=vscale, cmap='seismic')

plt.colorbar()

plt.title(r"Correlation of $w_p$ and $u_p$")
plt.xlabel(r"$Time\ [s]$")
plt.xlim([0, time[-1]])
plt.ylabel(r"$z\ [mm]$")
plt.ylim([evalZ[0], evalZ[-1]])

plt.gca().xaxis.set_major_locator(MultipleLocator(2))
plt.gca().xaxis.set_minor_locator(MultipleLocator(1))

imgname = imgdir + "wp-up"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

# Plot of wp and vp #
fig5 = plt.figure(figsize=fsize)
plt.imshow(wp_vp, origin="lower", aspect="auto", interpolation="none",
  extent=[time[0], time[-1], evalZ[0], evalZ[-1]],
  vmin=-vscale, vmax=vscale, cmap='seismic')

plt.colorbar()

plt.title(r"Correlation of $w_p$ and $u_p$")
plt.xlabel(r"$Time\ [s]$")
plt.xlim([0, time[-1]])
plt.ylabel(r"$z\ [mm]$")
plt.ylim([evalZ[0], evalZ[-1]])

plt.gca().xaxis.set_major_locator(MultipleLocator(2))
plt.gca().xaxis.set_minor_locator(MultipleLocator(1))

imgname = imgdir + "wp-vp"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

print "\n      ...Done!"
