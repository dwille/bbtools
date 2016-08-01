#!/usr/bin/env python2
from setup import *
os.system('clear')

print ""
print " ---- Histograms of 3D data ---- "
print ""

# Setup
simdir = simParams(sys)
(root, simdir, datadir, imgdir) = directoryStructureMarcc(simdir)

# Pull info
(binStart_vf, min_vf, dBin_vf, max_vf, binEnd_vf, mean_vf) = \
  np.genfromtxt(datadir + "info", skip_header=1,skip_footer=1,
    usecols=(1,2,3,4,5,6))
(binStart_wp, min_wp, dBin_wp, max_wp, binEnd_wp, mean_wp) = \
  np.genfromtxt(datadir + "info", skip_header=2,usecols=(1,2,3,4,5,6))

# Pull hists and normalize
hist_vfrac = np.genfromtxt(datadir + "volume_fraction_hist")
#hist_vfrac /= np.sum(hist_vfrac)
nBins_vf = np.size(hist_vfrac)

hist_wp = np.genfromtxt(datadir + "part_wp_hist")
#hist_wp /= np.sum(hist_wp)
nBins_wp = np.size(hist_wp)

# wp are rows, vf are columns
hist_vf_wp = np.genfromtxt(datadir + "vf_wp_hist")
#hist_vf_wp /= np.sum(hist_vf_wp)

# Set up bins
centers_vf = findCenters(binStart_vf, binEnd_vf, nBins_vf, dBin_vf)
centers_wp = findCenters(binStart_wp, binEnd_wp, nBins_wp, dBin_wp)

# Volume Fraction and Part Wp #
fig1 = plt.figure()

ax1 = fig1.add_subplot(211)
ax1.bar(centers_vf, hist_vfrac, dBin_vf)
ax1.plot([mean_vf, mean_vf], ax1.get_ylim(), 'k--')

ax1.set_xlabel(r'$\phi$')
ax1.set_ylabel(r'$P(\phi)$')

ax2 = fig1.add_subplot(212)
ax2.bar(centers_wp, hist_wp, dBin_wp)
ax2.plot([mean_wp, mean_wp], ax2.get_ylim(), 'k--')

ax2.set_xlabel(r'$w_p$')
ax2.set_ylabel(r'$P(w_p)$')


imgname = imgdir + "volume-fraction-hist"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

# Bivariate histogram of vf and wp #
# xlim is wp ylim is vf
fig2 = plt.figure()

ax1 = fig2.add_subplot(111)
ax1.imshow(hist_vf_wp, origin="lower", aspect="auto", interpolation="none",
  extent=[binStart_wp, binEnd_wp, binStart_vf, binEnd_vf], 
  vmin=-np.max(hist_vf_wp), cmap="seismic")
ax1.set_xlim([binStart_wp, binEnd_wp])
ax1.set_ylim([binStart_vf, binEnd_vf])
ax1.plot([mean_wp, mean_wp], ax1.get_ylim(), 'k--')
ax1.plot(ax1.get_xlim(), [mean_vf, mean_vf], 'k--')

imgname = imgdir + "bivariate-hist"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

print "      Done!"
