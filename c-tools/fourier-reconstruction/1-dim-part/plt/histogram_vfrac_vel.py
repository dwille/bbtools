#!/usr/bin/env python2
from setup import *
import scipy.stats as stats
os.system('clear')

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "              Histogram of VFrac and Vel"
print ""

# Setup simulation parameters
(partR, nparts, rho, vFracMean, simdir, tstart) = simParams(sys)

# Setup directory structures
(root, simdir, datadir, imgdir) = directoryStructure(simdir)

# Get time and z-position data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

# Find and load phase-averaged fluid velocity data
phaseVeldir = root + "simdata/phaseAveragedFluidVel"
phaseVelData = np.genfromtxt(phaseVeldir, skip_header=1)

for nn in np.arange(0,np.size(phaseVelData,0)):
  currN = phaseVelData[nn,0]
  currRho = phaseVelData[nn,1]
  if (nparts == currN) & (rho == currRho):
    phaseVel = phaseVelData[nn,2]

# Load output data -- each column is a different time
vFrac = np.genfromtxt(datadir + "volume-fraction").T[:,tsInd:]
up = np.genfromtxt(datadir + "part-u").T[:,tsInd:]
vp = np.genfromtxt(datadir + "part-v").T[:,tsInd:]
wp = np.genfromtxt(datadir + "part-w").T[:,tsInd:]

# Squash to 1-D array for statistics
vFrac = np.squeeze(np.reshape(vFrac,(-1,1)))
up = np.squeeze(np.reshape(up,(-1,1)))
vp = np.squeeze(np.reshape(vp,(-1,1)))
wp = np.squeeze(np.reshape(wp,(-1,1)))

#print stats.pearsonr(up, wp)[0]
#print stats.pearsonr(vp, wp)[0]
#print stats.pearsonr(up, vp)[0]
#sys.exit()

# Data statistics
vfMin = np.min(vFrac)
vfMax = np.max(vFrac)
upMin = np.min(up)
upMax = np.max(up)
vpMin = np.min(vp)
vpMax = np.max(vp)
wpMin = np.min(wp)
wpMax = np.max(wp)

recMeanVF = np.mean(vFrac)
recMeanWP = np.mean(wp)
recPhaseAvgWP = np.mean(vFrac * wp)

vfSDEV = np.std(vFrac)
wpSDEV = np.std(wp)

# print "\n      Mean (reconstructed) volume fraction = %.4f" % recMeanVF
# print "      Max/min volume fraction flucts = %.4f, %.4f" % (vfMax, vfMin)
# print "      Standard deviation of volume fraction flucts = %.4f\n" % vfSDEV
# 
# print "      Mean (reconstructed) velocity flucts = %.4f" % recMeanWP
# print "      Mean (reconstructed) phase averaged velocity flucts = %.4f" % recPhaseAvgWP
# print "      Max/min velocity flucts = %.4f, %.4f" % (wpMax, wpMin)
# print "      Standard deviation of velocity flucts = %.4f" % wpSDEV

# Histogram Setup
nBins = 30
vfEdges = np.linspace(vfMin, vfMax, nBins)
upEdges = np.linspace(upMin, upMax, nBins)
vpEdges = np.linspace(vpMin, vpMax, nBins)
wpEdges = np.linspace(wpMin, wpMax, nBins)

vfCenters = 0.5*(vfEdges[1:] + vfEdges[:-1])
upCenters = 0.5*(upEdges[1:] + upEdges[:-1])
vpCenters = 0.5*(vpEdges[1:] + vpEdges[:-1])
wpCenters = 0.5*(wpEdges[1:] + wpEdges[:-1])

# 1-D Histograms
Hvf,_ = np.histogram(vFrac, bins=vfEdges)
Hup,_ = np.histogram(up, bins=upEdges)
Hvp,_ = np.histogram(vp, bins=vpEdges)
Hwp,_ = np.histogram(wp, bins=wpEdges)

# 2-D Histogram -- phi increase left to right, {u,v,w}p increase top to bottom
H2_vfu,_,_ = np.histogram2d(vFrac, up, bins=(vfEdges, upEdges))
H2_vfv,_,_ = np.histogram2d(vFrac, vp, bins=(vfEdges, vpEdges))
H2_vfw,_,_ = np.histogram2d(vFrac, wp, bins=(vfEdges, wpEdges))

# Weighted least squares line fit for vfrac, wp
phiEval = np.linspace(vfMin, vfMax, 50)

vfX,wpY = np.meshgrid(vfCenters, wpCenters)   # use centers as x,y coords
vfX = np.squeeze(np.reshape(vfX,(-1,1)))
wpY = np.squeeze(np.reshape(wpY,(-1,1)))
H2temp = np.squeeze(np.reshape(H2_vfw,(-1,1)))    # use hist counts as weight

coeffs = np.polynomial.polynomial.polyfit(vfX,wpY,1,w=H2temp)
c = coeffs[0]
m = coeffs[1]

wpEval = m*phiEval + c

print "\n      wp' = %.3f phi + %.3f" % (m, c)
#print "      r_u,phi = %.3f" % (stats.pearsonr(vFrac, up)[0])
#print "      r_v,phi = %.3f" % (stats.pearsonr(vFrac, vp)[0])
#print "      r_w,phi = %.3f" % (stats.pearsonr(vFrac, wp)[0])

# Normalize to probablilities
inorm = 1./float(np.size(vFrac))
Hvf = Hvf.astype(float)*inorm
Hup = Hup.astype(float)*inorm
Hvp = Hvp.astype(float)*inorm
Hwp = Hwp.astype(float)*inorm
H2_vfu = H2_vfu.astype(float)*inorm
H2_vfv = H2_vfv.astype(float)*inorm
H2_vfw = H2_vfw.astype(float)*inorm

# Plot
histFig = plt.figure(figsize=(5,4))

# volume fraction pdf
ax1 = histFig.add_subplot(223)
ax1.plot(vfCenters, Hvf, 'o--')
ax1.plot([vFracMean, vFracMean], [0, np.max(Hvf)], 'k--')
ax1.text(vFracMean, 0.05, r"$\phi_{sdev}$ = %.4f" % vfSDEV)

ax1.set_xlabel(r'$\phi$')
ax1.set_ylabel(r'$P(\phi)$')
ax1.set_xticks(np.linspace(vfMin,vfMax,4))
ax1.xaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
ax1.set_yticks(np.linspace(0, np.max(Hvf), 4))
ax1.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))

# wp pdf
ax2 = histFig.add_subplot(222)
ax2.plot(Hwp, wpCenters, 'o--')
ax2.plot([0, np.max(Hwp)], [recMeanWP, recMeanWP], 'k--')
ax2.text(0.0, recMeanWP, r"$w_{p,sdev}$ = %.4f" % wpSDEV)

ax2.set_xlabel(r'$P(w_p^\prime)$')#/\langle w_f \rangle)$')
ax2.set_ylabel(r'$w_p^\prime$')#/\langle w_f \rangle$')
ax2.set_yticks(np.linspace(wpMin,wpMax,4))
ax2.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
ax2.set_xticks(np.linspace(0, np.max(Hwp),4))
ax2.xaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))

ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")

# histogram -- volume fraction, wp
ax3 = histFig.add_subplot(221, sharex=ax1, sharey=ax2)
im = ax3.imshow(H2_vfw, origin="lower", aspect="auto", interpolation="none",
  extent=[vfMin, vfMax, wpMin, wpMax], 
  vmin=-np.max(H2_vfw), vmax=np.max(H2_vfw), cmap = "seismic")
ax3.plot([vFracMean, vFracMean], [wpMin, wpMax], 'k--')
ax3.plot([vfMin, vfMax], [recMeanWP, recMeanWP], 'k--')
ax3.plot(vfCenters, wpCenters[np.argmax(H2_vfw, 0)], 'k-', alpha=0.6)
ax3.plot(phiEval, wpEval, 'w--')

ax3.text(vfMin, recMeanWP, "m = %.3f" % m)

plt.setp( ax3.get_xticklabels(), visible=False)
plt.setp( ax3.get_yticklabels(), visible=False)
ax3.set_xlim([vfMin, vfMax])
ax3.set_ylim([wpMin, wpMax])

imgname = imgdir + "hist-vfrac-wvel"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

# Comparision of vfrac histrograms
fig2 = plt.figure(figsize=(3,6))

ax3 = fig2.add_subplot(313)
im = ax3.imshow(H2_vfw, origin="lower", aspect="auto", interpolation="none",
  extent=[vfMin, vfMax, wpMin, wpMax], 
  vmin=-np.max(H2_vfw), vmax=np.max(H2_vfw), cmap = "seismic")
ax3.plot(vfCenters, wpCenters[np.argmax(H2_vfw, 0)], 'k-', alpha=0.6)
ax3.plot(phiEval, wpEval, 'w--')
ax3.plot([vFracMean, vFracMean], [wpMin, wpMax], 'k--')

ax3.set_xlabel(r'$\phi$')
ax3.set_xticks(np.linspace(vfMin,vfMax,4))
ax3.xaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
ax3.set_ylabel(r'$w$')
ax3.set_yticks(np.linspace(wpMin,wpMax,4))
ax3.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
ax3.set_xlim([vfMin, vfMax])
ax3.set_ylim([wpMin, wpMax])

ax3.text(0.3, -0.01, r"$r = %.3lf$" % stats.pearsonr(vFrac, wp)[0])
ax3.text(vfMin, recMeanWP, r"$m = %.3f$" % m)

ax1 = fig2.add_subplot(311)
im = ax1.imshow(H2_vfu, origin="lower", aspect="auto", interpolation="none",
  extent=[vfMin, vfMax, vpMin, vpMax], 
  vmin=-np.max(H2_vfu), vmax=np.max(H2_vfu), cmap = "seismic")

ax1.set_xticks(np.linspace(vfMin,vfMax,4))
ax1.xaxis.set_ticklabels([])
ax1.set_ylabel(r'$u$')
ax1.set_yticks(np.linspace(wpMin,wpMax,4))
ax1.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))

ax1.text(0.3, -0.01, r"$r = %.3lf$" % stats.pearsonr(vFrac, up)[0])

ax2 = fig2.add_subplot(312)
im = ax2.imshow(H2_vfv, origin="lower", aspect="auto", interpolation="none",
  extent=[vfMin, vfMax, upMin, upMax], 
  vmin=-np.max(H2_vfv), vmax=np.max(H2_vfv), cmap = "seismic")

ax2.set_xticks(np.linspace(vfMin,vfMax,4))
ax2.xaxis.set_ticklabels([])
ax2.set_ylabel(r'$v$')
ax2.set_yticks(np.linspace(wpMin,wpMax,4))
ax2.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))

ax2.text(0.3, -0.01, r"$r = %.3lf$" % stats.pearsonr(vFrac, vp)[0])

imgname = imgdir + "hist-vfrac-all_vels"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

print "\n      ...Done!"
