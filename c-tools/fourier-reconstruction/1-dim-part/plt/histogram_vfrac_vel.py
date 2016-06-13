#!/usr/bin/env python2
from setup import *
import matplotlib.ticker as tick
os.system('clear')

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "              Histogram of VFrac and Vel"
print ""

# Setup simulation parameters
(partR, simdir, tstart) = simParams(sys)
nparts = int(simdir.partition('/')[0])
rho = float(simdir.partition('/')[2][-4:-1])
vFracMean = nparts*(4./3.)*np.pi*(partR**3.)/(42.*42.*126.)
print "      Mean Volume Fraction = %.4f" % vFracMean

# Setup directory structures
#(root, simdir, datadir, imgdir) = directoryStructureDevel(simdir)
(root, simdir, datadir, imgdir) = directoryStructureMarcc(simdir)

# Get time and z data
(time, tsInd, nt, evalZ, nz) = initData(datadir, tstart)

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

# Find phase averaged data
phaseVeldir = root + "simdata/phaseAveragedFluidVel"
phaseVelData = np.genfromtxt(phaseVeldir, skip_header=1)

for nn in np.arange(0,np.size(phaseVelData,0)):
  currN = phaseVelData[nn,0]
  currRho = phaseVelData[nn,1]
  if (nparts == currN) & (rho == currRho):
    phaseVel = phaseVelData[nn,2]
print "\n      phaseVel = %.4f" % phaseVel

# Find output data -- each column is a different time
vFracFile = datadir + "volume-fraction"
wVelFile = datadir + "part-w"
vFrac = np.genfromtxt(vFracFile).T[:,tsInd:]
wp = np.genfromtxt(wVelFile).T[:,tsInd:] #/ phaseVel

vfTemp = np.squeeze(np.reshape(vFrac,(-1,1)))
wpTemp = np.squeeze(np.reshape(wp,(-1,1)))

# Data stats
vfMin = np.min(vFrac)
vfMax = np.max(vFrac)
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
# 
# Histogram
nBins = 30
vfEdges = np.linspace(vfMin, vfMax, nBins)
wpEdges = np.linspace(wpMin, wpMax, nBins)

vfCenters = 0.5*(vfEdges[1:] + vfEdges[:-1])
wpCenters = 0.5*(wpEdges[1:] + wpEdges[:-1])

Hvf,_ = np.histogram(vfTemp, bins=vfEdges)
Hwp,_ = np.histogram(wpTemp, bins=wpEdges)
# H2:   phi increasing left to right
# H2:   wp  increasing top to bottom
H2,_,_ = np.histogram2d(vfTemp,wpTemp, bins=(vfEdges, wpEdges))

# Weighted least squares fit
phiEval = np.linspace(vfMin, vfMax, 50)

vfX,wpY = np.meshgrid(vfCenters, wpCenters)   # use centers as x,y coords
vfX = np.squeeze(np.reshape(vfX,(-1,1)))
wpY = np.squeeze(np.reshape(wpY,(-1,1)))
H2temp = np.squeeze(np.reshape(H2,(-1,1)))    # use hist counts as weight

coeffs = np.polynomial.polynomial.polyfit(vfX,wpY,1,w=H2temp)
c = coeffs[0]
m = coeffs[1]

wpEval = m*phiEval + c

print "\n      wp' = %.3f phi + %.3f" % (m, c)

# normalize to probablilities
inorm = 1./float(np.size(vfTemp))
Hvf = Hvf.astype(float)*inorm
Hwp = Hwp.astype(float)*inorm
H2 = H2.astype(float)*inorm

# Plot
histFig = plt.figure(figsize=(5,4))

ax1 = histFig.add_subplot(223)
ax1.plot(vfCenters, Hvf, 'o--')
ax1.plot([vFracMean, vFracMean], [0, np.max(Hvf)], 'k--')

ax1.set_xlabel(r'$\phi$')
ax1.set_ylabel(r'$P(\phi)$')
ax1.set_xticks(np.linspace(vfMin,vfMax,4))
ax1.xaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
ax1.set_yticks(np.linspace(0, np.max(Hvf), 4))
ax1.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))

ax2 = histFig.add_subplot(222)
ax2.plot(Hwp, wpCenters, 'o--')
ax2.plot([0, np.max(Hwp)], [recMeanWP, recMeanWP], 'k--')

ax2.set_xlabel(r'$P(w_p^\prime)$')#/\langle w_f \rangle)$')
ax2.set_ylabel(r'$w_p^\prime$')#/\langle w_f \rangle$')
ax2.set_yticks(np.linspace(wpMin,wpMax,4))
ax2.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
ax2.set_xticks(np.linspace(0, np.max(Hwp),4))
ax2.xaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))

ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")

ax3 = histFig.add_subplot(221, sharex=ax1, sharey=ax2)
im = ax3.imshow(H2, origin="lower", aspect="auto", interpolation="none",
  extent=[vfMin, vfMax, wpMin, wpMax])
ax3.plot([vFracMean, vFracMean], [wpMin, wpMax], 'k--')
ax3.plot([vfMin, vfMax], [recMeanWP, recMeanWP], 'k--')
ax3.plot(vfCenters, wpCenters[np.argmax(H2, 0)], 'k-', alpha=0.6)
ax3.plot(phiEval, wpEval, 'w--')

plt.setp( ax3.get_xticklabels(), visible=False)
plt.setp( ax3.get_yticklabels(), visible=False)
ax3.set_xlim([vfMin, vfMax])
ax3.set_ylim([wpMin, wpMax])

imgname = imgdir + "hist-vfrac-wvel"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

print "\n      ...Done!"
