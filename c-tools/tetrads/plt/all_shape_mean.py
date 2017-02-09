#!/usr/bin/env python2
from setup import *

os.system('clear')

print ""
print " ---- Anisotropy Measures Plotting Utility ---- "
print "                     Mean"
print ""

# Setup directory structure and print
root = os.path.expanduser("~") + "/scratch/triply_per/"
imgdir = root + "simdata/img/tetrads/"
datadir = "/analysis/tetrads/data/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)
print " Root dir: " + root
print " Img dir:  " + imgdir

# Simulation information
partR = 2.1 # XXX mm
nu = 0.01715  # XXX mm^2/ms

# get phase averaged data
phaseVelDir = root + "simdata/phaseAveragedFluidVel"
phaseVelData = np.genfromtxt(phaseVelDir, skip_header=1)

# get terminal vel data
wTermDir = root + "simdata/singlePartSedi"
wTermData = np.genfromtxt(wTermDir, skip_header=1)

# pull variance ratios
varFile = root + "simdata/part_vel_variance"
vert_var = np.genfromtxt(varFile, skip_header=1, usecols=2)
hori_var = np.genfromtxt(varFile, skip_header=1, usecols=3)
var_ratio = vert_var / hori_var

# Parameter sweep -- include dead sims as '' so we get plot colors correct
simList = ['500/rho2.0', '500/rho3.3', '500/rho4.0', '500/rho5.0',
  '1000/rho2.0', '1000/rho3.3', '1000/rho4.0', '1000/rho5.0',
  '1500/rho2.0', '1500/rho3.3', '1500/rho4.0', '',
  '2000/rho2.0', '2000/rho3.3', '', '']
nSims = len(simList)

# Get colors for plots TODO this is terribly un-general
baseColors = ['r', 'g', 'b', 'k']
baseShades = [0.4, 0.57, 0.74, 0.9]
colors = ['']*nSims
shades = ['']*nSims
for ii in range(nSims):
   pp = ii / 4
   dd = ii % 4
   if simList[ii] != '':
     colors[ii] = baseColors[pp]
     shades[ii] = baseShades[dd]
colors = filter(lambda colors: colors != '', colors)     
shades = filter(lambda shades: shades != '', shades)     
 
# Now remove the dead sims for the rest
var_ratio = filter(lambda var_ratio: var_ratio > 0, var_ratio)
simList = filter(lambda simList: simList != '', simList)
nSims = len(simList)
legendText = ['']*nSims

# Initialize arrays
nTetrads = np.zeros(nSims)
data = [ structtype() for i in range(nSims) ]
RoG = [ structtype() for i in range(nSims) ]
EVar = [ structtype() for i in range(nSims) ]
Shape = [ structtype() for i in range(nSims) ]
I1 = [ structtype() for i in range(nSims) ]
I2 = [ structtype() for i in range(nSims) ]
I3 = [ structtype() for i in range(nSims) ]
S11 = [ structtype() for i in range(nSims) ]
S22 = [ structtype() for i in range(nSims) ]
S33 = [ structtype() for i in range(nSims) ]

# Loop over all directories
tau_p = np.zeros(nSims)
for ii in np.arange(nSims):
  simdir = root + simList[ii] + datadir
  infoFile = simdir + 'info.dat'
  nodeFile = simdir + 'regularNodes'
  statmean = simdir + 'stat.mean'

  nTetrads[ii] = np.genfromtxt(infoFile, skip_header=1, usecols=0)

  # Pull time [ms]
  tmpT = np.genfromtxt(statmean, skip_header=1, usecols=0)
  data[ii].time = np.zeros(np.size(tmpT))
  data[ii].time = tmpT - tmpT[0]
  data[ii].tau = np.zeros(np.size(tmpT))
  data[ii].tau = tmpT - tmpT[0]

  # Normalize time
  nparts = int(simList[ii].partition('/')[0])
  rho = float(simList[ii].partition('/')[2][-3:])
  # phase average fluid velocity
  for nn in np.arange(0,np.size(phaseVelData, 0)):
    currN = phaseVelData[nn,0]
    currRho = phaseVelData[nn,1]
    if (nparts == currN) & (rho == currRho):
      phaseVel = phaseVelData[nn,2]   # mm/ms
  # terminal velocity
  for nn in np.arange(0, np.size(wTermData, 0)):
    currRho = wTermData[nn,0]
    if (rho == currRho):
      termVel = wTermData[nn,1]
      Re_t = 2.*partR*wTermData[nn,1]/nu

  # normalize by wf/a
  data[ii].tau *= phaseVel/(2.*partR)      # ms*(mm/ms)/mm = []  
  normText = r"$t\langle w_f \rangle/a$"

  # calculate stokes relax time, (rho* d^2)/(18nu) / f(Ret)
  # (in terms of wf/2a)
  # mm^2 / (mm^s/ms) -> [ms]
  tau_p[ii] = (2.*rho + 1.)/18. * Re_t/(1. + 0.1935*Re_t**(0.6305))

  # mean
  RoG[ii].mean   = np.zeros(np.size(tmpT))
  EVar[ii].mean  = np.zeros(np.size(tmpT))
  Shape[ii].mean = np.zeros(np.size(tmpT))
  I1[ii].mean    = np.zeros(np.size(tmpT))
  I2[ii].mean    = np.zeros(np.size(tmpT))
  I3[ii].mean    = np.zeros(np.size(tmpT))
  S11[ii].mean    = np.zeros(np.size(tmpT))
  S22[ii].mean    = np.zeros(np.size(tmpT))
  S33[ii].mean    = np.zeros(np.size(tmpT))

  RoG[ii].mean   = np.genfromtxt(statmean, skip_header=1, usecols=1)
  EVar[ii].mean  = np.genfromtxt(statmean, skip_header=1, usecols=2)
  Shape[ii].mean = np.genfromtxt(statmean, skip_header=1, usecols=3)
  I1[ii].mean    = np.genfromtxt(statmean, skip_header=1, usecols=4)
  I2[ii].mean    = np.genfromtxt(statmean, skip_header=1, usecols=5)
  I3[ii].mean    = np.genfromtxt(statmean, skip_header=1, usecols=6)
  S11[ii].mean    = np.genfromtxt(statmean, skip_header=1, usecols=7)
  S22[ii].mean    = np.genfromtxt(statmean, skip_header=1, usecols=8)
  S33[ii].mean    = np.genfromtxt(statmean, skip_header=1, usecols=9)

  legendText[ii] = simList[ii] # + ': ' + str(nTetrads[ii])

## Radius of Gyration -- unscaled
rgFig = plt.figure()
rg_ax = rgFig.add_subplot(111)
for ii in np.arange(nSims):
  rg_ax.loglog(data[ii].time, RoG[ii].mean,
    color=colors[ii], alpha=shades[ii])

  # For estimating ensemble duration and lag
  #print simList[ii]
  #if ii == 0:
  #  #print "est dur = 6000"
  #else:
  #  #print data[ii].time[np.argwhere(RoG[ii].mean >= 42)[0]]
  #initRoG = RoG[ii].mean[0]
  #print data[ii].time[np.argwhere(RoG[ii].mean >= 2.*initRoG)[0]]

# Scaling
#xpnts = np.array([1e2, 1e4])
#ypnts = np.power(xpnts, 0.54) / 3.3
#rg_ax.loglog(xpnts, ypnts, 'k-.', linewidth=2, dashes=[4,2])
#rg_ax.text(400, 4, r'$t^{1/2}$')

# Domain extents
xl = 42 # XXX
#rg_ax.loglog(rg_ax.get_xlim(), [xl, xl], "k--")

#rg_ax.legend(legendText, ncol=2, loc='center left', bbox_to_anchor=(1.10,0.5))
rg_ax.set_xlabel(r"$t\ [ms]$")
rg_ax.set_ylabel(r"$\langle R_G \rangle \ [mm]$")

# Save
imgname = imgdir + "all_shape_mean_rog_unscaled"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## Radius of Gyration -- scaled
rgFig = plt.figure()
rg_ax = rgFig.add_subplot(111)
for ii in np.arange(nSims):
  alpha = 2./(2. + var_ratio[ii])
  #print alpha
  y_data = RoG[ii].mean / np.power(data[ii].tau, alpha)
  rg_ax.loglog(data[ii].tau, y_data/partR,
    color=colors[ii], alpha=shades[ii])

  # For estimating ensemble duration and lag
  #print simList[ii]
  #if ii == 0:
  #  #print "est dur = 6000"
  #else:
  #  #print data[ii].time[np.argwhere(RoG[ii].mean >= 42)[0]]
  #initRoG = RoG[ii].mean[0]
  #print data[ii].time[np.argwhere(RoG[ii].mean >= 2.*initRoG)[0]]

xpnts = np.array([100, 10000])
ypnts = np.power(xpnts, 0.6) / 5.
#rg_ax.loglog(xpnts, ypnts, 'k--', linewidth=2)
#rg_ax.text(25, 2.5, 'slope = 0.6')

xpnts = np.array([6e0, 5e2])
ypnts = np.power(xpnts, 0.5)  /1.25
#rg_ax.loglog(xpnts, ypnts, 'k-.', linewidth=2, dashes=[4,2])
rg_ax.text(15, 2, r'$t^{1/2}$')

# Domain extents
xl = 42./partR
#rg_ax.loglog(rg_ax.get_xlim(), [xl, xl], "k--")

#rg_ax.legend(legendText, ncol=2, loc='center left', bbox_to_anchor=(1.10,0.5))
rg_ax.set_xlabel(normText)
rg_ax.set_ylabel(r"$\langle R_G \rangle/a$")

rg_ax.set_xlim([4e0, 1e2])
rg_ax.set_ylim([1e0, 2e0])

# response time
for ii in np.arange(4):
  rg_ax.semilogx(tau_p[ii]*np.ones(2), [-1.,1.], color=colors[ii], alpha=shades[ii])

# Save
imgname = imgdir + "all_shape_mean_rog_scaled"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## Shape factor and variance ##
fig1 = plt.figure()
# EVar
v_ax = fig1.add_subplot(211)
for ii in np.arange(nSims):
  v_ax.semilogx(data[ii].tau, EVar[ii].mean, 
    color=colors[ii], alpha=shades[ii])

v_ax.set_xlim([1e-1, 6e2])
v_ax.xaxis.set_ticklabels([])
v_ax.set_ylim([0, 1])
v_ax.set_ylabel(r'$\langle \Delta \rangle$')
v_ax.yaxis.set_label_coords(-.15, 0.5)

v_ax.yaxis.set_major_locator(MultipleLocator(0.25))
v_ax.yaxis.set_minor_locator(MultipleLocator(0.125))

# Shape
s_ax = fig1.add_subplot(212)
for ii in np.arange(nSims):
  s_ax.semilogx(data[ii].tau, Shape[ii].mean, 
    color=colors[ii], alpha=shades[ii])

s_ax.set_xlim([1e-1, 6e2])
s_ax.set_xlabel(normText)
s_ax.set_ylim([-0.25, 2])
s_ax.set_ylabel(r'$\langle S \rangle$')
s_ax.yaxis.set_label_coords(-.15, 0.5)
s_ax.set_yticks([-0.25, 0, 0.5, 1, 1.5, 2])
s_ax.yaxis.set_minor_locator(MultipleLocator(0.25))

# response time
for ii in np.arange(4):
  s_ax.semilogx(tau_p[ii]*np.ones(2), [-1.,1.], color=colors[ii], alpha=shades[ii])

# Save
imgname = imgdir + "all_shape_mean_sf_var"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## I_j ##
iFig = plt.figure(figsize=(4,2))
ax1 = iFig.add_subplot(211)

for ii in np.arange(nSims):
  i_ratio = I1[ii].mean/I2[ii].mean
  i_ratio /= np.power(data[ii].tau, 1./3.)
  ax1.semilogx(data[ii].tau, i_ratio, color=colors[ii], alpha=shades[ii])

ax1.set_xlim([1e-1, 6e2])
ax1.set_xticklabels([])

ax1.set_ylabel(r"$(\langle I_1\rangle /\langle I_2\rangle)/t^{1/3} $", rotation=0)
ax1.yaxis.set_label_coords(-.25, 0.5)
#ax1.set_ylim([1, 10])

#xpnts = np.array([1e0, 2e1])
#ypnts = np.power(xpnts, 1./3.) * 2.
#ax1.semilogx(xpnts, ypnts, 'k--')

i_ax = iFig.add_subplot(212)
for ii in np.arange(nSims):
  i_ax.semilogx(data[ii].tau, I1[ii].mean, color=colors[ii], 
    alpha=shades[ii])
  i_ax.semilogx(data[ii].tau, I2[ii].mean, color=colors[ii], 
    alpha=shades[ii], dashes=[2,2])
  i_ax.semilogx(data[ii].tau, I3[ii].mean, color=colors[ii], 
    alpha=shades[ii], dashes=[5,2])

i_ax.text(5, 0.55, r'$I_1$')
i_ax.text(5, 0.35, r'$I_2$')
i_ax.text(5, 0.10, r'$I_3$')

i_ax.set_xlabel(normText)
i_ax.set_ylabel(r"$\langle I_k \rangle = \langle \lambda_k/R_G^2 \rangle$")
i_ax.set_ylim([0, 1])
i_ax.set_xlim([1e-1, 6e2])
#i_ax.legend(legendText, ncol=1, loc='center left', bbox_to_anchor=(1.10,0.5))

# response time
for ii in np.arange(4):
  i_ax.semilogx(tau_p[ii]*np.ones(2), [-1.,1.], color=colors[ii], alpha=shades[ii])

# Save
imgname = imgdir + "all_shape_mean_normi"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

# S11,S22,S33
sFig =  plt.figure()
ax1 = sFig.add_subplot(211)

for ii in np.arange(1):
  tr = S11[ii].mean + S22[ii].mean + S33[ii].mean
  ax1.semilogx(data[ii].tau, tr, 'k')
  ax1.semilogx(data[ii].tau, S11[ii].mean, color=colors[ii], alpha=shades[ii])# , #marker='.')
  ax1.semilogx(data[ii].tau, S22[ii].mean, color=colors[ii], alpha=shades[ii], #marker='s')
    dashes=[2,2])
  ax1.semilogx(data[ii].tau, S33[ii].mean, color=colors[ii], alpha=shades[ii], #marker='o')
    dashes=[5,2])

ax1.set_xlim([1e-1, 6e2])
ax1.set_xlabel(r"$\tau$")
ax1.set_ylim([-8e-2, 8e-2])
ax1.set_ylabel(r"$S_{k}$")

# Save
imgname = imgdir + "all_velgrad_eigs"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
