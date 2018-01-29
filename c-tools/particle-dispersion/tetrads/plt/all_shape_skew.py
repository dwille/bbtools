#!/usr/bin/env python2
from setup import *

os.system('clear')

print ""
print " ---- Anisotropy Measures Plotting Utility ---- "
print "                   Skewness"
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

# get phase averaged data
phaseVelDir = root + "simdata/phaseAveragedFluidVel"
phaseVelData = np.genfromtxt(phaseVelDir, skip_header=1)

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
 
# Loop over all directory's
for ii in np.arange(nSims):
  simdir = root + simList[ii] + datadir
  infoFile = simdir + 'info.dat'
  nodeFile = simdir + 'regularNodes'
  statskew = simdir + 'stat.skew'

  nTetrads[ii] = np.genfromtxt(infoFile, skip_header=1, usecols=0)

  # Pull time [ms]
  tmpT = np.genfromtxt(statskew, skip_header=1, usecols=0)
  data[ii].time = np.zeros(np.size(tmpT))
  data[ii].time = tmpT - tmpT[0]
  data[ii].tau = np.zeros(np.size(tmpT))
  data[ii].tau = tmpT - tmpT[0]

  # Normalize time
  nparts = int(simList[ii].partition('/')[0])
  rho = float(simList[ii].partition('/')[2][-3:])
  for nn in np.arange(0,np.size(phaseVelData, 0)):
    currN = phaseVelData[nn,0]
    currRho = phaseVelData[nn,1]
    if (nparts == currN) & (rho == currRho):
      phaseVel = phaseVelData[nn,2]   # mm/ms

  data[ii].tau *= phaseVel/partR      # ms*(mm/ms)/mm = []  

  # skew
  RoG[ii].skew   = np.zeros(np.size(tmpT))
  EVar[ii].skew  = np.zeros(np.size(tmpT))
  Shape[ii].skew = np.zeros(np.size(tmpT))
  I1[ii].skew    = np.zeros(np.size(tmpT))
  I2[ii].skew    = np.zeros(np.size(tmpT))
  I3[ii].skew    = np.zeros(np.size(tmpT))

  RoG[ii].skew   = np.genfromtxt(statskew, skip_header=1, usecols=1)
  EVar[ii].skew  = np.genfromtxt(statskew, skip_header=1, usecols=2)
  Shape[ii].skew = np.genfromtxt(statskew, skip_header=1, usecols=3)
  I1[ii].skew    = np.genfromtxt(statskew, skip_header=1, usecols=4)
  I2[ii].skew    = np.genfromtxt(statskew, skip_header=1, usecols=5)
  I3[ii].skew    = np.genfromtxt(statskew, skip_header=1, usecols=6)

  legendText[ii] = simList[ii] # + ': ' + str(nTetrads[ii])

## Radius of Gyration
rgFig = plt.figure()
rg_ax = rgFig.add_subplot(111)
for ii in np.arange(nSims):
  rg_ax.semilogx(data[ii].tau, RoG[ii].skew, color=colors[ii], alpha=shades[ii])

#rg_ax.legend(legendText, ncol=1, loc='center left', bbox_to_anchor=(1.10,0.5))
rg_ax.set_xlabel(r"$t\langle w_f \rangle/a$")
rg_ax.set_ylabel(r"$\mathrm{Skew}[R_g]$")

# Save
imgname = imgdir + "all_shape_skew_rog"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

# EVar
vFig = plt.figure()
v_ax = vFig.add_subplot(111)
for ii in np.arange(nSims):
  v_ax.semilogx(data[ii].tau, EVar[ii].skew, color=colors[ii], 
    alpha=shades[ii])

v_ax.set_xlabel(r"$t\langle w_f \rangle/a$")
v_ax.set_ylabel(r"$\mathrm{Skew}[\Delta]$")
#v_ax.legend(legendText, ncol=1, loc='center left', bbox_to_anchor=(1.10,0.5))

# Save
imgname = imgdir + "all_shape_skew_var"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

# Shape
sFig = plt.figure()
s_ax = sFig.add_subplot(111)
for ii in np.arange(nSims):
  s_ax.semilogx(data[ii].tau, Shape[ii].skew, color=colors[ii], 
    alpha=shades[ii])

s_ax.set_xlabel(r"$t\langle w_f \rangle/a$")
s_ax.set_ylabel(r"$\mathrm{Skew}[S]$")
#s_ax.legend(legendText, ncol=1, loc='center left', bbox_to_anchor=(1.10,0.5))

# Save
imgname = imgdir + "all_shape_skew_sf"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## I_j ##
iFig = plt.figure()
i_ax = iFig.add_subplot(111)
for ii in np.arange(nSims):
  i_ax.semilogx(data[ii].tau, I1[ii].skew, color=colors[ii],
    alpha=shades[ii])
  i_ax.semilogx(data[ii].tau, I2[ii].skew, color=colors[ii],
    alpha=shades[ii], label='_nolegend_')
  i_ax.semilogx(data[ii].tau, I3[ii].skew, color=colors[ii],
    alpha=shades[ii], label='_nolegend_')

i_ax.text(1.25, -0.50, 'I1')
i_ax.text(1.25, 0.00, 'I2')
i_ax.text(1.25, 0.50, 'I3')
i_ax.set_xlabel(r"$t\langle w_f \rangle/a$")
i_ax.set_ylabel(r"$\mathrm{Skew}[I_j]$")

#i_ax.legend(legendText, ncol=1, loc='center left', bbox_to_anchor=(1.10,0.5))

# Save
imgname = imgdir + "all_shape_skew_normi"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
