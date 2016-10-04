#!/usr/bin/env python2

from setup import *

os.system('clear')

print ""
print " ---- Shape/Strain Alignment Plotting Utility ---- "
print "                     Mean"
print ""
root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
print "      Sim root directory set to: " + root

# Setup directory structure and print
root = os.path.expanduser("~") + "/scratch/triply_per/"
imgdir = root + "simdata/img/tetrads/"
datadir = "/analysis/tetrads/data/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)
print " Root dir: " + root
print " Sim dir:  " + imgdir

# Simulation information
partR = 2.1   # XXX mm
nu = 0.01715  # XXX mm^2/ms

# get phase averaged data
phaseVelDir = root + "simdata/phaseAveragedFluidVel"
phaseVelData = np.genfromtxt(phaseVelDir, skip_header=1)

# get terminal vel data
wTermDir = root + "simdata/singlePartSedi"
wTermData = np.genfromtxt(wTermDir, skip_header=1)

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
g1_s1 = [ structtype() for i in range(nSims) ]
g2_s1 = [ structtype() for i in range(nSims) ]
g3_s1 = [ structtype() for i in range(nSims) ]
g1_s2 = [ structtype() for i in range(nSims) ]
g2_s2 = [ structtype() for i in range(nSims) ]
g3_s2 = [ structtype() for i in range(nSims) ]
g1_s3 = [ structtype() for i in range(nSims) ]
g2_s3 = [ structtype() for i in range(nSims) ]
g3_s3 = [ structtype() for i in range(nSims) ]
g1_z = [ structtype() for i in range(nSims) ]
g2_z = [ structtype() for i in range(nSims) ]
g3_z = [ structtype() for i in range(nSims) ]
s1_z = [ structtype() for i in range(nSims) ]
s2_z = [ structtype() for i in range(nSims) ]
s3_z = [ structtype() for i in range(nSims) ]
w_z = [ structtype() for i in range(nSims) ]
w_g1 = [ structtype() for i in range(nSims) ]
w_g2 = [ structtype() for i in range(nSims) ]
w_g3 = [ structtype() for i in range(nSims) ]
w_s1 = [ structtype() for i in range(nSims) ]
w_s2 = [ structtype() for i in range(nSims) ]
w_s3 = [ structtype() for i in range(nSims) ]
 
# Loop over all directory's
tau_p = np.zeros(nSims)
for ii in np.arange(nSims):
  simdir = root + simList[ii] + datadir
  infoFile = simdir + 'info.dat'
  nodeFile = simdir + 'regularNodes'
  alignMean = simdir + 'align.mean'

  nTetrads[ii] = np.genfromtxt(infoFile, skip_header=1, usecols=0)

  tmpT = np.genfromtxt(alignMean, skip_header=1, usecols=0)
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
  # terminal velocity
  for nn in np.arange(0, np.size(wTermData, 0)):
    currRho = wTermData[nn,0]
    if (rho == currRho):
      termVel = wTermData[nn,1]
      Re_t = 2.*partR*wTermData[nn,1]/nu

  # normalize by wf/a
  data[ii].tau *= phaseVel/(2.*partR)      # ms*(mm/ms)/mm = []  
  normText = r"$t\langle w_f \rangle/2a$"

  # calculate stokes relax time, (rho* d^2)/(18nu) / f(Ret)
  # mm^2 / (mm^s/ms) -> [ms]
  #tau_p[ii] = rho*(2.*partR)*(2.*partR)/(18.*nu)
  #tau_p[ii] *= termVel/phaseVel
  #tau_p[ii] /= (1. + 0.1935*Re_t**(0.6305))
  #tau_p[ii] *= phaseVel/(2.*partR)
  tau_p[ii] = rho/18. * Re_t/(1. + 0.1935*Re_t**(0.6305))
  #data[ii].tau /= taup
  #normText = r"$t 18\nu/(\rho^* (2a)^2)$"

  # MEAN
  g1_s1[ii].mean = np.zeros(np.size(tmpT))
  g2_s1[ii].mean = np.zeros(np.size(tmpT))
  g3_s1[ii].mean = np.zeros(np.size(tmpT))
  g1_s2[ii].mean = np.zeros(np.size(tmpT))
  g2_s2[ii].mean = np.zeros(np.size(tmpT))
  g3_s2[ii].mean = np.zeros(np.size(tmpT))
  g1_s3[ii].mean = np.zeros(np.size(tmpT))
  g2_s3[ii].mean = np.zeros(np.size(tmpT))
  g3_s3[ii].mean = np.zeros(np.size(tmpT))
  g1_z[ii].mean = np.zeros(np.size(tmpT))
  g2_z[ii].mean = np.zeros(np.size(tmpT))
  g3_z[ii].mean = np.zeros(np.size(tmpT))
  s1_z[ii].mean = np.zeros(np.size(tmpT))
  s2_z[ii].mean = np.zeros(np.size(tmpT))
  s3_z[ii].mean = np.zeros(np.size(tmpT))
  w_z[ii].mean = np.zeros(np.size(tmpT))
  w_g1[ii].mean = np.zeros(np.size(tmpT))
  w_g2[ii].mean = np.zeros(np.size(tmpT))
  w_g3[ii].mean = np.zeros(np.size(tmpT))
  w_s1[ii].mean = np.zeros(np.size(tmpT))
  w_s2[ii].mean = np.zeros(np.size(tmpT))
  w_s3[ii].mean = np.zeros(np.size(tmpT))

  g1_s1[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=1)
  g2_s1[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=2)
  g3_s1[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=3)
  g1_s2[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=4)
  g2_s2[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=5)
  g3_s2[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=6)
  g1_s3[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=7)
  g2_s3[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=8)
  g3_s3[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=9)

  g1_z[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=10)
  g2_z[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=11)
  g3_z[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=12)
  s1_z[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=13)
  s2_z[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=14)
  s3_z[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=15)

  w_z[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=16)

  w_g1[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=17)
  w_g2[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=18)
  w_g3[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=19)
  w_s1[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=20)
  w_s2[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=21)
  w_s3[ii].mean = np.genfromtxt(alignMean, skip_header=1, usecols=22)

  legendText[ii] = simList[ii] # + ': ' + str(nTetrads[ii])

print tau_p

# Plot things
labelx = -0.2
labely = 0.5

# Legend specs
rho2_spec = mlines.Line2D([],[], color='k', alpha=shades[0],
 label=r'$\rho^* = 2.0$')
rho3_spec = mlines.Line2D([],[], color='k', alpha=shades[1],
 label=r'$\rho^* = 3.3$')
rho4_spec = mlines.Line2D([],[], color='k', alpha=shades[2],
 label=r'$\rho^* = 4.0$')
rho5_spec = mlines.Line2D([],[], color='k', alpha=shades[3],
 label=r'$\rho^* = 5.0$')

align_spec_1 = mlines.Line2D([],[], color='k')
align_spec_2 = mlines.Line2D([],[], color='r')
align_spec_3 = mlines.Line2D([],[], color='b')

#### Alignment of principal axes with initial principal strain ####
gisj_fig = plt.figure()

## gi_s1 ##
gis1_ax = gisj_fig.add_subplot(311)
for ii in np.arange(nSims):
  gis1_ax.semilogx(data[ii].tau, g1_s1[ii].mean, 
    color=colors[ii], alpha=shades[ii])
  gis1_ax.semilogx(data[ii].tau, g2_s1[ii].mean, 
    color=colors[ii], alpha=shades[ii], dashes=[2,2])
  gis1_ax.semilogx(data[ii].tau, g3_s1[ii].mean, 
    color=colors[ii], alpha=shades[ii], dashes=[5,2])

#gis1_ax.set_xlim([1e-1, 6e2])
gis1_ax.xaxis.set_ticklabels([])

gis1_ax.set_ylim([0, 1])
gis1_ax.set_ylabel(r'$(\mathbf{v}_k,\ \mathbf{s}_1(0))$')
gis1_ax.yaxis.set_major_locator(MultipleLocator(0.5))
gis1_ax.yaxis.set_minor_locator(MultipleLocator(0.25))

gis1_ax.text(2e0, 0.82, r'$\mathbf{v}_1$', fontsize=9)
gis1_ax.text(2e0, 0.33, r'$\mathbf{v}_2$', fontsize=9)
gis1_ax.text(2e0, 0.10, r'$\mathbf{v}_3$', fontsize=9)

## gi_s2 ##
gis2_ax = gisj_fig.add_subplot(312)
for ii in np.arange(nSims):
  gis2_ax.semilogx(data[ii].tau, g1_s2[ii].mean, 
    color=colors[ii], alpha=shades[ii])
  gis2_ax.semilogx(data[ii].tau, g2_s2[ii].mean, 
    color=colors[ii], alpha=shades[ii], dashes=[2,2])
  gis2_ax.semilogx(data[ii].tau, g3_s2[ii].mean, 
    color=colors[ii], alpha=shades[ii], dashes=[5,2])

#gis2_ax.set_xlim([1e-1, 6e2])
gis2_ax.xaxis.set_ticklabels([])

gis2_ax.set_ylim([0, 1])
gis2_ax.set_ylabel(r'$(\mathbf{v}_k,\ \mathbf{s}_2(0))$')
gis2_ax.yaxis.set_major_locator(MultipleLocator(0.5))
gis2_ax.yaxis.set_minor_locator(MultipleLocator(0.25))

gis2_ax.text(2e0, 0.75, r'$\mathbf{v}_2$', fontsize=9)
gis2_ax.text(2e0, 0.25, r'$\mathbf{v}_1$', fontsize=9)
gis2_ax.text(4e0, 0.15, r'$\mathbf{v}_3$', fontsize=9)

## gi_s3 ##
gis3_ax = gisj_fig.add_subplot(313)
for ii in np.arange(nSims):
  gis3_ax.semilogx(data[ii].tau, g1_s3[ii].mean, 
    color=colors[ii], alpha=shades[ii])
  gis3_ax.semilogx(data[ii].tau, g2_s3[ii].mean, 
    color=colors[ii], alpha=shades[ii], dashes=[2,2])
  gis3_ax.semilogx(data[ii].tau, g3_s3[ii].mean, 
    color=colors[ii], alpha=shades[ii], dashes=[5,2])

#gis3_ax.set_xlim([1e-1, 6e2])
gis3_ax.set_xlabel(normText)

gis3_ax.set_ylim([0, 1])
gis3_ax.set_ylabel(r'$(\mathbf{v}_k,\ \mathbf{s}_3(0))$')
gis3_ax.yaxis.set_major_locator(MultipleLocator(0.5))
gis3_ax.yaxis.set_minor_locator(MultipleLocator(0.25))

gis3_ax.text(2e0, 0.83, r'$\mathbf{v}_3$', fontsize=9)
gis3_ax.text(2e0, 0.47, r'$\mathbf{v}_2$', fontsize=9)
gis3_ax.text(4e0, 0.05, r'$\mathbf{v}_1$', fontsize=9)

# stokes number
for ii in np.arange(nSims):
  gis1_ax.semilogx(tau_p[ii]*np.ones(2), [-1.,1.], color=colors[ii], alpha=shades[ii])

# Save
imgname = imgdir + "all_align_gisj_mean"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## gi_z ##
giz_fig = plt.figure(figsize=(4,2))
giz_ax = giz_fig.add_subplot(111)

for ii in np.arange(nSims):
  giz_ax.semilogx(data[ii].tau, g1_z[ii].mean, 
    color=colors[ii], alpha=shades[ii])
  giz_ax.semilogx(data[ii].tau, g2_z[ii].mean, 
    color=colors[ii], alpha=shades[ii], dashes=[2,2])
  giz_ax.semilogx(data[ii].tau, g3_z[ii].mean, 
    color=colors[ii], alpha=shades[ii], dashes=[5,2])

giz_ax.set_xlim([1e-1, 6e2])    
giz_ax.set_xlabel(normText)

giz_ax.set_ylim([0,1])
giz_ax.set_ylabel(r'$(\mathbf{v}_k,\ \mathbf{e}_z)$')
giz_ax.yaxis.set_major_locator(MultipleLocator(0.25))
giz_ax.yaxis.set_minor_locator(MultipleLocator(0.125))

giz_ax.text(3e1, 0.75, r'$\mathbf{v}_1$', fontsize=9)
giz_ax.text(7e0, 0.32, r'$\mathbf{v}_2$', fontsize=9)
giz_ax.text(2e0, 0.65, r'$\mathbf{v}_3$', fontsize=9)

# Save
imgname = imgdir + "all_align_giz_mean"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')


## si_z ##
siz_fig = plt.figure()
siz_ax = siz_fig.add_subplot(111)

for ii in np.arange(nSims):
  siz_ax.semilogx(data[ii].tau, s1_z[ii].mean, color=colors[ii], 
    alpha=shades[ii])
  siz_ax.semilogx(data[ii].tau, s2_z[ii].mean, color=colors[ii], 
    alpha=shades[ii], dashes=[2, 2])
  siz_ax.semilogx(data[ii].tau, s3_z[ii].mean, color=colors[ii], 
    alpha=shades[ii], dashes=[5, 2])

siz_ax.set_xlim([1e-1, 6e2])
siz_ax.set_xlabel(normText)

siz_ax.set_ylim([0, 1])
siz_ax.set_ylabel(r'$(\mathbf{s}_k,\ \mathbf{e}_z)$')

# Save
imgname = imgdir + "all_align_siz_mean"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## si_w ##
siw_fig = plt.figure()
siw_ax = siw_fig.add_subplot(111)

for ii in np.arange(nSims):
  siw_ax.semilogx(data[ii].tau, w_s1[ii].mean, color=colors[ii], 
    alpha=shades[ii])
  siw_ax.semilogx(data[ii].tau, w_s2[ii].mean, color=colors[ii], 
    alpha=shades[ii], dashes=[2, 2])
  siw_ax.semilogx(data[ii].tau, w_s3[ii].mean, color=colors[ii], 
    alpha=shades[ii], dashes=[5, 2])

siw_ax.set_xlim([1e-1, 6e2])
siw_ax.set_xlabel(normText)

#siw_ax.set_ylim([0, 1])
siw_ax.set_ylabel(r'$(\mathbf{s}_k,\ \mathbf{\omega})$')

# Save
imgname = imgdir + "all_align_siw_mean"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
