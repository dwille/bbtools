#!/usr/bin/env python2

import matplotlib.pyplot as plt
import numpy as np
import os,sys
os.system('clear')


print ""
print " ---- Phase-Averaged Fluid Velocity Plotting Utility"
print ""

# Simulation Params
d = 4.2
nu = .01715
rho = np.array([2.0, 3.3, 4.0, 5.0])
phi = np.array([0.087, 0.175, 0.262, 0.349])

# MARCC
home = os.path.expanduser("~")
root = home + "/scratch/triply_per/simdata/"
phaseFile = root + "phaseAveragedFluidVel"
termFile = root + "singlePartSedi"

print "      Root directory set to: " + root
print "      Using location: " + phaseFile
print "      Using location: " + termFile
print ""
print "      Using d = %.2f" % d
print "      Using nu = %.5f" % nu

# Check if datadir exists so we don't go creating extra dirs
if not os.path.exists(phaseFile):
  print "      " + phaseFile + " does not exist. Exiting..."
  print ""
  sys.exit()
elif not os.path.exists(termFile):
  print "      " + termFile + " does not exist. Exiting..."
  print ""
  sys.exit()

# Create imgdir if necessary
imgdir = root + "img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

imgdir = imgdir + "phaseAvgWf/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# Read data
flowVel = np.genfromtxt(phaseFile, skip_header=1, usecols=2)
termVel = np.genfromtxt(termFile, skip_header=1, usecols=1)

# Constant rho -- phi increasing
rho2_0 = flowVel[0:16:4]/termVel[0]
rho3_3 = flowVel[1:16:4]/termVel[1]
rho4_0 = flowVel[2:16:4]/termVel[2]
rho5_0 = flowVel[3:16:4]/termVel[3]

# Constant phi -- rho increasing
n0500 = flowVel[0:4]
n1000 = flowVel[4:8]
n1500 = flowVel[8:12]
n2000 = flowVel[12:16]

# R-Z relations
Ret = d*termVel/nu
# Richardson Zaki
#nRZ = 4.45*Ret**(-0.1)   

# Garside and Al-Dibouni
tmp = 0.1*Ret**0.9
nGAD = (5.1 + 2.7*tmp)/(1 + tmp)

# Determine average k -- (wf/wt) / (1-phi)^(n-1)
k = np.zeros(16)
for rr,_ in enumerate(rho):
  for aa,_ in enumerate(phi):
    # (wf/wt)
    if rr==0:     # rho=2.0
      lhs = rho2_0[aa]
    elif rr==1:   # rho=3.3
      lhs = rho3_3[aa]
    elif rr==2:   # rho=4.0
      lhs = rho4_0[aa]
    elif rr==3:   # rho=5.0
      lhs = rho5_0[aa]

    # (1-phi)^(n-1) -- from garside, al-dibouni
    rhs = (1 - phi[aa])**(nGAD[rr] - 1) 

    # k = lhs/rhs
    k[aa + 4*rr] = lhs/rhs

kMean = np.mean(k)
    
phiEval = np.linspace(0.05,0.4,81)
WfWtEval = np.zeros((np.size(phiEval),4))
for rr,_ in enumerate(rho):
  for pp,_ in enumerate(phiEval):
    WfWtEval[pp,rr] = kMean*(1 - phiEval[pp])**(nGAD[rr] - 1)

# Minimum Fluidization Velocity
# W_MF -- from Wen/Yu/Davidson/Harriston
wmf = np.array([0.0059785, 0.013575, 0.01766, 0.023])
wmf /= termVel
phimf = 0.6

## Plot Constant rho -- phi increasing ##
fig1 = plt.figure()

nAx = fig1.add_subplot(111)
nAx.plot(phi,rho2_0, 'b*', markersize=7, alpha=0.7)
nAx.plot(phi,rho3_3, 'gs', markersize=7, alpha=0.7)
nAx.plot(phi,rho4_0, 'ro', markersize=7, alpha=0.7)
nAx.plot(phi,rho5_0, 'c^', markersize=7, alpha=0.7)

nAx.set_xlabel(r'$\phi$', fontsize=14)
nAx.set_ylabel(r'$w_f / w_t$', fontsize=14)

nAx.set_xlim([0,0.45])
nAx.set_xticks([0, 0.15, 0.30, 0.45])
#nAx.set_ylim([0,0.7])
#nAx.set_yticks([0, 0.25, 0.5, 0.75])

nAx.grid(True)

lText = [r'$\rho^* = 2.0$', r'$\rho^* = 3.3$', 
         r'$\rho^* = 4.0$', r'$\rho^* = 5.0$']
nAx.legend(lText, bbox_to_anchor=(0,1.05,1,1),loc="lower left",mode="expand",
  ncol=2, borderaxespad=0)

nAx.plot(phiEval,WfWtEval[:,0], 'b--', zorder=1)
nAx.plot(phiEval,WfWtEval[:,1], 'g--', zorder=1)
nAx.plot(phiEval,WfWtEval[:,2], 'r--', zorder=1)
nAx.plot(phiEval,WfWtEval[:,3], 'c--', zorder=1)
labelText=r"$\frac{w_f}{w_t} = %.2f(1 - \phi)^{n - 1}$" % kMean
nAx.text(0.23,0.6,labelText)

#nAx.plot(phimf, wmf[0], 'bx')
#nAx.plot(phimf, wmf[1], 'gx')
#nAx.plot(phimf, wmf[2], 'rx')
#nAx.plot(phimf, wmf[3], 'cx')

#nAx.set_yscale("log")

# Subplot
a = plt.axes([0.25, 0.3, 0.2, 0.2])
plt.plot(phi,rho2_0, 'b*', markersize=7, alpha=0.7)
plt.plot(phi,rho3_3, 'gs', markersize=7, alpha=0.7)
plt.plot(phi,rho4_0, 'ro', markersize=7, alpha=0.7)
plt.plot(phi,rho5_0, 'c^', markersize=7, alpha=0.7)

plt.plot(phiEval,WfWtEval[:,0], 'b--', zorder=1)
plt.plot(phiEval,WfWtEval[:,1], 'g--', zorder=1)
plt.plot(phiEval,WfWtEval[:,2], 'r--', zorder=1)
plt.plot(phiEval,WfWtEval[:,3], 'c--', zorder=1)

plt.xlim([0.07, 0.11])
plt.ylim([0.69, 0.73])
plt.xticks([0.07, 0.09, 0.11])
plt.yticks([0.69, 0.71, 0.73])

imgname = imgdir + "phaseAverageWf-phi"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

# Constant phi -- rho increasing
fig2 = plt.figure()
rhoAx = fig2.add_subplot(111)
rhoAx.plot(rho, n0500, 'o')
rhoAx.plot(rho, n1000, 'o')
rhoAx.plot(rho, n1500, 'o')
rhoAx.plot(rho, n2000, 'o')

rhoAx.set_xlabel(r'$\rho^*$')
rhoAx.set_ylabel(r'$w_f\ [mm/ms]$')

rhoAx.set_xlim([1,6])
rhoAx.set_ylim([0,0.4])
rhoAx.set_yticks([0, 0.1, 0.2, 0.3, 0.4])

rhoAx.grid(True)

lText=[r'$\phi = 0.087$', r'$\phi = 0.175$', 
       r'$\phi = 0.262$',r'$\phi = 0.349$']
rhoAx.legend(lText, bbox_to_anchor=(0,1.1,1,1),loc="lower left",mode="expand",
  ncol=2, borderaxespad=0)

imgname = imgdir + "phaseAverageWf-rho"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
print "      ...Done!"
