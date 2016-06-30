#!/usr/bin/env python2

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy.polynomial.polynomial as poly
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

# Constant rho -- phi increasing -- left hand side
wfwt_rho20 = flowVel[0:16:4]/termVel[0]
wfwt_rho33 = flowVel[1:16:4]/termVel[1]
wfwt_rho40 = flowVel[2:16:4]/termVel[2]
wfwt_rho50 = flowVel[3:16:4]/termVel[3]

# Fit n,k for each density ratio
# Linear fit to log(wf/wt) = log(k) + (n-1)log(1-phi)
# Coefficients return as [A,B,C] in A + Bx + Cx^2...
[logk_20, nm1_20] = poly.polyfit(np.log(1. - phi), np.log(wfwt_rho20), 1)
[logk_33, nm1_33] = poly.polyfit(np.log(1. - phi), np.log(wfwt_rho33), 1)
[logk_40, nm1_40] = poly.polyfit(np.log(1. - phi), np.log(wfwt_rho40), 1)
[logk_50, nm1_50] = poly.polyfit(np.log(1. - phi), np.log(wfwt_rho50), 1)

k20 = np.exp(logk_20)
k33 = np.exp(logk_33)
k40 = np.exp(logk_40)
k50 = np.exp(logk_50)
n20 = nm1_20 + 1.
n33 = nm1_33 + 1.
n40 = nm1_40 + 1.
n50 = nm1_50 + 1.

#print "   2.0   3.3   4.0   5.0"
#print "n  %.2f  %.2f  %.2f  %.2f" % (n20,n33,n40,n50)
#print "k  %.2f  %.2f  %.2f  %.2f" % (k20,k33,k40,k50)

# R-Z relations -- n coefficient
Ret = d*termVel/nu

# Richardson Zaki
nRZ = 4.45*Ret**(-0.1)   
# Garside and Al-Dibouni
tmp = 0.1*Ret**0.9
nGAD = (5.1 + 2.7*tmp)/(1 + tmp)

# Determine average k -- (wf/wt) / (1-phi)^(n-1)
kGAD = np.zeros(16)
kRZ = np.zeros(16)
for rr,_ in enumerate(rho):
  for aa,_ in enumerate(phi):
    # (wf/wt)
    if rr==0:     # rho=2.0
      lhs = wfwt_rho20[aa]
    elif rr==1:   # rho=3.3
      lhs = wfwt_rho33[aa]
    elif rr==2:   # rho=4.0
      lhs = wfwt_rho40[aa]
    elif rr==3:   # rho=5.0
      lhs = wfwt_rho50[aa]

    # (1-phi)^(n-1) -- from garside, al-dibouni
    rhsGAD = (1 - phi[aa])**(nGAD[rr] - 1) 
    rhsRZ = (1 - phi[aa])**(nRZ[rr] - 1) 

    # k = lhs/rhs
    kGAD[aa + 4*rr] = lhs/rhsGAD
    kRZ[aa + 4*rr] = lhs/rhsRZ

kGAD_mean = np.mean(kGAD)
kRZ_mean = np.mean(kRZ)

## Plot k as a function of phi and rho ##
kFig = plt.figure(figsize=(2,2))

# Yin and Koch, 2007, data
Re_yk = [1.02, 1.83, 5.00, 9.98, 20.0]
k_yk = [0.92, 0.90, 0.86, 0.86, 0.88]

# as a function of rho 
ax2 = kFig.add_subplot(111)
ax2.plot(Ret, [k20, k33, k40, k50], 'k.', markersize=5)
ax2.plot(Re_yk, k_yk, 'bo', markersize=5, markerfacecolor='none')

ax2.set_xlim([0,130])
ax2.set_xticks([0,25,50,75,100,125])
ax2.set_xticklabels(["0","","50","","100",""])
ax2.set_xlabel(r"$Re_t$")
ax2.set_ylim([0.75, 0.95])
ax2.set_yticks([0.75, 0.8, 0.85, 0.9, 0.95])
ax2.set_yticklabels(["0.75", "", "0.85", "", "0.95"])
ax2.set_ylabel(r"$\kappa$")

ax2.grid(True)

#lText = [r'$\kappa_{WSP}$', '$\kappa_{YK}$']
#ax2.legend(lText, bbox_to_anchor=(0,1.05,1,1),loc="lower left",mode="expand",
#  ncol=2, borderaxespad=0)

imgname = imgdir + "k-relations"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## Plot n as a function of rho ##
Ret_eval = np.arange(0.1, 125, 0.1)
nRZ_eval = 4.45*Ret_eval**(-0.1)   
nGAD_eval = (5.1 + 2.7*(0.1*Ret_eval**0.9))/(1. + (0.1*Ret_eval**0.9))

# Yin and Koch, 2007, data
Re_yk = [1.02, 1.83, 5.00, 9.98, 20.0]
k_yk = [4.2, 4.1, 4.0, 3.8, 3.6]

nFig = plt.figure(figsize=(2,2))
ax1 = nFig.add_subplot(111)
## n from RZ
ax1.plot(Ret_eval, nRZ_eval, 'k-', markersize=7)
ax1.plot(Ret_eval, nGAD_eval, 'k--', markersize=7)
ax1.plot(Ret, [n20, n33, n40, n50], 'k.', markersize=5)
ax1.plot(Re_yk, k_yk, 'bo', markersize=5, markerfacecolor='none')

ax1.set_xlim([0,130])
ax1.set_xticks([0,25,50,75,100,125])
ax1.set_xticklabels(["0","","50","","100",""])
ax1.set_xlabel(r"$Re_t$")
ax1.set_ylim([2.5,5])
ax1.set_ylabel(r"$n$")

#lText = [r'$n_{RZ}$', r'$n_{GA}$', r'$n_{WSP}$', r'$n_{YK}$']
#ax1.legend(lText, bbox_to_anchor=(0,1.05,1,1),loc="lower left",mode="expand",
#  ncol=2, borderaxespad=0)

ax1.grid(True)

imgname = imgdir + "n-relations"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

# Evaluate relationship with RZ, GAD
phiEval = np.linspace(0.05,0.375,100)

## Plot as a function of phi ##
fig1 = plt.figure(figsize=(4,4))

# linear coodinates
ax1 = fig1.add_subplot(111)
ax1.plot(phi, wfwt_rho20, 'b*', markersize=7, alpha=0.7)
ax1.plot(phi, wfwt_rho33, 'gs', markersize=7, alpha=0.7)
ax1.plot(phi, wfwt_rho40, 'ro', markersize=7, alpha=0.7)
ax1.plot(phi, wfwt_rho50, 'c^', markersize=7, alpha=0.7)

ax1.plot(phiEval, k20*(1-phiEval)**(n20-1), 'b--', zorder=1)
ax1.plot(phiEval, k33*(1-phiEval)**(n33-1), 'g--', zorder=1)
ax1.plot(phiEval, k40*(1-phiEval)**(n40-1), 'r--', zorder=1)
ax1.plot(phiEval, k50*(1-phiEval)**(n50-1), 'c--', zorder=1)

ax1.set_xlim([0, 0.4])
#ax1.set_xticks([0, 0.1, 0.20, 0.30, 0.4])
ax1.xaxis.set_major_locator(MultipleLocator(0.1))
ax1.xaxis.set_minor_locator(MultipleLocator(0.05))
ax1.set_xlabel(r"$\phi$")
ax1.set_ylim([0.2, 0.8])
ax1.yaxis.set_major_locator(MultipleLocator(0.2))
ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
ax1.set_ylabel(r'$w_f / w_t$', fontsize=14)

lText = [r'$\rho^* = 2.0$', r'$\rho^* = 3.3$', 
         r'$\rho^* = 4.0$', r'$\rho^* = 5.0$']
ax1.legend(lText, bbox_to_anchor=(0,1.05,1,1),loc="lower left",mode="expand",
   ncol=2, borderaxespad=0)

ax1.grid(True)

# Subplot -- bottom left
a = plt.axes([0.25, 0.2, 0.25, 0.25])
plt.plot(phi,wfwt_rho20, 'b*', markersize=7, alpha=0.7)
plt.plot(phi,wfwt_rho33, 'gs', markersize=7, alpha=0.7)
plt.plot(phi,wfwt_rho40, 'ro', markersize=7, alpha=0.7)
plt.plot(phi,wfwt_rho50, 'c^', markersize=7, alpha=0.7)

plt.plot(phiEval,k20*(1-phiEval)**(n20-1), 'b--', zorder=1)
plt.plot(phiEval,k33*(1-phiEval)**(n33-1), 'g--', zorder=1)
plt.plot(phiEval,k40*(1-phiEval)**(n40-1), 'r--', zorder=1)
plt.plot(phiEval,k50*(1-phiEval)**(n50-1), 'c--', zorder=1)

plt.xlim([0.07, 0.11])
plt.ylim([0.69, 0.73])
plt.xticks([0.08, 0.10])
plt.yticks([0.7, 0.72])

# Subplot -- top right
a = plt.axes([0.625, 0.625, 0.25, 0.25])
plt.plot(phi,wfwt_rho20, 'b*', markersize=7, alpha=0.7)
plt.plot(phi,wfwt_rho33, 'gs', markersize=7, alpha=0.7)
plt.plot(phi,wfwt_rho40, 'ro', markersize=7, alpha=0.7)
plt.plot(phi,wfwt_rho50, 'c^', markersize=7, alpha=0.7)

plt.plot(phiEval,k20*(1-phiEval)**(n20-1), 'b--', zorder=1)
plt.plot(phiEval,k33*(1-phiEval)**(n33-1), 'g--', zorder=1)
plt.plot(phiEval,k40*(1-phiEval)**(n40-1), 'r--', zorder=1)
plt.plot(phiEval,k50*(1-phiEval)**(n50-1), 'c--', zorder=1)

plt.xlim([0.155, 0.195])
plt.ylim([0.55, 0.59])
plt.xticks([0.16, 0.19])
plt.yticks([0.56, 0.58])


## log log coords
#nAx = fig1.add_subplot(121)
#nAx.loglog(1 - phi, wfwt_rho20, 'b*', markersize=7, alpha=0.7)
#nAx.loglog(1 - phi, wfwt_rho33, 'gs', markersize=7, alpha=0.7)
#nAx.loglog(1 - phi, wfwt_rho40, 'ro', markersize=7, alpha=0.7)
#nAx.loglog(1 - phi, wfwt_rho50, 'c^', markersize=7, alpha=0.7)
#
#nAx.loglog(1 - phiEval, k20*(1- phiEval)**(n20-1), 'b--')
#nAx.loglog(1 - phiEval, k33*(1- phiEval)**(n33-1), 'g--')
#nAx.loglog(1 - phiEval, k40*(1- phiEval)**(n40-1), 'r--')
#nAx.loglog(1 - phiEval, k50*(1- phiEval)**(n50-1), 'c--')
#
#nAx.set_xlabel(r'$1 - \phi$', fontsize=14)
#nAx.set_ylabel(r'$w_f / w_t$', fontsize=14)
#
#nAx.set_xlim([0.6, 1])
#nAx.set_xticks(np.arange(0.6,1.001, 0.05))
#nAx.set_xticklabels(["0.6","","0.7","","0.8","","0.9","","1.0"])
#nAx.set_ylim([.3, 0.8])
#nAx.set_yticks(np.arange(0.3,.801, 0.1))
#nAx.set_yticklabels(["0.3","","0.5","","","0.8"])

imgname = imgdir + "phaseAverageWf-phi"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
