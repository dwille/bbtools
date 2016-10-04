#!/usr/bin/env python2

from confidence import fit_param

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy.polynomial.polynomial as poly
import numpy as np
import matplotlib.lines as mlines
import os,sys
os.system('clear')

# get axis limits (so 24125058)
def get_axis_limits(ax, scalex=0.8, scaley=0.85):
  xmin = ax.get_xlim()[0]
  xmax = ax.get_xlim()[1]
  ymin = ax.get_ylim()[0]
  ymax = ax.get_ylim()[1]

  dx = xmax - xmin
  dy = ymax - ymin

  return xmin + scalex*dx, ymin + scaley*dy

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

# Fit a line to the log log
(n20, k20, n20_err, k20_err) = fit_param(np.log(1.-phi), np.log(wfwt_rho20))
(n33, k33, n33_err, k33_err) = fit_param(np.log(1.-phi), np.log(wfwt_rho33))
(n40, k40, n40_err, k40_err) = fit_param(np.log(1.-phi), np.log(wfwt_rho40))
(n50, k50, n50_err, k50_err) = fit_param(np.log(1.-phi), np.log(wfwt_rho50))

n20 += 1
n33 += 1
n40 += 1
n50 += 1

# Gather reults and erros from our study
k_curr = [k20, k33, k40, k50]
kerrs = [k20_err, k33_err, k40_err, k50_err]
n_curr = [n20, n33, n40, n50]
nerrs = [n20_err, n33_err, n40_err, n50_err]

print "Error on coefficients"
print "   2.0   3.3   4.0   5.0"
print "n  %.2f  %.2f  %.2f  %.2f" % (n20,n33,n40,n50)
print "ne %.2f  %.2f  %.2f  %.2f" % (n20_err,n33_err,n40_err,n50_err)
print "k  %.2f  %.2f  %.2f  %.2f" % (k20,k33,k40,k50)
print "ke %.2f  %.2f  %.2f  %.2f" % (k20_err,k33_err,k40_err,k50_err)
print ""
print "Error on wavespeed"
print "k = k + k', n = n + n', c = c + c'"
print "Thus (c + c')/c = (k+k')/k (n+n')/n (1-phi)^n' "
print "              2.0   3.3   4.0   5.0"
phi = np.array([0.087, 0.175, 0.262, 0.349])
p20 = (1. + n20_err/n20)*(1. + k20_err/k20)*(1-phi)**(+n20_err)
p33 = (1. + n33_err/n33)*(1. + k33_err/k33)*(1-phi)**(+n33_err)
p40 = (1. + n40_err/n40)*(1. + k40_err/k40)*(1-phi)**(+n40_err)
p50 = (1. + n50_err/n50)*(1. + k50_err/k50)*(1-phi)**(+n50_err)
print "phi = 0.087   %.3f  %.3f  %.3f  %.3f" % (p20[0], p20[1], p20[2], p20[3])
print "phi = 0.175   %.3f  %.3f  %.3f  %.3f" % (p33[0], p33[1], p33[2], p33[3])
print "phi = 0.262   %.3f  %.3f  %.3f  %.3f" % (p40[0], p40[1], p40[2], p40[3])
print "phi = 0.349   %.3f  %.3f  %.3f  %.3f" % (p50[0], p50[1], p50[2], p50[3])

# R-Z relations -- n coefficient
Ret = d*termVel/nu

##
## Plot k as a function of phi and rho ##
##
kFig = plt.figure(figsize=(3.25,1.625))

# Yin and Koch, 2007, data
Re_yk = [1.02, 1.83, 5.00, 9.98, 20.0]
k_yk = [0.92, 0.90, 0.86, 0.86, 0.88]
kerr_yk = [0.03, 0.03, 0.06, 0.04, 0.02]
n_yk = [4.2, 4.1, 4.0, 3.8, 3.6]
nerr_yk = [0.1, 0.1, 0.3, 0.2, 0.1]

## k as a function of rho 
ax2 = kFig.add_subplot(121)
ax2.errorbar(Ret, k_curr, fmt='k.', markersize=5, 
  yerr=kerrs)
ax2.errorbar(Re_yk, k_yk, fmt='ko', markersize=5, markerfacecolor='none',
  yerr=kerr_yk)

ax2.set_xlabel(r"$Re_t$")
ax2.xaxis.set_major_locator(MultipleLocator(50))
ax2.xaxis.set_minor_locator(MultipleLocator(25))
ax2.set_xlim([0,130])

ax2.set_ylabel(r"$\kappa$", rotation=0)
ax2.yaxis.set_major_locator(MultipleLocator(.05))
ax2.yaxis.set_minor_locator(MultipleLocator(.025))
ax2.set_ylim([0.75, 1])
ax2.yaxis.set_label_coords(-.40, 0.5)

#ax2.grid(True)
ax2.annotate(r"$(a)$", xy=get_axis_limits(ax2))

imgname = imgdir + "k-relations"
#plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## Plot n as a function of rho ##
Ret_eval = np.arange(0.1, 150, 0.1)
nRZ_eval = 4.45*Ret_eval**(-0.1)   
nGAD_eval = (5.1 + 2.7*(0.1*Ret_eval**0.9))/(1. + (0.1*Ret_eval**0.9))

#nFig = plt.figure(figsize=(2,2))
ax1 = kFig.add_subplot(122)
## n from RZ
ax1.plot(Ret_eval, nRZ_eval, 'k-', markersize=7)
ax1.plot(Ret_eval, nGAD_eval, 'k--', markersize=7)
ax1.errorbar(Ret, n_curr, fmt='k.', markersize=5, yerr=nerrs)
ax1.errorbar(Re_yk, n_yk, fmt='ko', markersize=5, markerfacecolor='none',
  yerr=nerr_yk)

ax1.set_xlabel(r"$Re_t$")
ax1.xaxis.set_major_locator(MultipleLocator(50))
ax1.xaxis.set_minor_locator(MultipleLocator(25))
ax1.set_xlim([0,130])

ax1.set_ylabel(r"$n$", rotation=0)
ax1.yaxis.set_major_locator(MultipleLocator(0.5))
ax1.yaxis.set_minor_locator(MultipleLocator(0.25))
ax1.set_ylim([2.5,5])
ax1.yaxis.tick_right()
ax1.yaxis.set_ticks_position("both")
#ax1.yaxis.set_label_position("right")
ax1.yaxis.set_label_coords(1.35, 0.5)

#ax1.grid(True)
ax1.annotate(r"$(b)$", color="black", xy=get_axis_limits(ax1))

imgname = imgdir + "nk-relations"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
plt.savefig(imgname + ".eps", bbox_inches='tight', format='eps')

##
## Evaluate relationship with RZ, GAD ##
##
phiEval = np.linspace(0.05,0.375,100)

## Plot as a function of phi ##
fig1 = plt.figure(figsize=(3.25,3.25))

# linear coodinates (colors from seaborn)
# dark
bl= "#4C72B0"
gr= "#55A868"
re= "#C44E52"
pu= "#8172B2"
go= "#CCB974"
cy= "#64B5CD"

## muted
#bl= "#4878CF"
#gr= "#6ACC65"
#re= "#D65F5F"
#pu= "#B47CC7"
#go= "#C4AD66"
#cy= "#77BEDB"

ax1 = fig1.add_subplot(111)
ax1.plot(phi, wfwt_rho20, '*', color=bl, markersize=7)#, alpha=0.7)
ax1.plot(phi, wfwt_rho33, 's', color=gr, markersize=7)#, alpha=0.7)
ax1.plot(phi, wfwt_rho40, 'o', color=re, markersize=7)#, alpha=0.7)
ax1.plot(phi, wfwt_rho50, '^', color=cy, markersize=7)#, alpha=0.7)

ax1.plot(phiEval, k20*(1-phiEval)**(n20-1), ':', color=bl, zorder=1, 
  linewidth=2, dashes=[2, 2])
ax1.plot(phiEval, k33*(1-phiEval)**(n33-1),'-.', color=gr, zorder=1, 
  linewidth=2, dashes=[6, 2, 2, 2])
ax1.plot(phiEval, k40*(1-phiEval)**(n40-1), '-', color=re, zorder=1, linewidth=2)
ax1.plot(phiEval, k50*(1-phiEval)**(n50-1),'--', color=cy, zorder=1, linewidth=2)

ax1.set_xlim([0, 0.4])
#ax1.set_xticks([0, 0.1, 0.20, 0.30, 0.4])
ax1.xaxis.set_major_locator(MultipleLocator(0.1))
ax1.xaxis.set_minor_locator(MultipleLocator(0.05))
ax1.set_xlabel(r"$\phi$")
ax1.set_ylim([0.2, 0.8])
ax1.yaxis.set_major_locator(MultipleLocator(0.2))
ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
ax1.set_ylabel(r'$\frac{\langle w_f \rangle - \langle w_p \rangle}{w_t}$',
  fontsize=14)
ax1.yaxis.set_label_coords(-0.12, 0.5)

lText = [r'$\rho^* = 2.0$', r'$\rho^* = 3.3$', 
         r'$\rho^* = 4.0$', r'$\rho^* = 5.0$']

#h1 = mlines.Line2D([],[], linestyle=':', color=bl, marker='*', label=lText[0])
#h2 = mlines.Line2D([],[], linestyle='-.', color=gr, marker='s', label=lText[1])
#h3 = mlines.Line2D([],[], linestyle='-', color=re, marker='o', label=lText[2])
#h4 = mlines.Line2D([],[], linestyle='--', color=cy, marker='^', label=lText[3])
#ax1.legend(handles=[h1,h2,h3,h4], bbox_to_anchor=(0,1.05,1,1), loc="lower left",
#  mode="expand", ncol=2, borderaxespad=0)

ax1.grid(True)

# Subplot -- bottom left
a = plt.axes([0.25, 0.2, 0.25, 0.25])
plt.plot(phi,wfwt_rho20, '*', color=bl, markersize=7, alpha=0.7)
plt.plot(phi,wfwt_rho33, 's', color=gr, markersize=7, alpha=0.7)
plt.plot(phi,wfwt_rho40, 'o', color=re, markersize=7, alpha=0.7)
plt.plot(phi,wfwt_rho50, '^', color=cy, markersize=7, alpha=0.7)

plt.plot(phiEval,k20*(1-phiEval)**(n20-1), ':', color=bl, zorder=1,
  linewidth=2, dashes=[2, 2])
plt.plot(phiEval,k33*(1-phiEval)**(n33-1),'-.', color=gr, zorder=1,
  linewidth=2, dashes=[6, 2, 2, 2])
plt.plot(phiEval,k40*(1-phiEval)**(n40-1), '-', color=re, zorder=1, linewidth=2)
plt.plot(phiEval,k50*(1-phiEval)**(n50-1),'--', color=cy, zorder=1, linewidth=2)

ax1.arrow(0.09, 0.69, 0.033, -0.185, head_width=0.01, head_length=0.03, 
  fc='k', ec='k')

plt.xlim([0.07, 0.11])
plt.ylim([0.69, 0.73])
plt.xticks([0.08, 0.10])
plt.yticks([0.7, 0.72])

# Subplot -- top right
a = plt.axes([0.625, 0.625, 0.25, 0.25])
plt.plot(phi,wfwt_rho20, '*', color=bl, markersize=7, alpha=0.7)
plt.plot(phi,wfwt_rho33, 's', color=gr, markersize=7, alpha=0.7)
plt.plot(phi,wfwt_rho40, 'o', color=re, markersize=7, alpha=0.7)
plt.plot(phi,wfwt_rho50, '^', color=cy, markersize=7, alpha=0.7)

plt.plot(phiEval,k20*(1-phiEval)**(n20-1), ':',color=bl, zorder=1,
  linewidth=2, dashes=[2, 2])
plt.plot(phiEval,k33*(1-phiEval)**(n33-1),'-.',color=gr, zorder=1,
  linewidth=2, dashes=[6, 2, 2, 2])
plt.plot(phiEval,k40*(1-phiEval)**(n40-1), '-',color=re, zorder=1, linewidth=2)
plt.plot(phiEval,k50*(1-phiEval)**(n50-1),'--',color=cy, zorder=1, linewidth=2)

ax1.arrow(0.1825, 0.5825, 0.0425, 0.0125, head_width=0.01, head_length=0.02,
  fc='k', ec='k')

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

imgname = imgdir + "rz-comparison"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
plt.savefig(imgname + ".eps", bbox_inches='tight', format='eps')
