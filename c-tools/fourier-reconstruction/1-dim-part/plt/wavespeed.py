#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import os, sys
import matplotlib.lines as mlines
os.system('clear')

import matplotlib.ticker as tick
from matplotlib.ticker import MultipleLocator

print("\n ---- Wavespeed Plotting Utility ---- \n")

### Physical Parameters ### XXX should not be hardcoded
d = 4.2                                                 # mm
nu = .01715                                             # mm^2/ms
g = 0.00981                                             # mm/ms^2
rho_f = 8.75e-4                                         # g/mm^3
rho_p = np.array([0.001750, 0.00289, 0.0035, 0.004375]) # g/mm^3
Lx = 42                                                 # mm
Ly = 42                                                 # mm
Lz = 126                                                # mm

### Hardcoded directories ###
home = os.path.expanduser("~")
root = home + "/scratch/triply_per/simdata/"
imgdir = root + "img/f-rec-1D/"

phaseFile = root + "phaseAveragedFluidVel"
waveFile = root + "wave-data-with-sdev"
termFile = root + "singlePartSedi"
rzFile = root + "richardson_zaki_coeffs"

print("      Using location: " + waveFile)
print("")

if not os.path.exists(waveFile):
  print("      " + waveFile + " does not exist. Exiting...\n")
  sys.exit()
elif not os.path.exists(termFile):
  print("      " + termFile + " does not exist. Exiting...\n")
  sys.exit()
elif not os.path.exists(phaseFile):
  print("      " + phaseFile + " does not exist. Exiting...\n")
  sys.exit()
 
### Read and parse data ###
# wave_data -- contains mean and sdev wavespeeds
# flow_vel_data -- contains phase averaged fluid velocityh
# term_vel -- contains single particle terminal velocity
# n, kappa -- richardson-zaki coefficients
wave_data = np.genfromtxt(waveFile, skip_header=9) # wavespeed in mm/s
flow_vel_data = np.genfromtxt(phaseFile, skip_header=1) # mm/ms

nparts = np.unique(wave_data[:,0]).astype(int)
rho = np.unique(wave_data[:,1])
phi = 4./3.*np.pi*(0.5*d)**3.*nparts/(Lx * Ly * Lz)

wave_speed = np.zeros((rho.size, nparts.size))
wave_speed_sdev = np.zeros((rho.size, nparts.size))
phase_averaged_fluid_vel = np.zeros((rho.size, nparts.size))
for cc in np.arange(np.size(wave_data, 0)):
  curr_n = int(wave_data[cc, 0])
  curr_rho = wave_data[cc, 1]

  # Find indices to fill arrays
  pp = np.argwhere(curr_rho == rho)
  nn = np.argwhere(curr_n == nparts)

  wave_speed[pp,nn] = wave_data[cc,2]
  wave_speed_sdev[pp,nn] = wave_data[cc,3]
  phase_averaged_fluid_vel[pp,nn] = flow_vel_data[cc,2]

term_vel = np.genfromtxt(termFile, skip_header=2, usecols=1) # mm/ms
n = np.genfromtxt(rzFile, skip_header=1, usecols=1)
kappa = np.genfromtxt(rzFile, skip_header=1, usecols=2)

## Convert units to m/s ##
wave_speed /= 1000              ## mm/s -> m/s
wave_speed_sdev /= 1000         ## mm/s -> m/s
phase_averaged_fluid_vel *= 1.  ## mm/ms -> m/s
term_vel *= 1.                  ## mm/ms -> m/s

## Reynolds prefactor
Rep = d/nu

## Wallis relations ##
phi_eval = np.linspace(0.05,0.45,101)
wallis_speed = np.zeros((np.size(phi_eval),4))
for rr,_ in enumerate(rho):
  for pp,_ in enumerate(phi_eval):
    wallis_speed[pp,rr] = term_vel[rr]*kappa[rr]*n[rr]*phi_eval[pp]*(1-phi_eval[pp])**(n[rr] - 1)

######################
### plotting #########
######################
size = (3.25, 3.25)
colors = ["#4C72B0", "#55A868", "#C44E52", "#64B5CD"] # blue, green, red, cyan
markers = ['*', 's', 'o', '^']
lines = [':', '-.', '-', '--']

### c(phi), normalized by Rep = d/nu ###
fig1 = plt.figure(figsize=size)
ax1 = fig1.add_subplot(111)

for i in np.arange(rho.size):
  #ax1.plot(phi, Rep*wave_speed[i,:], markers[i], color=colors[i], markersize=7,
  #  markeredgecolor='k', markeredgewidth=0.5)
  ax1.errorbar(phi, Rep*wave_speed[i,:], yerr=(Rep*wave_speed_sdev[i,:]),
    marker=markers[i], linestyle='', color=colors[i], markersize=6,
    markeredgecolor='k', markeredgewidth=0.5,
    ecolor='k', capthick=2, capsize=2, elinewidth=1)

ax1.set_xlabel(r'$\phi$')
ax1.set_ylabel(r'$2ac/\nu$')
ax1.yaxis.set_label_coords(-0.12, 0.5)
ax1.set_ylim([0,45])
ax1.set_xlim([0,0.5])
ax1.grid(True)

ax1.xaxis.set_major_locator(MultipleLocator(0.1))
ax1.xaxis.set_minor_locator(MultipleLocator(0.05))
ax1.yaxis.set_major_locator(MultipleLocator(10))
ax1.yaxis.set_minor_locator(MultipleLocator(5))

# Wallis lines
for i in np.arange(rho.size):
  ax1.plot(phi_eval, Rep*wallis_speed[:,i], lines[i], color=colors[i], zorder=1,
    linewidth=2)

## XXX HACK FOR ERRORS
#s_wt = np.array([0.003, 0.007, 0.009, 0.011])
#s_n = np.array([0.15, 0.11, 0.11, 0.08])
#s_k = np.array([0.04, 0.03, 0.03, 0.02])
#theoretical_error = np.zeros((np.size(phi_eval),4))
#for rr,_ in enumerate(rho):
#  for pp,_ in enumerate(phi_eval):
#    theoretical_error[pp,rr] = 0 + \
#      (kappa[rr] * n[rr] * phi_eval[pp] * (1 - phi_eval[pp])**(n[rr] - 1))**2 * s_wt[rr]**2 + \
#      (term_vel[rr] * kappa[rr] * phi_eval[pp] * (1 - phi_eval[pp])**(n[rr] - 1)* \
#        (1 + n[rr] * (n[rr] - 1) * (1 - phi_eval[pp])**(-1)))**2 * s_n[rr]**2 + \
#      (term_vel[rr] * n[rr] * phi_eval[pp] * (1 - phi_eval[pp])**(n[rr] - 1))**2 * s_k[rr]**2
#    theoretical_error[pp,rr] = np.sqrt(theoretical_error[pp,rr])
#for i in np.arange(rho.size):
#  ax1.fill_between(phi_eval, Rep * (wallis_speed[:,i] + theoretical_error[:,i]),
#    Rep * (wallis_speed[:,i] - theoretical_error[:,i]), color=colors[i], alpha=0.25)


    


plt.savefig(imgdir + "wavespeed.png", bbox_inches='tight', format='png')
plt.savefig(imgdir + "wavespeed.eps", bbox_inches='tight', format='eps')
sys.exit()

### HERE BELOW IS OLD, UNMODIFIED. SOME MAY BE USEFUL BUT REQUIRE RENAMING VARIABLES,
### FOLLOWING THE ABOVE EXAMPLE


# wavespeed -- constant density, variable vfrac
rho20_c = waveSpeed[0:16:4] # / termVel[0]
rho33_c = waveSpeed[1:16:4] # / termVel[1]
rho40_c = waveSpeed[2:12:4] # / termVel[2]
rho50_c = waveSpeed[3:8:4]  # / termVel[3]

# fluid velocity -- constant density, variable vfrac
rho20_vel = flowVel[0:16:4]
rho33_vel = flowVel[1:16:4]
rho40_vel = flowVel[2:12:4]
rho50_vel = flowVel[3:8:4]


## Wave Speed Plot -- normed by 2a/nu
fig1 = plt.figure(figsize=(3.25,3.25))
ax1 = fig1.add_subplot(111)

ax1.plot(phi, Rep*rho20_c, '*', color=bl, markersize=7)
ax1.plot(phi, Rep*rho33_c, 's', color=gr, markersize=7)
ax1.plot(phi[0:3], Rep*rho40_c, 'o', color=re, markersize=7)
ax1.plot(phi[0:2], Rep*rho50_c, '^', color=cy, markersize=7)

ax1.set_xlabel(r'$\phi$')
#ax1.set_ylabel(r'$\frac{2ac}{\nu}$')
ax1.set_ylabel(r'$2ac/\nu$')
ax1.yaxis.set_label_coords(-0.12, 0.5)
ax1.set_ylim([0,45])
ax1.set_xlim([0,0.5])
ax1.grid(True)

ax1.xaxis.set_major_locator(MultipleLocator(0.1))
ax1.xaxis.set_minor_locator(MultipleLocator(0.05))
ax1.yaxis.set_major_locator(MultipleLocator(10))
ax1.yaxis.set_minor_locator(MultipleLocator(5))

lText = [r'$\rho^* = 2.0$', r'$\rho^* = 3.3$', 
         r'$\rho^* = 4.0$', r'$\rho^* = 5.0$']
h1 = mlines.Line2D([],[], linestyle=':', color=bl, marker='*', label=lText[0])
h2 = mlines.Line2D([],[], linestyle='-.', color=gr, marker='s', label=lText[1])
h3 = mlines.Line2D([],[], linestyle='-', color=re, marker='o', label=lText[2])
h4 = mlines.Line2D([],[], linestyle='--', color=cy, marker='^', label=lText[3])
#ax1.legend(handles=[h1,h2,h3,h4], bbox_to_anchor=(0,1.05,1,1), loc="lower left",
# mode="expand", ncol=2, borderaxespad=0)

# Wallis correlations
ax1.plot(phiEval, Rep*wallis_speed[:,0], ':', color=bl, zorder=1,
  linewidth=2, dashes=[2, 2])
ax1.plot(phiEval, Rep*wallis_speed[:,1], '-.', color=gr, zorder=1,
  linewidth=2, dashes=[6, 2, 2, 2])
ax1.plot(phiEval, Rep*wallis_speed[:,2], '-', color=re, zorder=1, linewidth=2)
ax1.plot(phiEval, Rep*wallis_speed[:,3], '--', color=cy, zorder=1, linewidth=2)

imgname = imgdir + "wavespeed_norm_2a_nu"
#imgname = imgdir + "wavespeed_normed"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
plt.savefig(imgname + ".eps", bbox_inches='tight', format='eps')

sys.exit("Exiting...")



## Wave Speed plot -- normed by w_t

## Wallis relations ##
phi_eval = np.linspace(0.05,0.45,101)
wallis_speed_norm_wt = np.zeros((np.size(phi_eval),4))
for rr,_ in enumerate(rho):
  for pp,_ in enumerate(phi_eval):
    wallis_speed_norm_wt[pp,rr] = kappa[rr]*n[rr]*phi_eval[pp]*(1-phi_eval[pp])**(n[rr] - 1)


fig2 = plt.figure(figsize=(3.25,3.25))
ax1 = fig2.add_subplot(111)

ax1.plot(phi, rho20_c/termVel[0], '*', color=bl, markersize=7)
ax1.plot(phi, rho33_c/termVel[1], 's', color=gr, markersize=7)
ax1.plot(phi[0:3], rho40_c/termVel[2], 'o', color=re, markersize=7)
ax1.plot(phi[0:2], rho50_c/termVel[3], '^', color=cy, markersize=7)

ax1.set_xlabel(r'$\phi$')
ax1.set_xlim([0,0.5])

ax1.set_ylabel(r'$c/w_t$')
ax1.yaxis.set_label_coords(-0.16, 0.5)
ax1.set_ylim([0,.45])
ax1.grid(True)

ax1.xaxis.set_major_locator(MultipleLocator(0.1))
ax1.xaxis.set_minor_locator(MultipleLocator(0.05))
ax1.yaxis.set_major_locator(MultipleLocator(.10))
ax1.yaxis.set_minor_locator(MultipleLocator(.05))

lText = [r'$\rho^* = 2.0$', r'$\rho^* = 3.3$', 
         r'$\rho^* = 4.0$', r'$\rho^* = 5.0$']
h1 = mlines.Line2D([],[], linestyle=':', color=bl, marker='*', label=lText[0])
h2 = mlines.Line2D([],[], linestyle='-.', color=gr, marker='s', label=lText[1])
h3 = mlines.Line2D([],[], linestyle='-', color=re, marker='o', label=lText[2])
h4 = mlines.Line2D([],[], linestyle='--', color=cy, marker='^', label=lText[3])
#ax1.legend(handles=[h1,h2,h3,h4], bbox_to_anchor=(0,1.05,1,1), loc="lower left",
# mode="expand", ncol=2, borderaxespad=0)

# Wallis correlations
ax1.plot(phiEval, wallis_speed_norm_wt[:,0], ':', color=bl, zorder=1,
  linewidth=2, dashes=[2, 2])
ax1.plot(phiEval, wallis_speed_norm_wt[:,1], '-.', color=gr, zorder=1,
  linewidth=2, dashes=[6, 2, 2, 2])
ax1.plot(phiEval, wallis_speed_norm_wt[:,2], '-', color=re, zorder=1, 
  linewidth=2)
ax1.plot(phiEval, wallis_speed_norm_wt[:,3], '--', color=cy, zorder=1, 
  linewidth=2)

imgname = imgdir + "wavespeed_norm_wt"
#imgname = imgdir + "wavespeed_normed"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
plt.savefig(imgname + ".eps", bbox_inches='tight', format='eps')


## Print Errors
rho20_c = waveSpeed[0:16:4] # / termVel[0]
rho33_c = waveSpeed[1:16:4] # / termVel[1]
rho40_c = waveSpeed[2:12:4] # / termVel[2]
rho50_c = waveSpeed[3:8:4]  # / termVel[3]
error = np.zeros((4,4))
for rr,_ in enumerate(rho):
  for pp,_ in enumerate(phi):
    wallis_speed[pp,rr] = termVel[rr]*kappa[rr]*n[rr]*phi[pp]*(1-phi[pp])**(n[rr] - 1)
    #wallis_speed[pp,rr] = kappa[rr]*n[rr]*phi[pp]*(1-phi[pp])**(n[rr] - 1)

    error[pp,0] = ((rho20_c[pp] - wallis_speed[pp,0]) / wallis_speed[pp,0])
    error[pp,1] = ((rho33_c[pp] - wallis_speed[pp,1]) / wallis_speed[pp,1])
    if pp < 3:
      error[pp,2] = ((rho40_c[pp] - wallis_speed[pp,2]) / wallis_speed[pp,2])
    if pp < 2:
      error[pp,3] = ((rho50_c[pp] - wallis_speed[pp,3]) / wallis_speed[pp,3])
                                                      
# 0500 #
print("n = 500")
print("  Error rho2.0 = %.3f" % ((rho20_c[0] - wallis_speed[0,0]) / wallis_speed[0,0]))
print("  Error rho3.3 = %.3f" % ((rho33_c[0] - wallis_speed[0,1]) / wallis_speed[0,1]))
print("  Error rho4.0 = %.3f" % ((rho40_c[0] - wallis_speed[0,2]) / wallis_speed[0,2]))
print("  Error rho5.0 = %.3f" % ((rho50_c[0] - wallis_speed[0,3]) / wallis_speed[0,3]))

# 1000 #
print("n = 1000")
print("  Error rho2.0 = %.3f" % ((rho20_c[1] - wallis_speed[1,0]) / wallis_speed[1,0]))
print("  Error rho3.3 = %.3f" % ((rho33_c[1] - wallis_speed[1,1]) / wallis_speed[1,1]))
print("  Error rho4.0 = %.3f" % ((rho40_c[1] - wallis_speed[1,2]) / wallis_speed[1,2]))
print("  Error rho5.0 = %.3f" % ((rho50_c[1] - wallis_speed[1,3]) / wallis_speed[1,3]))

# 1500 #
print("n = 1500")
print("  Error rho2.0 = %.3f" % ((rho20_c[2] - wallis_speed[2,0]) / wallis_speed[2,0]))
print("  Error rho3.3 = %.3f" % ((rho33_c[2] - wallis_speed[2,1]) / wallis_speed[2,1]))
print("  Error rho4.0 = %.3f" % ((rho40_c[2] - wallis_speed[2,2]) / wallis_speed[2,2]))

# 2000 #
print("n = 2000")
print("  Error rho2.0 = %.3f" % ((rho20_c[3] - wallis_speed[3,0]) / wallis_speed[3,0]))
print("  Error rho3.3 = %.3f" % ((rho33_c[3] - wallis_speed[3,1]) / wallis_speed[3,1]))

figErr = plt.figure(figsize=(3.25,1.625))
ax1 = figErr.add_subplot(111)
ax1.plot(phi, error[:,0])
ax1.plot(phi, error[:,1])
ax1.plot(phi[0:3], error[0:3,2])
ax1.plot(phi[0:2], error[0:2,3])

imgname = imgdir + "wavespeed_error"
#imgname = imgdir + "wavespeed_normed"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

