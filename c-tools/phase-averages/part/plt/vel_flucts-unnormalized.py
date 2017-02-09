#!/usr/bin/env python2

import matplotlib
matplotlib.use(u'Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np
import os,sys
os.system('clear')


print ""
print " ---- Phase-Averaged Particle Velocity Flucts Plotting Utility"
print ""

# Simulation Params
d = 4.2
nu = .01715
rhof = 8.75e-4
partR = 2.1
g = 0.00981
rhopf = np.array([2.0, 3.3, 4.0, 5.0])
rhop = np.array([0.001750, .00289, 0.0035, 0.004375])
phi = np.array([0.087, 0.175, 0.262, 0.349])

# MARCC
home = os.path.expanduser("~")
root = home + "/scratch/triply_per/simdata/"
fluctFile = root + "velFlucts_raw"
termFile = root + "singlePartSedi"

print "      Root directory set to: " + root
print "      Using location: " + fluctFile
print "      Using location: " + termFile
print ""
print "      Using d = %.2f" % d
print "      Using nu = %.5f" % nu

# Check if datadir exists so we don't go creating extra dirs
if not os.path.exists(fluctFile):
  print "      " + fluctFile + " does not exist. Exiting..."
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

imgdir = imgdir + "velocity_flucts/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# Read data
uflucts = np.genfromtxt(fluctFile, skip_header=1, usecols=2)
vflucts = np.genfromtxt(fluctFile, skip_header=1, usecols=3)
wflucts = np.genfromtxt(fluctFile, skip_header=1, usecols=4)
wrel = np.genfromtxt(fluctFile, skip_header=1, usecols=7)
termVel = np.genfromtxt(termFile, skip_header=1, usecols=1)
dwdphi = np.genfromtxt(fluctFile, skip_header=1, usecols=8)

Re_t = d*termVel/nu

# Normalize
#wflucts[0:4] /= termVel
#wflucts[4:8] /= termVel
#wflucts[8:12] /= termVel
#wflucts[12:16] /= termVel
#wflucts /= wrel

## Plot as a function of phi
fig1 = plt.figure(figsize=(6,1.5))

ax1 = fig1.add_subplot(121)
ax1.plot(phi, wflucts[0:16:4], 'b*', markersize=7, alpha=0.7)
ax1.plot(phi, wflucts[1:16:4], 'gs', markersize=7, alpha=0.7)
ax1.plot(phi, wflucts[2:16:4], 'ro', markersize=7, alpha=0.7)
ax1.plot(phi, wflucts[3:16:4], 'c^', markersize=7, alpha=0.7)


ax1.set_xlabel(r'$\phi$')
ax1.set_xlim([0,0.45])
ax1.set_xticks([0, 0.15, 0.30, 0.45])
#ax1.set_ylim([0.1, 0.25])
#ax1.set_yticks([0.1, 0.15, 0.2, 0.25])
ax1.set_ylabel(r'$w_p^\prime$', fontsize=14, rotation=0)
ax1.yaxis.set_label_coords(-0.25,0.5)

ax1.grid(True)

lText = [r'$\rho^* = 2.0$', r'$\rho^* = 3.3$', 
         r'$\rho^* = 4.0$', r'$\rho^* = 5.0$']
ax1.legend(lText, bbox_to_anchor=(0,1.05,1,1),loc="lower left",mode="expand",
  ncol=2, borderaxespad=0)


alphas = np.linspace(0.4, .9, 4)
rgba_colors_b = np.zeros((4,4))
rgba_colors_b[:,2] = 1.
rgba_colors_b[:,3] = alphas
rgba_colors_g = np.zeros((4,4))
rgba_colors_g[:,1] = 1.
rgba_colors_g[:,3] = alphas
rgba_colors_r = np.zeros((4,4))
rgba_colors_r[:,0] = 1.
rgba_colors_r[:,3] = alphas
rgba_colors_c = np.zeros((4,4))
rgba_colors_c[:,1] = .75
rgba_colors_c[:,2] = 1.
rgba_colors_c[:,3] = alphas
ax2 = fig1.add_subplot(122)

ax2.scatter(Re_t[0]*np.ones(4), wflucts[0:16:4], marker='*', color=rgba_colors_b)
ax2.scatter(Re_t[1]*np.ones(4), wflucts[1:16:4], marker='s', color=rgba_colors_g)
ax2.scatter(Re_t[2]*np.ones(4), wflucts[2:16:4], marker='o', color=rgba_colors_r)
ax2.scatter(Re_t[3]*np.ones(4), wflucts[3:16:4], marker='^', color=rgba_colors_c)
#ax2.plot(Re_t, wflucts[0:4], 'b*', markersize=7, alpha=0.7)
#ax2.plot(Re_t, wflucts[4:8], 'gs', markersize=7, alpha=0.7)
#ax2.plot(Re_t, wflucts[8:12], 'ro', markersize=7, alpha=0.7)
#ax2.plot(Re_t, wflucts[12:16], 'c^', markersize=7, alpha=0.7)

ax2.set_xlabel(r'$Re_t$')
ax2.set_xlim([0,130])
ax2.set_xticks([0,25,50,75,100,125])
#ax2.set_ylim([0.1, 0.25])
#ax2.set_yticks([0.1, 0.15, 0.2, 0.25])
ax2.set_yticklabels([])

ax2.grid(True)

lText = [r'$N = 500$', r'$N = 1000$', 
         r'$N = 1500$', r'$N = 2000$']
ax2.legend(lText, bbox_to_anchor=(0,1.05,1,1),loc="lower left",mode="expand",
  ncol=2, borderaxespad=0)


imgname = imgdir + "velflucts_phi_none"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## Plot as a function of rho (or Re_t)
fig2 = plt.figure(figsize=(4,3))
ax1 = fig2.add_subplot(211)
ax1.plot(phi, dwdphi[0:16:4], 'b*', markersize=7, alpha=0.7)
ax1.plot(phi, dwdphi[1:16:4], 'gs', markersize=7, alpha=0.7)
ax1.plot(phi, dwdphi[2:16:4], 'ro', markersize=7, alpha=0.7)
ax1.plot(phi, dwdphi[3:16:4], 'c^', markersize=7, alpha=0.7)

#ax1.set_ylim([0, 2])
#ax1.set_ylabel(r"$\frac{1}{w_t}\frac{dw}{d\phi}$", rotation=0)
ax1.set_ylabel(r"$\frac{dw'}{d\phi}$", rotation=0)
ax1.set_xticks([0, 0.15, 0.30, 0.45])
ax1.set_xlabel(r'$\phi$', fontsize=14)


ax1.grid(True)

ax2 = fig2.add_subplot(212)
ax2.plot(Re_t, dwdphi[0:4], 'b*', markersize=7, alpha=0.7)
ax2.plot(Re_t, dwdphi[4:8], 'gs', markersize=7, alpha=0.7)
ax2.plot(Re_t, dwdphi[8:12], 'ro', markersize=7, alpha=0.7)
ax2.plot(Re_t, dwdphi[12:16], 'c^', markersize=7, alpha=0.7)

#ax2.set_ylim([0, 2])
#ax2.set_ylabel(r"$\frac{1}{w_t}\frac{dw}{d\phi}$", rotation=0)
ax2.set_ylabel(r"$\frac{dw'}{d\phi}$", rotation=0)
#ax2.set_xticks([0, 0.15, 0.30, 0.45])
ax2.set_xlabel(r'$Re_t$', fontsize=14)


ax2.grid(True)

imgname = imgdir + "velflucts_rho"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
