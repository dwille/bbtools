#!/usr/bin/env python2

import sys, os
import glob
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import lines as mlines
import numpy as np
import re

## Get info
print ""
print " ---- Mean Measures plotting utility ---- "
print ""

# DEVEL
#root = "../sim"
#datadir = root + "data-tetrads/"

# MARCC
root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
print "      Sim root directory set to: " + root
simdir = raw_input("      Simulation directory: ")
if not simdir.endswith('/'):
  simdir = simdir + '/'
datadir = root + simdir + "data-tetrads/"

# Create imgdir if necessary
imgdir = root + simdir + "/img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

statMean = datadir + "stat.mean"
statSdev = datadir + "stat.sdev"
statSkew = datadir + "stat.skew"
statKurt = datadir + "stat.kurt"

# Find means and time from stat.dat
time = np.genfromtxt(statMean, skip_header=1, usecols = 0)
time = time - time[0]

meanR= np.genfromtxt(statMean, skip_header=1, usecols = 1)
meanEVar= np.genfromtxt(statMean, skip_header=1, usecols = 2)
meanShape = np.genfromtxt(statMean, skip_header=1, usecols = 3)
meanI1 = np.genfromtxt(statMean, skip_header=1, usecols = 4)
meanI2 = np.genfromtxt(statMean, skip_header=1, usecols = 5)
meanI3 = np.genfromtxt(statMean, skip_header=1, usecols = 6)
meanS11 = np.genfromtxt(statMean, skip_header=1, usecols = 7)
meanS22 = np.genfromtxt(statMean, skip_header=1, usecols = 8)
meanS33 = np.genfromtxt(statMean, skip_header=1, usecols = 9)

sdevR= np.genfromtxt(statSdev, skip_header=1, usecols = 1)
sdevEVar= np.genfromtxt(statSdev, skip_header=1, usecols = 2)
sdevShape = np.genfromtxt(statSdev, skip_header=1, usecols = 3)
sdevI1 = np.genfromtxt(statSdev, skip_header=1, usecols = 4)
sdevI2 = np.genfromtxt(statSdev, skip_header=1, usecols = 5)
sdevI3 = np.genfromtxt(statSdev, skip_header=1, usecols = 6)
sdevS11 = np.genfromtxt(statSdev, skip_header=1, usecols = 7)
sdevS22 = np.genfromtxt(statSdev, skip_header=1, usecols = 8)
sdevS33 = np.genfromtxt(statSdev, skip_header=1, usecols = 9)

skewR= np.genfromtxt(statSkew, skip_header=1, usecols = 1)
skewEVar= np.genfromtxt(statSkew, skip_header=1, usecols = 2)
skewShape = np.genfromtxt(statSkew, skip_header=1, usecols = 3)
skewI1 = np.genfromtxt(statSkew, skip_header=1, usecols = 4)
skewI2 = np.genfromtxt(statSkew, skip_header=1, usecols = 5)
skewI3 = np.genfromtxt(statSkew, skip_header=1, usecols = 6)
skewS11 = np.genfromtxt(statSkew, skip_header=1, usecols = 7)
skewS22 = np.genfromtxt(statSkew, skip_header=1, usecols = 8)
skewS33 = np.genfromtxt(statSkew, skip_header=1, usecols = 9)

kurtR= np.genfromtxt(statKurt, skip_header=1, usecols = 1)
kurtEVar= np.genfromtxt(statKurt, skip_header=1, usecols = 2)
kurtShape = np.genfromtxt(statKurt, skip_header=1, usecols = 3)
kurtI1 = np.genfromtxt(statKurt, skip_header=1, usecols = 4)
kurtI2 = np.genfromtxt(statKurt, skip_header=1, usecols = 5)
kurtI3 = np.genfromtxt(statKurt, skip_header=1, usecols = 6)
kurtS11 = np.genfromtxt(statKurt, skip_header=1, usecols = 7)
kurtS22 = np.genfromtxt(statKurt, skip_header=1, usecols = 8)
kurtS33 = np.genfromtxt(statKurt, skip_header=1, usecols = 9)

# Plot specs
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', labelsize=11)
plt.rc('figure', titlesize=14)
plt.rc('figure', figsize=(4,3))
plt.rc('legend', fontsize=11, numpoints=3)
plt.rc('lines', markersize=4)
labelx = -0.17

## SHAPE MEASURES ##
measures = plt.figure()

# Radius of Gyration
rg_ax = measures.add_subplot(311)
rg_ax.plot(time, meanR, 'ko', markevery=25)
rg_ax.plot(time, meanR + sdevR, 'k.', markevery=15, alpha=0.5)
rg_ax.plot(time, meanR - sdevR, 'k.', markevery=15, alpha=0.5)

rg_ax.set_ylabel(r"$\langle R \rangle$", rotation=0)
rg_ax.yaxis.set_label_coords(labelx, 0.5)
rg_ax.set_xlim([time[0], time[-1]])
rg_ax.tick_params(axis='x', labelbottom='off')

# Variance
var_ax = measures.add_subplot(312)
var_ax.plot(time, meanEVar, 'ko', markevery=25)
var_ax.plot(time, meanEVar + sdevEVar, 'k.', markevery=15, alpha=0.5)
var_ax.plot(time, meanEVar - sdevEVar, 'k.', markevery=15, alpha=0.5)
 
var_ax.set_ylabel(r"$\langle \Delta \rangle$", rotation=0)
var_ax.yaxis.set_label_coords(labelx, 0.5)
var_ax.set_xlim([time[0], time[-1]])
var_ax.set_ylim([0,1])
var_ax.set_yticks([0.00, 0.25, 0.50, 0.75, 1.00])
var_ax.tick_params(axis='x', labelbottom='off')

# Shape Factor
sf_ax = measures.add_subplot(313)
sf_ax.plot(time, meanShape, 'ko', markevery=25)
sf_ax.plot(time, meanShape + sdevShape, 'k.', markevery=15, alpha=0.5)
sf_ax.plot(time, meanShape - sdevShape, 'k.', markevery=15, alpha=0.5)

sf_ax.set_xlabel("Time [ms]")
sf_ax.set_ylabel(r"$\langle S \rangle$", rotation=0)
sf_ax.yaxis.set_label_coords(labelx, 0.5)
sf_ax.set_xlim([time[0], time[-1]])
sf_ax.set_ylim([-0.25,2])
sf_ax.set_yticks([-0.25, 0., 0.5, 1.0, 1.5, 2.0])

imgname = imgdir + "shape_mean_anisotropy"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## HIGHER ORDER MOMENTS ##
higher_measures = plt.figure()

# Radius of gyration
rg_ax2 = higher_measures.add_subplot(311)
rg_ax2.plot(time, skewR, 'ko', markevery=15)
rg_ax2.plot(time, kurtR, 'bo', markevery=15)

rg_ax2.set_ylabel(r"$\langle R \rangle$", rotation=0)
rg_ax2.yaxis.set_label_coords(labelx, 0.5)
rg_ax2.set_xlim([time[0], time[-1]])
rg_ax2.tick_params(axis='x', labelbottom='off')

legText = ['Skewness', 'Excess Kurtosis']
rg_ax2.legend(legText, bbox_to_anchor=(0, 1.20, 1, .1), loc=3, ncol=2,
  mode='expand', borderaxespad=0)
rg_ax2.plot([time[0], time[-1]], [0,0], 'k')

# Eigenvalue variance
var_ax2 = higher_measures.add_subplot(312)
var_ax2.plot(time, skewEVar, 'ko', markevery=15)
var_ax2.plot(time, kurtEVar, 'bo', markevery=15)

var_ax2.set_ylabel(r"$\langle \Delta \rangle$", rotation=0)
var_ax2.yaxis.set_label_coords(labelx, 0.5)
var_ax2.set_xlim([time[0], time[-1]])
var_ax2.set_ylim([-2, 2])
var_ax2.tick_params(axis='x', labelbottom='off')
var_ax2.set_yticks([-2, -1, 0, 1, 2])

var_ax2.plot([time[0], time[-1]], [0,0], 'k')

# Shape
s_ax2 = higher_measures.add_subplot(313)
s_ax2.plot(time, skewShape, 'ko', markevery=15)
s_ax2.plot(time, kurtShape, 'bo', markevery=15)

s_ax2.set_xlabel('Time [ms]')
s_ax2.set_ylabel(r"$\langle S \rangle$", rotation=0)
s_ax2.yaxis.set_label_coords(labelx, 0.5)
s_ax2.set_xlim([time[0], time[-1]])
s_ax2.set_ylim([-2, 6])
s_ax2.set_yticks([-2, 0, 2, 4, 6])

s_ax2.plot([time[0], time[-1]], [0,0], 'k')

imgname = imgdir + "shape_anisotropy_moments"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## PRINCIPAL AXIS ##
iFig = plt.figure()

# IMean and sdev
IMean = iFig.add_subplot(311)
IMean.plot(time, meanI1, 'ko', markevery=25)
IMean.plot(time, meanI1 + sdevI1, 'k.', markevery=15, alpha=0.5)
IMean.plot(time, meanI1 - sdevI1, 'k.', markevery=15, alpha=0.5)
IMean.plot(time, meanI2, 'ro', markevery=25)
IMean.plot(time, meanI2 + sdevI2, 'r.', markevery=15, alpha=0.5)
IMean.plot(time, meanI2 - sdevI2, 'r.', markevery=15, alpha=0.5)
IMean.plot(time, meanI3, 'bo', markevery=25)
IMean.plot(time, meanI3 + sdevI3, 'b.', markevery=15, alpha=0.5)
IMean.plot(time, meanI3 - sdevI3, 'b.', markevery=15, alpha=0.5)

IMean.set_ylabel(r'$\langle I_j \rangle$')
IMean.yaxis.set_label_coords(labelx, 0.5)
IMean.set_ylim([0, 1])
IMean.set_yticks([0, 0.25, 0.5, 0.75, 1])

IMean.set_xlim([time[0], time[-1]])
IMean.tick_params(axis='x', labelbottom='off')

I1s = mlines.Line2D([],[], color='k', marker='o', label=r'$I_1$')
I2s = mlines.Line2D([],[], color='r', marker='o', label=r'$I_2$')
I3s = mlines.Line2D([],[], color='b', marker='o', label=r'$I_3$')
IMean.legend(handles=[I1s, I2s, I3s], bbox_to_anchor=(0, 1.20, 1, 0.1), loc=3, 
  ncol=3, mode='expand', borderaxespad=0)

# Skew
ISkew = iFig.add_subplot(312)
ISkew.plot(time, skewI1, 'ko', markevery=25)
ISkew.plot(time, skewI2, 'ro', markevery=25)
ISkew.plot(time, skewI3, 'bo', markevery=25)

ISkew.set_ylabel('Skewness')
ISkew.yaxis.set_label_coords(labelx, 0.5)
ISkew.set_ylim([-1, 3])
ISkew.set_yticks([-1, 0, 1, 2, 3])
ISkew.set_xlim([time[0], time[-1]])
ISkew.tick_params(axis='x', labelbottom='off')

# kurt
IKurt = iFig.add_subplot(313)
IKurt.plot(time, kurtI1, 'ko', markevery=25)
IKurt.plot(time, kurtI2, 'ro', markevery=25)
IKurt.plot(time, kurtI3, 'bo', markevery=25)

IKurt.set_ylabel('Excess\nKurtosis')
IKurt.set_xlabel('Time [ms]')
IKurt.yaxis.set_label_coords(labelx, 0.5)
IKurt.set_ylim([-4, 12])
IKurt.set_yticks([-4,0,4,8,12])

IKurt.set_xlim([time[0], time[-1]])


imgname = imgdir + "shape_principal_axes_moments"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

# # Compressibility
# compFig = plt.figure(figsize=(12,8))
# 
# # mean
# meanAx = compFig.add_subplot(221)
# meanAx.set_xlabel('Time')
# meanAx.set_ylabel('Mean')
# meanAx.plot(time, meanS11, linewidth=2)
# meanAx.plot(time, meanS22, linewidth=2)
# meanAx.plot(time, meanS33, linewidth=2)
# meanAx.legend(['S11', 'S22', 'S33'])
# meanAx.set_xlim([500, 575])
# 
# # sdev
# sdevAx = compFig.add_subplot(222)
# sdevAx.set_xlabel('Time')
# sdevAx.set_ylabel('sdev')
# sdevAx.plot(time, sdevS11, linewidth=2)
# sdevAx.plot(time, sdevS22, linewidth=2)
# sdevAx.plot(time, sdevS33, linewidth=2)
# 
# # skew
# skewAx = compFig.add_subplot(223)
# skewAx.set_xlabel('Time')
# skewAx.set_ylabel('skew')
# skewAx.plot(time, skewS11, linewidth=2)
# skewAx.plot(time, skewS22, linewidth=2)
# skewAx.plot(time, skewS33, linewidth=2)
# 
# # kurt
# kurtAx = compFig.add_subplot(224)
# kurtAx.set_xlabel('Time')
# kurtAx.set_ylabel('kurt')
# kurtAx.plot(time, kurtS11, linewidth=2)
# kurtAx.plot(time, kurtS22, linewidth=2)
# kurtAx.plot(time, kurtS33, linewidth=2)



# # Surface plots
# rImg = rg_ax.imshow(rHist, origin='lower', aspect='auto', interpolation='none',
#          extent=[np.amin(time), np.amax(time), 0, rMax], cmap=cm.jet)
# pos1 = rg_ax.get_position()
# pos2 = [pos1.x0, pos1.y0 - 0.01, pos1.width*0.9, pos1.height]
# rg_ax.set_position(pos2)

#plt.show()
