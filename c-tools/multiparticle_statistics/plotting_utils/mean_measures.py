#!/usr/bin/env python2

import sys, os
import glob
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
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

# Plot over time
fs = 14
measures = plt.figure(figsize=(5,5))
measures.suptitle('Anisotropy Measures', fontsize=20)
labelx = -0.05

## SHAPE MEASURES ##
# Radius of Gyration
rg_ax = measures.add_subplot(311)
rg_ax.plot(time, meanR, 'ko-', linewidth=1.5, markevery=25)
rg_ax.plot(time, meanR + sdevR, 'k.', color='0', 
  linewidth=1, markevery=15, alpha=0.6)
rg_ax.plot(time, meanR - sdevR, 'k.', color='0', 
  linewidth=1, markevery=15, alpha=0.6)

rg_ax.set_ylabel("Radius of Gyration", fontsize=fs)
rg_ax.yaxis.set_label_coords(labelx, 0.5)
rg_ax.set_xlim([time[0], time[-1]])

# Variance
var_ax = measures.add_subplot(312)
var_ax.plot(time, meanEVar, 'ko-', markevery=25)
var_ax.plot(time, meanEVar + sdevEVar, 'k.', color='0', linewidth=1, 
  markevery=15, alpha=0.6)
var_ax.plot(time, meanEVar - sdevEVar, 'k.', color='0', linewidth=1, 
  markevery=15, alpha=0.6)
 
var_ax.set_ylabel("Eigenvalue Variance", fontsize=fs)
var_ax.yaxis.set_label_coords(labelx, 0.5)
var_ax.set_xlim([time[0], time[-1]])
var_ax.set_ylim([0,1])
var_ax.set_yticks(np.arange(0,1.1,0.2))

# Shape Factor
sf_ax = measures.add_subplot(313)
sf_ax.plot(time, meanShape, 'ko-', markevery=25)
sf_ax.plot(time, meanShape + sdevShape, 'k.', color='0', linewidth=1, 
  markevery=15, alpha=0.6)
sf_ax.plot(time, meanShape - sdevShape, 'k.', color='0', linewidth=1, 
  markevery=15, alpha=0.6)

sf_ax.set_xlabel("Time [ms]", fontsize=fs)
sf_ax.set_ylabel("Shape Factor", fontsize=fs)
sf_ax.yaxis.set_label_coords(labelx, 0.5)
sf_ax.set_xlim([time[0], time[-1]])
sf_ax.set_ylim([-0.25,2])
sf_ax.set_yticks([-0.25,0,0.5,1.0,1.5,2.0])

imgdir = root + simdir + "/img"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)
imgname = imgdir + "/test.png"
#plt.savefig(imgname, bbox_inches='tight')
plt.show()

## HIGHER ORDER MOMENTS ##
higher_measures = plt.figure(figsize=(12,8))

# Radius of gyration
rg_ax_2a = higher_measures.add_subplot(311)
rg_ax_2a.plot(time, skewR, 'ko-', linewidth=1.5, markevery=5)
rg_ax_2a.plot(time, kurtR, 'bo-', linewidth=1.5, markevery=5)
rg_ax_2a.set_xlabel('Time [ms]')
rg_ax_2a.set_ylabel('R', color='k')
rg_ax_2a.set_xlim([time[0], time[-1]])

rg_ax_2a.legend(['Skewness', 'Excess Kurtosis'], loc='upper left', ncol=2)
rg_ax_2a.xaxis.set_minor_locator(AutoMinorLocator())
rg_ax_2a.yaxis.set_minor_locator(AutoMinorLocator())
rg_ax_2a.tick_params(which='major', length=6)
rg_ax_2a.tick_params(which='minor', length=3)
rg_ax_2a.plot([time[0], time[-1]], [0,0], 'k')

#rg_ax_2b = rg_ax_2a.twinx()
#rg_ax_2b.plot(time, kurtR, 'bo-', linewidth=1.5, markevery=10)
#rg_ax_2b.set_xlabel('Time [ms]')
#rg_ax_2b.set_ylabel('Excess Kurtosis - R', color='b')
#rg_ax_2b.set_ylim([-1, 2.5])
#for t1 in rg_ax_2b.get_yticklabels():
#  t1.set_color('b')

# Eigenvalue variance
var_ax_2a = higher_measures.add_subplot(312)
var_ax_2a.plot(time, skewEVar, 'ko-', linewidth=1.5, markevery=10)
var_ax_2a.plot(time, kurtEVar, 'bo-', linewidth=1.5, markevery=10)
var_ax_2a.set_xlabel('Time [ms]')
var_ax_2a.set_ylabel('EVar', color='k')

var_ax_2a.xaxis.set_minor_locator(AutoMinorLocator())
var_ax_2a.yaxis.set_minor_locator(AutoMinorLocator())
var_ax_2a.tick_params(which='major', length=6)
var_ax_2a.tick_params(which='minor', length=3)
var_ax_2a.plot([time[0], time[-1]], [0,0], 'k')

# Shape
s_ax_2a = higher_measures.add_subplot(313)
s_ax_2a.plot(time, skewShape, 'ko-', linewidth=1.5, markevery=10)
s_ax_2a.plot(time, kurtShape, 'bo-', linewidth=1.5, markevery=10)
s_ax_2a.set_xlabel('Time [ms]')
s_ax_2a.set_ylabel('Shape', color='k')

s_ax_2a.xaxis.set_minor_locator(AutoMinorLocator())
s_ax_2a.yaxis.set_minor_locator(AutoMinorLocator())
s_ax_2a.tick_params(which='major', length=6)
s_ax_2a.tick_params(which='minor', length=3)
s_ax_2a.plot([time[0], time[-1]], [0,0], 'k')

## PRINCIPAL AXIS ##
iFig = plt.figure(figsize=(12,8))
iFig.suptitle('Size of Principal Axes', fontsize=20)

# IMean and sdev
IMean = iFig.add_subplot(311)
IMean.plot(time, meanI1, 'ko-', linewidth=1.5, markevery=10)
IMean.plot(time, meanI1 + sdevI1, 'k.', linewidth=1.5, markevery=2)
IMean.plot(time, meanI1 - sdevI1, 'k.', linewidth=1.5, markevery=2)
IMean.plot(time, meanI2, 'ro-', linewidth=1.5, markevery=10)
IMean.plot(time, meanI2 + sdevI2, 'r.', linewidth=1.5, markevery=2)
IMean.plot(time, meanI2 - sdevI2, 'r.', linewidth=1.5, markevery=2)
IMean.plot(time, meanI3, 'bo-', linewidth=1.5, markevery=10)
IMean.plot(time, meanI3 + sdevI3, 'b.', linewidth=1.5, markevery=2)
IMean.plot(time, meanI3 - sdevI3, 'b.', linewidth=1.5, markevery=2)
IMean.set_ylabel('I_j')

# Skew
ISkew = iFig.add_subplot(312)
ISkew.plot(time, skewI1, 'ko-', linewidth=1.5, markevery=10)
ISkew.plot(time, skewI2, 'ro-', linewidth=1.5, markevery=10)
ISkew.plot(time, skewI3, 'bo-', linewidth=1.5, markevery=10)
ISkew.set_ylabel('Skewness')

# kurt
IKurt = iFig.add_subplot(313)
IKurt.plot(time, kurtI1, 'ko-', linewidth=1.5, markevery=10)
IKurt.plot(time, kurtI2, 'ro-', linewidth=1.5, markevery=10)
IKurt.plot(time, kurtI3, 'bo-', linewidth=1.5, markevery=10)
IKurt.set_ylabel('excess kurtosis')
IKurt.legend(['I1', 'I2', 'I3'], loc='center right')

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
