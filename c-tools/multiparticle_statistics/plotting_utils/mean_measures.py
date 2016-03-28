#!/usr/bin/env python2

import sys
import glob
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import re

## Get info
print "Mean Measures plotting utility"
print ""
#root = raw_input("Simulation root: ")
#if not root.endswith('/'):
#  root = root + '/'
#ts = raw_input("Desired start time: ")
#te = raw_input("Desired end time: ")
root = "../sim/"
ts = "500"
te = "1000"
statFile = root + "data-tetrads/stat.dat"

# Find means and time from stat.dat
time = np.genfromtxt(statFile, skip_header=1, usecols = 0)
meanR= np.genfromtxt(statFile, skip_header=1, usecols = 1)
meanEVar= np.genfromtxt(statFile, skip_header=1, usecols = 2)
meanShape = np.genfromtxt(statFile, skip_header=1, usecols = 3)
meanS11 = np.genfromtxt(statFile, skip_header=1, usecols = 4)
meanS22 = np.genfromtxt(statFile, skip_header=1, usecols = 5)
meanS33 = np.genfromtxt(statFile, skip_header=1, usecols = 6)
sdevR = np.genfromtxt(statFile, skip_header=1, usecols = 7)
sdevEVar = np.genfromtxt(statFile, skip_header=1, usecols = 8)
sdevShape = np.genfromtxt(statFile, skip_header=1, usecols = 9)
sdevS11 = np.genfromtxt(statFile, skip_header=1, usecols = 10)
sdevS22 = np.genfromtxt(statFile, skip_header=1, usecols = 11)
sdevS33 = np.genfromtxt(statFile, skip_header=1, usecols = 12)
skewR = np.genfromtxt(statFile, skip_header=1, usecols = 13);
skewEVar = np.genfromtxt(statFile, skip_header=1, usecols = 14);
skewShape = np.genfromtxt(statFile, skip_header=1, usecols = 15);
skewS11 = np.genfromtxt(statFile, skip_header=1, usecols = 16)
skewS22 = np.genfromtxt(statFile, skip_header=1, usecols = 17)
skewS33 = np.genfromtxt(statFile, skip_header=1, usecols = 18)
kurtR = np.genfromtxt(statFile, skip_header=1, usecols = 19);
kurtEVar = np.genfromtxt(statFile, skip_header=1, usecols = 20);
kurtShape = np.genfromtxt(statFile, skip_header=1, usecols = 21);
kurtS11 = np.genfromtxt(statFile, skip_header=1, usecols = 22)
kurtS22 = np.genfromtxt(statFile, skip_header=1, usecols = 23)
kurtS33 = np.genfromtxt(statFile, skip_header=1, usecols = 24)

# Plot over time
fs = 14
measures = plt.figure(figsize=(12,8))
measures.suptitle('Anisotropy Measures', fontsize=20)
labelx = -0.05

# Radius of Gyration
rg_ax = measures.add_subplot(311)
rg_ax.plot(time, meanR, 'ko-', linewidth=1.5, markevery=10)
rg_ax.plot(time, meanR + sdevR, 'k.', color='0', 
  linewidth=1, markevery=2)
rg_ax.plot(time, meanR - sdevR, 'k.', color='0', 
  linewidth=1, markevery=2)

rg_ax.set_ylabel("Radius of Gyration", fontsize=fs)
rg_ax.yaxis.set_label_coords(labelx, 0.5)
rg_ax.set_xlim([time[0], time[-1]])

# Variance
var_ax = measures.add_subplot(312)
var_ax.plot(time, meanEVar, 'ko-', markevery=10)
var_ax.plot(time, meanEVar + sdevEVar, 'k.', color='0', linewidth=1, 
  markevery=2)
var_ax.plot(time, meanEVar - sdevEVar, 'k.', color='0', linewidth=1, 
  markevery=2)
 
var_ax.set_ylabel("Eigenvalue Variance", fontsize=fs)
var_ax.yaxis.set_label_coords(labelx, 0.5)
var_ax.set_xlim([time[0], time[-1]])
var_ax.set_ylim([0,1])
var_ax.set_yticks(np.arange(0,1.1,0.2))

# Shape Factor
sf_ax = measures.add_subplot(313)
sf_ax.plot(time, meanShape, 'ko-', markevery=10)
sf_ax.plot(time, meanShape + sdevShape, 'k.', color='0', linewidth=1, 
  markevery=2)
sf_ax.plot(time, meanShape - sdevShape, 'k.', color='0', linewidth=1, 
  markevery=2)

sf_ax.set_xlabel("Time [ms]", fontsize=fs)
sf_ax.set_ylabel("Shape Factor", fontsize=fs)
sf_ax.yaxis.set_label_coords(labelx, 0.5)
sf_ax.set_xlim([time[0], time[-1]])
sf_ax.set_ylim([-0.25,2])
sf_ax.set_yticks([-0.25,0,0.5,1.0,1.5,2.0])

## HIGHER ORDER MOMENTS ##
higher_measures = plt.figure(figsize=(12,8))

# Radius of gyration
rg_ax_2a = higher_measures.add_subplot(311)
rg_ax_2a.plot(time, skewR, 'ko-', linewidth=1.5, markevery=10)
rg_ax_2a.plot(time, kurtR, 'bo-', linewidth=1.5, markevery=10)
rg_ax_2a.set_xlabel('Time [ms]')
rg_ax_2a.set_ylabel('R', color='k')
rg_ax_2a.set_xlim([time[0], time[-1]])
#for t1 in rg_ax_2a.get_yticklabels():
#  t1.set_color('k')
rg_ax_2a.legend(['Skewness', 'Kurtosis'], loc='upper left', ncol=2)
rg_ax_2a.xaxis.set_minor_locator(AutoMinorLocator())
rg_ax_2a.yaxis.set_minor_locator(AutoMinorLocator())
rg_ax_2a.tick_params(which='major', length=6)
rg_ax_2a.tick_params(which='minor', length=3)

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
var_ax_2a.set_xlim([time[0], time[-1]])

var_ax_2a.xaxis.set_minor_locator(AutoMinorLocator())
var_ax_2a.yaxis.set_minor_locator(AutoMinorLocator())
var_ax_2a.tick_params(which='major', length=6)
var_ax_2a.tick_params(which='minor', length=3)

# Shape
s_ax_2a = higher_measures.add_subplot(313)
s_ax_2a.plot(time, skewShape, 'ko-', linewidth=1.5, markevery=10)
s_ax_2a.plot(time, kurtShape, 'bo-', linewidth=1.5, markevery=10)
s_ax_2a.set_xlabel('Time [ms]')
s_ax_2a.set_ylabel('Shape', color='k')
s_ax_2a.set_xlim([time[0], time[-1]])

s_ax_2a.xaxis.set_minor_locator(AutoMinorLocator())
s_ax_2a.yaxis.set_minor_locator(AutoMinorLocator())
s_ax_2a.tick_params(which='major', length=6)
s_ax_2a.tick_params(which='minor', length=3)

# Compressibility
compFig = plt.figure(figsize=(12,8))

# mean
meanAx = compFig.add_subplot(221)
meanAx.set_xlabel('Time')
meanAx.set_ylabel('Mean')
meanAx.plot(time, meanS11, linewidth=2)
meanAx.plot(time, meanS22, linewidth=2)
meanAx.plot(time, meanS33, linewidth=2)
meanAx.legend(['S11', 'S22', 'S33'])
meanAx.set_xlim([500, 575])

# sdev
sdevAx = compFig.add_subplot(222)
sdevAx.set_xlabel('Time')
sdevAx.set_ylabel('sdev')
sdevAx.plot(time, sdevS11, linewidth=2)
sdevAx.plot(time, sdevS22, linewidth=2)
sdevAx.plot(time, sdevS33, linewidth=2)

# skew
skewAx = compFig.add_subplot(223)
skewAx.set_xlabel('Time')
skewAx.set_ylabel('skew')
skewAx.plot(time, skewS11, linewidth=2)
skewAx.plot(time, skewS22, linewidth=2)
skewAx.plot(time, skewS33, linewidth=2)

# kurt
kurtAx = compFig.add_subplot(224)
kurtAx.set_xlabel('Time')
kurtAx.set_ylabel('kurt')
kurtAx.plot(time, kurtS11, linewidth=2)
kurtAx.plot(time, kurtS22, linewidth=2)
kurtAx.plot(time, kurtS33, linewidth=2)

plt.show()


# # Surface plots
# rImg = rg_ax.imshow(rHist, origin='lower', aspect='auto', interpolation='none',
#          extent=[np.amin(time), np.amax(time), 0, rMax], cmap=cm.jet)
# pos1 = rg_ax.get_position()
# pos2 = [pos1.x0, pos1.y0 - 0.01, pos1.width*0.9, pos1.height]
# rg_ax.set_position(pos2)
