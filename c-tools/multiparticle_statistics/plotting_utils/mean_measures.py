#!/usr/bin/env python2

import sys
import glob
import matplotlib.pyplot as plt
from matplotlib import cm, colors
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
time = np.genfromtxt(statFile, skip_header=1, usecols= 0)
R2Mean = np.genfromtxt(statFile, skip_header=1, usecols = 1)
varMean = np.genfromtxt(statFile, skip_header=1, usecols = 2)
shapeMean = np.genfromtxt(statFile, skip_header=1, usecols = 3)
R2Std = np.genfromtxt(statFile, skip_header=1, usecols = 4)
varStd = np.genfromtxt(statFile, skip_header=1, usecols = 5)
shapeStd = np.genfromtxt(statFile, skip_header=1, usecols = 6)

# Plot over time
fs = 14
measures = plt.figure(figsize=(12,8))
measures.suptitle('Anisotropy Measures', fontsize=20)
labelx = -0.05

# Radius of Gyration
rg_ax = measures.add_subplot(311)
rg_ax.plot(time, np.sqrt(R2Mean), 'ko-', linewidth=1.5, markevery=10)
rg_ax.plot(time, np.sqrt(R2Mean) + np.sqrt(R2Std), 'kx--', color='0', 
  linewidth=1, markevery=10)
rg_ax.plot(time, np.sqrt(R2Mean) - np.sqrt(R2Std), 'kx--', color='0', 
  linewidth=1, markevery=10)

rg_ax.set_ylabel("RMS Radius of Gyration", fontsize=fs)
rg_ax.yaxis.set_label_coords(labelx, 0.5)
rg_ax.set_xlim([time[0], time[-1]])

# Variance
var_ax = measures.add_subplot(312)
var_ax.plot(time, varMean, 'ko-', markevery=10)
var_ax.plot(time, varMean + varStd, 'kx--', color='0', linewidth=1, 
  markevery=10)
var_ax.plot(time, varMean - varStd, 'kx--', color='0', linewidth=1, 
  markevery=10)
 
var_ax.set_ylabel("Eigenvalue Variance", fontsize=fs)
var_ax.yaxis.set_label_coords(labelx, 0.5)
var_ax.set_xlim([time[0], time[-1]])
var_ax.set_ylim([0,1])
var_ax.set_yticks(np.arange(0,1.1,0.2))

# Shape Factor
sf_ax = measures.add_subplot(313)
sf_ax.plot(time, shapeMean, 'ko-', markevery=10)
sf_ax.plot(time, shapeMean + shapeStd, 'kx--', color='0', linewidth=1, 
  markevery=10)
sf_ax.plot(time, shapeMean - shapeStd, 'kx--', color='0', linewidth=1, 
  markevery=10)

sf_ax.set_xlabel("Time [ms]", fontsize=fs)
sf_ax.set_ylabel("Shape Factor", fontsize=fs)
sf_ax.yaxis.set_label_coords(labelx, 0.5)
sf_ax.set_xlim([time[0], time[-1]])
sf_ax.set_ylim([-0.25,2])
sf_ax.set_yticks([-0.25,0,0.5,1.0,1.5,2.0])


plt.show()

# # Surface plots
# rImg = rg_ax.imshow(rHist, origin='lower', aspect='auto', interpolation='none',
#          extent=[np.amin(time), np.amax(time), 0, rMax], cmap=cm.jet)
# pos1 = rg_ax.get_position()
# pos2 = [pos1.x0, pos1.y0 - 0.01, pos1.width*0.9, pos1.height]
# rg_ax.set_position(pos2)
