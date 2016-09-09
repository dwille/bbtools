#!/usr/bin/env python2

import sys, os
import glob
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import re

## Get info
print ""
print " ---- Alignment plotting utility ---- "
print ""
root = '/home-1/dwillen3@jhu.edu/scratch/triply_per/'
simdir = raw_input("      Simulation root: ")
if not simdir.endswith('/'):
  simdir = simdir + '/'
datadir = root + simdir + 'data-tetrads/'
print "      Data directory set to: " + datadir

# Check if datadir exists so we don't go creating extra dirs
if not os.path.exists(datadir):
  print "      " + datadir + " does not exist. Exiting..."
  print ""
  sys.exit()

# Create imgdir if necessary
imgdir = root + simdir + "/img/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

alignFile = datadir + 'align.mean'

time = np.genfromtxt(alignFile, skip_header=1, usecols=0)
time = time - time[0]

g1s1 = np.genfromtxt(alignFile, skip_header=1, usecols=1)
g2s1 = np.genfromtxt(alignFile, skip_header=1, usecols=2)
g3s1 = np.genfromtxt(alignFile, skip_header=1, usecols=3)

g1s2 = np.genfromtxt(alignFile, skip_header=1, usecols=4)
g2s2 = np.genfromtxt(alignFile, skip_header=1, usecols=5)
g3s2 = np.genfromtxt(alignFile, skip_header=1, usecols=6)

g1s3 = np.genfromtxt(alignFile, skip_header=1, usecols=7)
g2s3 = np.genfromtxt(alignFile, skip_header=1, usecols=8)
g3s3 = np.genfromtxt(alignFile, skip_header=1, usecols=9)

g1z = np.genfromtxt(alignFile, skip_header=1, usecols=10)
g2z = np.genfromtxt(alignFile, skip_header=1, usecols=11)
g3z = np.genfromtxt(alignFile, skip_header=1, usecols=12)

s1z = np.genfromtxt(alignFile, skip_header=1, usecols=13)
s2z = np.genfromtxt(alignFile, skip_header=1, usecols=14)
s3z = np.genfromtxt(alignFile, skip_header=1, usecols=15)

wz = np.genfromtxt(alignFile, skip_header=1, usecols=16)

wg1 = np.genfromtxt(alignFile, skip_header=1, usecols=17)
wg2 = np.genfromtxt(alignFile, skip_header=1, usecols=18)
wg3 = np.genfromtxt(alignFile, skip_header=1, usecols=19)

ws1 = np.genfromtxt(alignFile, skip_header=1, usecols=20)
ws2 = np.genfromtxt(alignFile, skip_header=1, usecols=21)
ws3 = np.genfromtxt(alignFile, skip_header=1, usecols=22)

#wnorm = np.genfromtxt(alignFile, skip_header=1, usecols=23)

## Principal axes / strain ##
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', labelsize=11)
plt.rc('figure', titlesize=14)
plt.rc('figure', figsize=(4,3))
plt.rc('legend', fontsize=11, numpoints=3)
plt.rc('lines', markersize=4)
plt.rc('savefig', dpi=250)
labelx = -0.30

g_s = plt.figure()

# major
g1s_Ax = g_s.add_subplot(311)
g1s_Ax.plot(time, g1s1, 'ko', markevery=2)
g1s_Ax.plot(time, g1s2, 'bo', markevery=2)
g1s_Ax.plot(time, g1s3, 'ro', markevery=2)

g1s_Ax.set_ylabel(r'$\cos \theta = (\mathbf{g}_1 \cdot \mathbf{s(0)}_i)$', 
  rotation=0)
g1s_Ax.yaxis.set_label_coords(labelx, 0.5)
g1s_Ax.tick_params(axis='x', labelbottom='off')

g1s_Ax.set_ylim([0, 1])
g1s_Ax.set_xlim([0, 300])

#legText = ['Major Strain Axis', 'Middle Strain Axis', 'Minor Strain Axis']
legText = [r'$\mathbf{s}_1$', r'$\mathbf{s}_2$', r'$\mathbf{s}_3$']
g1s_Ax.legend(legText, bbox_to_anchor=(0, 1.20, 1, .1), loc=3, ncol=3,
  mode='expand', borderaxespad=0)

# middle
g2s_Ax = g_s.add_subplot(312)
g2s_Ax.plot(time, g2s1, 'ko', markevery=2)
g2s_Ax.plot(time, g2s2, 'bo', markevery=2)
g2s_Ax.plot(time, g2s3, 'ro', markevery=2)


g2s_Ax.set_ylabel(r'$\cos \theta = (\mathbf{g}_2 \cdot \mathbf{s(0)}_i)$', 
  rotation=0)
g2s_Ax.yaxis.set_label_coords(labelx, 0.5)
g2s_Ax.tick_params(axis='x', labelbottom='off')

g2s_Ax.set_ylim([0, 1])
g2s_Ax.set_xlim([0, 300])

g3s_Ax = g_s.add_subplot(313)
g3s_Ax.plot(time, g3s1, 'ko', markevery=2)
g3s_Ax.plot(time, g3s2, 'bo', markevery=2)
g3s_Ax.plot(time, g3s3, 'ro', markevery=2)

# minor
g3s_Ax.set_ylabel(r'$(\mathbf{g}_3 \cdot \mathbf{s(0)}_i)$',
  rotation=0)
g3s_Ax.yaxis.set_label_coords(labelx, 0.5)

g3s_Ax.set_ylim([0, 1])
g3s_Ax.set_xlim([0, 300])
g3s_Ax.set_xlabel('Time [ms]')

# SAVE 
imgname = imgdir + "align_shape_strain"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## VORTICITY / PAXES, STRAIN ##
wFig = plt.figure()

# principal axes
gw_Ax = wFig.add_subplot(211)
gw_Ax.plot(time, wg1, 'ko', markevery=2)
gw_Ax.plot(time, wg2, 'bo', markevery=2)
gw_Ax.plot(time, wg3, 'ro', markevery=2)

gw_Ax.set_ylabel(r'$(\mathbf{g(0)}_i \cdot \mathbf{\omega})$',
  rotation=0)
gw_Ax.yaxis.set_label_coords(labelx, 0.5)
gw_Ax.tick_params(axis='x', labelbottom='off')

#gw_Ax.legend(['Major Shape Axis', 'Middle Shape Axis', 'Minor Shape Axis'])
legText = [r'$i=1$', r'$i=2$', r'$i=3$']
gw_Ax.legend(legText, bbox_to_anchor=(0, 1.20, 1, .1), loc=3, ncol=3,
  mode='expand', borderaxespad=0)

# strain
sw_Ax = wFig.add_subplot(212)
sw_Ax.plot(time, ws1, 'ko', markevery=2)
sw_Ax.plot(time, ws2, 'bo', markevery=2)
sw_Ax.plot(time, ws3, 'ro', markevery=2)

sw_Ax.set_ylabel(r'$(\mathbf{s(0)}_i \cdot \mathbf{\omega})$',
  rotation=0)
sw_Ax.yaxis.set_label_coords(labelx, 0.5)
sw_Ax.set_xlabel('Time [ms]')

#sw_Ax.legend(['Major Strain Axis', 'Middle Strain Axis', 'Minor Strain Axis'])

# SAVE 
imgname = imgdir + "align_shape-strain_vort"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## ALIGNMENT WITH Z ##
g_z = plt.figure()

gzAx = g_z.add_subplot(311)
gzAx.plot(time, g1z, 'ko', markevery=25)
gzAx.plot(time, g2z, 'bo', markevery=25)
gzAx.plot(time, g3z, 'ro', markevery=25)

gzAx.set_ylim([0, 1])
gzAx.set_ylabel(r'$\cos \theta = (\mathbf{g}_i \cdot \mathbf{z})$',
  rotation=0)
gzAx.yaxis.set_label_coords(labelx, 0.5)
gzAx.tick_params(axis='x', labelbottom='off')

#legText = ['Major Shape Axis', 'Middle Shape Axis', 'Minor Shape Axis']
legText = [r'$i=1$', r'$i=2$', r'$i=3$']
gzAx.legend(legText, bbox_to_anchor=(0, 1.20, 1, .1), loc=3, ncol=3,
  mode='expand', borderaxespad=0)

szAx = g_z.add_subplot(312)
szAx.plot(time, s1z, 'ko', markevery=25)
szAx.plot(time, s2z, 'bo', markevery=25)
szAx.plot(time, s3z, 'ro', markevery=25)

szAx.set_ylim([0, 1])
szAx.set_ylabel(r'$\cos \theta = (\mathbf{s}_i \cdot \mathbf{z})$',
  rotation=0)
szAx.yaxis.set_label_coords(labelx, 0.5)
szAx.tick_params(axis='x', labelbottom='off')

#szAx.legend(['Major Strain Axis', 'Middle Strain Axis', 'Minor Strain Axis'])

wzAx = g_z.add_subplot(313)
wzAx.plot(time, wz, 'ko', markevery=25)

wzAx.set_xlabel('Time [ms]')
wzAx.set_ylabel(r'$\cos \theta = (\mathbf{\omega} \cdot \mathbf{z})$',
  rotation=0)
wzAx.yaxis.set_label_coords(labelx, 0.5)
wzAx.set_ylim([-0.05,0.05])
wzAx.set_yticks([-0.05, -0.025, 0, 0.025, 0.05])

# SAVE 
imgname = imgdir + "align_shape-strain_gravity"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
