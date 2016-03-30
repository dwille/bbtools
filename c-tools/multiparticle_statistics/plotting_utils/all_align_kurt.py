#!/usr/bin/env python2

import sys, os
import glob
import matplotlib.pyplot as plt
from matplotlib import lines as mlines
import numpy as np
import re

os.system('clear')

## Define structure class
class structtype():
  pass

## Get file length
def file_len(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i    # bc don't count header line

print ""
print " ---- Shape/Strain Alignment Plotting Utility ---- "
print "                   Kurtosis"
print ""
root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
print "      Sim root directory set to: " + root

# Parameter sweep
partList = ['500', '1000', '1500', '2000']
densList = ['rho2.0', 'rho3.3', 'rho4.0', 'rho5.0']
legendText = ['']*16

# Initialize kurt/std arrays
nTetrads = np.zeros(16)
data = [ structtype() for i in range(16) ]
g1_s1 = [ structtype() for i in range(16) ]
g2_s1 = [ structtype() for i in range(16) ]
g3_s1 = [ structtype() for i in range(16) ]
g1_s2 = [ structtype() for i in range(16) ]
g2_s2 = [ structtype() for i in range(16) ]
g3_s2 = [ structtype() for i in range(16) ]
g1_s3 = [ structtype() for i in range(16) ]
g2_s3 = [ structtype() for i in range(16) ]
g3_s3 = [ structtype() for i in range(16) ]
g1_z = [ structtype() for i in range(16) ]
g2_z = [ structtype() for i in range(16) ]
g3_z = [ structtype() for i in range(16) ]
s1_z = [ structtype() for i in range(16) ]
s2_z = [ structtype() for i in range(16) ]
s3_z = [ structtype() for i in range(16) ]
w_z = [ structtype() for i in range(16) ]
w_g1 = [ structtype() for i in range(16) ]
w_g2 = [ structtype() for i in range(16) ]
w_g3 = [ structtype() for i in range(16) ]
w_s1 = [ structtype() for i in range(16) ]
w_s2 = [ structtype() for i in range(16) ]
w_s3 = [ structtype() for i in range(16) ]
 
# Loop over all directory's
for pp, part in enumerate(partList):
  for dd, dens in enumerate(densList):
    stride = 4*pp + dd
    caseDir = part + '/' + dens
    infoFile = root + caseDir + '/data-tetrads/info.dat'
    nodeFile = root + caseDir + '/data-tetrads/regularNodes'
    alignKurt = root + caseDir + '/data-tetrads/align.kurt'

    nTetrads[stride] = np.genfromtxt(infoFile, skip_header=1, usecols=0)

    tmpT = np.genfromtxt(alignKurt, skip_header=1, usecols=0)
    data[stride].time = np.zeros(np.size(tmpT))
    data[stride].time = tmpT - tmpT[0]

    # KURT
    g1_s1[stride].kurt = np.zeros(np.size(tmpT))
    g2_s1[stride].kurt = np.zeros(np.size(tmpT))
    g3_s1[stride].kurt = np.zeros(np.size(tmpT))
    g1_s2[stride].kurt = np.zeros(np.size(tmpT))
    g2_s2[stride].kurt = np.zeros(np.size(tmpT))
    g3_s2[stride].kurt = np.zeros(np.size(tmpT))
    g1_s3[stride].kurt = np.zeros(np.size(tmpT))
    g2_s3[stride].kurt = np.zeros(np.size(tmpT))
    g3_s3[stride].kurt = np.zeros(np.size(tmpT))
    g1_z[stride].kurt = np.zeros(np.size(tmpT))
    g2_z[stride].kurt = np.zeros(np.size(tmpT))
    g3_z[stride].kurt = np.zeros(np.size(tmpT))
    s1_z[stride].kurt = np.zeros(np.size(tmpT))
    s2_z[stride].kurt = np.zeros(np.size(tmpT))
    s3_z[stride].kurt = np.zeros(np.size(tmpT))
    w_z[stride].kurt = np.zeros(np.size(tmpT))
    w_g1[stride].kurt = np.zeros(np.size(tmpT))
    w_g2[stride].kurt = np.zeros(np.size(tmpT))
    w_g3[stride].kurt = np.zeros(np.size(tmpT))
    w_s1[stride].kurt = np.zeros(np.size(tmpT))
    w_s2[stride].kurt = np.zeros(np.size(tmpT))
    w_s3[stride].kurt = np.zeros(np.size(tmpT))

    g1_s1[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=1)
    g2_s1[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=2)
    g3_s1[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=3)
    g1_s2[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=4)
    g2_s2[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=5)
    g3_s2[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=6)
    g1_s3[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=7)
    g2_s3[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=8)
    g3_s3[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=9)

    g1_z[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=10)
    g2_z[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=11)
    g3_z[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=12)
    s1_z[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=13)
    s2_z[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=14)
    s3_z[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=15)

    w_z[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=16)

    w_g1[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=17)
    w_g2[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=18)
    w_g3[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=19)
    w_s1[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=20)
    w_s2[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=21)
    w_s3[stride].kurt = np.genfromtxt(alignKurt, skip_header=1, usecols=22)


    legendText[stride] = caseDir + ': ' + str(nTetrads[stride])

# Plotting specs
plt.rc('font', family='serif')
color = ['r', 'g', 'b', 'k']
shade = [0.4, 0.57, 0.74, 0.9]
fsAx = 14
fSize = (12,8)
lWidth = 2

# Legend specs
rho2_spec = mlines.Line2D([],[], color='k', alpha=shade[0], linewidth=lWidth,
 label=r'$\rho^* = 2.0$')
rho3_spec = mlines.Line2D([],[], color='k', alpha=shade[1], linewidth=lWidth,
 label=r'$\rho^* = 3.3$')
rho4_spec = mlines.Line2D([],[], color='k', alpha=shade[2], linewidth=lWidth,
 label=r'$\rho^* = 4.0$')
rho5_spec = mlines.Line2D([],[], color='k', alpha=shade[3], linewidth=lWidth,
 label=r'$\rho^* = 5.0$')

g1si_spec = mlines.Line2D([],[], color='k', linewidth=lWidth,
  label=r'$(g_1, s_i)$')
g2si_spec = mlines.Line2D([],[], color='r', linewidth=lWidth,
  label=r'$(g_2, s_i)$')
g3si_spec = mlines.Line2D([],[], color='b', linewidth=lWidth,
  label=r'$(g_3, s_i)$')

## Shape and Strain Alignment ##

###########
## gi_s1 ##
###########
gis1_Fig, gis1_Ax = plt.subplots(4,1, sharex=True, sharey=True, figsize=fSize)
gis1_Fig.suptitle(r'$(g_i(t), s_1(0)) $ -- kurtosis', fontsize=16)

for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    gis1_Ax[pp].plot(data[i].time, g1_s1[i].kurt, linewidth=lWidth, color='k',
      alpha=shade[dd])
    gis1_Ax[pp].plot(data[i].time, g2_s1[i].kurt, linewidth=lWidth, color='r',
      alpha=shade[dd])
    gis1_Ax[pp].plot(data[i].time, g3_s1[i].kurt, linewidth=lWidth, color='b',
      alpha=shade[dd])

gis1_Ax[0].set_ylim([-1.5,1])
gis1_Ax[0].set_xlim([0, 300])
gis1_Ax[3].set_xlabel(r'Time [ms]', fontsize=fsAx)
for ax in gis1_Ax:
  ax.grid(True)
for y, part in enumerate(partList):
  gis1_Ax[y].set_ylabel(r'n = ' + part + '\n\nkurt ' + r'$\cos\theta$', fontsize=fsAx)
 
gis1_Ax[2].legend(handles=[rho2_spec, rho3_spec, rho4_spec, rho5_spec], ncol=1,
  bbox_to_anchor=(1, 1.6), framealpha=0.9)
gis1_Ax[3].legend(handles=[g1si_spec, g2si_spec, g3si_spec], ncol=1,
  bbox_to_anchor=(1, 1.6), framealpha=0.9)

###########
## gi_s2 ##
###########
gis2_Fig, gis2_Ax = plt.subplots(4,1, sharex=True, sharey=True, figsize=fSize)
gis2_Fig.suptitle(r'$(g_i(t), s_2(0))$ -- kurtosis', fontsize=16)

for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    gis2_Ax[pp].plot(data[i].time, g1_s2[i].kurt, linewidth=lWidth, color='k',
      alpha=shade[dd])
    gis2_Ax[pp].plot(data[i].time, g2_s2[i].kurt, linewidth=lWidth, color='r',
      alpha=shade[dd])
    gis2_Ax[pp].plot(data[i].time, g3_s2[i].kurt, linewidth=lWidth, color='b',
      alpha=shade[dd])

gis2_Ax[0].set_ylim([-1.5,1])
gis2_Ax[0].set_xlim([0, 300])
gis2_Ax[3].set_xlabel('Time [ms]', fontsize=fsAx)
for ax in gis2_Ax:
  ax.grid(True)
for y, part in enumerate(partList):
  gis2_Ax[y].set_ylabel(r'n = ' + part + '\n\nkurt ' + r'$\cos\theta$', fontsize=fsAx)

gis2_Ax[2].legend(handles=[rho2_spec, rho3_spec, rho4_spec, rho5_spec], ncol=1,
  bbox_to_anchor=(1, 1.6), framealpha=0.9)
gis2_Ax[3].legend(handles=[g1si_spec, g2si_spec, g3si_spec], ncol=1,
  bbox_to_anchor=(1, 1.6), framealpha=0.9)

###########
## gi_s3 ##
###########
gis3_Fig, gis3_Ax = plt.subplots(4,1, sharex=True, sharey=True, figsize=fSize)
gis3_Fig.suptitle(r'$(g_i(t), s_3(0))$ -- kurtosis', fontsize=16)

for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    gis3_Ax[pp].plot(data[i].time, g1_s3[i].kurt, linewidth=lWidth, color='k',
      alpha=shade[dd])
    gis3_Ax[pp].plot(data[i].time, g2_s3[i].kurt, linewidth=lWidth, color='r',
      alpha=shade[dd])
    gis3_Ax[pp].plot(data[i].time, g3_s3[i].kurt, linewidth=lWidth, color='b',
      alpha=shade[dd])

gis3_Ax[0].set_ylim([-1.5,1])
gis3_Ax[0].set_xlim([0, 300])
gis3_Ax[3].set_xlabel('Time [ms]', fontsize=fsAx)
for ax in gis3_Ax:
  ax.grid(True)
for y, part in enumerate(partList):
  gis3_Ax[y].set_ylabel(r'n = ' + part + '\n\nkurt ' + r'$\cos\theta$', fontsize=fsAx)
 
gis3_Ax[2].legend(handles=[rho2_spec, rho3_spec, rho4_spec, rho5_spec], ncol=1,
  bbox_to_anchor=(1, 1.6), framealpha=0.9)
gis3_Ax[3].legend(handles=[g1si_spec, g2si_spec, g3si_spec], ncol=1,
  bbox_to_anchor=(1, 1.6), framealpha=0.9)

###########
## gi_z ##
###########
giz_Fig, giz_Ax = plt.subplots(4,1, sharex=True, sharey=True, figsize=fSize)
giz_Fig.suptitle(r'$(g_i(t), z)$ -- kurtosis', fontsize=16)

for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    giz_Ax[pp].plot(data[i].time, g1_z[i].kurt, linewidth=lWidth, color='k',
      alpha=shade[dd])
    giz_Ax[pp].plot(data[i].time, g2_z[i].kurt, linewidth=lWidth, color='r',
      alpha=shade[dd])
    giz_Ax[pp].plot(data[i].time, g3_z[i].kurt, linewidth=lWidth, color='b',
      alpha=shade[dd])

giz_Ax[0].set_ylim([-1.5,1])
giz_Ax[0].set_xlim([0, 1000])
giz_Ax[3].set_xlabel('Time [ms]', fontsize=fsAx)
for ax in giz_Ax:
  ax.grid(True)
for y, part in enumerate(partList):
  giz_Ax[y].set_ylabel(r'n = ' + part + '\n\nkurt ' + r'$\cos\theta$', fontsize=fsAx)
 
giz_Ax[2].legend(handles=[rho2_spec, rho3_spec, rho4_spec, rho5_spec], ncol=1,
  bbox_to_anchor=(1, 1.6), framealpha=0.9)
giz_Ax[3].legend(handles=[g1si_spec, g2si_spec, g3si_spec], ncol=1,
  bbox_to_anchor=(1, 1.6), framealpha=0.9)

###########
## si_z ##
###########
siz_Fig, siz_Ax = plt.subplots(4,1, sharex=True, sharey=True, figsize=fSize)
siz_Fig.suptitle(r'$(s_i(t), z)$ -- kurtosis', fontsize=16)

for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    siz_Ax[pp].plot(data[i].time, s1_z[i].kurt, linewidth=lWidth, color='k',
      alpha=shade[dd])
    siz_Ax[pp].plot(data[i].time, s2_z[i].kurt, linewidth=lWidth, color='r',
      alpha=shade[dd])
    siz_Ax[pp].plot(data[i].time, s3_z[i].kurt, linewidth=lWidth, color='b',
      alpha=shade[dd])

siz_Ax[0].set_ylim([-1.5,1])
siz_Ax[0].set_xlim([0, 1000])
siz_Ax[3].set_xlabel('Time [ms]', fontsize=fsAx)
for ax in siz_Ax:
  ax.grid(True)
for y, part in enumerate(partList):
  siz_Ax[y].set_ylabel(r'n = ' + part + '\n\nkurt ' + r'$\cos\theta$', fontsize=fsAx)
 
siz_Ax[2].legend(handles=[rho2_spec, rho3_spec, rho4_spec, rho5_spec], ncol=1,
  bbox_to_anchor=(1, 1.6), framealpha=0.9)
siz_Ax[3].legend(handles=[g1si_spec, g2si_spec, g3si_spec], ncol=1,
  bbox_to_anchor=(1, 1.6), framealpha=0.9)


plt.show()
