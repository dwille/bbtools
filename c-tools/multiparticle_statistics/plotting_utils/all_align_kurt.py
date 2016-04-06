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

# Create imgdir if necessary
imgdir = root + "img/align/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# Parameter sweep
partList = ['500', '1000', '1500', '2000']
densList = ['rho2.0', 'rho3.3', 'rho4.0', 'rho5.0']
legendText = ['']*16

# Initialize arrays
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

# Plot specs
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('axes', labelsize=11)
plt.rc('figure', titlesize=14)
plt.rc('figure', figsize=(4,4))
plt.rc('legend', fontsize=10, numpoints=3)
plt.rc('lines', markersize=4, linewidth=1.5)
plt.rc('savefig', dpi=250)
labelx = -0.175
labely = 0.5
colors = ['r', 'g', 'b', 'k']
shades = [0.4, 0.57, 0.74, 0.9]

# Legend specs
rho2_spec = mlines.Line2D([],[], color='k', alpha=shades[0],
 label=r'$\rho^* = 2.0$')
rho3_spec = mlines.Line2D([],[], color='k', alpha=shades[1],
 label=r'$\rho^* = 3.3$')
rho4_spec = mlines.Line2D([],[], color='k', alpha=shades[2],
 label=r'$\rho^* = 4.0$')
rho5_spec = mlines.Line2D([],[], color='k', alpha=shades[3],
 label=r'$\rho^* = 5.0$')

align_spec_1 = mlines.Line2D([],[], color='k')
align_spec_2 = mlines.Line2D([],[], color='r')
align_spec_3 = mlines.Line2D([],[], color='b')

## gi_s1 ##
gis1_Fig, gis1_Ax = plt.subplots(4,1, sharex=True, sharey=True)

for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    gis1_Ax[pp].plot(data[i].time, g1_s1[i].kurt, color='k', alpha=shades[dd])
    gis1_Ax[pp].plot(data[i].time, g2_s1[i].kurt, color='r', alpha=shades[dd])
    gis1_Ax[pp].plot(data[i].time, g3_s1[i].kurt, color='b', alpha=shades[dd])

gis1_Ax[0].set_ylim([-2,4])
gis1_Ax[0].set_yticks([-2, 0, 2, 4])
gis1_Ax[0].set_xlim([0, 300])
gis1_Ax[3].set_xlabel(r'$Time\ [ms]$')
for ax in gis1_Ax:
  ax.grid(True)
for y, part in enumerate(partList):
  labelText = r'$n = ' + part + "$"
  gis1_Ax[y].set_ylabel(labelText, rotation=0)
  gis1_Ax[y].yaxis.set_label_coords(labelx, labely)

align_spec_1.set_label(r'$\mathrm{Kurt} [(g_1(t), s_1(0))] - 3$')
align_spec_2.set_label(r'$\mathrm{Kurt} [(g_2(t), s_1(0))] - 3$')
align_spec_3.set_label(r'$\mathrm{Kurt} [(g_3(t), s_1(0))] - 3$')
 
gis1_Ax[1].legend(handles=[rho2_spec, rho3_spec, rho4_spec, rho5_spec], ncol=1,
  bbox_to_anchor=(1, 0.65), loc='center')
gis1_Ax[2].legend(handles=[align_spec_1, align_spec_2, align_spec_3], ncol=1,
  bbox_to_anchor=(1, 0.1), loc='center')
gis1_Ax[2].set_zorder(1)

# Save
imgname = imgdir + "all_align_gis1_kurt"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## gi_s2 ##
gis2_Fig, gis2_Ax = plt.subplots(4,1, sharex=True, sharey=True)

for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    gis2_Ax[pp].plot(data[i].time, g1_s2[i].kurt, color='k', alpha=shades[dd])
    gis2_Ax[pp].plot(data[i].time, g2_s2[i].kurt, color='r', alpha=shades[dd])
    gis2_Ax[pp].plot(data[i].time, g3_s2[i].kurt, color='b', alpha=shades[dd])

gis2_Ax[0].set_ylim([-1.5, 0])
gis2_Ax[0].set_yticks([-1.5, -1, -0.5, 0])
gis2_Ax[0].set_xlim([0, 300])
gis2_Ax[3].set_xlabel(r'$Time\ [ms]$')
for ax in gis2_Ax:
  ax.grid(True)
for y, part in enumerate(partList):
  labelText = r'$n = ' + part + "$"
  gis2_Ax[y].set_ylabel(labelText, rotation=0)
  gis2_Ax[y].yaxis.set_label_coords(labelx - 0.05, labely)

align_spec_1.set_label(r'$\mathrm{Kurt} [(g_1(t), s_2(0))] - 3$')
align_spec_2.set_label(r'$\mathrm{Kurt} [(g_2(t), s_2(0))] - 3$')
align_spec_3.set_label(r'$\mathrm{Kurt} [(g_3(t), s_2(0))] - 3$')

gis2_Ax[1].legend(handles=[rho2_spec, rho3_spec, rho4_spec, rho5_spec], ncol=1,
  bbox_to_anchor=(1, 0.65), loc='center')
gis2_Ax[2].legend(handles=[align_spec_1, align_spec_2, align_spec_3], ncol=1,
  bbox_to_anchor=(1, 0.35), loc='center')
gis2_Ax[2].set_zorder(1)

# Save
imgname = imgdir + "all_align_gis2_kurt"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## gi_s3 ##
gis3_Fig, gis3_Ax = plt.subplots(4,1, sharex=True, sharey=True)

for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    gis3_Ax[pp].plot(data[i].time, g1_s3[i].kurt, color='k', alpha=shades[dd])
    gis3_Ax[pp].plot(data[i].time, g2_s3[i].kurt, color='r', alpha=shades[dd])
    gis3_Ax[pp].plot(data[i].time, g3_s3[i].kurt, color='b', alpha=shades[dd])

gis3_Ax[0].set_ylim([-2, 4])
gis3_Ax[0].set_yticks([-2, 0, 2, 4])
gis3_Ax[0].set_xlim([0, 300])
gis3_Ax[3].set_xlabel(r'$Time\ [ms]$')
for ax in gis3_Ax:
  ax.grid(True)
for y, part in enumerate(partList):
  labelText = r'$n = ' + part + "$"
  gis3_Ax[y].set_ylabel(labelText, rotation=0)
  gis3_Ax[y].yaxis.set_label_coords(labelx, labely)

align_spec_1.set_label(r'$\mathrm{Kurt} [(g_1(t), s_3(0))] - 3$')
align_spec_2.set_label(r'$\mathrm{Kurt} [(g_2(t), s_3(0))] - 3$')
align_spec_3.set_label(r'$\mathrm{Kurt} [(g_3(t), s_3(0))] - 3$')
 
gis3_Ax[1].legend(handles=[rho2_spec, rho3_spec, rho4_spec, rho5_spec], ncol=1,
  bbox_to_anchor=(1, .65), loc='center')
gis3_Ax[2].legend(handles=[align_spec_1, align_spec_2, align_spec_3], ncol=1,
  bbox_to_anchor=(1, .35), loc='center')
gis3_Ax[2].set_zorder(1)

# Save
imgname = imgdir + "all_align_gis3_kurt"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## gi_z ##
giz_Fig, giz_Ax = plt.subplots(4,1, sharex=True, sharey=True)

for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    giz_Ax[pp].plot(data[i].time, g1_z[i].kurt, color='k', alpha=shades[dd])
    giz_Ax[pp].plot(data[i].time, g2_z[i].kurt, color='r', alpha=shades[dd])
    giz_Ax[pp].plot(data[i].time, g3_z[i].kurt, color='b', alpha=shades[dd])

giz_Ax[0].set_ylim([-5, 15])
giz_Ax[0].set_yticks([-5, 0, 5, 10, 15])
giz_Ax[3].set_xlabel(r'$Time\ [ms]$')
for ax in giz_Ax:
  ax.grid(True)
for y, part in enumerate(partList):
  labelText = r'$n = ' + part + "$"
  giz_Ax[y].set_ylabel(labelText, rotation=0)
  giz_Ax[y].yaxis.set_label_coords(labelx - 0.025, labely)

align_spec_1.set_label(r'$\mathrm{Kurt} [(g_1(t), z)] - 3$')
align_spec_2.set_label(r'$\mathrm{Kurt} [(g_2(t), z)] - 3$')
align_spec_3.set_label(r'$\mathrm{Kurt} [(g_3(t), z)] - 3$')
 
giz_Ax[1].legend(handles=[rho2_spec, rho3_spec, rho4_spec, rho5_spec], ncol=1,
  bbox_to_anchor=(1, .6), loc='center')
giz_Ax[2].legend(handles=[align_spec_1, align_spec_2, align_spec_3], ncol=1,
  bbox_to_anchor=(1, .35), loc='center')
giz_Ax[2].set_zorder(1)

# Save
imgname = imgdir + "all_align_giz_kurt"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

## si_z ##
siz_Fig, siz_Ax = plt.subplots(4,1, sharex=True, sharey=True)

for pp in range(4):
  for dd in range(4):
    i = 4*pp + dd
    siz_Ax[pp].plot(data[i].time, s1_z[i].kurt, color='k', alpha=shades[dd])
    siz_Ax[pp].plot(data[i].time, s2_z[i].kurt, color='r', alpha=shades[dd])
    siz_Ax[pp].plot(data[i].time, s3_z[i].kurt, color='b', alpha=shades[dd])

siz_Ax[0].set_ylim([-1.25, -0.75])
siz_Ax[0].set_yticks([-1.25, -1, -0.75])
siz_Ax[0].set_xlim([0, 1000])
siz_Ax[3].set_xlabel(r'$Time\ [ms]$')
for ax in siz_Ax:
  ax.grid(True)
for y, part in enumerate(partList):
  labelText = r'$n = ' + part + "$"
  siz_Ax[y].set_ylabel(labelText, rotation=0)
  siz_Ax[y].yaxis.set_label_coords(labelx - 0.075, labely)

align_spec_1.set_label(r'$\mathrm{Kurt} [(s_1(t), z)] - 3$')
align_spec_2.set_label(r'$\mathrm{Kurt} [(s_2(t), z)] - 3$')
align_spec_3.set_label(r'$\mathrm{Kurt} [(s_3(t), z)] - 3$')
 
siz_Ax[1].legend(handles=[rho2_spec, rho3_spec, rho4_spec, rho5_spec], ncol=1,
  bbox_to_anchor=(1, .6), loc='center')
siz_Ax[2].legend(handles=[align_spec_1, align_spec_2, align_spec_3], ncol=1,
  bbox_to_anchor=(1, 0.35), loc='center')
siz_Ax[2].set_zorder(1)

# Save
imgname = imgdir + "all_align_siz_kurt"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
