#!/usr/bin/env python2

import sys
import glob
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import re

## Nicely sorted function
def sorted_nicely( l ):
  """ Sorts the given iterable in a natural way

  Required arguments:
  l -- The iterable to be sroted.

  courtesy stackoverflow/2669059
  """
  convert = lambda text: int(text) if text.isdigit() else text
  alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
  return sorted(l, key = alphanum_key)

## Number of lines in files function (not includinng header line)
def file_len(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i    # Does not include header

## Alignment of two vectors -- <(a.b)^2>
def alignVects(vec1, vec2):
  tmp = np.multiply(vec1, vec2)
  tmp = np.sum(tmp, axis=1)
  #tmp = np.square(tmp)    # why square/L2norm?
  tmp = np.mean(tmp)
  return tmp
def erralign(vec1, vec2):
  tmp = np.multiply(vec1, vec2)
  tmp = np.sum(tmp, axis=1)
  #tmp = np.square(tmp)
  tmp = np.std(tmp)
  return tmp

## Define structure class
class structtype():
  pass

## Get info
print "Anisotropy Measures plotting utility"
print ""
#root = raw_input("Simulation root: ")
#if not root.endswith('/'):
#  root = root + '/'
#ts = raw_input("Desired start time: ")
#te = raw_input("Desired end time: ")
root = "../sim/"
ts = "500"
te = "1000"

## Sort output files and find number within time range
files = sorted_nicely(glob.glob(root + "data-tetrads/raw-data-*"))
inRange = np.zeros(len(files))
n = 0
for fname in files:
  ftime = fname.split('/')[-1]

  if ftime.startswith('raw-data-'):
    time = ftime[9:]

  if float(time) >= float(ts) and float(time) <= float(te):
    inRange[n] = 1
    n += 1

nFiles = n

print "Found " + str(nFiles) + " files in time range"
print ""

# Find number of tetrads
nTetrads = file_len(files[0])

# Initialize numpy arrays
gEigVals = [ structtype() for i in range(nFiles) ]
gEigVecs = [ structtype() for i in range(nFiles) ]
sEigVals = [ structtype() for i in range(nFiles) ]
sEigVecs = [ structtype() for i in range(nFiles) ]
vorticity = [ structtype() for i in range(nFiles) ]
alignment = [ structtype() for i in range(1) ]

for i in range(nFiles):
  gEigVals[i].vals = np.zeros((nTetrads, 3))
  gEigVecs[i].vals = np.zeros((nTetrads, 9))
  sEigVals[i].vals = np.zeros((nTetrads, 3))
  sEigVecs[i].vals = np.zeros((nTetrads, 9))
  vorticity[i].vals = np.zeros((nTetrads, 3))

alignment[0].g1_s1 = np.zeros(nFiles)
alignment[0].g1_s2 = np.zeros(nFiles)
alignment[0].g1_s3 = np.zeros(nFiles)
alignment[0].g2_s1 = np.zeros(nFiles)
alignment[0].g2_s2 = np.zeros(nFiles)
alignment[0].g2_s3 = np.zeros(nFiles)
alignment[0].g3_s1 = np.zeros(nFiles)
alignment[0].g3_s2 = np.zeros(nFiles)
alignment[0].g3_s3 = np.zeros(nFiles)
alignment[0].g1_w = np.zeros(nFiles)
alignment[0].g2_w = np.zeros(nFiles)
alignment[0].g3_w = np.zeros(nFiles)
alignment[0].s1_w = np.zeros(nFiles)
alignment[0].s2_w = np.zeros(nFiles)
alignment[0].s3_w = np.zeros(nFiles)

alignment[0].g1_z = np.zeros(nFiles)
alignment[0].g2_z = np.zeros(nFiles)
alignment[0].g3_z = np.zeros(nFiles)
alignment[0].s1_z = np.zeros(nFiles)
alignment[0].s2_z = np.zeros(nFiles)
alignment[0].s3_z = np.zeros(nFiles)
alignment[0].w_z = np.zeros(nFiles)

alignment[0].err_g1_s1 = np.zeros(nFiles)
alignment[0].err_g1_s2 = np.zeros(nFiles)
alignment[0].err_g1_s3 = np.zeros(nFiles)
alignment[0].err_g2_s1 = np.zeros(nFiles)
alignment[0].err_g2_s2 = np.zeros(nFiles)
alignment[0].err_g2_s3 = np.zeros(nFiles)
alignment[0].err_g3_s1 = np.zeros(nFiles)
alignment[0].err_g3_s2 = np.zeros(nFiles)
alignment[0].err_g3_s3 = np.zeros(nFiles)
alignment[0].err_g1_w = np.zeros(nFiles)
alignment[0].err_g2_w = np.zeros(nFiles)
alignment[0].err_g3_w = np.zeros(nFiles)
alignment[0].err_s1_w = np.zeros(nFiles)
alignment[0].err_s2_w = np.zeros(nFiles)
alignment[0].err_s3_w = np.zeros(nFiles)

alignment[0].err_g1_z = np.zeros(nFiles)
alignment[0].err_g2_z = np.zeros(nFiles)
alignment[0].err_g3_z = np.zeros(nFiles)
alignment[0].err_s1_z = np.zeros(nFiles)
alignment[0].err_s2_z = np.zeros(nFiles)
alignment[0].err_s3_z = np.zeros(nFiles)
alignment[0].err_w_z = np.zeros(nFiles)

z = np.zeros((nTetrads, 3))
z[:,2] = 1;
vortNorm = np.zeros((nTetrads, nFiles))
time = np.zeros(nFiles)

# Enumerate the columns
gEigValsCols = (3,4,5)
gEigVecsCols = (6,7,8,9,10,11,12,13,14)
sEigValsCols = (15,16,17)
sEigVecsCols = (18,19,20,21,22,23,24,25,26)
vorticityCols = (27,28,29)

maxAx = (0,1,2)
medAx = (3,4,5)
minAx = (6,7,8)

# Initialize histogram bins
nBins = float(100)
g1s1 = np.zeros([nBins, nFiles])
edges = np.linspace(0,0.3,nBins + 1)
centers = edges[0:nBins] + 0.5*0.3/nBins

# Loop over all output files, pull data to structures
i = 0;
for fname in files:
  ftime = fname.split('/')[-1]

  if ftime.startswith('raw-data-'):
    ftime = ftime[9:]

  if float(ftime) >= float(ts) and float(ftime) <= float(te):
    time[i] = float(ftime)
    ifile = open(fname)
    ifile.readline()

    # If this gets slow, see
    # -- stackoverlfow 8956832 -- python out of memory on large csv file numpy
    # -- "  " 18259393 -- numpy loading csv too slow compared to matlab
    # -- "  " 26482209 -- fastest way to load huge dat into array

    gEigVals[i].vals = np.genfromtxt(fname, skip_header=1, usecols = 
      gEigValsCols)
    gEigVecs[i].vals = np.genfromtxt(fname, skip_header=1, usecols = 
      gEigVecsCols)
    sEigVals[i].vals = np.genfromtxt(fname, skip_header=1, usecols = 
      sEigValsCols)
    sEigVecs[i].vals = np.genfromtxt(fname, skip_header=1, usecols = 
      sEigVecsCols)
    vorticity[i].vals = np.genfromtxt(fname, skip_header=1, usecols = 
      vorticityCols)

    vortNorm[:,i] = np.sqrt(np.square(vorticity[i].vals[:,0]) 
                          + np.square(vorticity[i].vals[:,1])
                          + np.square(vorticity[i].vals[:,2]))

    vorticity[i].vals[:,0] = np.divide(vorticity[i].vals[:,0], vortNorm[:,i])
    vorticity[i].vals[:,1] = np.divide(vorticity[i].vals[:,1], vortNorm[:,i])
    vorticity[i].vals[:,2] = np.divide(vorticity[i].vals[:,2], vortNorm[:,i])

    # Alignment of principal axes of tetrads with prncpl axes of strain
    alignment[0].g1_s1[i] = alignVects(gEigVecs[i].vals[:,maxAx],
                                       sEigVecs[0].vals[:,maxAx])
    alignment[0].g1_s2[i] = alignVects(gEigVecs[i].vals[:,maxAx],
                                       sEigVecs[0].vals[:,medAx])
    alignment[0].g1_s3[i] = alignVects(gEigVecs[i].vals[:,maxAx],
                                       sEigVecs[0].vals[:,minAx])

    alignment[0].g2_s1[i] = alignVects(gEigVecs[i].vals[:,medAx],
                                       sEigVecs[0].vals[:,maxAx])
    alignment[0].g2_s2[i] = alignVects(gEigVecs[i].vals[:,medAx],
                                       sEigVecs[0].vals[:,medAx])
    alignment[0].g2_s3[i] = alignVects(gEigVecs[i].vals[:,medAx],
                                       sEigVecs[0].vals[:,minAx])

    alignment[0].g3_s1[i] = alignVects(gEigVecs[i].vals[:,minAx],
                                       sEigVecs[0].vals[:,maxAx])
    alignment[0].g3_s2[i] = alignVects(gEigVecs[i].vals[:,minAx],
                                       sEigVecs[0].vals[:,medAx])
    alignment[0].g3_s3[i] = alignVects(gEigVecs[i].vals[:,minAx],
                                       sEigVecs[0].vals[:,minAx])

    alignment[0].g1_w[i] = alignVects(vorticity[i].vals,
                                      gEigVecs[0].vals[:,maxAx])

    alignment[0].g2_w[i] = alignVects(gEigVecs[0].vals[:,medAx],
                                      vorticity[i].vals)
    alignment[0].g3_w[i] = alignVects(gEigVecs[0].vals[:,minAx],
                                      vorticity[i].vals)

    alignment[0].s1_w[i] = alignVects(sEigVecs[0].vals[:,maxAx],
                                      vorticity[i].vals)
    alignment[0].s2_w[i] = alignVects(sEigVecs[0].vals[:,medAx],
                                      vorticity[i].vals)
    alignment[0].s3_w[i] = alignVects(sEigVecs[0].vals[:,minAx],
                                      vorticity[i].vals)
    
    alignment[0].g1_z[i] = alignVects(gEigVecs[i].vals[:,maxAx], z)
    alignment[0].g2_z[i] = alignVects(gEigVecs[i].vals[:,medAx], z)
    alignment[0].g3_z[i] = alignVects(gEigVecs[i].vals[:,minAx], z)
    alignment[0].s1_z[i] = alignVects(sEigVecs[i].vals[:,maxAx], z)
    alignment[0].s2_z[i] = alignVects(sEigVecs[i].vals[:,medAx], z)
    alignment[0].s3_z[i] = alignVects(sEigVecs[i].vals[:,minAx], z)
    alignment[0].w_z[i] = alignVects(vorticity[i].vals, z)

    # STD of alignment
    alignment[0].err_g1_s1[i] = alignVects(gEigVecs[i].vals[:,maxAx],
                                           sEigVecs[0].vals[:,maxAx])
    alignment[0].err_g1_s2[i] = alignVects(gEigVecs[i].vals[:,maxAx],
                                           sEigVecs[0].vals[:,medAx])
    alignment[0].err_g1_s3[i] = alignVects(gEigVecs[i].vals[:,maxAx],
                                           sEigVecs[0].vals[:,minAx])

    alignment[0].err_g2_s1[i] = alignVects(gEigVecs[i].vals[:,medAx],
                                           sEigVecs[0].vals[:,maxAx])
    alignment[0].err_g2_s2[i] = alignVects(gEigVecs[i].vals[:,medAx],
                                           sEigVecs[0].vals[:,medAx])
    alignment[0].err_g2_s3[i] = alignVects(gEigVecs[i].vals[:,medAx],
                                           sEigVecs[0].vals[:,minAx])

    alignment[0].err_g3_s1[i] = alignVects(gEigVecs[i].vals[:,minAx],
                                           sEigVecs[0].vals[:,maxAx])
    alignment[0].err_g3_s2[i] = alignVects(gEigVecs[i].vals[:,minAx],
                                           sEigVecs[0].vals[:,medAx])
    alignment[0].err_g3_s3[i] = alignVects(gEigVecs[i].vals[:,minAx],
                                           sEigVecs[0].vals[:,minAx])

    alignment[0].err_g1_w[i] = alignVects(gEigVecs[i].vals[:,maxAx],
                                          vorticity[0].vals)
    alignment[0].err_g2_w[i] = alignVects(gEigVecs[i].vals[:,medAx],
                                          vorticity[0].vals)
    alignment[0].err_g3_w[i] = alignVects(gEigVecs[i].vals[:,minAx],
                                          vorticity[0].vals)

    alignment[0].err_s1_w[i] = alignVects(sEigVecs[i].vals[:,maxAx],
                                          vorticity[0].vals)
    alignment[0].err_s2_w[i] = alignVects(sEigVecs[i].vals[:,medAx],
                                          vorticity[0].vals)
    alignment[0].err_s3_w[i] = alignVects(sEigVecs[i].vals[:,minAx],
                                          vorticity[0].vals)
    
    alignment[0].err_g1_z[i] = alignVects(gEigVecs[i].vals[:,maxAx], z)
    alignment[0].err_g2_z[i] = alignVects(gEigVecs[i].vals[:,medAx], z)
    alignment[0].err_g3_z[i] = alignVects(gEigVecs[i].vals[:,minAx], z)
    alignment[0].err_s1_z[i] = alignVects(sEigVecs[i].vals[:,maxAx], z)
    alignment[0].err_s2_z[i] = alignVects(sEigVecs[i].vals[:,medAx], z)
    alignment[0].err_s3_z[i] = alignVects(sEigVecs[i].vals[:,minAx], z)
    alignment[0].err_w_z[i] = alignVects(vorticity[i].vals, z)

    i += 1
    ifile.close()

plt.rc('font', family='serif')


## Principal axes / strain ##
g_s = plt.figure(figsize=(12,8))

g_s.suptitle('Alignment of Shape Principal Axes with Strain', fontsize=16)

# major
g1s_Ax = g_s.add_subplot(311)
g1s_Ax.plot(time, alignment[0].g1_s1[:], 'ko-', linewidth=1.5)
g1s_Ax.plot(time, alignment[0].g1_s2[:], 'bo-', linewidth=1.5)
g1s_Ax.plot(time, alignment[0].g1_s3[:], 'ro-', linewidth=1.5)

g1s_Ax.set_ylabel(r'$\cos \theta = (\mathbf{g}_1 \cdot \mathbf{s(0)}_i)$', 
  fontsize=15)
g1s_Ax.set_ylim([0,0.6])
g1s_Ax.xaxis.set_minor_locator(AutoMinorLocator())
g1s_Ax.yaxis.set_minor_locator(AutoMinorLocator())
g1s_Ax.tick_params(which='major', length=6)
g1s_Ax.tick_params(which='minor', length=3)

g2s_Ax = g_s.add_subplot(312)
g2s_Ax.plot(time, alignment[0].g2_s1[:], 'ko-', linewidth=1.5)
g2s_Ax.plot(time, alignment[0].g2_s2[:], 'bo-', linewidth=1.5)
g2s_Ax.plot(time, alignment[0].g2_s3[:], 'ro-', linewidth=1.5)

# middle
g2s_Ax.set_ylabel(r'$\cos \theta = (\mathbf{g}_2 \cdot \mathbf{s(0)}_i)$', 
  fontsize=15)
g2s_Ax.set_ylim([0,0.6])
g2s_Ax.xaxis.set_minor_locator(AutoMinorLocator())
g2s_Ax.yaxis.set_minor_locator(AutoMinorLocator())
g2s_Ax.tick_params(which='major', length=6)
g2s_Ax.tick_params(which='minor', length=3)

g3s_Ax = g_s.add_subplot(313)
g3s_Ax.plot(time, alignment[0].g3_s1[:], 'ko-', linewidth=1.5)
g3s_Ax.plot(time, alignment[0].g3_s2[:], 'bo-', linewidth=1.5)
g3s_Ax.plot(time, alignment[0].g3_s3[:], 'ro-', linewidth=1.5)

# minor
g3s_Ax.set_ylabel(r'$\cos \theta = (\mathbf{g}_3 \cdot \mathbf{s(0)}_i)$', 
  fontsize=15)
g3s_Ax.set_ylim([0,0.6])
g3s_Ax.set_xlabel('Time [ms]')
g3s_Ax.legend(['Major Strain Axis', 'Middle Strain Axis', 'Minor Strain Axis'])
g3s_Ax.xaxis.set_minor_locator(AutoMinorLocator())
g3s_Ax.yaxis.set_minor_locator(AutoMinorLocator())
g3s_Ax.tick_params(which='major', length=6)
g3s_Ax.tick_params(which='minor', length=3)

## vorticity / paxes, strain ##
wFig = plt.figure(figsize=(12,8))
wFig.suptitle('Alignment of shape and strain with vorticity', fontsize=16)

# principal axes
gw_Ax = wFig.add_subplot(311)
gw_Ax.plot(time, alignment[0].g1_w[:], 'ko-', linewidth=1.5)
gw_Ax.plot(time, alignment[0].g2_w[:], 'bo-', linewidth=1.5)
gw_Ax.plot(time, alignment[0].g3_w[:], 'ro-', linewidth=1.5)

gw_Ax.set_ylabel(r'$\cos \theta = (\mathbf{g(0)}_i \cdot \mathbf{\omega})$', 
  fontsize=15)
gw_Ax.legend(['Major Shape Axis', 'Middle Shape Axis', 'Minor Shape Axis'])
gw_Ax.xaxis.set_minor_locator(AutoMinorLocator())
gw_Ax.yaxis.set_minor_locator(AutoMinorLocator())
gw_Ax.tick_params(which='major', length=6)
gw_Ax.tick_params(which='minor', length=3)

# strain
sw_Ax = wFig.add_subplot(312)
sw_Ax.plot(time, alignment[0].s1_w[:], 'ko-', linewidth=1.5)
sw_Ax.plot(time, alignment[0].s2_w[:], 'bo-', linewidth=1.5)
sw_Ax.plot(time, alignment[0].s3_w[:], 'ro-', linewidth=1.5)
sw_Ax.legend(['Major Strain Axis', 'Middle Strain Axis', 'Minor Strain Axis'])

sw_Ax.set_ylabel(r'$\cos \theta = (\mathbf{s(0)}_i \cdot \mathbf{\omega})$', 
  fontsize=15)
sw_Ax.xaxis.set_minor_locator(AutoMinorLocator())
sw_Ax.yaxis.set_minor_locator(AutoMinorLocator())
sw_Ax.tick_params(which='major', length=6)
sw_Ax.tick_params(which='minor', length=3)

# vorticity magnitude
w_Ax = wFig.add_subplot(313)

w_Ax.plot(time, np.mean(vortNorm,axis=0), 'ko-', linewidth=1.5)

w_Ax.set_xlabel('Time [ms]')
w_Ax.set_ylabel('Vorticity')
w_Ax.xaxis.set_minor_locator(AutoMinorLocator())
w_Ax.yaxis.set_minor_locator(AutoMinorLocator())
w_Ax.tick_params(which='major', length=6)
w_Ax.tick_params(which='minor', length=3)

## alignment with z ##
g_z = plt.figure(figsize=(12,8))
g_z.suptitle('Alignment of shape, strain, and vorticity with gravity', 
  fontsize=16)

gzAx = g_z.add_subplot(311)
gzAx.errorbar(time, alignment[0].g1_z[:], yerr=alignment[0].err_g1_z[:], 
  fmt='ko-', linewidth=1.5)
gzAx.errorbar(time, alignment[0].g2_z[:], yerr=alignment[0].err_g2_z[:], 
  fmt='bo-', linewidth=1.5)
gzAx.errorbar(time, alignment[0].g3_z[:], yerr=alignment[0].err_g3_z[:], 
  fmt='ro-', linewidth=1.5)

gzAx.set_ylabel(r'$\cos \theta = (\mathbf{g}_i \cdot \mathbf{z})$', fontsize=15)
gzAx.legend(['Major Shape Axis', 'Middle Shape Axis', 'Minor Shape Axis'])
gzAx.xaxis.set_minor_locator(AutoMinorLocator())
gzAx.yaxis.set_minor_locator(AutoMinorLocator())
gzAx.tick_params(which='major', length=6)
gzAx.tick_params(which='minor', length=3)

szAx = g_z.add_subplot(312)
szAx.plot(time, alignment[0].s1_z[:], 'ko-', linewidth=1.5)
szAx.plot(time, alignment[0].s2_z[:], 'bo-', linewidth=1.5)
szAx.plot(time, alignment[0].s3_z[:], 'ro-', linewidth=1.5)

szAx.set_ylabel(r'$\cos \theta = (\mathbf{s}_i \cdot \mathbf{z})$', fontsize=15)
szAx.legend(['Major Strain Axis', 'Middle Strain Axis', 'Minor Strain Axis'])
szAx.xaxis.set_minor_locator(AutoMinorLocator())
szAx.yaxis.set_minor_locator(AutoMinorLocator())
szAx.tick_params(which='major', length=6)
szAx.tick_params(which='minor', length=3)

wzAx = g_z.add_subplot(313)
wzAx.plot(time, alignment[0].w_z[:], 'ko-', linewidth=1.5)
wzAx.set_xlabel('Time')
wzAx.set_ylabel(r'$\cos \theta = (\mathbf{\omega} \cdot \mathbf{z})$', fontsize=15)
wzAx.xaxis.set_minor_locator(AutoMinorLocator())
wzAx.yaxis.set_minor_locator(AutoMinorLocator())
wzAx.tick_params(which='major', length=6)
wzAx.tick_params(which='minor', length=3)

plt.show()

