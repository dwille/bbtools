#!/usr/bin/env python2
# Combine various runs of a tetrad analysis into master dat files

import sys, os, re
import glob
import csv
import numpy as np
import matplotlib.pyplot as plt

## Define structure class
class structtype():
  pass

## Get file length
def file_len(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i    # bc don't count header line

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

## Interpolate
def interp(time, allTime, array, nRuns):
  for rr in np.arange(0,nRuns):
    array[rr].mean = np.interp(time, allTime[:,rr], array[rr].mean)
    array[rr].sdev = np.interp(time, allTime[:,rr], array[rr].sdev)
    array[rr].skew = np.interp(time, allTime[:,rr], array[rr].skew)
    array[rr].kurt = np.interp(time, allTime[:,rr], array[rr].kurt)

  return array

## Overall moments
def stats(data):
  mean = np.zeros(maxTsteps)
  var = np.zeros(maxTsteps)
  skew = np.zeros(maxTsteps)
  kurt = np.zeros(maxTsteps)

  # Mean
  for rr in np.arange(0,nRuns):
    mean += nTetrads[rr] * data[rr].mean
  mean /= np.sum(nTetrads)

  # sdev
  for rr in np.arange(0,nRuns):
    Nr = nTetrads[rr]
    diff = data[rr].mean - mean
    var += (Nr - 1.) * np.square(data[rr].sdev) + Nr * diff * diff
  var /= (np.sum(nTetrads) - 1.)
  sdev = np.sqrt(var)

  # skew
  for rr in np.arange(0,nRuns):
    Nr = nTetrads[rr]
    diff = data[rr].mean - mean
    skew += Nr*np.power(data[rr].sdev, 3.)*data[rr].skew 
    + 3.*(Nr - 1.)*diff 
    + diff*diff*diff
  skew /= np.sum(nTetrads)
  skew /= np.power(sdev, 3.)

  # kurt
  kurt = np.zeros(maxTsteps)
  for rr in np.arange(0,nRuns):
    Nr = nTetrads[rr]
    diff = data[rr].mean - mean
    kurt += Nr*np.power(data[rr].sdev, 4.)*data[rr].kurt 
    + 4.*Nr*np.power(data[rr].sdev, 3.)*data[rr].skew*diff
    + 6.*Nr*np.power(data[rr].sdev, 2.)*diff*diff
    + diff*diff*diff*diff
  kurt /= np.sum(nTetrads)
  kurt /= np.power(sdev, 4.)

  moments = np.zeros((maxTsteps,4))
  moments[:,0] = mean
  moments[:,1] = sdev
  moments[:,2] = skew
  moments[:,3] = kurt
  return moments

print " ---- Combine-Runs Utility ----"
print ""

## Set root dir, simdir, and datadir
root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
simdir = raw_input("  Simulation directory: ")
#root = "/home/dwille/bbtools/c-tools/multiparticle_statistics/"
#simdir =  "sim/"

if not simdir.endswith('/'):
  simdir = simdir + '/'

print "      Sim root directory set to: " + root

## Get simdir and create datadir
datadir = root + simdir + "data-tetrads/"

# Get list of directories in the datadir -- corresponds to different runs
runs = sorted_nicely(glob.glob(datadir + "ts_*"))

# Initalize data structs and arrays
nRuns = int(len(runs))
global nTetrads 
nTetrads = np.zeros(nRuns)
global nTsteps 
nTsteps = np.zeros(nRuns)

RoG = [ structtype() for i in range(nRuns) ]
EVar = [ structtype() for i in range(nRuns) ]
Shape = [ structtype() for i in range(nRuns) ]
I1 = [ structtype() for i in range(nRuns) ]
I2 = [ structtype() for i in range(nRuns) ]
I3 = [ structtype() for i in range(nRuns) ]

# Loop over all runs, find nTsteps
for rr, run in enumerate(runs):
  infoFile = run + '/info.dat'

  nTetrads[rr] = np.genfromtxt(infoFile, skip_header=1, usecols=0)
  nTsteps[rr] = np.genfromtxt(infoFile, skip_header=1, usecols=1)

# find maxTsteps
global maxTsteps
maxTsteps = np.max(nTsteps)
allTime = np.zeros((maxTsteps, nRuns))

for rr, run in enumerate(runs):
  meanFile = run + '/stat.mean'
  sdevFile = run + '/stat.sdev'
  skewFile = run + '/stat.skew'
  kurtFile = run + '/stat.kurt'

  # init data to maxTstep size
  RoG[rr].mean = np.zeros(maxTsteps)
  RoG[rr].sdev = np.zeros(maxTsteps)
  RoG[rr].skew = np.zeros(maxTsteps)
  RoG[rr].kurt = np.zeros(maxTsteps)
  EVar[rr].mean = np.zeros(maxTsteps)
  EVar[rr].sdev = np.zeros(maxTsteps)
  EVar[rr].skew = np.zeros(maxTsteps)
  EVar[rr].kurt = np.zeros(maxTsteps)
  Shape[rr].mean = np.zeros(maxTsteps)
  Shape[rr].sdev = np.zeros(maxTsteps)
  Shape[rr].skew = np.zeros(maxTsteps)
  Shape[rr].kurt = np.zeros(maxTsteps)
  I1[rr].mean = np.zeros(maxTsteps)
  I1[rr].sdev = np.zeros(maxTsteps)
  I1[rr].skew = np.zeros(maxTsteps)
  I1[rr].kurt = np.zeros(maxTsteps)
  I2[rr].mean = np.zeros(maxTsteps)
  I2[rr].sdev = np.zeros(maxTsteps)
  I2[rr].skew = np.zeros(maxTsteps)
  I2[rr].kurt = np.zeros(maxTsteps)
  I3[rr].mean = np.zeros(maxTsteps)
  I3[rr].sdev = np.zeros(maxTsteps)
  I3[rr].skew = np.zeros(maxTsteps)
  I3[rr].kurt = np.zeros(maxTsteps)

  # only fill up to nTsteps[rr] -- this pads the end of the array with 0s
  nt = nTsteps[rr]
  tmpT = np.genfromtxt(meanFile, skip_header=1, usecols=0)
  allTime[0:nt, rr] = tmpT - tmpT[0]

  RoG[rr].mean[0:nt] = np.genfromtxt(meanFile, skip_header=1, usecols=1)
  EVar[rr].mean[0:nt] = np.genfromtxt(meanFile, skip_header=1, usecols=2)
  Shape[rr].mean[0:nt] = np.genfromtxt(meanFile, skip_header=1, usecols=3)
  I1[rr].mean[0:nt] = np.genfromtxt(meanFile, skip_header=1, usecols=4)
  I2[rr].mean[0:nt] = np.genfromtxt(meanFile, skip_header=1, usecols=5)
  I3[rr].mean[0:nt] = np.genfromtxt(meanFile, skip_header=1, usecols=6)

  RoG[rr].sdev[0:nt] = np.genfromtxt(sdevFile, skip_header=1, usecols=1)
  EVar[rr].sdev[0:nt] = np.genfromtxt(sdevFile, skip_header=1, usecols=2)
  Shape[rr].sdev[0:nt] = np.genfromtxt(sdevFile, skip_header=1, usecols=3)
  I1[rr].sdev[0:nt] = np.genfromtxt(sdevFile, skip_header=1, usecols=4)
  I2[rr].sdev[0:nt] = np.genfromtxt(sdevFile, skip_header=1, usecols=5)
  I3[rr].sdev[0:nt] = np.genfromtxt(sdevFile, skip_header=1, usecols=6)

  RoG[rr].skew[0:nt] = np.genfromtxt(skewFile, skip_header=1, usecols=1)
  EVar[rr].skew[0:nt] = np.genfromtxt(skewFile, skip_header=1, usecols=2)
  Shape[rr].skew[0:nt] = np.genfromtxt(skewFile, skip_header=1, usecols=3)
  I1[rr].skew[0:nt] = np.genfromtxt(skewFile, skip_header=1, usecols=4)
  I2[rr].skew[0:nt] = np.genfromtxt(skewFile, skip_header=1, usecols=5)
  I3[rr].skew[0:nt] = np.genfromtxt(skewFile, skip_header=1, usecols=6)

  RoG[rr].kurt[0:nt] = np.genfromtxt(kurtFile, skip_header=1, usecols=1)
  EVar[rr].kurt[0:nt] = np.genfromtxt(kurtFile, skip_header=1, usecols=2)
  Shape[rr].kurt[0:nt] = np.genfromtxt(kurtFile, skip_header=1, usecols=3)
  I1[rr].kurt[0:nt] = np.genfromtxt(kurtFile, skip_header=1, usecols=4)
  I2[rr].kurt[0:nt] = np.genfromtxt(kurtFile, skip_header=1, usecols=5)
  I3[rr].kurt[0:nt] = np.genfromtxt(kurtFile, skip_header=1, usecols=6)

totalTetrads = np.sum(nTetrads)
print '      Total tetrads tracked: ' + str(totalTetrads)

# Average timestep and corresponding time to interplate to
avgTime = np.sum(allTime, 1)
avgTime /= nRuns
avgDt =  np.round(np.mean(np.diff(avgTime)))
time = np.arange(0, avgDt*np.size(avgTime), avgDt) 

## Interpolate data
RoG = interp(time, allTime, RoG, nRuns)
EVar = interp(time, allTime, EVar, nRuns)
Shape = interp(time, allTime, Shape, nRuns)
I1 = interp(time, allTime, I1, nRuns)
I2 = interp(time, allTime, I2, nRuns)
I3 = interp(time, allTime, I3, nRuns)

## Find Moments
m_RoG = stats(RoG)
m_EVar = stats(EVar)
m_Shape = stats(Shape)
m_I1 = stats(I1)
m_I2 = stats(I2)
m_I3 = stats(I3)

# Print info to infofile
infofile = datadir + 'info.dat'
with open(infofile, 'wb') as outfile:
  a = csv.writer(outfile, delimiter=' ')
  headers = [['nTetrads', 'nTsteps']]
  a.writerows(headers)
  data = [[totalTetrads, maxTsteps]]
  a.writerows(data)
  #TODO: maybe should be a minTsteps as well?

# Print info to datadir -- mean
allMeanFile = datadir + 'stat.mean'
with open(allMeanFile, 'wb') as outfile:
  a = csv.writer(outfile, delimiter=' ')
  headers = [['time', 'RoG', 'EVar', 'Shape', 'I1', 'I2', 'I3']]
  a.writerows(headers)
  for tt in np.arange(0, maxTsteps):
    data = [[time[tt], m_RoG[tt,0], m_EVar[tt,0], m_Shape[tt,0],
      m_I1[tt,0], m_I1[tt,0], m_I1[tt,0]]]
    a.writerows(data)

# Print info to datadir -- sdev
allSdevFile = datadir + 'stat.sdev'
with open(allSdevFile, 'wb') as outfile:
  a = csv.writer(outfile, delimiter=' ')
  headers = [['time', 'RoG', 'EVar', 'Shape', 'I1', 'I2', 'I3']]
  a.writerows(headers)
  for tt in np.arange(0, maxTsteps):
    data = [[time[tt], m_RoG[tt,1], m_EVar[tt,1], m_Shape[tt,1],
      m_I1[tt,1], m_I1[tt,1], m_I1[tt,1]]]
    a.writerows(data)

# Print info to datadir -- skew
allSkewFile = datadir + 'stat.skew'
with open(allSkewFile, 'wb') as outfile:
  a = csv.writer(outfile, delimiter=' ')
  headers = [['time', 'RoG', 'EVar', 'Shape', 'I1', 'I2', 'I3']]
  a.writerows(headers)
  for tt in np.arange(0, maxTsteps):
    data = [[time[tt], m_RoG[tt,2], m_EVar[tt,2], m_Shape[tt,2],
      m_I1[tt,2], m_I1[tt,2], m_I1[tt,2]]]
    a.writerows(data)

# Print info to datadir -- kurt
allKurtFile = datadir + 'stat.kurt'
with open(allKurtFile, 'wb') as outfile:
  a = csv.writer(outfile, delimiter=' ')
  headers = [['time', 'RoG', 'EVar', 'Shape', 'I1', 'I2', 'I3']]
  a.writerows(headers)
  for tt in np.arange(0, maxTsteps):
    data = [[time[tt], m_RoG[tt,3], m_EVar[tt,3], m_Shape[tt,3],
      m_I1[tt,3], m_I1[tt,3], m_I1[tt,3]]]
    a.writerows(data)

# # R
# R = plt.figure(figsize=(12,8))
# m1_R = R.add_subplot(221)
# plt.plot(time, m_RoG[:,0], 'k', linewidth=2)
# for rr in np.arange(0,nRuns):
#   plt.plot(allTime[:,rr], RoG[rr].mean)
# m1_R.set_title('Mean R')
# 
# m2_R = R.add_subplot(222)
# plt.plot(time, m_RoG[:,1], 'k', linewidth=2)
# for rr in np.arange(0,nRuns):
#   plt.plot(allTime[:,rr], RoG[rr].sdev)
# m2_R.set_title('sdev R')
# 
# m3_R = R.add_subplot(223)
# plt.plot(time, m_RoG[:,2], 'k', linewidth=2)
# for rr in np.arange(0,nRuns):
#   plt.plot(allTime[:,rr], RoG[rr].skew)
# m3_R.set_title('skew R')
# 
# m4_R = R.add_subplot(224)
# plt.plot(time, m_RoG[:,3], 'k', linewidth=2)
# for rr in np.arange(0,nRuns):
#   plt.plot(allTime[:,rr], RoG[rr].kurt)
# m4_R.set_title('kurt R')
# 
# # EVar
# EVar_fig = plt.figure(figsize=(12,8))
# m1_EVar = EVar_fig.add_subplot(221)
# plt.plot(time, m_EVar[:,0], 'k', linewidth=2)
# for rr in np.arange(0,nRuns):
#   plt.plot(allTime[:,rr], EVar[rr].mean)
# m1_EVar.set_title('Mean EVar')
# 
# m2_EVar = EVar_fig.add_subplot(222)
# plt.plot(time, m_EVar[:,1], 'k', linewidth=2)
# for rr in np.arange(0,nRuns):
#   plt.plot(allTime[:,rr], EVar[rr].sdev)
# m2_EVar.set_title('sdev EVar')
# 
# m3_EVar = EVar_fig.add_subplot(223)
# plt.plot(time, m_EVar[:,2], 'k', linewidth=2)
# for rr in np.arange(0,nRuns):
#   plt.plot(allTime[:,rr], EVar[rr].skew)
# m3_EVar.set_title('skew EVar')
# 
# m4_EVar = EVar_fig.add_subplot(224)
# plt.plot(time, m_EVar[:,3], 'k', linewidth=2)
# for rr in np.arange(0,nRuns):
#   plt.plot(allTime[:,rr], EVar[rr].kurt)
# m4_EVar.set_title('kurt EVar')
# 
# # Shape
# Shape_fig = plt.figure(figsize=(12,8))
# m1_Shape = Shape_fig.add_subplot(221)
# plt.plot(time, m_Shape[:,0], 'k', linewidth=2)
# for rr in np.arange(0,nRuns):
#   plt.plot(allTime[:,rr], Shape[rr].mean)
# m1_Shape.set_title('Mean Shape')
# 
# m2_Shape = Shape_fig.add_subplot(222)
# plt.plot(time, m_Shape[:,1], 'k', linewidth=2)
# for rr in np.arange(0,nRuns):
#   plt.plot(allTime[:,rr], Shape[rr].sdev)
# m2_Shape.set_title('sdev Shape')
# 
# m3_Shape = Shape_fig.add_subplot(223)
# plt.plot(time, m_Shape[:,2], 'k', linewidth=2)
# for rr in np.arange(0,nRuns):
#   plt.plot(allTime[:,rr], Shape[rr].skew)
# m3_Shape.set_title('skew Shape')
# 
# m4_Shape = Shape_fig.add_subplot(224)
# plt.plot(time, m_Shape[:,3], 'k', linewidth=2)
# for rr in np.arange(0,nRuns):
#   plt.plot(allTime[:,rr], Shape[rr].kurt)
# m4_Shape.set_title('kurt Shape')

# I factors
# iFig = plt.figure(figsize=(12,8))
# m1_I = iFig.add_subplot(221)
# plt.plot(time, m_I1[:,0], 'k', linewidth=2)
# plt.plot(time, m_I2[:,0], 'r', linewidth=2)
# plt.plot(time, m_I3[:,0], 'b', linewidth=2)
# for rr in np.arange(0,nRuns):
#   plt.plot(allTime[:,rr], I1[rr].mean, 'k--')
#   plt.plot(allTime[:,rr], I2[rr].mean, 'r--')
#   plt.plot(allTime[:,rr], I3[rr].mean, 'b--')
# m1_I.set_title('MeanI')
# 
# m2_I = iFig.add_subplot(222)
# plt.plot(time, m_I1[:,1], 'k', linewidth=2)
# plt.plot(time, m_I2[:,1], 'r', linewidth=2)
# plt.plot(time, m_I3[:,1], 'b', linewidth=2)
# for rr in np.arange(0,nRuns):
#   plt.plot(allTime[:,rr], I1[rr].sdev, 'k--')
#   plt.plot(allTime[:,rr], I2[rr].sdev, 'r--')
#   plt.plot(allTime[:,rr], I3[rr].sdev, 'b--')
# m2_I.set_title('sdevI')
# 
# m3_I = iFig.add_subplot(223)
# plt.plot(time, m_I1[:,2], 'k', linewidth=2)
# plt.plot(time, m_I2[:,2], 'r', linewidth=2)
# plt.plot(time, m_I3[:,2], 'b', linewidth=2)
# for rr in np.arange(0,nRuns):
#   plt.plot(allTime[:,rr], I1[rr].skew, 'k--')
#   plt.plot(allTime[:,rr], I2[rr].skew, 'r--')
#   plt.plot(allTime[:,rr], I3[rr].skew, 'b--')
# m3_I.set_title('skewI')
# 
# m4_I = iFig.add_subplot(224)
# plt.plot(time, m_I1[:,3], 'k', linewidth=2)
# plt.plot(time, m_I2[:,3], 'r', linewidth=2)
# plt.plot(time, m_I3[:,3], 'b', linewidth=2)
# for rr in np.arange(0,nRuns):
#   plt.plot(allTime[:,rr], I1[rr].kurt, 'k--')
#   plt.plot(allTime[:,rr], I2[rr].kurt, 'r--')
#   plt.plot(allTime[:,rr], I3[rr].kurt, 'b--')
# m4_I.set_title('kurtI')
# 
# plt.show()
