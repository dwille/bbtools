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
  l -- The iterable to be sorted.

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

print ""
print "---- Combine-Runs Utility -- Alignment ----"
print ""

# Set root dir
#root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
root = "/home/dwille/bbtools/c-tools/multiparticle_statistics/"
print "    Sim root directory set to: " + root

# Get simdir and create datadir
#simdir = raw_input("  Simulation directory: ")
simdir =  "sim/"
#if not simdir.endswith('/'):
#  simdir = simdir + '/'
datadir = root + simdir + "data-tetrads/"

# Get list of directories in the datadir -- corresponds to different runs
runs = sorted_nicely(glob.glob(datadir + "ts_*"))

# Initalize data structs and arrays
nRuns = int(len(runs))
global nTetrads 
nTetrads = np.zeros(nRuns)
global nTsteps 
nTsteps = np.zeros(nRuns)
g1s1 = [ structtype() for i in range(nRuns) ]

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
  meanFile = run + '/align.mean'
  sdevFile = run + '/align.sdev'
  skewFile = run + '/align.skew'
  kurtFile = run + '/align.kurt'

  # init data to maxTstep size
  g1s1[rr].mean = np.zeros(maxTsteps)
  g1s1[rr].sdev = np.zeros(maxTsteps)
  g1s1[rr].skew = np.zeros(maxTsteps)
  g1s1[rr].kurt = np.zeros(maxTsteps)

  # only fill up to nTsteps[rr] -- this pads the end of the array with 0s
  nt = nTsteps[rr]
  tmpT = np.genfromtxt(meanFile, skip_header=1, usecols=0)
  allTime[0:nt, rr] = tmpT - tmpT[0]

  g1s1[rr].mean[0:nt] = np.genfromtxt(meanFile, skip_header=1, usecols=1)
  g1s1[rr].sdev[0:nt] = np.genfromtxt(sdevFile, skip_header=1, usecols=1)
  g1s1[rr].skew[0:nt] = np.genfromtxt(skewFile, skip_header=1, usecols=1)
  g1s1[rr].kurt[0:nt] = np.genfromtxt(kurtFile, skip_header=1, usecols=1)

totalTetrads = np.sum(nTetrads)
print '      Total tetrads tracked: ' + str(totalTetrads)

# Average timestep and corresponding time to interplate to
avgTime = np.sum(allTime, 1)
avgTime /= nRuns
avgDt =  np.round(np.mean(np.diff(avgTime)))
time = np.arange(0, avgDt*np.size(avgTime), avgDt) 

## Interpolate data
g1s1 = interp(time, allTime, g1s1, nRuns)

## Find moments
m_g1s1 = stats(g1s1)

# R
R = plt.figure(figsize=(12,8))
m1_R = R.add_subplot(221)
plt.plot(time, m_g1s1[:,0], 'k', linewidth=2)
for rr in np.arange(0,nRuns):
  plt.plot(allTime[:,rr], g1s1[rr].mean)
m1_R.set_title('Mean R')

m2_R = R.add_subplot(222)
plt.plot(time, m_g1s1[:,1], 'k', linewidth=2)
for rr in np.arange(0,nRuns):
  plt.plot(allTime[:,rr], g1s1[rr].sdev)
m2_R.set_title('sdev R')

m3_R = R.add_subplot(223)
plt.plot(time, m_g1s1[:,2], 'k', linewidth=2)
for rr in np.arange(0,nRuns):
  plt.plot(allTime[:,rr], g1s1[rr].skew)
m3_R.set_title('skew R')

m4_R = R.add_subplot(224)
plt.plot(time, m_g1s1[:,3], 'k', linewidth=2)
for rr in np.arange(0,nRuns):
  plt.plot(allTime[:,rr], g1s1[rr].kurt)
m4_R.set_title('kurt R')

plt.show()
