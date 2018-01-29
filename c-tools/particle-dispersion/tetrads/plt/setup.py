#!/usr/bin/env python2

# TODO: add partR, root, rho, dom.{X,Y,Z} to a config file

# Load common modules
import sys, os, glob, re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.ticker import AutoMinorLocator
from matplotlib import lines as mlines
import matplotlib.ticker as tick
from matplotlib.ticker import MultipleLocator

# Function Definition
# (partR, simdir, tstart) = simParams(sys)
#   Input:
#     sys: Command-line input args
#   Output:
#     partR: particle radius.
#     simdir: simulation directory in-between "root" and "analysis"
#     tstart: desired starting time to run data
def simParams(sys):
  # Parameters
  partR = 2.1 # XXX

  # Parse command line args
  if len(sys.argv) > 2:
    simdir = sys.argv[1]
    tstart = float(sys.argv[2]) / 1000
  else:
    simdir = raw_input("      Simulation directory: ")
    tstart = float(raw_input("      Starting time [ms]: ")) / 1000
    # TODO if tstart is -1 or empty, choose statsimtime

  if not simdir.endswith('/'):
    simdir = simdir + '/'

  return (partR, simdir, tstart)

# Function Definition
# (root, simdir, datadir, imgdir) = directoryStructure(simdir)
#   Input:
#     simdir: simulation directory as defined above
#   Output:
#     root:
#     simdir: simulation directory with analysis appended to it
#     datadir: directory within simdir that data is stored in 
#     imgdir: directory within simdir to store images in
def directoryStructure(simdir):
  home = os.path.expanduser("~")
  root = home + "/scratch/triply_per/"
  simdir = simdir + 'analysis/tetrads/'
  datadir = root + simdir + "data/"

  # Check if datadir exists so we don't go creating extra dirs
  if not os.path.exists(datadir):
    print "      " + datadir + " does not exist. Exiting..."
    print ""
    sys.exit()

  # Create imgdir if necessary
  imgdir = root + simdir + "/img/"
  if not os.path.exists(imgdir):
    os.makedirs(imgdir)

  return (root, simdir, datadir, imgdir)

# Initialize time and z data from simulation
# info.dat contains nTetrads and nTsteps
def initData(datadir, tstart):
  infoFile = datadir + "info.dat"
  # Time -- convert to secs
  # time = np.genfromtxt(infoFile, skip_footer=1)[1:] / 1000
  # tsInd = np.squeeze(np.argwhere(time >= tstart)[0])
  # print "      Starting time set to: %.3f [s]" % time[tsInd]
  # time = time[tsInd:] - time[tsInd]
  # nt = np.size(time)

  # Z-Locations
  # evalZ = np.genfromtxt(infoFile, skip_header=1)[1:]
  # nz = np.size(evalZ)

  # return (time, tsInd, nt, evalZ, nz)

# Print simulation data
def printSimulationData(partR, root, simdir, datadir):
  print "      Particle Radius set to: %.2f [mm]\n" % partR
  print "      Sim root directory set to: " + root
  print "      Sim directory set to: " + simdir
  print "      Data directory set to: " + datadir

## Define structure class
class structtype():
  pass

## Get file length
def file_len(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i    # bc don't count header line


## Interpolate
def interp(time, allTime, array, nRuns):
  for rr in np.arange(0,nRuns):
    array[rr].mean = np.interp(time, allTime[:,rr], array[rr].mean)
    array[rr].sdev = np.interp(time, allTime[:,rr], array[rr].sdev)
    array[rr].skew = np.interp(time, allTime[:,rr], array[rr].skew)
    array[rr].kurt = np.interp(time, allTime[:,rr], array[rr].kurt)

  return array

## Overall moments
def stats(data, minTsteps, nTetrads, nRuns):
  mean = np.zeros(minTsteps)
  var = np.zeros(minTsteps)
  skew = np.zeros(minTsteps)
  kurt = np.zeros(minTsteps)

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
  kurt = np.zeros(minTsteps)
  for rr in np.arange(0,nRuns):
    Nr = nTetrads[rr]
    diff = data[rr].mean - mean
    kurt += Nr*np.power(data[rr].sdev, 4.)*data[rr].kurt 
    + 4.*Nr*np.power(data[rr].sdev, 3.)*data[rr].skew*diff
    + 6.*Nr*np.power(data[rr].sdev, 2.)*diff*diff
    + diff*diff*diff*diff
  kurt /= np.sum(nTetrads)
  kurt /= np.power(sdev, 4.)

  moments = np.zeros((minTsteps,4))
  moments[:,0] = mean
  moments[:,1] = sdev
  moments[:,2] = skew
  moments[:,3] = kurt
  return moments

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
