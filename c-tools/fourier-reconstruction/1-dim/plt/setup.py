#!/usr/bin/env python2
import sys, os, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import scipy.fftpack as scifft

# Define simulation parameters, get simulation directory and startin time
def simParams(sys):
  # Parameters
  partR = 2.1

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

# Setup up directory paths for lucan
def directoryStructureDevel(simdir):
  root = "/home/dwille/bbtools/c-tools/fourier-reconstruction/1-dim/"
  simdir = "sim/"
  datadir = root + simdir + "data/reconstruct-1D/"

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


# Setup up directory paths for marcc
def directoryStructureMarcc(simdir):
  root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
  simdir = simdir + 'f-rec-1D/'
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
def initData(datadir, tstart):
  infoFile = datadir + "info"
  # Time -- convert to secs
  time = np.genfromtxt(infoFile, skip_footer=1)[1:] / 1000
  tsInd = np.argwhere(time >= tstart)[0]
  print "      Starting time set to: %.3f [s]" % time[tsInd]
  time = time[tsInd:] - time[tsInd]
  nt = np.size(time)

  # Z-Locations
  evalZ = np.genfromtxt(infoFile, skip_header=1)[1:]
  nz = np.size(evalZ)

  return (time, tsInd, nt, evalZ, nz)

# Print simulation data
def printSimulationData(partR, root, simdir, datadir):
  print "      Particle Radius set to: %.2f [mm]\n" % partR
  print "      Sim root directory set to: " + root
  print "      Sim directory set to: " + simdir
  print "      Data directory set to: " + datadir

# FFT Autocorrelation
def AutoCorrelationFFT(x1):
  x1 = np.asarray(x1)
  y1 = x1 - x1.mean()
  result = signal.fftconvolve(y1[::-1],y1,mode="full")
  # Reflip array
  result = result[::-1]
  result = result[len(result)/2:]
  result /= result[0]
  return result

def CrossCorrelationFFT(x1,x2):
  #x1 = np.asarray(x1)
  y1 = x1 - x1.mean()
  #x2 = np.asarray(x2[::-1])
  y2 = x2 - x2.mean()
  result = signal.fftconvolve(y2[::-1],y1,mode="full")
  # Re-flip result array 
  result = result[::-1]
  result = result[np.size(result)/2:]
  return result

# stackoverflow 16044491 statistical scaling of autocorrelation using numpy.fft
def AutoCorrelation(x):
  x = np.asarray(x)
  y = x - x.mean()
  result = np.correlate(y,y,mode="full")
  result = result[len(result)/2:]
  result /= result[0]
  return result

# stackoverflow 16044491 statistical scaling of autocorrelation using numpy.fft
def CrossCorrelation(x1,x2):
  x1 = np.asarray(x1)
  y1 = x1 - x1.mean()
  x2 = np.asarray(x2)
  y2 = x2 - x2.mean()
  result = np.correlate(y2,y1,mode="full")
  result = result[len(result)/2:]
  return result


