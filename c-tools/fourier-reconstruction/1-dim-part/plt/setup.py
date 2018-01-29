#!/usr/bin/env python
import sys, os, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from scipy import signal
import scipy.fftpack as scifft
import matplotlib.ticker as tick

# Define simulation parameters
#  Input: sys -- commandline args
#  Output: partR -- particle radius
#          nparts -- particles in simulation
#          rho_star -- particle-to-fluid density ratio
#          vFracMean -- mean volume fraction
#          simdir -- input simulation directory
#          tstart -- starting time for analysis (useful for transients)
def simParams(sys):
  # TODO better job of which parameters needs to be set
  # TODO better way to return these values, i.e. call global -- maybe settings?
  # User Input Parameters
  partR = 2.1 

  # Parse command line arguments
  if len(sys.argv) > 2:
    simdir = sys.argv[1]
    tstart = float(sys.argv[2]) / 1000 # CONVERT FROM MS TO S!!!
  else:
    simdir = raw_input("      Simulation directory: ")
    tstart = float(raw_input("      Starting time [ms]: ")) / 1000
    # TODO if tstart is -1 or empty, choose statsimtime

  if not simdir.endswith('/'):
    simdir = simdir + '/'

  nparts = int(simdir.split('/')[-3])
  rho = float(simdir.split('/')[-2][-3:])
  vFracMean = nparts*(4./3.)*np.pi*(partR**3.)/(42.*42.*126.)

  return (partR, nparts, rho, vFracMean, simdir, tstart)

# Setup up directory paths
def directoryStructure(simdir):
  home = os.path.expanduser("~")
  root = home + "/scratch/"
  simdir = simdir + 'analysis/fourier-reconstruction/1-dim-part/'
  datadir = root + simdir + "data/"

  # Check if datadir exists so we don't go creating extra dirs
  if not os.path.exists(datadir):
    print("      " + datadir + " does not exist. Exiting...")
    print("")
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
  tsInd = np.squeeze(np.argwhere(time >= tstart)[0])
  print("      Starting time set to: %.3f [s]" % time[tsInd])
  time = time[tsInd:] - time[tsInd]
  nt = np.size(time)

  # Z-Locations
  evalZ = np.genfromtxt(infoFile, skip_header=1)[1:]
  nz = np.size(evalZ)

  return (time, tsInd, nt, evalZ, nz)

# Print simulation data
def printSimulationData(partR, root, simdir, datadir):
  print("      Particle Radius set to: %.2f [mm]\n" % partR)
  print("      Sim root directory set to: " + root)
  print("      Sim directory set to: " + simdir)
  print("      Data directory set to: " + datadir)

# FFT Autocorrelation
def AutoCorrelationFFT(x1):
  x1 = np.asarray(x1)
  y1 = x1 - x1.mean()
  result = signal.fftconvolve(y1[::-1],y1,mode="full")
  # Reflip array
  result = result[::-1]
  result = result[int(np.size(result)/2):]
  result /= result[0]
  return result

def AutoCorrelationSpaceFFT(arrayIn):
  arrayFlucts = arrayIn - np.mean(arrayIn)
  arrayFFT1 = np.fft.fft(arrayFlucts)
  arrayFFT2 = np.fft.fft(arrayFlucts[::-1])
  arrayOut = np.fft.ifft(arrayFFT1*arrayFFT2)

  arrayOut /= arrayOut[0]
  powerSpec = np.absolute(arrayFFT1)**2

  return (arrayOut, powerSpec)

  #arrayOutReal = np.real(arrayOut)
  #arrayOutReal /= arrayOutReal[0]
  #arrayOutImag = np.imag(arrayOut)

  ##powerSpec = np.absolute(np.fft.fft(arrayIn))**2
  #powerSpec = np.absolute(arrayFFT1)**2

  #return (arrayOutReal, arrayOutImag, powerSpec)


def CrossCorrelationFFT(x1,x2):
  #x1 = np.asarray(x1)
  y1 = x1 - x1.mean()
  #x2 = np.asarray(x2[::-1])
  y2 = x2 - x2.mean()
  result = signal.fftconvolve(y2[::-1],y1,mode="full")
  # Re-flip result array 
  result = result[::-1]
  print("len(result) = %d" % len(result))
  result = result[np.size(result)/2:]
  return result

# Same as CrossCorrelationFFT except with logic for odd/even
def xcorr_fft(x1,x2):
  y1 = x1 - x1.mean()
  y2 = x2 - x2.mean()
  result = signal.fftconvolve(y2[::-1],y1,mode="full")

  # Re-flip result array 
  result = result[::-1]
  print("len(result) = %d" % len(result))
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


# get axis limits (so 24125058)
# ax2.annotate(r"$(a)$",xy=get_axis_limits(ax2))
def get_axis_limits(ax, scale=0.85):
  xmin = ax.get_xlim()[0]
  xmax = ax.get_xlim()[1]
  ymin = ax.get_ylim()[0]
  ymax = ax.get_ylim()[1]

  dx = xmax - xmin
  dy = ymax - ymin

  return xmin + scale*dx, ymin + scale*dy

