#!/usr/bin/env python2

import sys, os, csv
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as tick
from matplotlib.ticker import MultipleLocator
# from scipy import signal
# import scipy.fftpack as scifft
# matplotlib.use('Agg')      ## marcc savefig
# matplotlib.use(u'Qt4Agg')  ## marcc plt.show
# matplotlib.use(u'GTKAgg')  ## lucan plt.show

# Define simulation parameters, get simulation directory and startin time
def simParams(sys):
  # Parameters
  partR = 2.1 # XXX

  # Parse command line args
  if len(sys.argv) > 2:
    simdir = sys.argv[1]
    tstart = float(sys.argv[2]) / 1000. ## XXX convert from [ms] to [s]  
  else:
    simdir = raw_input("      Simulation directory: ")
    tstart = float(raw_input("      Starting time [ms]: ")) / 1000
    # TODO if tstart is -1 or empty, choose statsimtime

  if not simdir.endswith('/'):
    simdir = simdir + '/'

  return (partR, simdir, tstart)

# Setup up directory paths
def directoryStructure(simdir):
  home = os.path.expanduser("~")
  root = home + "/scratch/triply_per/"
  simdir = simdir + 'analysis/phase-averaged-velocity/'
  pdatadir = root + simdir + "part/data/"
  fdatadir = root + simdir + "flow/data/"

  # Check if datadir exists so we don't go creating extra dirs
  if not os.path.exists(pdatadir):
    print "      " + pdatadir + " does not exist. Exiting..."
    print ""
    sys.exit()
  elif not os.path.exists(fdatadir):
    print "      " + fdatadir + " does not exist. Exiting..."
    print ""
    sys.exit()

  # Create imgdir if necessary
  imgdir = root + simdir + "img/"
  if not os.path.exists(imgdir):
    os.makedirs(imgdir)

  return (root, simdir, pdatadir, fdatadir, imgdir)

# Print simulation data
def printSimulationData(partR, root, simdir, pdatadir, fdatadir):
  print "      Particle Radius set to: %.2f [mm]\n" % partR
  print "      Root directory set to: " + root
  print "      Sim directory set to: " + simdir
  print "      Flow data directory set to: " + fdatadir
  print "      Part data directory set to: " + pdatadir

