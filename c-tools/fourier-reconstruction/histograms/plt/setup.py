#!/usr/bin/env python2
import sys, os, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Process command line arguments
def simParams(sys):
  # 0 = program name
  # 1 = simulation directory

  # Parse command line args
  if len(sys.argv) > 1:
    simdir = sys.argv[1]
  else:
    simdir = raw_input("      Simulation directory: ")

  if not simdir.endswith('/'):
    simdir = simdir + '/'

  return simdir

# Setup up directory paths
def directoryStructureMarcc(simdir):
  home = os.path.expanduser("~")
  root = home + "/scratch/triply_per/"
  simdir = simdir + 'analysis/fourier-reconstruction/histograms/'
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

def findCenters(binStart, binEnd, nBins, dBin):
  centers = np.zeros(nBins)
  for i in np.arange(0, nBins):
    centers[i] = binStart + (i+0.5)*dBin
  return centers
