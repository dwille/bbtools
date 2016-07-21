#!/usr/bin/env python2
import sys, os, csv, glob, re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Define simulation parameters, get simulation directory and startin time
def simParams(sys):
  # Parameters
  partR = 2.1

  # Parse command line args
  if len(sys.argv) == 2:
    simdir = sys.argv[1]
  else:
    simdir = raw_input("      Simulation directory: ")

  if not simdir.endswith('/'):
    simdir = simdir + '/'

  return (partR, simdir)

# Setup up directory paths
def directoryStructure(simdir):
  home = os.path.expanduser("~")
  root = home + "/scratch/triply_per/"
  simdir = simdir + 'f-rec-part-phase-3D/'
  datadir = root + simdir + "track_data/"

  # Check if datadir exists so we don't go creating extra dirs
  if not os.path.exists(datadir):
    print "      " + datadir + " does not exist. Exiting..."
    print ""
    sys.exit()

  # Create imgdir if necessary
  imgdir = root + simdir + "/track_img/"
  if not os.path.exists(imgdir):
    os.makedirs(imgdir)

  return (root, simdir, datadir, imgdir)

# Initialize time and z data from simulation
def initData(datadir):
  infoFile = datadir + "info"
  # Time -- convert to secs
  time = np.genfromtxt(infoFile, skip_header=1, usecols=0) / 1000
  print "      Starting time set to: %.3f [s]" % time[0]
  nt = np.size(time)

  # Z-Locations
  evalZ = np.genfromtxt(infoFile, skip_header=1, usecols=1)
  nz = np.size(evalZ)

  return (time, nt, evalZ, nz)

# Print simulation data
def printSimulationData(partR, root, simdir, datadir):
  print "      Particle Radius set to: %.2f [mm]\n" % partR
  print "      Sim root directory set to: " + root
  print "      Sim directory set to: " + simdir
  print "      Data directory set to: " + datadir
