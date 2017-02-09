#!/usr/bin/env python2
import sys, os, glob, re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib.ticker import MultipleLocator

# Simulation Parameters
diam = 4.2    # mm
nu = .01715   # mm^2/ms
g = 0.00981   # mm/ms^2
rho_f = 8.75e-4   # g/mm^3
rho_p = np.array([0.001750, 0.00289, 0.0035, 0.004375]) # g/mm^3
rho_star = np.array([2., 3.3, 4., 5.])
Lx = 42.
Ly = 42.
Lz = 126.
vol = Lx * Ly * Lz;
nparts = np.array([500, 1000, 1500, 2000])
vfrac = (4./3. * np.pi * (0.5*diam)**3. * nparts) / (vol)

# Directory structure
home = os.path.expanduser("~")
root = home + "/scratch/collision/"
datadir = root + "/*/rho*/analysis/phase-averages/collision/data/collision_stats"
imgdir = root + "/simdata/img/colls/"

if not os.path.exists(imgdir):
  os.makedirs(imgdir)

## Struct class
class structtype():
  pass

