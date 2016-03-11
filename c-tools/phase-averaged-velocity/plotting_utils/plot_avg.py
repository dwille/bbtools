#!/usr/bin/env python2

import matplotlib.pyplot as plt
import numpy as np

# From stackoverflow 8956832, User Joe Kington
def iter_loadtxt(filename, delimiter=' ', skiprows=1, dtype=float):
  def iter_func():
    with open(filename, 'r') as infile:
      for _ in range(skiprows):
        next(infile)
      for line in infile:
        line = line.rstrip().split(delimiter)

        # Could possibly only take certain columns here?
        for item in line:
          yield dtype(item)
    iter_loadtxt.rowlength = len(line)

  data = np.fromiter(iter_func(), dtype=dtype)
  data = data.reshape((-1, iter_loadtxt.rowlength))
  return data

# Test -- time, u, v, w
data = iter_loadtxt("../sim/data/phaseAveragedVel")

# Plot
phaseVel = plt.figure()
phaseVel.suptitle('Phase Averaged Fluid Velocity', fontsize=20)

maxY = 1.1*np.max(data[:,3])
minY = np.min(data[:,3]) - 0.1*maxY

# u
uAx = phaseVel.add_subplot(311)
uAx.plot(data[:,0], data[:,1], 'ko-')
uAx.set_ylim([minY, maxY])
uAx.set_ylabel('u [m/s]')

# v
vAx = phaseVel.add_subplot(312)
vAx.plot(data[:,0], data[:,2], 'ko-')
vAx.set_ylim([minY, maxY])
vAx.set_ylabel('u [m/s]')

# w
wAx = phaseVel.add_subplot(313)
wAx.plot(data[:,0], data[:,3], 'ko-')
wAx.set_ylim([minY, maxY])
wAx.set_ylabel('u [m/s]')
wAx.set_xlabel('Time [ms]')

plt.show()
