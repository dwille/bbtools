#!/usr/bin/env python

# PURPOSE
#   Calculate the autocorrelation of particle velocity from the supplied time
#   series.
# USAGE
#   ./single-particle-vel-autocorrelation.py </path/to/sim/output> <start_time>
#
# OUTPUT
#   </path/to/sim/analysis/...XXX

# Imports:
import sys
import numpy as np
from scipy import signal
import bluebottle_particle_reader as bbparts

##########

# Parse output directory from command line
if len(sys.argv) >= 2:    # output directory given
  data_dir = sys.argv[1]

  if len(sys.argv) >= 3:  # start time given
    t_start = float(sys.argv[2])

else:                     # nothing given
  print("Error: Invalid commandline arguments.")
  print("Usage: ")
  print("   ./single-particle-vel-autocorrelation.py <./path/to/sim/output> <start_time>")
  print(" or")
  print("   ./single-particle-vel-autocorrelation.py <./path/to/sim/output>")
  sys.exit()

# Initialize the reader
times = bbparts.init(data_dir)
time_flt = np.array([tt for tt in times], dtype='float')
ind = np.argwhere(time_flt >= t_start)[0][0]
times = times[ind:]
time_flt = time_flt[ind:] - time_flt[ind]

print("Running from t = %.2f to t = %.2f" % (time_flt[0], time_flt[-1]))

# Get number of particles
bbparts.open(times[0])
nparts = bbparts.read_nparts()
bbparts.close()

# Initialize data arrays for particle data
u_total = np.zeros((len(times), nparts))
u_verti = np.zeros((len(times), nparts))
u_horiz = np.zeros((len(times), nparts))
t = np.zeros((len(times), nparts))

# Loop over and pull data
print("Reading in cgns files...")
for tt,time in enumerate(times):
  print(" file %d of %d\r" % (tt+1, len(times)), end="")
  bbparts.open(time)

  t[tt] = bbparts.read_time()
  (u, v, w) = bbparts.read_part_velocity()
  u_total[tt,:] = np.sqrt(u*u + v*v + w*w)
  u_verti[tt,:] = w
  u_horiz[tt,:] = np.sqrt(u*u + v*v)

  bbparts.close()

print() # for pretty output

# Perform autocorrelation, put result back in array
print("Performing autocorrelation...")
for n in np.arange(nparts):
  print(" part %d of %d\r" % (n+1, nparts), end="")
  
  ## u_total ##
  # Subtract mean from velocity
  temp = u_total[:,n] - np.mean(u_total[:,n]) 

  # Convolve with itself
  result = signal.fftconvolve(temp[::-1], temp, mode="full")

  # Place back into array
  result = result[::-1]                 # Flip result

  # Take positive freqs
  pstart = int(np.floor(len(result)/2.))
  u_total[:,n] = result[pstart:]

  # Normalize
  u_total[:,n] /= u_total[0,n]

  ## u_verti ##
  temp = u_verti[:,n] - np.mean(u_verti[:,n]) 
  result = signal.fftconvolve(temp[::-1], temp, mode="full")
  result = result[::-1]                 # Flip result
  pstart = int(np.floor(len(result)/2.))
  u_verti[:,n] = result[pstart:] # Take positive freq
  u_verti[:,n] /= u_verti[0,n]

  ## u_horiz ##
  temp = u_horiz[:,n] - np.mean(u_horiz[:,n]) 
  result = signal.fftconvolve(temp[::-1], temp, mode="full")
  result = result[::-1]                 # Flip result
  pstart = int(np.floor(len(result)/2.))
  u_horiz[:,n] = result[pstart:] # Take positive freq
  u_horiz[:,n] /= u_horiz[0,n]

print() # for pretty output

# Average results over particles
u_total = np.mean(u_total, axis=1)
u_verti = np.mean(u_verti, axis=1)
u_horiz = np.mean(u_horiz, axis=1)

# Save results to file
#import csv
# Want:
# time,utotal,uverti,uhoriz
#savefile = "./tmp-save-file.csv"
#with open(savefile, 'wb') as outfile:
#  out = csv.writer(outfile, delimiter=' ')
#  out.write

# plot temporarily here
import matplotlib.pyplot as plt

a = 2.1 # mm
rho = 3.3
nu = 0.01715
tau_p = (2.*a)**2 * rho / (18. * nu)

fig = plt.figure()
plt.semilogx(time_flt/tau_p, u_total)
plt.semilogx(time_flt/tau_p, u_verti)
plt.semilogx(time_flt/tau_p, u_horiz)
plt.legend(["total", "verti", "horiz"])
plt.xlabel(r"$t/\tau_p$")
plt.ylabel("xcorr")
#plt.xlim(xmin=0)
plt.tight_layout()
plt.savefig("tmp-xcorr.png", bbox_inches='tight', format='png')
   





