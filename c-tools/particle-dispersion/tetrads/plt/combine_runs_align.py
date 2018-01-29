#!/usr/bin/env python2
# Combine various runs of a tetrad analysis into master dat files

from setup import *
import csv

#os.system('clear')

print ""
print "---- Combine-Runs Utility -- Alignment ----"
print ""

# Setup directory structure and print
(_, simdir, _) = simParams(sys)
(root, simdir, datadir, _) = directoryStructure(simdir)

print " Root dir: " + root
print " Sim dir:  " + simdir
print " Data dir:  " + datadir

# Get list of directories in the datadir -- corresponds to different runs
runs = sorted_nicely(glob.glob(datadir + "ts_*"))

# Initalize data structs and arrays
nRuns = int(len(runs))
nTetrads = np.zeros(nRuns)
global nTsteps 
nTsteps = np.zeros(nRuns)

# Structs
g1s1 = [ structtype() for i in range(nRuns) ]
g2s1 = [ structtype() for i in range(nRuns) ]
g3s1 = [ structtype() for i in range(nRuns) ]
g1s2 = [ structtype() for i in range(nRuns) ]
g2s2 = [ structtype() for i in range(nRuns) ]
g3s2 = [ structtype() for i in range(nRuns) ]
g1s3 = [ structtype() for i in range(nRuns) ]
g2s3 = [ structtype() for i in range(nRuns) ]
g3s3 = [ structtype() for i in range(nRuns) ]

g1_z = [ structtype() for i in range(nRuns) ]
g2_z = [ structtype() for i in range(nRuns) ]
g3_z = [ structtype() for i in range(nRuns) ]

s1_z = [ structtype() for i in range(nRuns) ]
s2_z = [ structtype() for i in range(nRuns) ]
s3_z = [ structtype() for i in range(nRuns) ]

w_z = [ structtype() for i in range(nRuns) ]

w_g1 = [ structtype() for i in range(nRuns) ]
w_g2 = [ structtype() for i in range(nRuns) ]
w_g3 = [ structtype() for i in range(nRuns) ]

w_s1 = [ structtype() for i in range(nRuns) ]
w_s2 = [ structtype() for i in range(nRuns) ]
w_s3 = [ structtype() for i in range(nRuns) ]

# Loop over all runs, find nTsteps
for rr, run in enumerate(runs):
  infoFile = run + '/info.dat'

  nTetrads[rr] = np.genfromtxt(infoFile, skip_header=1, usecols=0)
  nTsteps[rr] = np.genfromtxt(infoFile, skip_header=1, usecols=1)

# find minTsteps
minTsteps = np.min(nTsteps)
allTime = np.zeros((minTsteps, nRuns))

for rr, run in enumerate(runs):
  meanFile = run + '/align.mean'
  sdevFile = run + '/align.sdev'
  skewFile = run + '/align.skew'
  kurtFile = run + '/align.kurt'

  # init data to minTstep size
  g1s1[rr].mean = np.zeros(minTsteps)
  g1s1[rr].sdev = np.zeros(minTsteps)
  g1s1[rr].skew = np.zeros(minTsteps)
  g1s1[rr].kurt = np.zeros(minTsteps)
  g2s1[rr].mean = np.zeros(minTsteps) 
  g2s1[rr].sdev = np.zeros(minTsteps) 
  g2s1[rr].skew = np.zeros(minTsteps)
  g2s1[rr].kurt = np.zeros(minTsteps)
  g3s1[rr].mean = np.zeros(minTsteps)
  g3s1[rr].sdev = np.zeros(minTsteps)
  g3s1[rr].skew = np.zeros(minTsteps)
  g3s1[rr].kurt = np.zeros(minTsteps)

  g1s2[rr].mean = np.zeros(minTsteps)
  g1s2[rr].sdev = np.zeros(minTsteps)
  g1s2[rr].skew = np.zeros(minTsteps)
  g1s2[rr].kurt = np.zeros(minTsteps)
  g2s2[rr].mean = np.zeros(minTsteps)
  g2s2[rr].sdev = np.zeros(minTsteps)
  g2s2[rr].skew = np.zeros(minTsteps)
  g2s2[rr].kurt = np.zeros(minTsteps)
  g3s2[rr].mean = np.zeros(minTsteps)
  g3s2[rr].sdev = np.zeros(minTsteps)
  g3s2[rr].skew = np.zeros(minTsteps)
  g3s2[rr].kurt = np.zeros(minTsteps)

  g1s3[rr].mean = np.zeros(minTsteps)
  g1s3[rr].sdev = np.zeros(minTsteps)
  g1s3[rr].skew = np.zeros(minTsteps)
  g1s3[rr].kurt = np.zeros(minTsteps)
  g2s3[rr].mean = np.zeros(minTsteps)
  g2s3[rr].sdev = np.zeros(minTsteps)
  g2s3[rr].skew = np.zeros(minTsteps)
  g2s3[rr].kurt = np.zeros(minTsteps)
  g3s3[rr].mean = np.zeros(minTsteps)
  g3s3[rr].sdev = np.zeros(minTsteps)
  g3s3[rr].skew = np.zeros(minTsteps)
  g3s3[rr].kurt = np.zeros(minTsteps)

  g1_z[rr].mean = np.zeros(minTsteps)
  g1_z[rr].sdev = np.zeros(minTsteps)
  g1_z[rr].skew = np.zeros(minTsteps)
  g1_z[rr].kurt = np.zeros(minTsteps)
  g2_z[rr].mean = np.zeros(minTsteps)
  g2_z[rr].sdev = np.zeros(minTsteps)
  g2_z[rr].skew = np.zeros(minTsteps)
  g2_z[rr].kurt = np.zeros(minTsteps)
  g3_z[rr].mean = np.zeros(minTsteps)
  g3_z[rr].sdev = np.zeros(minTsteps)
  g3_z[rr].skew = np.zeros(minTsteps)
  g3_z[rr].kurt = np.zeros(minTsteps)

  s1_z[rr].mean = np.zeros(minTsteps)
  s1_z[rr].sdev = np.zeros(minTsteps)
  s1_z[rr].skew = np.zeros(minTsteps)
  s1_z[rr].kurt = np.zeros(minTsteps)
  s2_z[rr].mean = np.zeros(minTsteps)
  s2_z[rr].sdev = np.zeros(minTsteps)
  s2_z[rr].skew = np.zeros(minTsteps)
  s2_z[rr].kurt = np.zeros(minTsteps)
  s3_z[rr].mean = np.zeros(minTsteps)
  s3_z[rr].sdev = np.zeros(minTsteps)
  s3_z[rr].skew = np.zeros(minTsteps)
  s3_z[rr].kurt = np.zeros(minTsteps)

  w_z[rr].mean  = np.zeros(minTsteps)
  w_z[rr].sdev  = np.zeros(minTsteps)
  w_z[rr].skew  = np.zeros(minTsteps)
  w_z[rr].kurt  = np.zeros(minTsteps)

  w_g1[rr].mean = np.zeros(minTsteps)
  w_g1[rr].sdev = np.zeros(minTsteps)
  w_g1[rr].skew = np.zeros(minTsteps)
  w_g1[rr].kurt = np.zeros(minTsteps)
  w_g2[rr].mean = np.zeros(minTsteps)
  w_g2[rr].sdev = np.zeros(minTsteps)
  w_g2[rr].skew = np.zeros(minTsteps)
  w_g2[rr].kurt = np.zeros(minTsteps)
  w_g3[rr].mean = np.zeros(minTsteps)
  w_g3[rr].sdev = np.zeros(minTsteps)
  w_g3[rr].skew = np.zeros(minTsteps)
  w_g3[rr].kurt = np.zeros(minTsteps)

  w_s1[rr].mean = np.zeros(minTsteps)
  w_s1[rr].sdev = np.zeros(minTsteps)
  w_s1[rr].skew = np.zeros(minTsteps)
  w_s1[rr].kurt = np.zeros(minTsteps)
  w_s2[rr].mean = np.zeros(minTsteps)
  w_s2[rr].sdev = np.zeros(minTsteps)
  w_s2[rr].skew = np.zeros(minTsteps)
  w_s2[rr].kurt = np.zeros(minTsteps)
  w_s3[rr].mean = np.zeros(minTsteps)
  w_s3[rr].sdev = np.zeros(minTsteps)
  w_s3[rr].skew = np.zeros(minTsteps)
  w_s3[rr].kurt = np.zeros(minTsteps)

  # only pull minTsteps 
  # -- skip_header=1, skip_footer = nTsteps[rr] - minTsteps
  skip = int(nTsteps[rr] - minTsteps)

  tmpT = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, usecols=0)
  allTime[:,rr] = tmpT - tmpT[0]

  g1s1[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip,
    usecols=1)
  g1s1[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=1)
  g1s1[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=1)
  g1s1[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=1)
  g2s1[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=2)
  g2s1[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=2)
  g2s1[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=2)
  g2s1[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=2)
  g3s1[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=3)
  g3s1[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=3)
  g3s1[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=3)
  g3s1[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=3)

  g1s2[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=4)
  g1s2[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=4)
  g1s2[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=4)
  g1s2[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=4)
  g2s2[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=5)
  g2s2[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=5)
  g2s2[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=5)
  g2s2[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=5)
  g3s2[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=6)
  g3s2[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=6)
  g3s2[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=6)
  g3s2[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=6)

  g1s3[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=7)
  g1s3[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=7)
  g1s3[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=7)
  g1s3[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=7)
  g2s3[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=8)
  g2s3[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=8)
  g2s3[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=8)
  g2s3[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=8)
  g3s3[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=9)
  g3s3[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=9)
  g3s3[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=9)
  g3s3[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=9)

  g1_z[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=10) 
  g1_z[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=10)
  g1_z[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=10)
  g1_z[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=10)
  g2_z[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=11)
  g2_z[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=11)
  g2_z[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=11)
  g2_z[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=11)
  g3_z[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=12)
  g3_z[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=12)
  g3_z[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=12)
  g3_z[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=12)

  s1_z[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=13)
  s1_z[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=13)
  s1_z[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=13)
  s1_z[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=13)
  s2_z[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=14)
  s2_z[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=14)
  s2_z[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=14)
  s2_z[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=14)
  s3_z[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=15)
  s3_z[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=15)
  s3_z[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=15)
  s3_z[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=15)

  w_z[rr].mean[:]  = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=16) 
  w_z[rr].sdev[:]  = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=16) 
  w_z[rr].skew[:]  = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=16) 
  w_z[rr].kurt[:]  = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=16) 

  w_g1[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=17)
  w_g1[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=17)
  w_g1[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=17)
  w_g1[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=17)
  w_g2[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=18)
  w_g2[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=18)
  w_g2[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=18)
  w_g2[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=18)
  w_g3[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=19)
  w_g3[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=19)
  w_g3[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=19)
  w_g3[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=19)

  w_s1[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=20)
  w_s1[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=20)
  w_s1[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=20)
  w_s1[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=20)
  w_s2[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=21)
  w_s2[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=21)
  w_s2[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=21)
  w_s2[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=21)
  w_s3[rr].mean[:] = np.genfromtxt(meanFile, skip_header=1, skip_footer=skip, 
    usecols=22)
  w_s3[rr].sdev[:] = np.genfromtxt(sdevFile, skip_header=1, skip_footer=skip, 
    usecols=22)
  w_s3[rr].skew[:] = np.genfromtxt(skewFile, skip_header=1, skip_footer=skip, 
    usecols=22)
  w_s3[rr].kurt[:] = np.genfromtxt(kurtFile, skip_header=1, skip_footer=skip, 
    usecols=22)

totalTetrads = np.sum(nTetrads)
print '    Total tetrads tracked: ' + str(totalTetrads)

# Average timestep and corresponding time to interplate to
avgTime = np.sum(allTime, 1)
avgTime /= nRuns
avgDt =  np.round(np.mean(np.diff(avgTime)))
time = np.arange(0, avgDt*np.size(avgTime), avgDt) 

## Interpolate data
g1s1 = interp(time, allTime, g1s1, nRuns)
g2s1 = interp(time, allTime, g2s1, nRuns) 
g3s1 = interp(time, allTime, g3s1, nRuns) 
g1s2 = interp(time, allTime, g1s2, nRuns) 
g2s2 = interp(time, allTime, g2s2, nRuns) 
g3s2 = interp(time, allTime, g3s2, nRuns) 
g1s3 = interp(time, allTime, g1s3, nRuns) 
g2s3 = interp(time, allTime, g2s3, nRuns) 
g3s3 = interp(time, allTime, g3s3, nRuns) 
                             
g1_z = interp(time, allTime, g1_z, nRuns) 
g2_z = interp(time, allTime, g2_z, nRuns) 
g3_z = interp(time, allTime, g3_z, nRuns) 
                             
s1_z = interp(time, allTime, s1_z, nRuns) 
s2_z = interp(time, allTime, s2_z, nRuns) 
s3_z = interp(time, allTime, s3_z, nRuns) 
                             
w_z =  interp(time, allTime, w_z, nRuns)
                             
w_g1 = interp(time, allTime, w_g1, nRuns) 
w_g2 = interp(time, allTime, w_g2, nRuns) 
w_g3 = interp(time, allTime, w_g3, nRuns) 
                             
w_s1 = interp(time, allTime, w_s1, nRuns) 
w_s2 = interp(time, allTime, w_s2, nRuns) 
w_s3 = interp(time, allTime, w_s3, nRuns) 

## Find moments
m_g1s1 = stats(g1s1, minTsteps, nTetrads, nRuns)
m_g2s1 = stats(g2s1, minTsteps, nTetrads, nRuns)
m_g3s1 = stats(g3s1, minTsteps, nTetrads, nRuns)
m_g1s2 = stats(g1s2, minTsteps, nTetrads, nRuns)
m_g2s2 = stats(g2s2, minTsteps, nTetrads, nRuns)
m_g3s2 = stats(g3s2, minTsteps, nTetrads, nRuns)
m_g1s3 = stats(g1s3, minTsteps, nTetrads, nRuns)
m_g2s3 = stats(g2s3, minTsteps, nTetrads, nRuns)
m_g3s3 = stats(g3s3, minTsteps, nTetrads, nRuns)

m_g1_z = stats(g1_z, minTsteps, nTetrads, nRuns)
m_g2_z = stats(g2_z, minTsteps, nTetrads, nRuns)
m_g3_z = stats(g3_z, minTsteps, nTetrads, nRuns)

m_s1_z = stats(s1_z, minTsteps, nTetrads, nRuns)
m_s2_z = stats(s2_z, minTsteps, nTetrads, nRuns)
m_s3_z = stats(s3_z, minTsteps, nTetrads, nRuns)

m_w_z =  stats(w_z, minTsteps, nTetrads, nRuns)

m_w_g1 = stats(w_g1, minTsteps, nTetrads, nRuns)
m_w_g2 = stats(w_g2, minTsteps, nTetrads, nRuns)
m_w_g3 = stats(w_g3, minTsteps, nTetrads, nRuns)

m_w_s1 = stats(w_s1, minTsteps, nTetrads, nRuns)
m_w_s2 = stats(w_s2, minTsteps, nTetrads, nRuns)
m_w_s3 = stats(w_s3, minTsteps, nTetrads, nRuns)

# Print to file -- mean
allFile = datadir + 'align.mean'
with open(allFile, 'wb') as outfile:
  a = csv.writer(outfile, delimiter=' ')
  headers = [['time', 'g1s1', 'g2s1', 'g3s1', 'g1s2', 'g2s2', 'g3s2',
    'g1s3', 'g2s3', 'g3s3',
    'g1_z', 'g2_z', 'g3_z',
    's1_z', 's2_z', 's3_z',
    'w_z', 
    'w_g1', 'w_g2', 'w_g3',
    'w_s1', 'w_s2', 'w_s3']]
  a.writerows(headers)
  for tt in np.arange(0, minTsteps):
    data = [[time[tt], m_g1s1[tt,0], m_g2s1[tt,0], m_g3s1[tt,0], m_g1s2[tt,0],
      m_g2s2[tt,0], m_g3s2[tt,0],
      m_g1s3[tt,0], m_g2s3[tt,0], m_g3s3[tt,0],
      m_g1_z[tt,0], m_g2_z[tt,0], m_g3_z[tt,0],
      m_s1_z[tt,0], m_s2_z[tt,0], m_s3_z[tt,0],
      m_w_z[tt,0],
      m_w_g1[tt,0], m_w_g2[tt,0], m_w_g3[tt,0],
      m_w_s1[tt,0], m_w_s2[tt,0], m_w_s3[tt,0]]]
    a.writerows(data)


# Print to file -- sdev
allFile = datadir + 'align.sdev'
with open(allFile, 'wb') as outfile:
  a = csv.writer(outfile, delimiter=' ')
  headers = [['time', 'g1s1', 'g2s1', 'g3s1', 'g1s2', 'g2s2', 'g3s2',
    'g1s3', 'g2s3', 'g3s3',
    'g1_z', 'g2_z', 'g3_z',
    's1_z', 's2_z', 's3_z',
    'w_z', 
    'w_g1', 'w_g2', 'w_g3',
    'w_s1', 'w_s2', 'w_s3']]
  a.writerows(headers)
  for tt in np.arange(0, minTsteps):
    data = [[time[tt], m_g1s1[tt,1], m_g2s1[tt,1], m_g3s1[tt,1], m_g1s2[tt,1],
      m_g2s2[tt,1], m_g3s2[tt,1],
      m_g1s3[tt,1], m_g2s3[tt,1], m_g3s3[tt,1],
      m_g1_z[tt,1], m_g2_z[tt,1], m_g3_z[tt,1],
      m_s1_z[tt,1], m_s2_z[tt,1], m_s3_z[tt,1],
      m_w_z[tt,1],
      m_w_g1[tt,1], m_w_g2[tt,1], m_w_g3[tt,1],
      m_w_s1[tt,1], m_w_s2[tt,1], m_w_s3[tt,1]]]
    a.writerows(data)

# Print to file -- skew
allFile = datadir + 'align.skew'
with open(allFile, 'wb') as outfile:
  a = csv.writer(outfile, delimiter=' ')
  headers = [['time', 'g1s1', 'g2s1', 'g3s1', 'g1s2', 'g2s2', 'g3s2',
    'g1s3', 'g2s3', 'g3s3',
    'g1_z', 'g2_z', 'g3_z',
    's1_z', 's2_z', 's3_z',
    'w_z', 
    'w_g1', 'w_g2', 'w_g3',
    'w_s1', 'w_s2', 'w_s3']]
  a.writerows(headers)
  for tt in np.arange(0, minTsteps):
    data = [[time[tt], m_g1s1[tt,2], m_g2s1[tt,2], m_g3s1[tt,2], m_g1s2[tt,2],
      m_g2s2[tt,2], m_g3s2[tt,2],
      m_g1s3[tt,2], m_g2s3[tt,2], m_g3s3[tt,2],
      m_g1_z[tt,2], m_g2_z[tt,2], m_g3_z[tt,2],
      m_s1_z[tt,2], m_s2_z[tt,2], m_s3_z[tt,2],
      m_w_z[tt,2],
      m_w_g1[tt,2], m_w_g2[tt,2], m_w_g3[tt,2],
      m_w_s1[tt,2], m_w_s2[tt,2], m_w_s3[tt,2]]]
    a.writerows(data)

# Print to file -- kurt
allFile = datadir + 'align.kurt'
with open(allFile, 'wb') as outfile:
  a = csv.writer(outfile, delimiter=' ')
  headers = [['time', 'g1s1', 'g2s1', 'g3s1', 'g1s2', 'g2s2', 'g3s2',
    'g1s3', 'g2s3', 'g3s3',
    'g1_z', 'g2_z', 'g3_z',
    's1_z', 's2_z', 's3_z',
    'w_z', 
    'w_g1', 'w_g2', 'w_g3',
    'w_s1', 'w_s2', 'w_s3']]
  a.writerows(headers)
  for tt in np.arange(0, minTsteps):
    data = [[time[tt], m_g1s1[tt,3], m_g2s1[tt,3], m_g3s1[tt,3], m_g1s2[tt,3],
      m_g2s2[tt,3], m_g3s2[tt,3],
      m_g1s3[tt,3], m_g2s3[tt,3], m_g3s3[tt,3],
      m_g1_z[tt,3], m_g2_z[tt,3], m_g3_z[tt,3],
      m_s1_z[tt,3], m_s2_z[tt,3], m_s3_z[tt,3],
      m_w_z[tt,3],
      m_w_g1[tt,3], m_w_g2[tt,3], m_w_g3[tt,3],
      m_w_s1[tt,3], m_w_s2[tt,3], m_w_s3[tt,3]]]
    a.writerows(data)

## g1_z
#fig = plt.figure(figsize=(12,8))
#m1 = fig.add_subplot(221)
#plt.plot(time, m_g1_z[:,0], 'k', linewidth=2)
#for rr in np.arange(0,nRuns):
#  plt.plot(allTime[:,rr], g1_z[rr].mean)
#m1.set_title('Mean g1_z')
#
#m2 = fig.add_subplot(222)
#plt.plot(time, m_g1_z[:,1], 'k', linewidth=2)
#for rr in np.arange(0,nRuns):
#  plt.plot(allTime[:,rr], g1_z[rr].sdev)
#m2.set_title('sdev g1_z')
#
#m3 = fig.add_subplot(223)
#plt.plot(time, m_g1_z[:,2], 'k', linewidth=2)
#for rr in np.arange(0,nRuns):
#  plt.plot(allTime[:,rr], g1_z[rr].skew)
#m3.set_title('skew g1_z')
#
#m4 = fig.add_subplot(224)
#plt.plot(time, m_g1_z[:,3], 'k', linewidth=2)
#for rr in np.arange(0,nRuns):
#  plt.plot(allTime[:,rr], g1_z[rr].kurt)
#m4.set_title('kurt g1_z')
#
#plt.show()
