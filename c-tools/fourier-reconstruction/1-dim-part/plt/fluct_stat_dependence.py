#!/usr/bin/env python2
from setup import *
os.system('clear')

print 
print " ---- 1D Fourier Reconstruction ---- "
print "  Fluct Dependence on Reconstruction"
print

# Setup simulation parameters
(partR, simdir, tstart) = simParams(sys)
nparts = int(simdir.partition('/')[0])
domX = 42
domY = 42
domZ = 126
avgVolumeFraction = 4./3.*np.pi*partR*partR*partR*nparts/(domX*domY*domZ)

# Setup directory structures
(root, simdir, datadir, imgdir) = directoryStructure(simdir)

# Print simulation data
printSimulationData(partR, root, simdir, datadir)

# Load dependendence
cols = (1,2,3,4,5,6,7,8,9)
modes = np.genfromtxt(root + simdir + "stat_dependence", skip_footer=10,
  usecols=cols)
data = np.genfromtxt(root + simdir + "stat_dependence", skip_header=1, 
  usecols=cols)
mins = data[0,:]
p01 =  data[1,:]
p05 =  data[2,:]
mean = data[4,:]
sdev = data[5,:]
p95 =  data[7,:]
p99 =  data[8,:]
maxs = data[9,:]

fig1 = plt.figure()

plt.plot(modes, mean, 'ko-', linewidth=2)
plt.plot(modes, mins, 'k--o', linewidth=2)
plt.plot(modes, (mean+sdev), 'ro', linewidth=1.5)
plt.plot(modes, p01, 'bo')
plt.plot(modes, p05, 'go')

lText = [r'$\langle \phi \rangle $', 
  r'$\langle \phi \rangle_{min/max}$', 
  r'$\langle \phi \rangle \pm \sigma_{\langle \phi \rangle}$', 
  r'$\pm 45\%$ intervals', 
  r'$\pm 49\%$ intervals'] 
plt.legend(lText, loc="lower right", framealpha=.5, fontsize=8)

plt.plot(modes, p95, 'go')
plt.plot(modes, p99, 'bo')
plt.plot(modes, (mean-sdev), 'ro', linewidth=1.5)
plt.plot(modes, maxs, 'k--o', linewidth=2)

plt.xlabel('Fourier Modes in Reconstruction')
plt.ylabel(r'$\phi$', rotation=0)

# Save
imgname = imgdir + "vfrac-flucts-dependence"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')
print "      Done!"
