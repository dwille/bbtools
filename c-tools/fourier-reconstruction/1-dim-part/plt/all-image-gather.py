#!/usr/bin/env python
from setup import *
os.system("clear")

from scipy.ndimage import imread

#print ""
#print " ---- Fourier Reconstruction Plotting Utility ---- "
#print "              Image Gather"
#print ""

root = os.path.expanduser("~") +  "/scratch/triply_per/"
analysis_dir = "analysis/fourier-reconstruction/1-dim-part/img/"
imgdir = root + "simdata/img/f-rec-1D/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

#File = "volume-fraction"
#File = "avg-power-spectrum-time-vf"
#File = "avg-autocorr-time-vf"
File = "avg-power-spectrum-space-vf"
#File = "avg-autocorr-space-vf"
#File = "vf-coeffs"

#File = "autocorr-time-vf"
#File = "autocorr-space-vf"
#File = "avg-crosscorr-spacetime-vf"
#File = "hist-vfrac-wvel"
print("     Using " + File)

# file type
ftype = ".png"

# Load images
i0500_20 = imread(root + "500/rho2.0/" + analysis_dir + File + ftype)
i0500_33 = imread(root + "500/rho3.3/" + analysis_dir + File + ftype)
i0500_40 = imread(root + "500/rho4.0/" + analysis_dir + File + ftype)
i0500_50 = imread(root + "500/rho5.0/" + analysis_dir + File + ftype)

i1000_20 = imread(root + "1000/rho2.0/" + analysis_dir + File + ftype)
i1000_33 = imread(root + "1000/rho3.3/" + analysis_dir + File + ftype)
i1000_40 = imread(root + "1000/rho4.0/" + analysis_dir + File + ftype)
i1000_50 = imread(root + "1000/rho5.0/" + analysis_dir + File + ftype)

i1500_20 = imread(root + "1500/rho2.0/" + analysis_dir + File + ftype)
i1500_33 = imread(root + "1500/rho3.3/" + analysis_dir + File + ftype)
i1500_40 = imread(root + "1500/rho4.0/" + analysis_dir + File + ftype)

i2000_20 = imread(root + "2000/rho2.0/" + analysis_dir + File + ftype)
i2000_33 = imread(root + "2000/rho3.3/" + analysis_dir + File + ftype)

# Plot
fig = plt.figure(figsize=(6,6))
## 16, 12

a11 = fig.add_subplot(4,4,1)
a11.imshow(i0500_20,cmap='gray')
a11.xaxis.set_ticks([])
a11.yaxis.set_ticks([])
a11.set_title(r'$\phi = 0.087$')
a11.set_ylabel(r'$\rho^* = 2.0$')

a12 = fig.add_subplot(4,4,2)
a12.imshow(i1000_20,cmap='gray')
a12.xaxis.set_ticks([])
a12.yaxis.set_ticks([])
a12.set_title(r'$\phi = 0.175$')

a13 = fig.add_subplot(4,4,3)
a13.imshow(i1500_20,cmap='gray')
a13.xaxis.set_ticks([])
a13.yaxis.set_ticks([])
a13.set_title(r'$\phi = 0.262$')

a14 = fig.add_subplot(4,4,4)
a14.imshow(i2000_20,cmap='gray')
a14.xaxis.set_ticks([])
a14.yaxis.set_ticks([])
a14.set_title(r'$\phi = 0.349$')

a21 = fig.add_subplot(4,4,5)
a21.imshow(i0500_33,cmap='gray')
a21.xaxis.set_ticks([])
a21.yaxis.set_ticks([])
a21.set_ylabel(r'$\rho^* = 3.3$')

a22 = fig.add_subplot(4,4,6)
a22.imshow(i1000_33,cmap='gray')
a22.xaxis.set_ticks([])
a22.yaxis.set_ticks([])

a23 = fig.add_subplot(4,4,7)
a23.imshow(i1500_33,cmap='gray')
a23.xaxis.set_ticks([])
a23.yaxis.set_ticks([])

a24 = fig.add_subplot(4,4,8)
a24.imshow(i2000_33,cmap='gray')
a24.xaxis.set_ticks([])
a24.yaxis.set_ticks([])

a31 = fig.add_subplot(4,4,9)
a31.imshow(i0500_40,cmap='gray')
a31.xaxis.set_ticks([])
a31.yaxis.set_ticks([])
a31.set_ylabel(r'$\rho^* = 4.0$')

a32 = fig.add_subplot(4,4,10)
a32.imshow(i1000_40,cmap='gray')
a32.xaxis.set_ticks([])
a32.yaxis.set_ticks([])

a33 = fig.add_subplot(4,4,11)
a33.imshow(i1500_40,cmap='gray')
a33.xaxis.set_ticks([])
a33.yaxis.set_ticks([])

a41 = fig.add_subplot(4,4,13)
a41.imshow(i0500_50,cmap='gray')
a41.xaxis.set_ticks([])
a41.yaxis.set_ticks([])
a41.set_ylabel(r'$\rho^* = 5.0$')

a42 = fig.add_subplot(4,4,14)
a42.imshow(i1000_50,cmap='gray')
a42.xaxis.set_ticks([])
a42.yaxis.set_ticks([])

# SAVE
imgname = imgdir + "all-" + File
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
plt.savefig(imgname + ".eps", bbox_inches='tight', format='eps')

print("\n      ...Done!")
