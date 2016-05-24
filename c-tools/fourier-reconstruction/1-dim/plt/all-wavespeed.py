#!/usr/bin/env python2
from setup import *
os.system("clear")

from scipy.ndimage import imread

print ""
print " ---- Fourier Reconstruction Plotting Utility ---- "
print "              All Time Autocorrs"
print ""

root = "/home-1/dwillen3@jhu.edu/scratch/triply_per/"
imgdir = root + "simdata/img/f-rec-1D/"
if not os.path.exists(imgdir):
  os.makedirs(imgdir)

# Load images
i0500_20 = imread(root + "500/rho2.0/f-rec-1D/img/autocorr-time-vf.png")
i0500_33 = imread(root + "500/rho3.3/f-rec-1D/img/autocorr-time-vf.png")
i0500_40 = imread(root + "500/rho4.0/f-rec-1D/img/autocorr-time-vf.png")
i0500_50 = imread(root + "500/rho5.0/f-rec-1D/img/autocorr-time-vf.png")

i1000_20 = imread(root + "1000/rho2.0/f-rec-1D/img/autocorr-time-vf.png")
i1000_33 = imread(root + "1000/rho3.3/f-rec-1D/img/autocorr-time-vf.png")
i1000_40 = imread(root + "1000/rho4.0/f-rec-1D/img/autocorr-time-vf.png")
i1000_50 = imread(root + "1000/rho5.0/f-rec-1D/img/autocorr-time-vf.png")

i1500_20 = imread(root + "1500/rho2.0/f-rec-1D/img/autocorr-time-vf.png")
i1500_33 = imread(root + "1500/rho3.3/f-rec-1D/img/autocorr-time-vf.png")
i1500_40 = imread(root + "1500/rho4.0/f-rec-1D/img/autocorr-time-vf.png")

i2000_20 = imread(root + "2000/rho2.0/f-rec-1D/img/autocorr-time-vf.png")
i2000_33 = imread(root + "2000/rho3.3/f-rec-1D/img/autocorr-time-vf.png")

# Plot
autocorrTimeFig = plt.figure(figsize=(16,12))

a11 = autocorrTimeFig.add_subplot(4,4,1)
a11.imshow(i0500_20,cmap='gray')
a11.xaxis.set_ticks([])
a11.yaxis.set_ticks([])
a11.set_title(r'$\phi = 0.087$')
a11.set_ylabel(r'$\rho^* = 2.0$')

a12 = autocorrTimeFig.add_subplot(4,4,2)
a12.imshow(i1000_20,cmap='gray')
a12.xaxis.set_ticks([])
a12.yaxis.set_ticks([])
a12.set_title(r'$\phi = 0.175$')

a13 = autocorrTimeFig.add_subplot(4,4,3)
a13.imshow(i1500_20,cmap='gray')
a13.xaxis.set_ticks([])
a13.yaxis.set_ticks([])
a13.set_title(r'$\phi = 0.262$')

a14 = autocorrTimeFig.add_subplot(4,4,4)
a14.imshow(i2000_20,cmap='gray')
a14.xaxis.set_ticks([])
a14.yaxis.set_ticks([])
a14.set_title(r'$\phi = 0.349$')

a21 = autocorrTimeFig.add_subplot(4,4,5)
a21.imshow(i0500_33,cmap='gray')
a21.xaxis.set_ticks([])
a21.yaxis.set_ticks([])
a21.set_ylabel(r'$\rho^* = 3.3$')

a22 = autocorrTimeFig.add_subplot(4,4,6)
a22.imshow(i1000_33,cmap='gray')
a22.xaxis.set_ticks([])
a22.yaxis.set_ticks([])

a23 = autocorrTimeFig.add_subplot(4,4,7)
a23.imshow(i1500_33,cmap='gray')
a23.xaxis.set_ticks([])
a23.yaxis.set_ticks([])

a24 = autocorrTimeFig.add_subplot(4,4,8)
a24.imshow(i2000_33,cmap='gray')
a24.xaxis.set_ticks([])
a24.yaxis.set_ticks([])

a31 = autocorrTimeFig.add_subplot(4,4,9)
a31.imshow(i0500_40,cmap='gray')
a31.xaxis.set_ticks([])
a31.yaxis.set_ticks([])
a31.set_ylabel(r'$\rho^* = 4.0$')

a32 = autocorrTimeFig.add_subplot(4,4,10)
a32.imshow(i1000_40,cmap='gray')
a32.xaxis.set_ticks([])
a32.yaxis.set_ticks([])

a33 = autocorrTimeFig.add_subplot(4,4,11)
a33.imshow(i1500_40,cmap='gray')
a33.xaxis.set_ticks([])
a33.yaxis.set_ticks([])

a41 = autocorrTimeFig.add_subplot(4,4,13)
a41.imshow(i0500_50,cmap='gray')
a41.xaxis.set_ticks([])
a41.yaxis.set_ticks([])
a41.set_ylabel(r'$\rho^* = 5.0$')

a42 = autocorrTimeFig.add_subplot(4,4,14)
a42.imshow(i1000_50,cmap='gray')
a42.xaxis.set_ticks([])
a42.yaxis.set_ticks([])

# SAVE
imgname = imgdir + "power-spectrum-time-all"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')

print "\n      ...Done!"
