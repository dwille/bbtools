#!/usr/bin/env python2
from setup import *
os.system('clear')

print ""
print " ---- Histograms of 3D data ---- "
print ""

# Setup
simdir = simParams(sys)
rhopf = float(simdir.partition('/')[2][-4:-1])
(root, simdir, datadir, imgdir) = directoryStructure(simdir)

# Pull info
cols = (1,2,3,4)
(U_min, U_dBin, U_max, U_mean) = np.genfromtxt(datadir + "info",
    skip_header=1, skip_footer=9, usecols=cols)

(ke_min, ke_dBin, ke_max, ke_mean) = np.genfromtxt(datadir + "info",
    skip_header=2, skip_footer=8, usecols=cols)
(ker_min, ker_dBin, ker_max, ker_mean) = np.genfromtxt(datadir + "info",
    skip_header=3, skip_footer=7, usecols=cols)
(ket_min, ket_dBin, ket_max, ket_mean) = np.genfromtxt(datadir + "info",
    skip_header=4, skip_footer=6, usecols=cols)

(Tz_min, Tz_dBin, Tz_max, Tz_mean) = np.genfromtxt(datadir + "info",
    skip_header=5, skip_footer=5, usecols=cols)
(Tperp_min, Tperp_dBin, Tperp_max, Tperp_mean) = np.genfromtxt(datadir + "info",
    skip_header=6, skip_footer=4, usecols=cols)
(T_min, T_dBin, T_max, T_mean) = np.genfromtxt(datadir + "info",
    skip_header=7, skip_footer=3, usecols=cols)

(up_min, up_dBin, up_max, up_mean) = np.genfromtxt(datadir + "info",
    skip_header=8, skip_footer=2, usecols=cols)
(vp_min, vp_dBin, vp_max, vp_mean) = np.genfromtxt(datadir + "info",
    skip_header=9, skip_footer=1, usecols=cols)
(wp_min, wp_dBin, wp_max, wp_mean) = np.genfromtxt(datadir + "info",
    skip_header=10, skip_footer=0, usecols=cols)

# ke -- in g m^2/s^2, convert by /= 1000
# U  -- in mm/ms = m/s
# T  -- in g m^2/s^2, convert by /= 1000

# Pull hists and normalize
hist_U = np.genfromtxt(datadir + "U_hist")
nBins = np.size(hist_U)
hist_U /= (U_dBin * np.sum(hist_U))

hist_ke = np.genfromtxt(datadir + "ke_hist")
hist_ke /= (ke_dBin * np.sum(hist_ke))
hist_ker = np.genfromtxt(datadir + "ke_hist_rot")
hist_ker /= (ker_dBin * np.sum(hist_ker))
hist_ket = np.genfromtxt(datadir + "ke_hist_trans")
hist_ket /= (ket_dBin * np.sum(hist_ket))

hist_Tz = np.genfromtxt(datadir + "T_hist_z")
hist_Tz /= (Tz_dBin * np.sum(hist_Tz))
hist_Tperp = np.genfromtxt(datadir + "T_hist_perp")
hist_Tperp /= (Tperp_dBin * np.sum(hist_Tperp))
hist_T = np.genfromtxt(datadir + "T_hist")
hist_T /= (T_dBin * np.sum(hist_T))

hist_up = np.genfromtxt(datadir + "ux_hist")
hist_up /= (up_dBin * np.sum(hist_up))
hist_vp = np.genfromtxt(datadir + "vy_hist")
hist_vp /= (vp_dBin * np.sum(hist_vp))
hist_wp = np.genfromtxt(datadir + "wz_hist")
hist_wp /= (wp_dBin * np.sum(hist_wp))

# Set up bins
centers_U = findCenters(U_min, U_max, nBins, U_dBin)
centers_ke = findCenters(ke_min, ke_max, nBins, ke_dBin)
centers_ket = findCenters(ket_min, ket_max, nBins, ket_dBin)
centers_ker = findCenters(ker_min, ker_max, nBins, ker_dBin)
centers_Tz = findCenters(Tz_min, Tz_max, nBins, Tz_dBin)
centers_Tperp = findCenters(Tperp_min, Tperp_max, nBins, Tperp_dBin)
centers_T = findCenters(T_min, T_max, nBins, T_dBin)
centers_up = findCenters(up_min, up_max, nBins, up_dBin)
centers_vp = findCenters(vp_min, vp_max, nBins, vp_dBin)
centers_wp = findCenters(wp_min, wp_max, nBins, wp_dBin)

# constants
rhof = 8.75e-4
rhop = rhopf * rhof
partR = 2.1
mass = 4./3.*np.pi*(partR**3.) * rhop # g -- particle mass

# Speed probability relations from maxwell distributions
print "Maxwell estimates vs simulations estimates"
vp = centers_U[np.argmax(hist_U[0:-1])] # mm/ms or m/s (last index is trash)
a = vp/np.sqrt(2.)                      # scale parameter

# Mean speed
meanV = 2.*vp/np.sqrt(np.pi)
err = (meanV - U_mean)/meanV
print "  (meanV_maxwell - meanV_sim)/meanV_maxwell = %.3f" % err

# RMS speed
rmsV_sim = np.sqrt(2.*(ket_mean/1000.)/(mass/1000.))
rmsV_max = np.sqrt(3./2.)*vp
err = (rmsV_max - rmsV_sim)/rmsV_max
print "  (rmsV_maxwell - rmsV_sim)/rmsV_maxwell    = %.3f" % err

# kT measurements
print "\nEstimates of kT"
# from vp = sqrt(2kT/m)
kt_1 = 0.5*(vp**2.)*(mass/1000.)
print "  kT from most probable speed = %e" % kt_1

# from meanU = sqrt(8kT/pi*m)
kt_2 = np.pi*(mass/1000.)*(U_mean**2.)/8.
print "  kT from mean speed          = %e" % kt_2

# from mean KE_t = 3./2. kbT
kt_3 = 2./3.*(ket_mean/1000.)
print "  kT from rms speed (mean KE) = %e" % kt_3


# fit data to maxwell distr; params = (shift, scale)
xval = np.linspace(0, 0.2, 1000)
yval = stats.maxwell.pdf(xval, loc=0, scale=a)

# velocity pdfs #
fig1 = plt.figure(figsize=(5.5,5.5))

ax1 = fig1.add_subplot(222)
ax1.plot(centers_U[0:-1], hist_U[0:-1], 'b')
ax1.plot([U_mean, U_mean], ax1.get_ylim(), 'b--')
ax1.plot(xval, yval, 'k')

ax1.set_xlabel(r'$\sqrt{(u_i \cdot u_i)}$')
ax1.xaxis.set_major_locator(MultipleLocator(0.05))
ax1.xaxis.set_minor_locator(MultipleLocator(0.025))
ax1.set_ylabel(r'$pdf$')

ax2 = fig1.add_subplot(223)
ax2.plot(centers_Tz[0:-1], hist_Tz[0:-1], 'b')
ax2.plot(centers_Tperp[0:-1], hist_Tperp[0:-1], 'g')
ax2.plot(centers_T[0:-1], hist_T[0:-1], 'k')

ax2.set_xlabel(r'$T_i$')
#ax2.xaxis.set_major_locator(MultipleLocator(0.0005))
#ax2.xaxis.set_minor_locator(MultipleLocator(0.00025))
ax2.set_ylabel(r'$pdf$')
ax2.legend([r"$T_z$",r"$T_\perp$",r"$T = \frac{1}{3}(T_z + T_\perp)$"],
  loc="lower left", framealpha=0.5)
#ax2.set_xlim([0, 0.001])
#ax2.set_ylim([0, 10000])

xval = np.array((2e-4, 2e-3))
yval = xval**(-3.2)*10e-9
ax2.loglog(xval, yval, 'r--')

perp_mean_idx = (np.abs(centers_Tperp - Tperp_mean)).argmin()
z_mean_idx = (np.abs(centers_Tz - Tz_mean)).argmin()
mean_idx = (np.abs(centers_T - T_mean)).argmin()

ax2.loglog(Tperp_mean, hist_Tperp[perp_mean_idx], 'go')
ax2.loglog(Tz_mean, hist_Tz[z_mean_idx], 'bo')
ax2.loglog(T_mean, hist_T[mean_idx], 'ko')

ax3 = fig1.add_subplot(221)
ax3.semilogy(centers_up[1:-2], hist_up[1:-2], 'g')
ax3.semilogy(centers_vp[1:-2], hist_vp[1:-2], 'b--')
ax3.semilogy(centers_wp[1:-2], hist_wp[1:-2], 'k')
ax3.set_xlabel(r'$u_i$')
ax3.xaxis.set_major_locator(MultipleLocator(0.10))
ax3.xaxis.set_minor_locator(MultipleLocator(0.05))
ax3.set_ylabel(r'$pdf$')
ax3.legend([r'$u$',r'$v$',r'$w$'])

ax3.plot([up_mean, up_mean], ax1.get_ylim(), 'g')
ax3.plot([vp_mean, vp_mean], ax1.get_ylim(), 'b--')
ax3.plot([wp_mean, wp_mean], ax1.get_ylim(), 'k')

ax4 = fig1.add_subplot(224)
centers_uperp = np.sqrt(centers_Tperp/mass)
centers_uz = np.sqrt(centers_Tz/mass)
ax4.loglog(centers_uperp[0:-1], hist_Tperp[0:-1], 'g')
ax4.loglog(centers_uz[0:-1], hist_Tz[0:-1], 'b')

ax4.set_xlabel(r'$|u_i|$')

xval = np.array((2e-2, 2e-1))
yval = xval**(-5.)*10e-2
ax4.loglog(xval, yval, 'r--')

imgname = imgdir + "hist"
plt.savefig(imgname + ".png", bbox_inches='tight', format='png')
#plt.savefig(imgname + ".pdf", bbox_inches='tight', format='pdf')

print "      Done!"
