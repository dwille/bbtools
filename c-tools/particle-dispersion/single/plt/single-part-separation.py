#!/usr/bin/env python

# Temporary plotting utility for particle separation

import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("./data/data.csv", skip_header=1, delimiter=",")

time = data[:,0]      #  ms
time -= time[0]
r_total = data[:,1]   #  mm
r_verti = data[:,2]   #  mm
r_horiz = data[:,3]   #  mm

# predefined constants
a = 2.1           # mm
nu = 0.01715      # mm^2/ms
rho_f = 8.75e-4   # g/mm^3
rho_p = 0.00289   # g/mm^3
rho = rho_p / rho_f
Lx=42.            # mm
Ly=42.            # mm
Lz=126.           # mm
wf=0.22204        # from phaseAveragedFluidVel, mm/ms

tau_p = (2.*a)**2 * rho / (18. * nu)


fig = plt.figure()
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')

plt.plot(time / tau_p, r_total / a, '.', markersize=2)
plt.plot(time / tau_p, r_verti / a, '.', markersize=2)
plt.plot(time / tau_p, r_horiz / a, '.', markersize=2)

plt.axhline(Lx / a, color='k', linestyle='--')
plt.axhline(Lz / a, color='k', linestyle='--')
plt.text(2e-2, 1e1, r"$L_x$")
plt.text(2e-2, 4e1, r"$L_z$")
#plt.axhline(1., color='k', linestyle='--')

# r ~ t^1
xpts = [1e-2, 1.5e0]
ypts = np.power(xpts, 1) * 10.
plt.plot(xpts, ypts, 'k--')
plt.text(2e-2, 3e-1, r"$t^1$")

# r ~ t^0.5
xpts = [1e0, 1e2]
ypts = np.power(xpts, 0.5) * 2.
plt.plot(xpts, ypts, 'k--')
plt.text(2e0, 1.5e0, r"$t^{1/2}$")

plt.xlabel(r"$t / \tau_p$")
plt.ylabel(r"$r/a$")

plt.legend(["total", "vertical", "horizontal"])

plt.tight_layout()
plt.savefig("separation.png", bbox_inches='tight', format='png')
