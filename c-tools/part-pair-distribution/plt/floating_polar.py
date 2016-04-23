from __future__ import division
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.transforms import Affine2D
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist import angle_helper
from mpl_toolkits.axisartist.grid_finder import MaxNLocator, DictFormatter
from mpl_toolkits.axisartist.floating_axes import GridHelperCurveLinear, FloatingSubplot

# From stackoverflow 32232285
def fractional_polar_axes(f, thlim=(0, 180), rlim=(0, 1), step=(30, 0.25),
  thlabel='theta', rlabel='r', ticklabels=True, theta_offset=0, rlabels = None):
  '''Return polar axes that adhere to desired theta (in deg) and r limits. steps for theta
  and r are really just hints for the locators.'''
  th0, th1 = thlim # deg
  r0, r1 = rlim
  thstep, rstep = step

  tr_rotate = Affine2D().translate(theta_offset, 0)
  # scale degrees to radians:
  tr_scale = Affine2D().scale(np.pi/180., 1.)
  #pa = axes(polar="true") # create a polar axis
  pa = PolarAxes
  tr = tr_rotate + tr_scale + pa.PolarTransform()

  #theta_grid_locator = angle_helper.LocatorDMS((th1-th0)//thstep)
  angle_ticks = [(0, r"$0$"),
                 (.25*np.pi, r"$\frac{1}{4}\pi$"),
                 (.5*np.pi, r"$\frac{1}{2}\pi$")]
  theta_grid_locator1 = FixedLocator([v for v, s in angle_ticks])

  r_grid_locator = MaxNLocator((r1-r0)//rstep)
  #theta_tick_formatter = angle_helper.FormatterDMS()
  theta_tick_formatter = DictFormatter(dict(angle_ticks))
  if rlabels:
    rlabels = DictFormatter(rlabels)

  grid_helper = GridHelperCurveLinear(tr,
                                      extremes=(0, 90, r0, r1),
                                      grid_locator1=theta_grid_locator,
                                      grid_locator2=r_grid_locator,
                                      tick_formatter1=theta_tick_formatter,
                                      tick_formatter2=rlabels)

  a = FloatingSubplot(f, 111, grid_helper=grid_helper)
  f.add_subplot(a)
  # adjust x axis (theta):
  print(a)
  a.axis["bottom"].set_visible(False)
  a.axis["top"].set_axis_direction("bottom") # tick direction
  a.axis["top"].toggle(ticklabels=ticklabels, label=bool(thlabel))
  a.axis["top"].major_ticklabels.set_axis_direction("top")
  a.axis["top"].label.set_axis_direction("top")
  a.axis["top"].major_ticklabels.set_pad(10)

  # adjust y axis (r):
  a.axis["left"].set_axis_direction("bottom") # tick direction
  a.axis["right"].set_axis_direction("top") # tick direction
  a.axis["left"].toggle(ticklabels=True, label=bool(rlabel))
  # add labels:
  a.axis["top"].label.set_text(thlabel)
  a.axis["left"].label.set_text(rlabel)
  # create a parasite axes whose transdata is theta, r:
  auxa = a.get_aux_axes(tr)
  print(auxa)
  # make aux_ax to have a clip path as in a?:
  auxa.patch = a.patch 
  # this has a side effect that the patch is drawn twice, and possibly over some other
  # artists. so, we decrease the zorder a bit to prevent this:
  a.patch.zorder = -2

  # add sector lines for both dimensions:
  thticks = grid_helper.grid_info['lon_info'][0]
  rticks = grid_helper.grid_info['lat_info'][0]
  print(grid_helper.grid_info['lat_info'])
  for th in thticks[1:-1]: # all but the first and last
    auxa.plot([th, th], [r0, r1], ':', c='grey', zorder=-1, lw=0.5)
  for ri, r in enumerate(rticks):
    # plot first r line as axes border in solid black only if it  isn't at r=0
    if ri == 0 and r != 0:
      ls, lw, color = 'solid', 1, 'black'
    else:
      ls, lw, color = 'dashed', 0.5, 'grey'
    # from http://stackoverflow.com/a/19828753/2020363
    auxa.add_artist(plt.Circle([0, 0], radius=r, ls=ls, lw=lw, color=color, fill=False,
    transform=auxa.transData._b, zorder=-1))

  return auxa
