#!/usr/bin/env python
import sys
import argparse
import numpy as np
from scipy.optimize import leastsq
import matplotlib as mpl
mpl.use('GTKCairo')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from settings import *
from kandiLib import *

'''
Plot 15-minute averaged horizontal profiles.
'''

#==========================================================#

parser = argparse.ArgumentParser(
    prog='plotProfiles.py', description='''Plot time-averaged profiles.''')
parser.add_argument("-f", "--file", type=str,default=None,
                    help="Name of the input netCDF4 file.")
parser.add_argument("-v", "--variable", type=str,
                    default="u", help="Variable to be plotted.")
parser.add_argument("-cm", "--colormap", type=str, default="viridis",
                    help="Colormap style. Default viridis.")
parser.add_argument("-z","--height", type=float, help="Height of the cross-sectional slice.")
parser.add_argument("-stat", "--statistics", type=str, default="avg", help="Statistics to be displayed. Options: avg, max, min")
parser.add_argument("-d", "--discretize", type=int, default=15,
                    help="Discretize color map to number of colors. Default: 15")
parser.add_argument("--lims", type=float, nargs=2, help="Colorbar limits.")
parser.add_argument("-t","--tstep", type=float, help="Time-averaging interval.")
parser.add_argument("-fs", "--figsize", type=int, nargs=2, help="Subplot figures grid NxM", default=[2,4])
args = parser.parse_args()

#==========================================================#

clevels = np.linspace(args.lims[0], args.lims[1], args.discretize + 1)
cmap = mpl.cm.get_cmap(args.colormap)


ds = openDataSet(args.file)

var, x_dims, y_dims, z_dims = readVariableFromMask(ds, [0, 7200.], args.variable)
h_index, = np.where(z_dims == args.height)[0]
var = var[:,h_index,:,:]

t_set=np.arange(0,7201,args.tstep)

fig = plt.figure()
for i in xrange(len(t_set)-1):
  t_inds, = np.where(np.logical_and(ds.variables['time'][:] >= t_set[i], ds.variables['time'][:] <= t_set[i+1]))
  var_temp= calculateTemporalStatistics(var[t_inds,:,:], args.statistics)
  ax = fig.add_subplot(args.figsize[0],args.figsize[1],i+1)
  imgplot = ax.contourf(x_dims, y_dims, var_temp, levels=clevels, cmap=cmap, origin='lower',zorder=0)
  ax.title.set_text("{:.0f}-{:.0f} seconds".format(t_set[i],t_set[i+1]))
  fig.colorbar(imgplot, ax=ax)

plt.show()
