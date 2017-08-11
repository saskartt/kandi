#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from PIL import Image
from kandiLib import *
from settings import *


'''
Plot time-averaged horizontal cross-sections from mask data (following M. Hefny Salim et al / J. Wind. Eng. Ind. Aerodyn. 144 (2015) 84-95)
'''

#==========================================================#

parser = argparse.ArgumentParser(prog='plotHSlice.py',
                                 description='''Plot time-averaged horizontal cross sections from output mask data. Some settings are defined in settings.py''')
parser.add_argument("-f", "--file", type=str, help="Name of the input netCDF4 file.")
parser.add_argument("-var", "--variable", type=str, help="Variable to be plotted")
parser.add_argument("-z","--height", type=float, help="Height of the cross-sectional slice.")
parser.add_argument("-stat", "--statistics", type=str, default="avg", help="Statistics to be displayed. Options: avg, max, min")
parser.add_argument("-cm", "--colormap", type=str, default="RdBu_r",
                    help="Colormap style. Default RdBu_r.")
parser.add_argument("-d", "--discretize", type=int, default=15,
                    help="Discretize color map to number of colors. Default: 15")
parser.add_argument("--lims", type=float, nargs=2, help="Colorbar limits.")
parser.add_argument("-t", "--topography", type=str, help="Topography mask file.")
parser.add_argument("-q", "--quiver", type=str, help="File to read quiver vectors from.")
parser.add_argument("-qs", "--quiverscale", type=float, default=1, help="Scaling factor for the quiver plot arrows.")
parser.add_argument("--canopy", type=str, help="Read single tree image file. The script locates this raster tree file into the plot according to parameters in settings file.")
parser.add_argument("-ca", "--canopyalpha", default=0.5, type=float, help="Alpha (transparency) used for canopy layer with 0 being completely transparent and 1 being opaque.")
parser.add_argument("-s", "--save", type=str,help="Save the resulting figure as an PDF file.")

args = parser.parse_args()

#==========================================================#

# Required to correctly save as vector graphics
mpl.rcParams['image.composite_image'] = False

# Read data from the main dataset
ds = openDataSet(args.file)
var, x_dims, y_dims, z_dims = readVariableFromMask(ds, timespan, args.variable)
# Cut out a slice
h_index, = np.where(z_dims == args.height)[0]
var = var[:,h_index,:,:]
var = calculateTemporalStatistics(var, args.statistics)

# Read a slice of topography
if (args.topography):
  tds = openDataSet(args.topography)
  topo = tds.variables["buildings_0"][:]
  topo = topo[:, :, :]
  topo_mask = np.ma.masked_where(topo == 0, topo)
  # Offset topography by 0.5 m in plot x direction
  topo_x_dims = tds.variables["x"][:]
  topo_y_dims = tds.variables["y"][:]
  tds.close()

# Interpolate variable to topography resolution
x_dims_i = np.arange(x_dims[0], x_dims[-1], 1)  # dimensions to interpolate to
y_dims_i = np.arange(y_dims[0], y_dims[-1], 1)
var = interpolateScalarField(var, x_dims, y_dims, x_dims_i, y_dims_i)

# Set custom axes
hlen = bld_height
xlen = x_dims[-1]-x_dims[0]
ylen = y_dims[-1]-y_dims[0]
y_mp = np.mean([y_dims[0], y_dims[-1]])

xticks = np.arange(0.0, xlen, hlen / 2)
xticklabels = np.arange(0.0, xlen/hlen, 0.5)

yticklabels = np.arange(-1.5, 2.0, 0.5)
yticks = np.arange(y_mp - 1.5 * hlen, y_mp + 2. * hlen, hlen / 2.)

#Set colormap
cmap = mpl.cm.get_cmap(args.colormap)

# Construct the figure
fig = plt.figure(figsize=(9, 4), dpi=120)
ax = plt.gca()
plt.xticks(xticks, xticklabels)
plt.yticks(yticks, yticklabels)
# ax.set_ylim([y_mp - 1.7 * hlen, y_mp + 1.7 * hlen ])
# ax.set_ylim([1, args.ymax * hlen])
ax.set_xlim([topo_x_dims[0], topo_x_dims[-1]])
# NOTE: Should not be hardcoded
plt.axes().set_aspect('equal')


# Color levels
clevels = np.linspace(args.lims[0], args.lims[1], args.discretize + 1)

# Plot filled contour velocity field
imgplot = plt.contourf(x_dims_i-0.5, y_dims_i, var, levels=clevels, cmap=cmap, origin='lower',zorder=0)


# if (args.topography):
  # topoplot = plt.contourf(topo_x_dims, topo_y_dims, topo_mask,
                          # cmap=mpl.cm.binary, vmin=0, vmax=1, interpolation=None, zorder=1)

# Add colormap
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.15)
cbar = fig.colorbar(imgplot, cax=cax, extend="both", extendfrac='auto')

# Set colorbar ticks to cover a range instead of the color bounds
tick_len = (np.abs(args.lims[0]) + np.abs(args.lims[1]))/args.discretize
cbar.set_ticks(np.linspace(args.lims[0]+tick_len/2.0, args.lims[1]-tick_len/2.0, args.discretize))
ticklabellst = [str(round(args.lims[0]+i*tick_len,3)) + u" {} ".format(unichr(8211)) + str(round(args.lims[0]+i*tick_len+tick_len,3)) for i in xrange(args.discretize)]
cbar.set_ticklabels(ticklabellst)
# if (args.compare):
#   cbar.ax.set_title('$\Delta |\mathbf{'+args.variable+'}| \; \mathrm{[m/s]}$', loc="left")
# else:
cbar.ax.set_title('$\mathbf{'+args.variable+'} \; \mathrm{[m/s]}$', loc="left")
cbar.ax.tick_params(labelsize=10)

plt.tight_layout(w_pad=1.25, h_pad=1.25)

plt.show()
