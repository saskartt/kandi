#!/usr/bin/env python
import sys
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from kandiLib import *
from settings import *


'''
Plot time-averaged cross-sections from mask data.
'''

#==========================================================#

parser = argparse.ArgumentParser(prog='createCanopyRaster.py',
                                 description='''Plot time-averaged masks from output mask data.''')
parser.add_argument("-f", "--file", type=str, help="Name of the input netCDF4 file.")
parser.add_argument("-c", "--compare", type=str,
                    help="Compare to velocity field to another netCDF4 file.")
parser.add_argument("-var", "--variable", type=str, help="Variable to be plotted")
parser.add_argument("-cm", "--colormap", type=str, default="RdBu_r",
                    help="Colormap style. Default RdBu_r.")
parser.add_argument("-d", "--discretize", type=int, default=15,
                    help="Discretize color map to number of colors. Default: 15")
parser.add_argument("-l", "--lims", type=float, nargs=2, help="Colorbar limits.")
parser.add_argument("-y", "--ymax", type=float, help="Maximum (relative) height.")
parser.add_argument("-t", "--topography", type=str, help="Topography mask file.")
parser.add_argument("-q", "--quiver", type=str, help="File to read quiver vectors from.")
args = parser.parse_args()

#==========================================================#


# Read data from the main dataset
ds = openDataSet(args.file)

t_inds, = np.where(ds.variables['time'][:] >= skip_time_avg)
var = ds.variables[args.variable][t_inds, :, :, :]
# Set Dimensions
x_dims = ds.variables["y"][:]  # plot x dimension
y_dims = ds.variables["zw_3d"][:]  # plot y dimension

ds.close()

# Read a slice of topography from topography dataset
if (args.topography):
  tds = openDataSet(args.topography)
  topo = tds.variables["buildings_0"][:]
  topo = topo[:, :, 64]  # Should not be hardcoded
  topo_mask = np.ma.masked_where(topo == 0, topo)
  # Offset topography by 0.5 m in plot x direction
  topo_x_dims = tds.variables["y"][:] + 0.5
  topo_y_dims = tds.variables["z"][:]
  tds.close()

# Read interpolated vector data field from another file (output of groupVectorData.py)
if (args.quiver):
  qds = openDataSet(args.quiver)
  v = qds.variables["v"][t_inds, :, :, :]
  v = np.mean(v, 0)
  v = np.mean(v, 2)
  w = qds.variables["w"][t_inds, :, :, :]
  w = np.mean(w, 0)
  w = np.mean(w, 2)
  quiver_x_dims = qds.variables["y"][:]
  quiver_y_dims = qds.variables["z"][:]
  qds.close()

# Time averaging and averaging along x
var = np.mean(var, 0)
var = np.mean(var, 2)

# Interpolate variable to topography resolution
x_dims_i = np.arange(x_dims[0], x_dims[-1], 1)  # dimensions to interpolate to
y_dims_i = np.arange(y_dims[0], y_dims[-1], 1)
var = interpolateScalarField(var, x_dims, y_dims, x_dims_i, y_dims_i)

# Calculate comparsion and do the same averaging
if (args.compare):
  cds = openDataSet(args.compare)
  cvar = cds.variables[args.variable][t_inds, :, :, :]
  # Time averaging and averaging along x
  cvar = np.mean(cvar, 0)
  cvar = np.mean(cvar, 2)
  # Dimensions needed for interpolation
  cx_dims = ds.variables["y"][:]
  cy_dims = ds.variables["zw_3d"][:]
  cvar = interpolateScalarField(cvar, cx_dims, cy_dims, x_dims_i, y_dims_i)
  var = np.abs(cvar) - np.abs(var)
  cds.close()

# Set custom axes (ticks relative to building height)
# NOTE: These should not be hardcoded
hlen = bld_height  # This could be detected from the topography
xlen = x_dims[-1]
ylen = y_dims[-1]

xticklabels = np.arange(-1.5, 2.0, 0.5)
xticks = np.arange(xlen / 2 - 1.5 * hlen - 0.5, xlen / 2 + 2 * hlen - 0.5, hlen / 2)
yticks = np.arange(0.0, 3.0 * hlen + hlen / 2, hlen / 2)
yticklabels = np.arange(0.0, 3.0, 0.5)

# Set colormap
cmap = mpl.cm.get_cmap(args.colormap)

# Construct the figure
fig = plt.figure()
ax = plt.gca()
plt.xticks(xticks, xticklabels)
plt.yticks(yticks, yticklabels)
ax.set_ylim([1, args.ymax * hlen])
# NOTE: Should not be hardcoded
ax.set_xlim([xlen / 2 - 1.7 * hlen - 0.5, xlen / 2 + 1.7 * hlen - 0.5])
plt.axes().set_aspect('equal', 'datalim')

# Color levels
clevels = np.linspace(args.lims[0], args.lims[1], args.discretize + 1)

# Plot filled contour velocity field
imgplot = plt.contourf(x_dims_i, y_dims_i, var, levels=clevels, cmap=cmap, origin='lower')

# Plot topohraphy
if (args.topography):
  topoplot = plt.contourf(topo_x_dims, topo_y_dims, topo_mask,
                          cmap=mpl.cm.binary, vmin=0, vmax=1, interpolation=None)

# Plot vector arrows
if (args.quiver):
  quiverplot = plt.quiver(quiver_x_dims, quiver_y_dims, v, w,
                          color="k", scale=1, units='xy', minshaft=3)
  # quiverplot = plt.streamplot(quiver_x_dims, quiver_y_dims, v, w, color="k", linewidth=1.5, arrowsize=2)

# Add colormap
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(imgplot, cax=cax, extend="both", extendfrac='auto')
# fig.savefig( "dundun.eps", format='eps', dpi=600, bbox_inches='tight')
plt.show()
