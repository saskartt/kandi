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
Plot time-averaged vertical cross-sections from mask data (following M. Hefny Salim et al / J. Wind. Eng. Ind. Aerodyn. 144 (2015) 84-95)
'''

#==========================================================#

parser = argparse.ArgumentParser(prog='plotVSlice.py',
                                 description='''Plot time-averaged masks from output mask data. Some settings are defined in settings.py''')
parser.add_argument("-f", "--file", nargs='+', type=str, help="Name of the input netCDF4 file(s).")
parser.add_argument("-c", "--compare", nargs='+', type=str,
                    help="Compare to velocity field to another netCDF4 file.")
parser.add_argument("-var", "--variable", type=str, help="Variable to be plotted")
parser.add_argument("-stat", "--statistics", type=str, default="avg", help="Statistics to be displayed. Options: avg, max, min")
parser.add_argument("-cm", "--colormap", type=str, default="RdBu_r",
                    help="Colormap style. Default RdBu_r.")
parser.add_argument("-d", "--discretize", type=int, default=15,
                    help="Discretize color map to number of colors. Default: 15")
parser.add_argument("--lims", type=float, nargs=2, help="Colorbar limits.")
parser.add_argument("--ymax", type=float, help="Maximum (relative) height.")
parser.add_argument("-t", "--topography", type=str, help="Topography mask file.")
parser.add_argument("-q", "--quiver", type=str, help="File to read quiver vectors from.")
parser.add_argument("-qt", "--quivertop", type=float, help="Draw quiver plot only to given level.")
parser.add_argument("-qs", "--quiverscale", type=float, default=1, help="Scaling factor for the quiver plot arrows.")
parser.add_argument("--canopy", type=str, help="Read single tree image file. The script locates this raster tree file into the plot according to parameters in settings file.")
parser.add_argument("-ca", "--canopyalpha", default=0.5, type=float, help="Alpha (transparency) used for canopy layer with 0 being completely transparent and 1 being opaque.")
parser.add_argument("-s", "--save", type=str,help="Save the resulting figure as an PDF file.")

args = parser.parse_args()

#==========================================================#

# Required to correctly save as vector graphics
mpl.rcParams['image.composite_image'] = False


# Read data from the main datasets
ds = openDataSet(args.file[0])

# Read variable and dimensions, dims correspond to plot dimensions
var, __, x_dims, y_dims = readVariableFromMask(ds, timespan, args.variable)
# y_dims = y_dims+0.5
# h_index, = np.where(y_dims <= bld_height)[0]
# var = var[:,h_index,:,:]
var = calculateTemporalStatistics(var, args.statistics)
var = np.mean(var, 2)
ds.close()

# Read in second dataset and average with first one
if (len(args.file)>1):
  for filename in args.file[1:]:
    ds = openDataSet(filename)
    dsvar, __, __, __ = readVariableFromMask(ds, timespan, args.variable)
    # Cut out a slice
    # dsvar = dsvar[:,h_index,:,:]
    dsvar = calculateTemporalStatistics(dsvar, args.statistics)
    dsvar = np.mean(dsvar, 2)
    var=(var+dsvar)/2

# Read a slice of topography from topography dataset
if (args.topography):
  tds = openDataSet(args.topography)
  topo = tds.variables["buildings_0"][:]
  topo = topo[:, :, 64]  # Should not be hardcoded
  topo_mask = np.ma.masked_where(topo == 0, topo)
  # Offset topography by 0.5 m in plot x direction
  topo_x_dims = tds.variables["y"][:] +0.5
  topo_y_dims = tds.variables["z"][:]
  tds.close()

# Read interpolated vector data field from another file (output of groupVectorData.py)
if (args.quiver):
  qds = openDataSet(args.quiver)
  # There might be differences in dt_domask between runs
  qt_inds, = np.where(np.logical_and(qds.variables['time'][:] >= timespan[0], qds.variables['time'][:] <= timespan[1]))
  v = qds.variables["v"][qt_inds, :, :, :]
  v = np.mean(v, 0)
  v = np.mean(v[:,:,0:30], 2)
  w = qds.variables["w"][qt_inds, :, :, :]
  w = np.mean(w, 0)
  w = np.mean(w[:,:,0:30], 2)
  quiver_x_dims = qds.variables["y"][:]

  # Crop to user-defined top level
  if (args.quivertop):
    q_inds, = np.where(qds.variables["z"][:] <= args.quivertop)
    quiver_y_dims = qds.variables["z"][q_inds]
    v = v[q_inds,:]
    w = w[q_inds,:]
  else:
    quiver_y_dims = qds.variables["z"][:]
  # Calculate maximum arrow length for matplotlib quiver minshaft length (scales vector arrows)
  q_lengths = np.hypot(v,w)
  # minshaft = np.amax(q_lengths)
  minshaft = 3.5 # Arbitrary for now
  # Remove the smallest vector arrows because for some reason matplotlib shows them as big dots
  v[q_lengths <= minshaft/20.*args.quiverscale] = np.nan
  w[q_lengths <= minshaft/20.*args.quiverscale] = np.nan
  qds.close()

# Interpolate variable to topography resolution
x_dims_i = np.arange(x_dims[0], x_dims[-1], 1)  # dimensions to interpolate to
y_dims_i = np.arange(y_dims[0], y_dims[-1], 1)
var = interpolateScalarField(var, x_dims, y_dims, x_dims_i, y_dims_i)

# Calculate comparsion and do the same averaging
if (args.compare):
  cds = openDataSet(args.compare[0])
  cvar, __, cx_dims, cy_dims = readVariableFromMask(cds, timespan, args.variable)
  cvar = calculateTemporalStatistics(cvar, args.statistics)
  cvar = np.mean(cvar, 2)
  cds.close()

  # Read in second dataset and average with first one
  if (len(args.compare)>1):
    for filename in args.file[1:]:
      cds = openDataSet(filename)
      cdsvar, __, __, __ = readVariableFromMask(cds, timespan, args.variable)
      cdsvar = calculateTemporalStatistics(cdsvar, args.statistics)
      cdsvar = np.mean(cdsvar, 2)
      cvar=(cvar+cdsvar)/2

  cvar = interpolateScalarField(cvar, cx_dims, cy_dims, x_dims_i, y_dims_i)
  var = (var-cvar)/u_ref
  cds.close()

# Set custom axes (ticks relative to building height)
# NOTE: These should not be hardcoded
hlen = bld_height  # This could be detected from the topography
xlen = 128.
ylen = 128.

xticklabels = np.arange(0, 4.5, 0.5)
xticks = np.array([0.5])
xticks = np.append(xticks,np.arange(hlen/2., xlen/2.0, hlen / 2))
xticks = np.append(xticks,xlen/2.-0.5) # Offset last element by -0.5 m in order it to show
yticks = np.arange(0.0, args.ymax * hlen + hlen / 2, hlen / 2)
yticklabels = np.arange(0.0, args.ymax+1.0, 0.5)

# Set colormap
cmap = mpl.cm.get_cmap(args.colormap)

# Construct the figure
fig = plt.figure(figsize=(10,4))
ax = plt.gca()
plt.xticks(xticks, xticklabels)
plt.yticks(yticks, yticklabels)
ax.set_ylim([0, args.ymax * hlen])
# NOTE: Should not be hardcoded
ax.set_xlim([7.5, xlen/2. -7.5 ])
plt.axes().set_aspect('equal', 'datalim')

ax.set_xlabel('$x/H$ [-]')
ax.set_ylabel('$z/H$ [-]')

# Color levels
clevels = np.linspace(args.lims[0], args.lims[1], args.discretize + 1)

# Plot filled contour velocity field
imgplot = plt.contourf(x_dims_i, y_dims_i, var, levels=clevels, cmap=cmap, origin='lower',zorder=0)

# Plot topography
if (args.topography):
  topoplot = plt.contourf(topo_x_dims, topo_y_dims, topo_mask,
                          cmap=mpl.cm.binary, vmin=0, vmax=1, interpolation=None, zorder=1)

# Plot vector arrows
if (args.quiver):
  quiverplot = plt.quiver(quiver_x_dims, quiver_y_dims, v, w,
                          color="k", scale=args.quiverscale, units='xy', minshaft=minshaft, pivot="tail", zorder=3)
  # quiverplot = plt.streamplot(quiver_x_dims, quiver_y_dims, v, w, color="k", linewidth=1.5, arrowsize=2)

# Add raster tree images to plot
if (args.canopy):
  for ntree in xrange(len(tree_location)):
    try:
      img = Image.open(args.canopy)
    except RuntimeError:
      raise IOError("Input file {} not found!".format(args.canopy))
    imgarr = np.array(img)
    t_left = tree_location[ntree] - tree_radius[ntree]
    t_right =  tree_location[ntree] + tree_radius[ntree]
    t_top = tree_height[ntree]
    plt.imshow(imgarr,extent=[t_left, t_right, 0.0, t_top],alpha=args.canopyalpha,interpolation=None, zorder=2)

# Add colormap
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.15)
cbar = fig.colorbar(imgplot, cax=cax, extend="both", extendfrac='auto')

# Set colorbar ticks to cover a range instead of the color bounds
tick_len = (np.abs(args.lims[0]) + np.abs(args.lims[1]))/args.discretize
cbar.set_ticks(np.linspace(args.lims[0]+tick_len/2.0, args.lims[1]-tick_len/2.0, args.discretize))
ticklabellst = [str(round(args.lims[0]+i*tick_len,3)) + u" {} ".format(unichr(8211)) + str(round(args.lims[0]+i*tick_len+tick_len,3)) for i in xrange(args.discretize)]
cbar.set_ticklabels(ticklabellst)
if (args.compare):
  cbar.ax.set_title('$\Delta '+args.variable+' /u^{*}_{ref} $', loc="left")
  plt.suptitle(args.file[0].split("/")[0] + " vs. " + args.compare[0].split("/")[0])
else:
  cbar.ax.set_title('$\mathbf{'+args.variable+'} \; \mathrm{[m/s]}$', loc="left")
cbar.ax.tick_params(labelsize=10)

plt.tight_layout(w_pad=1.25, h_pad=1.25)

if (args.save):
  fig.savefig( args.save, format='pdf', dpi=600, bbox_inches='tight')
  print("Figure {} saved.".format(args.save))
# plt.show()
