#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse
import numpy as np
from scipy.ndimage import interpolation
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from PIL import Image
from kandiLib import *
from settings import *


'''
Plot Hovmöller diagram.
'''

#==========================================================#

parser = argparse.ArgumentParser(prog='plotHovmollerDiagram.py',
                                 description='''Plot Hovmöller diagram.''')
parser.add_argument("-f", "--file", type=str, help="Name of the input netCDF4 file(s).")
parser.add_argument("-var", "--variable", type=str, help="Variable to be plotted")
args = parser.parse_args()

#==========================================================#


ds = openDataSet(args.file)

var, xdims, __, __ = readVariableFromMask(ds,[0., 7200.], args.variable)
t_inds, = np.where(np.logical_and(ds.variables['time'][:] >= 0., ds.variables['time'][:] <= 7200.))
var=var[:,25,0,:]
tdims=ds.variables['time'][t_inds]
tdims = interpolation.zoom(tdims,0.25)
xdims = interpolation.zoom(xdims,0.25)
var = interpolation.zoom(var,0.25)
X,Y = np.meshgrid(xdims,tdims)
cmap = mpl.cm.get_cmap("viridis")
bounds = np.linspace(0,11,10)
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
print(np.amax(var))
cp=plt.contourf(X,Y,var,8,vmin=0., vmax=11.,cmap=cmap)
plt.ylabel("Time [s]")
plt.xlabel("x [m]")
cbar = plt.colorbar(cp,ticks=bounds)
plt.title(u"x-directional diagram of variable {}".format(args.variable))
plt.show()
