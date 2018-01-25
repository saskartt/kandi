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
Plot 1-dimensional profile of one variable at given point in the domain grid.
'''

#==========================================================#

parser = argparse.ArgumentParser(
    prog='gridPointPlot1D.py', description='''Plot 1-dimensional profile of one variable at given point in the domain grid.''')
parser.add_argument("-f", "--file", type=str,default=None,
                    help="Name of the input netCDF4 file.")
parser.add_argument("-v", "--variable", type=str,
                    default="u", help="Variable to be plotted.")
parser.add_argument("-p", "--gridpoint", type=float, nargs=3, help="Grid point position [x y z]")
args = parser.parse_args()

#==========================================================#


ds = openDataSet(args.file)

var, x_dims, y_dims, z_dims = readVariableFromMask(ds, [0, 18000.], args.variable)
h_index, = np.where(z_dims == args.gridpoint[2])[0]
x_index, = np.where(x_dims == args.gridpoint[0])[0]
y_index, = np.where(y_dims == args.gridpoint[1])[0]
var = var[:,h_index,y_index,x_index]
print(np.shape(var))
plt.plot(var)

plt.show()

