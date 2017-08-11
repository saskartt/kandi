#!/usr/bin/env python

import sys
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from kandiLib import *
from settings import *

'''
Calculate and plot velocity reduction coefficients (following M. Hefny Salim et al / J. Wind. Eng. Ind. Aerodyn. 144 (2015) 84-95)
'''

#==========================================================#

parser = argparse.ArgumentParser(prog='plotReductionCoefficients.py',
                                 description='''Plot time-averaged reduction coefficients for every grid point and ''')
parser.add_argument("-f", "--file", type=str, help="Name of the input netCDF4 file.")
parser.add_argument("-c", "--compare", type=str, help="A netCDF4 file to compare to.")
parser.add_argument("-var", "--variable", type=str, default="w", help="A variable to be processed.")
parser.add_argument("-x", "--xlims", type=float, nargs=2, help="Clip the mask into specified limits in x-direction.")
parser.add_argument("-y", "--ylims", type=float, nargs=2, help="Clip the mask into specified limits in y-direction.")
parser.add_argument("-z", "--zlims", type=float, nargs=2, help="Clip the mask into specified limits in z-direction.")
args = parser.parse_args()

#==========================================================#

# Read in the tree-free case
cds = openDataSet(args.compare)
carr, cx_dims, cy_dims, cz_dims = readVariableFromMask(cds, skip_time_avg, args.variable)
print(np.shape(carr))
# Clip desired domain
carr = clipMask(carr, cx_dims, cy_dims, cz_dims, args.xlims, args.ylims, args.zlims)
carr = calculateTemporalStatistics(carr, "avg")
carr = carr.flatten()

# Read dataset(s)
ds = openDataSet(args.file)
arr, x_dims, y_dims, z_dims = readVariableFromMask(ds, skip_time_avg, args.variable)
arr = clipMask(arr, x_dims, y_dims, z_dims, args.xlims, args.ylims, args.zlims)
arr = calculateTemporalStatistics(arr, "avg")
arr = arr.flatten()

fig = plt.figure()
ax = plt.gca()
plt.scatter(carr,arr)
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
ax.set_xlim([-0.4, 0.4])
ax.set_ylim([-0.4, 0.4])
plt.plot(np.linspace(-0.4,0.4),np.linspace(-0.4,0.4))
plt.show()
