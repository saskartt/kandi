#!/usr/bin/env python
import sys
import argparse
import numpy as np
from scipy.optimize import curve_fit
import matplotlib as mpl
mpl.use('GTKCairo')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib import rc
from itertools import cycle
from settings import *
from kandiLib import *

'''
Plot z_0 sensitivity.
'''

#==========================================================#

parser = argparse.ArgumentParser(
    prog='plotProfiles.py', description='''Plot time-averaged profiles.''')
parser.add_argument("-f", "--file", type=str, nargs='+',default=None,
                    help="Name of the input netCDF4 file.")
parser.add_argument("-i", "--interval", type=int, nargs=2, default=[15,96], help="Maximum and minimum fit values")
parser.add_argument("-s", "--stride", type=int, help="Stride of fit.")
args = parser.parse_args()

#==========================================================#

interval=args.interval
stride=args.stride

ds = openDataSet(args.files)
t_inds, = np.where(np.logical_and(ds.variables['time'][:] >= timespan[0], ds.variables['time'][:] <= timespan[1]))
plt.figure(1)
plt.grid()
axes = plt.gca()

flux1=averageProfilesMomentumFluxes(0, tpList[ds], pr_heights_plot, ds)[0]
flux2=averageProfilesMomentumFluxes(0, tpList[ds], pr_heights_plot, ds)[1]
flux = (flux1**2.0 + flux2**2.0)**0.25
hwind=averageProfilesWS(domain, tpList[ds], pr_heights_plot, ds)[3]

z0=np.array([])

for start in xrange(interval[0],interval[1]-stride):
  fric_vel= np.mean(flux[start:start+stride])
  z=pr_heights_plot[start:start+stride]
  u_profile = hwind[start:start+stride]
  funcLogProfile = lambda val,z : (fric_vel/0.4)*np.log((z-val[1])/val[0])
  ErrorFunc = lambda val,z,pr:  np.abs(funcLogProfile(val,z)-u_profile)
  valInitial=(1.0,0.0)
  valFinal,success = leastsq(ErrorFunc,valInitial[:],args=(z,u_profile))
  z0.append(valFinal[0])


plt.plot(np.arange(interval[0],interval[1]-stride),z0)
