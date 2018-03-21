#!/usr/bin/env python
import sys
import argparse
import numpy as np
from scipy.interpolate import interp2d as ip2d
from scipy.interpolate import RegularGridInterpolator as rgi
from scipy.optimize import curve_fit
import matplotlib as mpl
import matplotlib.pyplot as plt
from settings import *
from kandiLib import *

'''
Plot 2D contours of logarithmic fits.
'''

#==========================================================#

parser = argparse.ArgumentParser(
    prog='plotProfiles.py', description='''Plot time-averaged profiles.''')
parser.add_argument("-f", "--file", type=str, default=None,
                    help="Name of the input netCDF4 file containing horizontal slices of wind data.")
parser.add_argument("-p", "--variable", type=str, default="residual", help="Variable to be plotted. Options: residual, z0, zd")
parser.add_argument("-u*", "--fricvel", type=float, help="Friction velocity of the run.")
parser.add_argument("-s", "--save", type=str, help="Save resulting figure as.")
parser.add_argument("-ft", "--fit", type=int, nargs=2, default=[30,60], help="Vertical range of the fit in meters.")
parser.add_argument("-gp", "--gryningProfile", help="Use Gryning profile instead of logprofile.",
                    action="store_true", default=False)
args = parser.parse_args()

#=========================================================#

# Read in the dataset
print("Reading the dataset...")
ds = openDataSet(args.file)
t_inds, = np.where(np.logical_and(ds.variables['time'][:] >= timespan[0], ds.variables['time'][:] <= timespan[1]))

# Read 2D wind vectors and average temporally
u_field, xu_dims, yu_dims, zu_dims = readVariableFromXYSlice(ds, timespan, 'u_xy')
u_field = calculateTemporalStatistics(u_field, 'avg')
u_field = np.swapaxes(u_field,0,2) # [z,y,x] -> [x,y,z]
v_field, xv_dims, yv_dims, zv_dims = readVariableFromXYSlice(ds, timespan, 'v_xy')
v_field = calculateTemporalStatistics(v_field, 'avg')
v_field = np.swapaxes(v_field,0,2) # [z,y,x] -> [x,y,z]

print("Interpolating wind vectors...")
# 2-dimensional interpolation of the two vector fields into same grid
xGridPoints = np.arange(np.ceil(np.amin(xv_dims)),np.amax(xu_dims)+1.0,1)
yGridPoints = np.arange(np.ceil(np.amin(yu_dims)),np.amax(yv_dims)+1.0,1)
zGridPoints = np.arange(np.ceil(np.amin(zu_dims)),np.floor(np.amax(zu_dims))+1.0,1)


# Create a grid to interpolate to
ipGridX, ipGridY, ipGridZ = np.meshgrid(xGridPoints,yGridPoints,zGridPoints)

# bounds_error=False extrapolates outside the original grid (max 0.5m extrapolation)
ip3DFuncU = rgi(points=[xu_dims,yu_dims,zu_dims],values=u_field,bounds_error=False)
u = ip3DFuncU((ipGridX,ipGridY,ipGridZ))
u_field = None # Clear out old field
ip3DFuncV = rgi(points=[xv_dims,yv_dims,zv_dims],values=v_field,bounds_error=False)
v = ip3DFuncV((ipGridX,ipGridY,ipGridZ))
v_field = None
print("...done.")

print("Calculating norms of wind vectors...")
# Euclidean norm to calculate total horizontal wind
hWind = np.sqrt(u**2+v**2)
print("...done.")

'''
Calculate logarithmic profile column-by-column
'''
print("Calculating logarithmic profiles column-by-column...")
fricVel = args.fricvel
fitIndices = np.where(np.logical_and(zGridPoints >= args.fit[0], zGridPoints <= args.fit[1]))
zFitRange = zGridPoints[fitIndices]
hWindFit=hWind[:,:,fitIndices]
z0 = np.zeros((len(xGridPoints),len(yGridPoints)))
z0StdDev = np.zeros((len(xGridPoints),len(yGridPoints)))
z0Residual = np.zeros((len(xGridPoints),len(yGridPoints)))
zd = np.zeros((len(xGridPoints),len(yGridPoints)))
np.seterr(invalid='ignore') # Ignore invalid logarithm values

if (not(args.gryningProfile)):
  funcLogProfile = lambda z,a,b : (fricVel/0.4)*np.log((z-b)/a)

  for ix in np.arange(len(xGridPoints)):
    for iy in np.arange(len(yGridPoints)):
        hWindColumn, = hWindFit[ix,iy,:]
        fitSolution, pcov = curve_fit(funcLogProfile,zFitRange,hWindColumn)
        z0[ix,iy]=fitSolution[0]
        z0StdDev[ix,iy]=np.sqrt(np.diag(pcov))[0]/fitSolution[0]
        z0Residual[ix,iy]=np.sum((hWindColumn-funcLogProfile(zFitRange,fitSolution[0],fitSolution[1]))**2)
        zd[ix,iy]=fitSolution[1]
else:
    h=320.
    f=1e-4
    funcLogProfile = lambda z,a,b : (fricVel/0.4)*(np.log((z-b)/a)+(z-b)/(h/(2*(((np.log(fricVel/(1e-4*a))-1.9)**2+4.9**2)**0.5-np.log(h/a))))-((z-b)/h)*((z-b)/(h/(((np.log(fricVel/(1e-4*a))-1.9)**2+4.9**2)**0.5-np.log(h/a)))))

    for ix in np.arange(len(xGridPoints)):
       for iy in np.arange(len(yGridPoints)):
           hWindColumn, = hWindFit[ix,iy,:]
           fitSolution, pcov = curve_fit(funcLogProfile,zFitRange,hWindColumn)
           z0[ix,iy]=fitSolution[0]
           z0StdDev[ix,iy]=np.sqrt(np.diag(pcov))[0]/fitSolution[0]
           z0Residual[ix,iy]=np.sum((hWindColumn-funcLogProfile(zFitRange,fitSolution[0],fitSolution[1]))**2)
           zd[ix,iy]=fitSolution[1]


print("...done.")
print("Plotting results...")
plt.figure()
plt.contourf(z0,10)
plt.title("$z_0$")
plt.colorbar()
print("z_0 avg: {}".format(np.mean(z0)))

hWindColumn, = np.mean(hWindFit,(0,1))
fitSolution, pcov = curve_fit(funcLogProfile,zFitRange,hWindColumn)
pr=funcLogProfile(zGridPoints,fitSolution[0],fitSolution[1])
plt.figure()
plt.plot(np.mean(hWind,(0,1)),zGridPoints)
plt.plot(pr,zGridPoints)


plt.figure()
plt.contourf(z0StdDev,10)
plt.title("$\sigma_{z_0} /z_0$")
plt.colorbar()
print("Mean residual: {}".format(np.mean(z0Residual)))
plt.figure()
plt.contourf(z0Residual,levels=np.arange(0.0,0.2,0.01))
plt.title("Least squares fit residual")
plt.colorbar()

plt.figure()
plt.contourf(zd,10)
plt.title("$z_d$")
plt.colorbar()
plt.show()
