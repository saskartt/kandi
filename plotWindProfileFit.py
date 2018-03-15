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
Plot fit for wind profiles
'''

#==========================================================#

parser = argparse.ArgumentParser(
    prog='plotWindProfileFit.py', description='''Plot time-averaged profiles.''')
parser.add_argument("-f", "--files", type=str, nargs='+',default=None,
                    help="Name of the input netCDF4 files.")
parser.add_argument("-d", "--domains", type=str, nargs='+',
                    default=['0'], help="Statistical domains to process. Default: 00")
parser.add_argument("-s", "--save", type=str, help="Save resulting figure as.")
parser.add_argument("-ft", "--fit", type=int, nargs=2, default=[30,60], help="Range of vertical grid points to fit to.")
parser.add_argument("-x", "--xlims", type=float, nargs=2, help="Set x axis limits manually.")
parser.add_argument("-y", "--ylims", type=float, nargs=2, help="Set y axis limits manually.")
parser.add_argument("-pr", "--profile", help="Profile to be used", type=str, default="log")
parser.add_argument("-blh", "--blh", type=float, help="Boundary layer height for the Gryning profile, default is fricVel/(12*1e-4)")
args = parser.parse_args()

#==========================================================#

mpl.rcParams["mathtext.fontset"] ="cm"

# Read all datasets into a list
dsList = []; tpList = {}; nameList = {}
for fname in args.files:
  ds = openDataSet(fname)
  nameList[ds] = fname
  dsList.append(ds)
  t_inds, = np.where(np.logical_and(ds.variables['time'][:] >= timespan[0], ds.variables['time'][:] <= timespan[1]))
  tpList[ds] = t_inds

plt.figure(1)
plt.grid()
axes = plt.gca()
if (args.xlims):
  axes.set_xlim(args.xlims)
if (args.ylims):
  axes.set_ylim(args.ylims)
else:
  axes.set_ylim([0, 128])
plt.ylabel("$z\/\mathrm{(m)}$",fontsize=14)
color_cycle = ['b', 'g', 'r', 'c', 'm', 'y', 'fuchsia', 'gold', 'orange', 'lightcoral', 'lightslategrey','tan']
i=0

for ds in dsList:
  for domain in args.domains:
    flux1=averageProfilesMomentumFluxes(domain, tpList[ds], pr_heights_plot, ds)[0]
    flux2=averageProfilesMomentumFluxes(domain, tpList[ds], pr_heights_plot, ds)[1]
    flux = (flux1**2.0 + flux2**2.0)**0.25
    fricVel= np.mean(flux[args.fit[0]:args.fit[1]])

    datalist=averageProfilesWS(domain, tpList[ds], pr_heights_plot, ds)
    hwind = datalist[3]
    plt.xlabel("$\mathbf{u}/\mathbf{u_*}\/\mathrm{(m/s)}$",fontsize=14)
    plt.plot(hwind,pr_heights_plot, label=r'Run: {}, simulated'.format(nameList[ds][4], domain),color=color_cycle[i])

    z=pr_heights_plot[args.fit[0]:args.fit[1]]
    u_profile = hwind[args.fit[0]:args.fit[1]]

    if (args.profile=="log"):
      funcLogProfile = lambda z,a,b : (fricVel/0.4)*np.log((z-b)/a)
    elif (args.profile=="blackadar"):
      zi=1000. # Wind maximum
      eta=63e-4*fricVel/1e-4 # eta=kb1*u*/f_c
      funcLogProfile = lambda z,a,b : (fricVel/0.4)*(np.log((z-b)/a)+(0.4*(z-b)/eta)-(z/zi)*((0.4*z)/(2*eta)+1))
    elif (args.profile=="gryning"):
      if(args.blh):
        h=args.blh
      else:
        h=fricVel/(12*1e-4)
      funcLogProfile = lambda z,a,b : (fricVel/0.4)*(np.log((z-b)/a)+(z-b)/(fricVel/(2*1e-4*np.log(fricVel/(1e-4*a))+55*1e-4))-(z-b)/(h)*((z-b)/(fricVel/(1e-4*np.log(fricVel/(1e-4*a))+55*1e-4))))

    fitSolution, pcov = curve_fit(funcLogProfile,z,u_profile)
    pr=funcLogProfile(pr_heights_plot,fitSolution[0],fitSolution[1])
    print("Least squrares fit: {}".format(fitSolution))
    np.seterr(invalid='ignore')
    # plt.plot(pr,pr_heights_plot, label=r'Run: {}, logprofile'.format(nameList[ds][4:], domain))

    if (args.ylims):
      axes.set_ylim([args.ylims[0],args.ylims[1]])
    if (args.xlims):
      axes.set_xlim([args.xlims[0],args.xlims[1]])

    plt.plot(pr,pr_heights_plot, label=r'Run: {}, log profile'.format(nameList[ds][4]), linestyle='--', color=color_cycle[i])
    i=i+1

#axes.fill_between(np.linspace(0,12.0), 16, 32, facecolor='yellow', alpha=0.3,
#                label='Roof level < h < 0.3*BLH')
leg = plt.legend(loc=0, fontsize=9)
for legobj in leg.legendHandles:
    legobj.set_linewidth(2.0)

if (args.save):
  plt.savefig(args.save)
  print("Figure {} saved.".format(args.save))
plt.show()
