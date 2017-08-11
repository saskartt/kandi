#!/usr/bin/env python
import sys
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('GTKCairo')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib import rc
from settings import *
from kandiLib import *

'''
Plot time-averaged profiles.
'''

#==========================================================#

parser = argparse.ArgumentParser(
    prog='plotProfiles.py', description='''Plot time-averaged profiles.''')
parser.add_argument("-f", "--files", type=str, nargs='+',default=None,
                    help="Name of the input netCDF4 files.")
parser.add_argument("-d", "--domains", type=str, nargs='+',
                    default=['0'], help="Statistical domains to process. Default: 00")
parser.add_argument("-v", "--variable", type=str,
                    default="u", help="Variable to be plotted.")
parser.add_argument("-s", "--save", type=str, help="Save resulting figure as.")
parser.add_argument("-x", "--xlims", type=float, nargs=2, help="Set x axis limits manually.")
parser.add_argument("-y", "--ylims", type=float, nargs=2, help="Set y axis limits manually.")
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
for ds in dsList:
  for domain in args.domains:

    if (args.variable == "u"):
      datalist=averageProfilesWS(domain, tpList[ds], pr_heights_plot, ds)
      pr = datalist[0]
      plt.xlabel("$\mathbf{u}\/\mathrm{(m/s)}$",fontsize=14)

    elif (args.variable == "v"):
      datalist=averageProfilesWS(domain, tpList[ds], pr_heights_plot, ds)
      pr = datalist[1]
      plt.xlabel("$\mathbf{v}\/\mathrm{(m/s)}$",fontsize=14)

    elif (args.variable == "w"):
      datalist=averageProfilesWS(domain, tpList[ds], pr_heights_plot, ds)
      pr = datalist[2]
      plt.xlabel("$\mathbf{w}\/\mathrm{(m/s)}$",fontsize=14)

    elif (args.variable == "U"):
      datalist=averageProfilesWS(domain, tpList[ds], pr_heights_plot, ds)
      pr = datalist[3]
      plt.xlabel("$\mathbf{U}\/\mathrm{(m/s)}$",fontsize=14)

    elif (args.variable == "u*2"):
      datalist=averageProfilesVariances(domain, tpList[ds], pr_heights_plot, ds)
      pr = datalist[0]
      plt.xlabel("$\sigma^2_{u}\/\mathrm{(m^2/s^2)}$",fontsize=14)

    elif (args.variable == "v*2"):
      datalist=averageProfilesVariances(domain, tpList[ds], pr_heights_plot, ds)
      pr = datalist[1]
      plt.xlabel("$\sigma^2_{v}\/\mathrm{(m^2/s^2)}$",fontsize=14)

    elif (args.variable == "w*2"):
      datalist=averageProfilesVariances(domain, tpList[ds], pr_heights_plot, ds)
      pr = datalist[2]
      plt.xlabel("$\sigma^2_{w}\/\mathrm{(m^2/s^2)}$",fontsize=14)

    elif (args.variable == "wu"):
      datalist=averageProfilesMomentumFluxes(domain, tpList[ds], pr_heights_plot, ds)
      pr = datalist[0]
      plt.xlabel("$\overline{u'w'}\/\mathrm{(m^2/s^2)}$",fontsize=14)

    elif (args.variable == "wv"):
      datalist=averageProfilesMomentumFluxes(domain, tpList[ds], pr_heights_plot, ds)
      pr = datalist[1]
      plt.xlabel("$\overline{v'w'}\/\mathrm{(m^2/s^2)}$",fontsize=14)

    elif (args.variable == "flux"):
      datalist=averageProfilesMomentumFluxes(domain, tpList[ds], pr_heights_plot, ds)
      pr = datalist[2]
      plt.xlabel(r"$\overline{u'w'}+ \overline{v'w'} \/\mathrm{(m^2/s^2)}$",fontsize=14)

    elif (args.variable == "tke"):
      datalist=averageProfilesTKE(domain, tpList[ds], pr_heights_plot, ds)
      pr = datalist[0]
      plt.xlabel("$TKE\/\mathrm{(m^2/s^2)}$",fontsize=14)

    elif (args.variable == "u*"):
      datalist=averageProfilesMomentumFluxes(domain, tpList[ds], pr_heights_plot, ds)
      pr1 = np.power(datalist[0],2.0)
      pr2 = np.power(datalist[1],2.0)
      pr = np.power(pr1+pr2,0.25)

    elif (args.variable == "z0"):
      fluxes=averageProfilesMomentumFluxes(domain, tpList[ds], pr_heights_plot, ds)
      flx1 = np.power(fluxes[0],2.0)
      flx2 = np.power(fluxes[1],2.0)
      fric_vel = np.power(flx1+flx2,0.25)
      print(fric_vel)

      hwind = datalist=averageProfilesWS(domain, tpList[ds], pr_heights_plot, ds)
      hwind = hwind[3]

      disp_height = 28.-((34.-28.)/(np.exp(0.41*((hwind[34]/fric_vel[34])-(hwind[28]/fric_vel[28])))-1))
      print(disp_height)
      uu = (hwind/fric_vel)
      pr = (pr_heights_plot-disp_height)*np.exp(-uu*0.41)
      pr_heights_plot = pr_heights_plot
      axes.set_xlim([0.0,0.5])
      axes.set_ylim([disp_height,60.0])
    else:
      raise NameError("Unknown variable "+args.variable)

    plt.plot(pr,pr_heights_plot, label=r'Run: {}, stat. domain: {}'.format(nameList[ds][4:], domain))

leg = plt.legend(loc=0, fontsize=9)
for legobj in leg.legendHandles:
    legobj.set_linewidth(2.0)

if (args.save):
  plt.savefig(args.save)
  print("Figure {} saved.".format(args.save))
plt.show()
