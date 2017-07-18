#!/usr/bin/env python
import sys
import argparse
import numpy as np
import netCDF4 as nc
import tabulate as tb
from settings import *
from kandiLib import *

'''
Calculate various time-averaged quatntities from profiles.
'''

#==========================================================#

parser = argparse.ArgumentParser(prog='analyzeProfiles.py', description='''Calcluates various time-averaged quantities \
                                 from PALM output NetCDF4 profiles. More settings available in settings.py file.''')
parser.add_argument("-f", "--file", type=str, default=None, help="Name of the input NetCDF4 profile file.")
parser.add_argument("-d", "--domains", type=str, nargs='+', default=['00'], help="Statistical domains to process.\
                    Default: 00")
parser.add_argument("-c", "--compare", type=str, default=None, help="Print out a comparsion of the results to the\
                    results of given NetCDF4 profile file.")
args = parser.parse_args()

#==========================================================#

# Read in the dataset
try:
  ds = nc.Dataset(args.file)
except RuntimeError:
  raise IOError("Input file {} not found!".format(args.file))

# Skip first hour of the simulation
timepoints, = np.where(ds.variables['time'][:]>=skip_time_avg)

sep = [[''] * (len(pr_heights)+1)] # empty line string


if (not(args.compare)):
  print("#=============== Time-averaged variables for dataset {} ===============#".format(args.file))
  for domain in args.domains:
    print("***** Statistical domain {} *****".format(domain))
    pr_01 = averageProfilesWS(domain, timepoints, pr_heights, ds)
    pr_02 = averageProfilesVariances(domain, timepoints, pr_heights,ds)
    pr_03 = averageProfilesMomentumFluxes(domain, timepoints, pr_heights, ds)
    pr_04 = averageProfilesTKE(domain, timepoints, pr_heights, ds)
    qty_list = np.concatenate((pr_01, sep, pr_02, sep, pr_03, sep, pr_04), axis=0)
    print(tb.tabulate(qty_list, headers=pr_header)+"\n")

else:
  # Compare to another dataset
  try:
    cds = nc.Dataset(args.compare)
  except RuntimeError:
    raise IOError("Input file {} not found!".format(args.compare))

  ctimepoints, = np.where(cds.variables['time'][:]>=skip_time_avg)
  print("#=============== Comparing dataset {} against {} ===============#".format(args.file, args.compare))
  for domain in args.domains:
    print("***** Statistical domain {} *****".format(domain))
    pr_01 = compareProfilesWS(domain, timepoints, ctimepoints, pr_heights, ds, cds)
    pr_02 = compareProfilesVariances(domain, timepoints, ctimepoints, pr_heights, ds, cds)
    pr_03 = compareProfilesMomentumFluxes(domain, timepoints, ctimepoints, pr_heights, ds, cds)
    pr_04 = compareProfilesTKE(domain, timepoints, ctimepoints, pr_heights, ds, cds)
    qty_list = np.concatenate((pr_01, sep, pr_02, sep, pr_03, sep, pr_04), axis=0)
    print(tb.tabulate(qty_list, headers=pr_header)+"\n")
