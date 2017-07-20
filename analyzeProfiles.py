#!/usr/bin/env python
import sys
import argparse
import numpy as np
import tabulate as tb
from settings import *
from kandiLib import *

'''
Calculate various time-averaged quatntities from profiles.
'''

#==========================================================#

parser = argparse.ArgumentParser(
    prog='analyzeProfiles.py', description='''Calcluates various time-averaged quantities from PALM output NetCDF4 profiles. More settings available in settings.py file.''')
parser.add_argument("-f", "--file", type=str, default=None,
                    help="Name of the input NetCDF4 profile file.")
parser.add_argument("-d", "--domains", type=str, nargs='+',
                    default=['00'], help="Statistical domains to process. Default: 00")
parser.add_argument("-c", "--compare", type=str, default=None,
                    help="Print out a comparsion of the results to the results of given NetCDF4 profile file.")
parser.add_argument("-v", "--variables", type=str, nargs='+',
                    default=['ws', 'var', 'flux', 'tke'], help="Variables to be plotted.")
args = parser.parse_args()

#==========================================================#

# Read in the dataset
ds = openDataSet(args.file)

# Skip time from beginning
t_inds, = np.where(ds.variables['time'][:] >= skip_time_avg)

if (not(args.compare)):
  print("#=============== Time-averaged variables for dataset {} ===============#".format(args.file))
  for domain in args.domains:
    print("***** Statistical domain {} *****".format(domain))
    qty_list = compileDataListAverages(domain, t_inds, pr_heights_table, ds, args.variables)
    print(tb.tabulate(qty_list, headers=pr_table_header) + "\n")

else:
  # Compare to another dataset
  cds = openDataSet(args.compare)
  ct_inds, = np.where(cds.variables['time'][:] >= skip_time_avg)
  print("#=============== Comparing dataset {} against {} ===============#".format(
      args.file, args.compare))
  for domain in args.domains:
    print("***** Statistical domain {} *****".format(domain))
    qty_list = compileDataListCompare(domain, t_inds, ct_inds, pr_heights_table, ds, cds, args.variables)
    print(tb.tabulate(qty_list, headers=pr_table_header) + "\n")
