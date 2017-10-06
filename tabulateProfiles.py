#!/usr/bin/env python
import sys
import argparse
import numpy as np
import tabulate as tb
from settings import *
from kandiLib import *

'''
Calculate various time-averaged quantities from profiles.
'''

#==========================================================#

parser = argparse.ArgumentParser(
    prog='tabulateProfiles.py', description='''Calcluates various time-averaged quantities from PALM output NetCDF4 profiles. More settings available in settings.py file.''')
parser.add_argument("-f", "--file", type=str, default=None,
                    help="Name of the input NetCDF4 profile file.")
parser.add_argument("-c", "--compare", type=str, default=None,
                    help="Print out a comparsion of the results to the results of given NetCDF4 profile file.")
parser.add_argument("-v", "--variables", type=str, nargs='+',
                    default=['ws', 'var', 'flux', 'tke'], help="Variables to be plotted.")
args = parser.parse_args()

#==========================================================#

# Read in the dataset
ds = openDataSet(args.file)

# Skip time from beginning
t_inds, = np.where(np.logical_and(ds.variables['time'][:] >= timespan[0], ds.variables['time'][:] <= timespan[1]))

# Generate table header
pr_table_header = str(pr_heights_table)
pr_table_header = [`s` + ' m' for s in pr_heights_table]
pr_table_header.insert(0, "Variable")
print(pr_table_header)

if (not(args.compare)):
  print("#=============== Time-averaged variables for dataset {} ===============#".format(args.file))
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
